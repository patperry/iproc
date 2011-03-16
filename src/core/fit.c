
#include <assert.h>
#include <math.h>
#include <string.h>
#include "config.h"
#include "memory.h"
#include "fit.h"

int dcsrch_(double *stp, double *f, double *g, 
            double *ftol, double *gtol, double *xtol, char *task, 
            double *stpmin, double *stpmax, f77int *isave, double *
            dsave, f77int task_len);



static void
eval_objective (iproc_loglik *loglik,
                double        penalty,
                double       *valuep,
                iproc_vector *grad)
{
    iproc_model *model = loglik->model;
    iproc_vector *coefs = model->coefs;
    
    double ll_value = iproc_loglik_value(loglik);
    iproc_vector *ll_grad = iproc_loglik_grad(loglik);
    
    double norm = iproc_vector_norm(coefs);
    double norm2 = norm * norm;
    double value = -ll_value + 0.5 * penalty * norm2;
    *valuep = value;
    
    iproc_vector_copy(grad, coefs);
    iproc_vector_scale(grad, penalty);
    iproc_vector_acc(grad, -1.0, ll_grad);
}


static int
iproc_fit_init (iproc_fit *fit)
{
    int64_t dim = iproc_model_dim(fit->model);

    fit->inv_hess = NULL;
    fit->grad = iproc_vector_new(dim);
    fit->grad0 = iproc_vector_new(dim);
    fit->search_dir = iproc_vector_new(dim);
    fit->loglik = iproc_loglik_new(fit->model, fit->messages);
    
    if (!(fit->inv_hess && fit->grad && fit->grad0 && fit->search_dir
          && fit->loglik))
        return 0;
    
    eval_objective(fit->loglik, fit->penalty, &fit->value, fit->grad);
    fit->value0 = NAN;
    fit->step = 1.0;
    
    return 1;
}

iproc_fit *
iproc_fit_new (iproc_model    *model0,
               iproc_messages *messages,
               double          penalty)
{
    assert(model0);
    assert(messages);
    assert(penalty >= 0.0);
    assert(iproc_messages_max_from(messages) < iproc_model_nsender(model0));
    assert(iproc_messages_max_to(messages) < iproc_model_nreceiver(model0));
    
    iproc_fit *fit = iproc_malloc(sizeof(*fit));
    if (!fit)
        return NULL;


    fit->model = iproc_model_ref(model0);
    fit->messages = iproc_messages_ref(messages);
    fit->penalty = penalty;
    fit->loglik = NULL;


    if (!iproc_fit_init(fit)) {
        iproc_fit_free(fit);
        fit = NULL;
    }
    
    return fit;
}

void
iproc_fit_free (iproc_fit *fit)
{
    if (fit) {
        iproc_loglik_unref(fit->loglik);
        iproc_vector_unref(fit->search_dir);
        iproc_vector_unref(fit->grad0);
        iproc_vector_unref(fit->grad);
        iproc_matrix_unref(fit->inv_hess);
        iproc_messages_unref(fit->messages);        
        iproc_model_unref(fit->model);
        iproc_free(fit);
    }
}

static void
linesearch (iproc_fit *fit)
{
    assert(fit);

    iproc_model *model = fit->model;
    iproc_messages *messages = fit->messages;
    double penalty = fit->penalty;
    iproc_design *design = model->design;
    int has_loops = model->has_loops;
    iproc_loglik *loglik = fit->loglik;
    
    iproc_vector *x0 = iproc_vector_ref(model->coefs);
    iproc_vector *x = iproc_vector_new_copy(x0);
    iproc_vector *search_dir = fit->search_dir;
    iproc_vector *grad0 = fit->grad; /* swap grad0 and grad */
    iproc_vector *grad = fit->grad0;
    double value0 = fit->value;
    double value;
    
    double f = value0;
    double g = iproc_vector_dot(grad0, search_dir);
    double ftol = 1e-4;
    double gtol = 0.9;
    double xtol = 0.1;
    f77int task_len = 60;
    char task[task_len + 1];
    double stpmin = 0;
    double stpmax = f / (gtol * g);
    double stp = g == 0 ? 1.0 : IPROC_MIN(1.0, stpmax);
    f77int isave[2];
    double dsave[13];
    
    if (g == 0)
        goto cleanup;
    
    strcpy(task, "START");
    dcsrch_(&stp, &f, &g, &ftol, &gtol, &xtol, task, &stpmin, &stpmax, isave,
            dsave, task_len);
    
    while (strncmp(task, "FG", strlen("FG")) == 0) {
        iproc_loglik_unref(loglik);        
        iproc_model_unref(model);

        /* Take a step and create a new model */
        iproc_vector_copy(x, x0);
        iproc_vector_acc(x, stp, search_dir);
        model = iproc_model_new(design, x, has_loops);
        
        /* Update the loglik, value, and gradient */
        loglik = iproc_loglik_new(model, messages);
        eval_objective(loglik, penalty, &value, grad);
        
        f = value;
        g = iproc_vector_dot(grad, search_dir);
        
        dcsrch_(&stp, &f, &g, &ftol, &gtol, &xtol, task, &stpmin, &stpmax,
                isave, dsave, task_len);

    }

    const char * xtol_msg = "WARNING: XTOL TEST SATISFIED";
    assert(!(strncmp(task, xtol_msg, strlen(xtol_msg))
             && strncmp(task, "CONV", strlen("CONV"))));

cleanup:
    
    iproc_vector_unref(x);
    iproc_vector_unref(x0);
    
    fit->step = stp;
    fit->value = value;
    fit->value0 = value0;
    fit->grad = grad;
    fit->grad0 = grad0;
    fit->model = model;
    fit->loglik = loglik;
}

static void
update_hess (iproc_fit *fit)
{
    iproc_vector *s = iproc_vector_new_copy(fit->search_dir);
    iproc_vector_scale(s, fit->step);

    iproc_vector *y = iproc_vector_new_copy(fit->grad);
    iproc_vector_acc(y, -1.0, fit->grad0);

    double s_y = iproc_vector_dot(s, y);
    
    iproc_matrix *H = fit->inv_hess;
    
    if (H == NULL) {
        double scale = 1.0;
        
        if (s_y > 0) {
            double s_s = iproc_vector_dot(s, s);
            scale = s_y / s_s;
        }

        int64_t i, n = iproc_vector_dim(y);
        H = iproc_matrix_new(n, n);
        iproc_matrix_set_all(H, 0.0);
        
        for (i = 0; i < n; i++) {
            iproc_matrix_set(H, i, i, scale);
        }
        
        fit->inv_hess = H;
    } else {
        iproc_vector *H_y = iproc_vector_new(iproc_vector_dim(y));
        iproc_matrix_mul(1.0, IPROC_TRANS_NOTRANS, H, y, 0.0, H_y);
    
        double y_H_y = iproc_vector_dot(H_y, y);
        double scale1 = (1.0 + (y_H_y / s_y)) / s_y;
        double rho = 1.0 / s_y;
    
        iproc_matrix_update1(H, scale1, s, s);
        iproc_matrix_update1(H, -rho, H_y, s);
        iproc_matrix_update1(H, -rho, s, H_y);
    }
    
    iproc_vector_unref(y);
    iproc_vector_unref(s);
}

static void
update_searchdir (iproc_fit *fit)
{
    assert(fit);
    
    iproc_vector *s = fit->search_dir;
    
    if (fit->inv_hess) {
        iproc_matrix_mul(-1.0, IPROC_TRANS_NOTRANS, fit->inv_hess, fit->grad,
                         +0.0, s);
    } else {
        iproc_vector_copy(s, fit->grad);
        double scale = iproc_vector_norm(s);
        int64_t i, n = iproc_vector_dim(s);

        if (scale != 0) {
            for (i = 0; i < n; i++) {
                double s_i = iproc_vector_get(s, i);
                iproc_vector_set(s, i, s_i / scale);
            }
        }
    }
}

void
iproc_fit_step (iproc_fit *fit)
{
    update_searchdir(fit);    
    linesearch(fit);
    update_hess(fit);
}

bool
iproc_fit_converged (iproc_fit *fit,
                     double     abs_tol,
                     double     rel_tol)
{
    double f = fit->value;
    double f0 = fit->value0;
    double step = fit->step;
    double g0 = iproc_vector_dot(fit->search_dir, fit->grad);
    bool converged = false;
    
    if (fabs(f - f0) <= abs_tol
        && step * fabs(g0) <= abs_tol) {
        /* abs_tol test satisfied */
        converged = true;
    } else if (fabs(f - f0) <= rel_tol * fabs(f0)
               && step * fabs(g0) <= rel_tol * fabs(f0)) {
        /* rel_tol test satisfied */
        converged = true;
    }
    
    return converged;
}
