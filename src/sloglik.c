#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sloglik.h"
#include "util.h"


static void mean_init(struct recv_sloglik_mean *mean,
		      const struct model *model,
		      ssize_t isend)
{
	assert(mean);

	const struct recv_model *rm = model_recv_model(model, isend);
	const struct vector *mean0 = recv_model_mean0(rm);
	const struct design *design = model_design(model);
	ssize_t dyn_dim = design_recv_dyn_dim(design);
	
	mean->mean0 = mean0;
	mean->gamma = 1.0;
	vector_init(&mean->dp, 0);
	vector_init(&mean->mean_dx, dyn_dim);
}

static void imat_init(struct recv_sloglik_imat *imat,
		      const struct model *model,
		      ssize_t isend)
{
	assert(imat);
	
	const struct recv_model *rm = model_recv_model(model, isend);
	const struct matrix *imat0 = recv_model_imat0(rm);
	const struct design *design = model_design(model);
	ssize_t dyn_dim = design_recv_dyn_dim(design);
	
	imat->imat0 = imat0;
	imat->gamma2 = 0.0;
	vector_init(&imat->gamma_dp, 0);
	vector_init(&imat->gamma_mean_dx, dyn_dim);
	matrix_init(&imat->dx_p, dyn_dim, 0);
	matrix_init(&imat->mean_dx_dp, dyn_dim, 0);
	matrix_init(&imat->dp2, 0, 0);
	matrix_init(&imat->var_dx, dyn_dim, dyn_dim);	
}

static void mean_deinit(struct recv_sloglik_mean *mean)
{
	assert(mean);

	vector_deinit(&mean->mean_dx);
	vector_deinit(&mean->dp);
}

static void imat_deinit(struct recv_sloglik_imat *imat)
{
	assert(imat);

	matrix_deinit(&imat->var_dx);
	matrix_deinit(&imat->dp2);
	matrix_deinit(&imat->mean_dx_dp);
	matrix_deinit(&imat->dx_p);
	vector_deinit(&imat->gamma_mean_dx);
	vector_deinit(&imat->gamma_dp);	
}

static void mean_clear(struct recv_sloglik_mean *mean)
{
	assert(mean);
	mean->gamma = 1.0;
	vector_reinit(&mean->dp, 0);
	vector_fill(&mean->mean_dx, 0.0);
}

static void imat_clear(struct recv_sloglik_imat *imat)
{
	assert(imat);
	imat->gamma2 = 0.0;
	vector_reinit(&imat->gamma_dp, 0);
	vector_fill(&imat->gamma_mean_dx, 0.0);
	matrix_reinit(&imat->dx_p, matrix_nrow(&imat->dx_p), 0);
	matrix_reinit(&imat->mean_dx_dp, matrix_nrow(&imat->mean_dx_dp), 0);
	matrix_reinit(&imat->dp2, 0, 0);
	matrix_fill(&imat->var_dx, 0.0);
}

static ssize_t mean_active_count(const struct recv_sloglik_mean *mean)
{
	assert(mean);
	return vector_dim(&mean->dp);
}

static ssize_t imat_active_count(const struct recv_sloglik_imat *imat)
{
	assert(imat);
	return vector_dim(&imat->gamma_dp);
}

static void mean_insert_active(struct recv_sloglik_mean *mean, ssize_t i)
{
	assert(mean);
	assert(0 <= i && i <= mean_active_count(mean));
	
	ssize_t n0 = mean_active_count(mean);
	ssize_t n1 = n0 + 1;
	struct vector dp, src, dst;
	vector_init(&dp, n1);
	
	src = vector_slice(&mean->dp, 0, i);
	dst = vector_slice(&dp, 0, i);
	vector_assign_copy(&dst, &src);
	
	src = vector_slice(&mean->dp, i, n0 - i);
	dst = vector_slice(&dp, i + 1, n0 - i);
	vector_assign_copy(&dst, &src);
	
	vector_deinit(&mean->dp);
	mean->dp = dp;
	
	assert(mean_active_count(mean) == n1);
}

static void imat_insert_active(struct recv_sloglik_imat *imat, ssize_t i)
{
	assert(imat);
	assert(0 <= i && i <= imat_active_count(imat));
	
	ssize_t n0 = imat_active_count(imat);
	ssize_t n1 = n0 + 1;
	struct vector vdst, vsrc;
	struct matrix dst, src;

	/* gamma_dp */
	struct vector gamma_dp;
	vector_init(&gamma_dp, n1);
	
	vsrc = vector_slice(&imat->gamma_dp, 0, i);
	vdst = vector_slice(&gamma_dp, 0, i);
	vector_assign_copy(&vdst, &vsrc);

	vsrc = vector_slice(&imat->gamma_dp, i, n0 - i);
	vdst = vector_slice(&gamma_dp, i + 1, n0 - i);
	vector_assign_copy(&vdst, &vsrc);
	
	vector_deinit(&imat->gamma_dp);
	imat->gamma_dp = gamma_dp;

	
	/* dx_p */
	struct matrix dx_p;
	matrix_init(&dx_p, matrix_nrow(&imat->dx_p), n1);
	
	src = matrix_slice_cols(&imat->dx_p, 0, i);
	dst = matrix_slice_cols(&dx_p, 0, i);
	matrix_assign_copy(&dst, &src);
	
	src = matrix_slice_cols(&imat->dx_p, i, n0 - i);
	dst = matrix_slice_cols(&dx_p, i + 1, n0 - i);
	matrix_assign_copy(&dst, &src);
	
	matrix_deinit(&imat->dx_p);
	imat->dx_p = dx_p;
	

	/* mean_dx_dp */
	struct matrix mean_dx_dp;
	matrix_init(&mean_dx_dp, matrix_nrow(&imat->mean_dx_dp), n1);
	
	src = matrix_slice_cols(&imat->mean_dx_dp, 0, i);
	dst = matrix_slice_cols(&mean_dx_dp, 0, i);
	matrix_assign_copy(&dst, &src);
	
	src = matrix_slice_cols(&imat->mean_dx_dp, i, n0 - i);
	dst = matrix_slice_cols(&mean_dx_dp, i + 1, n0 - i);
	matrix_assign_copy(&dst, &src);
	
	matrix_deinit(&imat->mean_dx_dp);
	imat->mean_dx_dp = mean_dx_dp;

	
	/* dp2 */
	struct matrix dp2;
	matrix_init(&dp2, n1, n1);

	dst = matrix_slice(&dp2, 0, 0, i, i);
	src = matrix_slice(&imat->dp2, 0, 0, i, i);
	matrix_assign_copy(&dst, &src);
	
	dst = matrix_slice(&dp2, i + 1, 0, n0 - i, i);
	src = matrix_slice(&imat->dp2, i, 0, n0 - i, i);
	matrix_assign_copy(&dst, &src);
	
	dst = matrix_slice(&dp2, 0, i + 1, i, n0 - i);
	src = matrix_slice(&imat->dp2, 0, i, i, n0 - i);
	matrix_assign_copy(&dst, &src);
	
	dst = matrix_slice(&dp2, i + 1, i + 1, n0 - i, n0 - i);
	src = matrix_slice(&imat->dp2, i, i, n0 - i, n0 - i);
	matrix_assign_copy(&dst, &src);
	matrix_deinit(&imat->dp2);
	imat->dp2 = dp2;
}

static void mean_set(struct recv_sloglik_mean *mean,
		     const struct recv_model *model,
		     const struct frame *f)
{
	assert(mean);
	assert(model);
	assert(f);
	
	const ssize_t isend = model->isend;
	const struct array *active = &model->active;
	const ssize_t n = array_count(active);
	ssize_t i;

	mean->mean0 = recv_model_mean0(model);
	
	const double gamma = model->gamma;
	mean->gamma = gamma;
	
	if (vector_dim(&mean->dp) != n) {
		vector_reinit(&mean->dp, n);
	}
	vector_fill(&mean->mean_dx, 0.0);

	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);
		
		/* mean */		
		double p = recv_model_prob(model, jrecv);
		double p0 = recv_model_prob0(model, jrecv);
		double dp = p - gamma * p0;
		vector_set_item(&mean->dp, i, dp);
		
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		vector_axpy(p, dx, &mean->mean_dx);
	}
}

static void imat_set(struct recv_sloglik_imat *imat,
		     const struct recv_model *model,
		     const struct frame *f,
		     const struct recv_sloglik_mean *mean)
{
	assert(imat);
	assert(model);
	assert(f);
	
	const ssize_t isend = model->isend;
	const struct array *active = &model->active;
	const ssize_t n = array_count(active);
	ssize_t i;
	
	imat->imat0 = recv_model_imat0(model);
	
	const double gamma = model->gamma;
	
	
	if (vector_dim(&imat->gamma_dp) != n) {
		vector_reinit(&imat->gamma_dp, n);
		matrix_reinit(&imat->dx_p, matrix_nrow(&imat->dx_p), n);
		matrix_reinit(&imat->mean_dx_dp, matrix_nrow(&imat->mean_dx_dp), n);		
		matrix_reinit(&imat->dp2, n, n);
	}
	matrix_fill(&imat->dp2, 0.0);	
	matrix_fill(&imat->var_dx, 0.0);	
	struct vector y;
	double ptot = 0.0;
	
	/* gamma2 */
	imat->gamma2 = gamma * (1 - gamma);
	
	/* gamma_mean_dx */
	vector_assign_copy(&imat->gamma_mean_dx, &mean->mean_dx);
	vector_scale(&imat->gamma_mean_dx, gamma);
	
	/* dp2 */
	matrix_update1(&imat->dp2, -1.0, &mean->dp, &mean->dp);
	
	vector_init(&y, vector_dim(&mean->mean_dx));
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);
		double dp = vector_item(&mean->dp, i);
		double p = recv_model_prob(model, jrecv);
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		
		/* gamma_dp */
		vector_set_item(&imat->gamma_dp, i, gamma * dp);
		
		/* mean_dx_dp */
		struct vector mean_dx_dp_j = matrix_col(&imat->mean_dx_dp, i);
		vector_assign_copy(&mean_dx_dp_j, &mean->mean_dx);
		vector_scale(&mean_dx_dp_j, dp);
		
		/* dp2 */
		*matrix_item_ptr(&imat->dp2, i, i) += dp;
		
		/* dx_p */
		struct vector dx_p_j = matrix_col(&imat->dx_p, i);
		vector_assign_copy(&dx_p_j, dx);
		vector_scale(&dx_p_j, p);

		/* var_dx */
		vector_assign_copy(&y, dx);
		vector_sub(&y, &mean->mean_dx);
		matrix_update1(&imat->var_dx, p, &y, &y);
		ptot += p;
	}
	/* var_dx */
	matrix_update1(&imat->var_dx, 1 - ptot, &mean->mean_dx, &mean->mean_dx);
	vector_deinit(&y);
}


static void vector_mean_update(double scale, const struct vector *val,
			       struct vector *mean, struct vector *work)
{
	struct vector diff = vector_slice(work, 0, vector_dim(mean));
	vector_assign_copy(&diff, val);
	vector_sub(&diff, mean);
	vector_axpy(scale, &diff, mean);
}

static void matrix_mean_update(double scale, const struct matrix *val,
			       struct matrix *mean, struct vector *work)
{
	ssize_t m = matrix_nrow(mean);
	ssize_t n = matrix_ncol(mean);
	struct vector vwork = vector_slice(work, 0, m * n);
	struct matrix diff = matrix_make(&vwork, m, n);
	
	matrix_assign_copy(&diff, val);
	matrix_sub(&diff, mean);
	matrix_axpy(scale, &diff, mean);
}


static void mean_update(const struct array *active,
			ssize_t n0,
			const struct recv_sloglik_mean *mean0,
			ssize_t n1,
			struct recv_sloglik_mean *mean1)
{
	assert(active);
	assert(n0 >= 0);
	assert(mean0);
	assert(n1 >= 0);
	assert(mean1);
	assert(mean_active_count(mean0) == array_count(active));
	assert(mean_active_count(mean0) == mean_active_count(mean1));
	assert(mean0->mean0 == mean1->mean0);

	ssize_t ntot = n0 + n1;
	double scale = ((double)n0) / ntot;
	struct vector work;
	vector_init(&work, MAX(vector_dim(&mean0->dp), vector_dim(&mean0->mean_dx)));
	
	mean1->gamma += scale * (mean0->gamma - mean1->gamma);

	vector_mean_update(scale, &mean0->dp, &mean1->dp, &work);
	vector_mean_update(scale, &mean0->mean_dx, &mean1->mean_dx, &work);
	
	vector_deinit(&work);
}

static void imat_update(const struct array *active,
			ssize_t n0,
			const struct recv_sloglik_imat *imat0,
			ssize_t n1,
			struct recv_sloglik_imat *imat1)
{
	assert(active);
	assert(n0 >= 0);
	assert(imat0);
	assert(n1 >= 0);
	assert(imat1);
	assert(imat_active_count(imat0) == array_count(active));
	assert(imat_active_count(imat0) == imat_active_count(imat1));
	assert(imat0->imat0 == imat1->imat0);
	
	ssize_t ntot = n0 + n1;
	double scale = ((double)n0) / ntot;
	
	ssize_t nactive = vector_dim(&imat0->gamma_dp);
	ssize_t dyn_dim = vector_dim(&imat0->gamma_mean_dx);
	ssize_t nwork = MAX(nactive * nactive, dyn_dim * dyn_dim);
	
	struct vector work;
	vector_init(&work, nwork);
	
	imat1->gamma2 += scale * (imat0->gamma2 - imat1->gamma2);
	
	vector_mean_update(scale, &imat0->gamma_dp, &imat1->gamma_dp, &work);
	vector_mean_update(scale, &imat0->gamma_mean_dx, &imat1->gamma_mean_dx, &work);

	matrix_mean_update(scale, &imat0->dx_p, &imat1->dx_p, &work);
	matrix_mean_update(scale, &imat0->mean_dx_dp, &imat1->mean_dx_dp, &work);
	matrix_mean_update(scale, &imat0->dp2, &imat1->dp2, &work);
	matrix_mean_update(scale, &imat0->var_dx, &imat1->var_dx, &work);	
	
	vector_deinit(&work);
}


static void mean_axpy(double alpha,
		      const struct recv_sloglik_mean *mean,
		      const struct array *active,
		      const struct design *design,
		      ssize_t isend,
		      struct vector *y)
{
	const struct vector *mean0 = mean->mean0;
	double gamma = mean->gamma;
	vector_axpy(alpha * gamma, mean0, y);
		
	ssize_t i, n = array_count(active);
	struct svector dp;
	svector_init(&dp, design_recv_count(design));
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);
		double val = vector_item(&mean->dp, i);
		svector_set_item(&dp, jrecv, val);
	}
	
	design_recv_muls0(alpha, TRANS_TRANS, design, isend, &dp, 1.0, y);
	svector_deinit(&dp);
	
	const struct vector *mean_dx = &mean->mean_dx;
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);
	struct vector ysub = vector_slice(y, off, dim);
	vector_axpy(alpha, mean_dx, &ysub);
}

static void imat_axpy(double alpha,
		      const struct recv_sloglik_imat *imat,
		      const struct recv_sloglik_mean *mean,
		      const struct array *active,
		      const struct design *design,
		      ssize_t isend,
		      struct matrix *y)
{
	const ssize_t dim = design_recv_dim(design);
	const ssize_t nrecv = design_recv_count(design);
	const struct vector *mean0 = mean->mean0;
	const struct matrix *var0 = imat->imat0;
	const double gamma = mean->gamma;
	const double gamma2 = imat->gamma2;
	
	const ssize_t n = array_count(active);
	ssize_t i, j, k;
	struct svector gamma_dp;
	svector_init(&gamma_dp, nrecv);
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);
		double val = vector_item(&imat->gamma_dp, i);
		svector_set_item(&gamma_dp, jrecv, val);
	}
	
	struct vector gamma_x0_dp;
	vector_init(&gamma_x0_dp, dim);
	design_recv_muls0(1.0, TRANS_TRANS, design, isend, &gamma_dp, 0.0, &gamma_x0_dp);
	
	struct matrix x0_dp2;
	struct svector dp2_j;
	
	matrix_init(&x0_dp2, dim, n);
	svector_init(&dp2_j, nrecv);
	for (j = 0; j < n; j++) {
		svector_clear(&dp2_j);
		for (i = 0; i < n; i++) {
			ssize_t jrecv = *(ssize_t *)array_item(active, i);
			double val = matrix_item(&imat->dp2, i, j);
			svector_set_item(&dp2_j, jrecv, val);
		}
		struct vector dst = matrix_col(&x0_dp2, j);
		design_recv_muls0(1.0, TRANS_TRANS, design, isend, &dp2_j, 0.0, &dst);
	}
	
	matrix_axpy(alpha * gamma, var0, y);
	matrix_update1(y, alpha * gamma2, mean0, mean0);
	matrix_update1(y, -alpha, mean0, &gamma_x0_dp);
	matrix_update1(y, -alpha, &gamma_x0_dp, mean0);
	
	struct svector x0_dp2_k;
	svector_init(&x0_dp2_k, nrecv);
	for (k = 0; k < dim; k++) {
		svector_clear(&x0_dp2_k);
		for (j = 0; j < n; j++) {
			ssize_t jrecv = *(ssize_t *)array_item(active, j);
			double val = matrix_item(&x0_dp2, k, j);
			svector_set_item(&x0_dp2_k, jrecv, val);
		}
		struct vector dst = matrix_col(y, k);
		design_recv_muls0(alpha, TRANS_TRANS, design, isend, &x0_dp2_k, 1.0, &dst);
	}
	
	ssize_t dyn_off = design_recv_dyn_index(design);
	ssize_t dyn_dim = design_recv_dyn_dim(design);	
	struct matrix y_1 = matrix_slice(y, 0, dyn_off, dim, dyn_dim);
	struct matrix y1_ = matrix_slice(y, dyn_off, 0, dyn_dim, dim);
	struct matrix y11 = matrix_slice(y, dyn_off, dyn_off, dyn_dim, dyn_dim);
	
	matrix_axpy(alpha, &imat->var_dx, &y11);
	
	struct vector x0_j;	
	struct svector e_j;
	
	vector_init(&x0_j, dim);
	svector_init(&e_j, design_recv_count(design));
	
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);
		
		svector_set_basis(&e_j, jrecv);
		design_recv_muls0(1.0, TRANS_TRANS, design, isend, &e_j, 0.0, &x0_j);
		
		const struct vector dx_p_j = matrix_col(&imat->dx_p, i);
		matrix_update1(&y_1, alpha, &x0_j, &dx_p_j);
		matrix_update1(&y1_, alpha, &dx_p_j, &x0_j);
		
		struct vector mean_dx_dp_j = matrix_col(&imat->mean_dx_dp, i);
		matrix_update1(&y_1, -alpha, &x0_j, &mean_dx_dp_j);
		matrix_update1(&y1_, -alpha, &mean_dx_dp_j, &x0_j);

	}
	
	matrix_update1(&y_1, -alpha, mean0, &imat->gamma_mean_dx);
	matrix_update1(&y1_, -alpha, &imat->gamma_mean_dx, mean0);
	
	vector_deinit(&x0_j);
	svector_deinit(&e_j);
	svector_deinit(&x0_dp2_k);
	svector_deinit(&dp2_j);
	matrix_deinit(&x0_dp2);
	vector_deinit(&gamma_x0_dp);
	svector_deinit(&gamma_dp);
}
		      

void recv_sloglik_init(struct recv_sloglik *ll, const struct model *model, ssize_t isend)
{
	assert(ll);
	assert(model);
	assert(0 <= isend && isend < design_send_count(model_design(model)));

	ssize_t n = model_receiver_count(model);
	ssize_t p = model_dim(model);
	const struct design *design = model_design(model);
	ssize_t dyn_dim = design_recv_dyn_dim(design);
	
	ll->model = (struct model *)model;
	ll->isend = isend;
	
	ll->n = 0;
	ll->n_last = 0;
	ll->dev_avg = 0.0;
	ll->dev_last = 0.0;
	array_init(&ll->active, sizeof(ssize_t));
	mean_init(&ll->mean_last, model, isend);
	mean_init(&ll->mean_avg, model, isend);
	imat_init(&ll->imat_last, model, isend);
	imat_init(&ll->imat_avg, model, isend);
	
	/* deprecated */
	ll->f = 0.0;
	vector_init(&ll->grad, p);
	ll->grad_cached = false;
	
	ll->nsend = 0;
	svector_init(&ll->nrecv, n);
	vector_init(&ll->dxobs, dyn_dim);
	

	ll->gamma = 0.0;
	svector_init(&ll->dp, n);
	vector_init(&ll->dxbar, dyn_dim);
}

void recv_sloglik_deinit(struct recv_sloglik *ll)
{
	assert(ll);

	imat_deinit(&ll->imat_avg);
	imat_deinit(&ll->imat_last);
	mean_deinit(&ll->mean_avg);
	mean_deinit(&ll->mean_last);
	array_deinit(&ll->active);
	vector_deinit(&ll->dxbar);
	svector_deinit(&ll->dp);
	vector_deinit(&ll->dxobs);
	svector_deinit(&ll->nrecv);
	vector_deinit(&ll->grad);
}

static void recv_sloglik_update_active(struct recv_sloglik *ll,
				       const struct array *active)
{
	assert(ll);
	assert(active);
	
	if (array_count(active) > array_count(&ll->active)) {
		const ssize_t n0 = array_count(&ll->active);
		const ssize_t n1 = array_count(active);
		const ssize_t *begin0 = array_to_ptr(&ll->active);
		const ssize_t *end0 = begin0 + n0;
		const ssize_t *begin1 = array_to_ptr(active);
		const ssize_t *end1 = begin1 + n1;
		const ssize_t *i0, *i1;
		
		for (i0 = begin0, i1 = begin1; i1 < end1; i1++) {
			if (i0 < end0 && *i0 == *i1) {
				i0++;
			} else {
				assert(i0 == end0 || *i1 < *i0);
				mean_insert_active(&ll->mean_avg, i1 - begin1);
				imat_insert_active(&ll->imat_avg, i1 - begin1);				
			}
		}
		assert(i0 == end0);
		assert(mean_active_count(&ll->mean_avg) == n1);
		assert(imat_active_count(&ll->imat_avg) == n1);		

		array_assign_copy(&ll->active, active);
	}
}

void recv_sloglik_add(struct recv_sloglik *ll,
		     const struct frame *f, ssize_t *jrecv, ssize_t n)
{
	ssize_t isend = ll->isend;
	const struct recv_model *model = model_recv_model(ll->model, isend);
	ssize_t nreceiver = recv_model_count(model);

	double ntot = ll->n + n;
	double scale1 = n / ntot;
	//double scale0 = 1 - scale1;
	ssize_t i;

	recv_sloglik_update_active(ll, &model->active);
	
	ll->dev_last = 0.0;
	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);
		
		double lp = recv_model_logprob(model, jrecv[i]);
		ll->dev_last += -2 * lp;
	}
	ll->dev_avg += scale1 * (ll->dev_last / n - ll->dev_avg);

	mean_set(&ll->mean_last, model, f);
	mean_update(&model->active, n, &ll->mean_last, ll->n, &ll->mean_avg);
	
	imat_set(&ll->imat_last, model, f, &ll->mean_last);
	imat_update(&model->active, n, &ll->imat_last, ll->n, &ll->imat_avg);	
		    
	ll->n_last = n;
	ll->n += n;
		    

		    
	/*
	 struct svector *wt = svector_alloc(nreceiver);
	 double lpbar = 0.0;


	
	ll->grad_cached = false;

	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);

		double lp = recv_model_logprob(model, jrecv[i]);
		lpbar += (lp - lpbar) / (i + 1);

		// *svector_item_ptr(wt, jrecv[i]) += 1.0;

		// update number of receives
		*svector_item_ptr(&ll->nrecv, jrecv[i]) += 1.0;
	}

	// update log likelihood
	ll->f += scale1 * ((-lpbar) - ll->f);

	// update observed variable diffs
	frame_recv_dmuls(scale1 / n, TRANS_TRANS, f, model->isend,
			 wt, scale0, &ll->dxobs);
	ll->gamma += scale1 * (model->gamma - ll->gamma);

	svector_scale(&ll->dp, scale0);
	
	ssize_t nactive = array_count(&model->active);
	for (i = 0; i < nactive; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(&model->active, i);
		double dp = vector_item(&model->dp, i);
		double *dst = svector_item_ptr(&ll->dp, jrecv);
		*dst += scale1 * dp;
	}
	//svector_axpys(scale1, &model->dp, &ll->dp);

	vector_scale(&ll->dxbar, scale0);
	vector_axpy(scale1, &model->mean_dx, &ll->dxbar);

	// update number of sends
	ll->nsend += n;

	svector_free(wt);*/
}

double recv_sloglik_value(const struct recv_sloglik *ll)
{
	assert(ll);
	return (-ll->nsend) * ll->f;
}

/*
 *         g = [ sum{gamma[t,i]} * xbar[0,i]
 *               + ( X[0,i])^T * sum{dP[t,i]}
 *               + sum{dxbar[t,i]} ]
 *             -
 *             [ (X[0,i])^T n[i] + sum{dx[t,i,j]} ]
 */
static void
recv_sloglik_axpy_grad_nocache(double alpha, const struct recv_sloglik *ll, struct vector *y)
{
	double scale = (-ll->nsend) * alpha;
	const struct recv_model *rm = model_recv_model(ll->model, ll->isend);
	const struct vector *xbar0 = recv_model_mean0(rm);

	const struct design *design = model_design(ll->model);
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);
	struct vector ysub = vector_slice(y, off, dim);


	// sum{gamma[t,i]} * xbar[0,i]
	vector_axpy(scale * ll->gamma, xbar0, y);

	// (X[0,i])^T * sum{dP[t,i]}
	design_recv_muls0(scale, TRANS_TRANS,
			  ll->model->design, ll->isend, &ll->dp, 1.0,
			  y);

	// sum{dxbar[t,i]}
	vector_axpy(scale, &ll->dxbar, &ysub);

	// - (X[0,i])^T n[i]
	design_recv_muls0(-scale / ll->nsend, TRANS_TRANS,
			  ll->model->design, ll->isend, &ll->nrecv,
			  1.0, y);

	// -sum{dx[t,i,j]}
	vector_axpy(-scale, &ll->dxobs, &ysub);
}

static void
recv_sloglik_cache_grad(struct recv_sloglik *ll)
{
	vector_fill(&ll->grad, 0.0);
	recv_sloglik_axpy_grad_nocache(1.0, ll, &ll->grad);
	ll->grad_cached = true;
}

void recv_sloglik_axpy_grad(double alpha, const struct recv_sloglik *ll, struct vector *y)
{
	assert(ll);
	assert(y);

	if (!ll->grad_cached) {
		recv_sloglik_cache_grad((struct recv_sloglik *)ll);
	}
	
	vector_axpy(alpha, &ll->grad, y);
}

ssize_t recv_sloglik_count(const struct recv_sloglik *sll)
{
	assert(sll);
	return sll->n;
}

double recv_sloglik_avg_dev(const struct recv_sloglik *sll)
{
	assert(sll);
	return sll->dev_avg;
}

double recv_sloglik_last_dev(const struct recv_sloglik *sll)
{
	assert(sll);
	return sll->dev_last;
}

void recv_sloglik_axpy_avg_mean(double alpha, const struct recv_sloglik *sll, struct vector *y)
{
	assert(sll);
	assert(y);
	
	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;
	
	mean_axpy(alpha, &sll->mean_avg, &sll->active, design, isend, y);
}

void recv_sloglik_axpy_last_mean(double alpha, const struct recv_sloglik *sll, struct vector *y)
{
	assert(sll);
	assert(y);
	
	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;
	
	mean_axpy(sll->n_last * alpha, &sll->mean_last, &sll->active, design, isend, y);
}

void recv_sloglik_axpy_avg_imat(double alpha, const struct recv_sloglik *sll, struct matrix *y)
{
	assert(sll);
	assert(y);

	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;
	
	imat_axpy(alpha, &sll->imat_avg, &sll->mean_avg, &sll->active, design, isend, y);
}

void recv_sloglik_axpy_last_imat(double alpha, const struct recv_sloglik *sll, struct matrix *y)
{
	assert(sll);
	assert(y);

	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;
	
	imat_axpy(sll->n_last * alpha, &sll->imat_last, &sll->mean_last, &sll->active, design, isend, y);
}

