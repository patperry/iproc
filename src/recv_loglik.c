#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "xalloc.h"
#include "recv_loglik.h"

static void cohort_init(struct recv_loglik_cohort *cll,
			const struct recv_model *m, ssize_t c);
static void cohort_deinit(struct recv_loglik_cohort *cll);
static void cohort_clear(struct recv_loglik_cohort *cll);

static void sender_init(struct recv_loglik_sender *ll,
			const struct recv_model *model, ssize_t isend);
static void sender_deinit(struct recv_loglik_sender *ll);

static void sender_add(struct recv_loglik_sender *ll,
		       const struct frame *f, const ssize_t *jrecv, ssize_t n);
static void sender_clear(struct recv_loglik_sender *ll);

static ssize_t sender_count(const struct recv_loglik_sender *sll);
static double sender_avg_dev(const struct recv_loglik_sender *sll);
static void sender_axpy_avg_mean(double alpha,
				 const struct recv_loglik_sender *sll,
				 struct vector *y);
static void sender_axpy_avg_score(double alpha,
				  const struct recv_loglik_sender *sll,
				  struct vector *y);
static void sender_axpy_avg_imat(double alpha,
				 const struct recv_loglik_sender *sll,
				 struct matrix *y);

static ssize_t sender_last_count(const struct recv_loglik_sender *sll);
static double sender_last_dev(const struct recv_loglik_sender *sll);
static void sender_axpy_last_mean(double alpha,
				  const struct recv_loglik_sender *sll,
				  struct vector *y);
static void sender_axpy_last_score(double alpha,
				   const struct recv_loglik_sender *sll,
				   struct vector *y);
static void sender_axpy_last_imat(double alpha,
				  const struct recv_loglik_sender *sll,
				  struct matrix *y);

static void score_init(struct recv_loglik_sender_score *score,
		       const struct recv_model *model, ssize_t isend)
{
	assert(score);

	ssize_t ic = recv_model_cohort(model, isend);
	const struct vector *mean0 = recv_model_mean0(model, ic);
	const struct design *design = recv_model_design(model);
	ssize_t nrecv = design_recv_count(design);
	ssize_t dyn_dim = design_recv_dyn_dim(design);

	score->mean0 = mean0;
	svector_init(&score->nrecv, nrecv);
	vector_init(&score->mean_obs_dx, dyn_dim);
	score->gamma = 1.0;
	vector_init(&score->dp, 0);
	vector_init(&score->mean_dx, dyn_dim);
}

static void imat_init(struct recv_loglik_sender_imat *imat,
		      const struct recv_model *model, ssize_t isend)
{
	assert(imat);

	ssize_t ic = recv_model_cohort(model, isend);
	const struct matrix *imat0 = recv_model_imat0(model, ic);
	const struct design *design = recv_model_design(model);
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

static void score_deinit(struct recv_loglik_sender_score *score)
{
	assert(score);

	vector_deinit(&score->mean_dx);
	vector_deinit(&score->dp);
	vector_deinit(&score->mean_obs_dx);
	svector_deinit(&score->nrecv);
}

static void imat_deinit(struct recv_loglik_sender_imat *imat)
{
	assert(imat);

	matrix_deinit(&imat->var_dx);
	matrix_deinit(&imat->dp2);
	matrix_deinit(&imat->mean_dx_dp);
	matrix_deinit(&imat->dx_p);
	vector_deinit(&imat->gamma_mean_dx);
	vector_deinit(&imat->gamma_dp);
}

static void score_clear(struct recv_loglik_sender_score *score)
{
	assert(score);
	svector_clear(&score->nrecv);
	vector_fill(&score->mean_obs_dx, 0.0);
	score->gamma = 1.0;
	vector_reinit(&score->dp, 0);
	vector_fill(&score->mean_dx, 0.0);
}

static void imat_clear(struct recv_loglik_sender_imat *imat)
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

#ifndef NDEBUG
static ssize_t score_active_count(const struct recv_loglik_sender_score *score)
{
	assert(score);
	return vector_dim(&score->dp);
}
#endif

#ifndef NDEBUG
static ssize_t imat_active_count(const struct recv_loglik_sender_imat *imat)
{
	assert(imat);
	return vector_dim(&imat->gamma_dp);
}
#endif

static void score_insert_active(struct recv_loglik_sender_score *score,
				ssize_t i)
{
	assert(score);
	assert(0 <= i && i <= score_active_count(score));

	vector_insert(&score->dp, i);
}

static void imat_insert_active(struct recv_loglik_sender_imat *imat, ssize_t i)
{
	assert(imat);
	assert(0 <= i && i <= imat_active_count(imat));

	vector_insert(&imat->gamma_dp, i);
	matrix_insert_col(&imat->dx_p, i);
	matrix_insert_col(&imat->mean_dx_dp, i);
	matrix_insert_row_col(&imat->dp2, i, i);
}

static void score_set_obs(struct recv_loglik_sender_score *score,
			  const struct frame *f,
			  ssize_t isend, const ssize_t *jrecv, const ssize_t n)
{
	ssize_t i;
	svector_clear(&score->nrecv);

	for (i = 0; i < n; i++) {
		*svector_item_ptr(&score->nrecv, jrecv[i]) += 1.0;
	}

	frame_recv_dmuls(1.0 / n, BLAS_TRANS, f, isend,
			 &score->nrecv, 0.0, &score->mean_obs_dx);
}

static void score_set_mean(struct recv_loglik_sender_score *score,
			   const struct frame *f,
			   const struct recv_model *model, ssize_t isend)
{
	assert(score);
	assert(model);
	assert(f);

	ssize_t *active, i, n;
	ssize_t ic = recv_model_cohort(model, isend);
	recv_model_get_active(model, isend, &active, &n);

	score->mean0 = recv_model_mean0(model, ic);

	const double gamma = recv_model_invgrow(model, isend);
	score->gamma = gamma;

	if (vector_dim(&score->dp) != n) {
		vector_reinit(&score->dp, n);
	}
	vector_fill(&score->mean_dx, 0.0);

	for (i = 0; i < n; i++) {
		ssize_t jrecv = active[i];

		/* mean */
		double p = recv_model_prob(model, isend, jrecv);
		double p0 = recv_model_prob0(model, ic, jrecv);
		double dp = p - gamma * p0;
		vector_set_item(&score->dp, i, dp);

		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		vector_axpy(p, dx, &score->mean_dx);
	}
}

static void score_set(struct recv_loglik_sender_score *score,
		      const struct frame *f,
		      ssize_t isend, const ssize_t *jrecv, ssize_t n,
		      const struct recv_model *model)
{
	score_set_obs(score, f, isend, jrecv, n);
	score_set_mean(score, f, model, isend);
}

static void imat_set(struct recv_loglik_sender_imat *imat,
		     const struct frame *f,
		     const struct recv_model *model,
		     ssize_t isend,
		     const struct recv_loglik_sender_score *score)
{
	assert(imat);
	assert(model);
	assert(f);

	ssize_t ic = recv_model_cohort(model, isend);
	ssize_t *active, i, n;
	recv_model_get_active(model, isend, &active, &n);

	imat->imat0 = recv_model_imat0(model, ic);

	const double gamma = recv_model_invgrow(model, isend);

	if (vector_dim(&imat->gamma_dp) != n) {
		vector_reinit(&imat->gamma_dp, n);
		matrix_reinit(&imat->dx_p, matrix_nrow(&imat->dx_p), n);
		matrix_reinit(&imat->mean_dx_dp, matrix_nrow(&imat->mean_dx_dp),
			      n);
		matrix_reinit(&imat->dp2, n, n);
	}
	matrix_fill(&imat->dp2, 0.0);
	matrix_fill(&imat->var_dx, 0.0);
	struct vector y;
	double ptot = 0.0;

	/* gamma2 */
	imat->gamma2 = gamma * (1 - gamma);

	/* gamma_mean_dx */
	vector_assign_copy(&imat->gamma_mean_dx, &score->mean_dx);
	vector_scale(&imat->gamma_mean_dx, gamma);

	/* dp2 */
	matrix_update1(&imat->dp2, -1.0, &score->dp, &score->dp);

	vector_init(&y, vector_dim(&score->mean_dx));
	for (i = 0; i < n; i++) {
		ssize_t jrecv = active[i];
		double dp = vector_item(&score->dp, i);
		double p = recv_model_prob(model, isend, jrecv);
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);

		/* gamma_dp */
		vector_set_item(&imat->gamma_dp, i, gamma * dp);

		/* mean_dx_dp */
		struct vector mean_dx_dp_j = matrix_col(&imat->mean_dx_dp, i);
		vector_assign_copy(&mean_dx_dp_j, &score->mean_dx);
		vector_scale(&mean_dx_dp_j, dp);

		/* dp2 */
		*matrix_item_ptr(&imat->dp2, i, i) += dp;

		/* dx_p */
		struct vector dx_p_j = matrix_col(&imat->dx_p, i);
		vector_assign_copy(&dx_p_j, dx);
		vector_scale(&dx_p_j, p);

		/* var_dx */
		vector_assign_copy(&y, dx);
		vector_sub(&y, &score->mean_dx);
		matrix_update1(&imat->var_dx, p, &y, &y);
		ptot += p;
	}
	/* var_dx */
	matrix_update1(&imat->var_dx, 1 - ptot, &score->mean_dx,
		       &score->mean_dx);
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

	matrix_assign_copy(&diff, BLAS_NOTRANS, val);
	matrix_sub(&diff, mean);
	matrix_axpy(scale, &diff, mean);
}

static void score_update(ssize_t n0,
			 const struct recv_loglik_sender_score *score0,
			 ssize_t n1, struct recv_loglik_sender_score *score1)
{
	assert(n0 >= 0);
	assert(score0);
	assert(n1 >= 0);
	assert(score1);
	assert(score_active_count(score0) == score_active_count(score1));
	assert(score0->mean0 == score1->mean0);

	ssize_t ntot = n0 + n1;
	double scale = ((double)n0) / ntot;
	struct vector work;
	vector_init(&work,
		    MAX(vector_dim(&score0->dp), vector_dim(&score0->mean_dx)));

	svector_axpys(1.0, &score0->nrecv, &score1->nrecv);
	vector_mean_update(scale, &score0->mean_obs_dx, &score1->mean_obs_dx,
			   &work);

	score1->gamma += scale * (score0->gamma - score1->gamma);

	vector_mean_update(scale, &score0->dp, &score1->dp, &work);
	vector_mean_update(scale, &score0->mean_dx, &score1->mean_dx, &work);

	vector_deinit(&work);
}

static void imat_update(ssize_t n0, const struct recv_loglik_sender_imat *imat0,
			ssize_t n1, struct recv_loglik_sender_imat *imat1)
{
	assert(n0 >= 0);
	assert(imat0);
	assert(n1 >= 0);
	assert(imat1);
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
	vector_mean_update(scale, &imat0->gamma_mean_dx, &imat1->gamma_mean_dx,
			   &work);

	matrix_mean_update(scale, &imat0->dx_p, &imat1->dx_p, &work);
	matrix_mean_update(scale, &imat0->mean_dx_dp, &imat1->mean_dx_dp,
			   &work);
	matrix_mean_update(scale, &imat0->dp2, &imat1->dp2, &work);
	matrix_mean_update(scale, &imat0->var_dx, &imat1->var_dx, &work);

	vector_deinit(&work);
}

static void score_axpy_obs(double alpha,
			   const struct recv_loglik_sender_score *score,
			   ssize_t nsend,
			   const struct design *design,
			   ssize_t isend, struct vector *y)
{
	(void)isend; // unused
	
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);
	struct vector ysub = vector_slice(y, off, dim);

	// (X[0,i])^T n[i]
	design_recv_muls0(alpha / nsend, BLAS_TRANS,
			  design, &score->nrecv, 1.0, y);

	// sum{dx[t,i,j]}
	vector_axpy(alpha, &score->mean_obs_dx, &ysub);
}

static void score_axpy_mean(double alpha,
			    const struct recv_loglik_sender_score *score,
			    const struct array *active,
			    const struct design *design,
			    ssize_t isend, struct vector *y)
{
	(void)isend; // unused

	const struct vector *mean0 = score->mean0;
	double gamma = score->gamma;
	vector_axpy(alpha * gamma, mean0, y);

	ssize_t i, n = array_count(active);
	struct svector dp;
	svector_init(&dp, design_recv_count(design));
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);
		double val = vector_item(&score->dp, i);
		svector_set_item(&dp, jrecv, val);
	}

	design_recv_muls0(alpha, BLAS_TRANS, design, &dp, 1.0, y);
	svector_deinit(&dp);

	const struct vector *mean_dx = &score->mean_dx;
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);
	struct vector ysub = vector_slice(y, off, dim);
	vector_axpy(alpha, mean_dx, &ysub);
}

static void score_axpy(double alpha,
		       const struct recv_loglik_sender_score *score,
		       ssize_t nsend,
		       const struct array *active,
		       const struct design *design,
		       ssize_t isend, struct vector *y)
{
	score_axpy_obs(alpha, score, nsend, design, isend, y);
	score_axpy_mean(-alpha, score, active, design, isend, y);
}

static void imat_axpy(double alpha,
		      const struct recv_loglik_sender_imat *imat,
		      const struct recv_loglik_sender_score *score,
		      const struct array *active,
		      const struct design *design,
		      ssize_t isend, struct matrix *y)
{
	(void)isend; // unused

	const ssize_t dim = design_recv_dim(design);
	const ssize_t nrecv = design_recv_count(design);
	const struct vector *mean0 = score->mean0;
	const struct matrix *var0 = imat->imat0;
	const double gamma = score->gamma;
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
	design_recv_muls0(1.0, BLAS_TRANS, design, &gamma_dp, 0.0,
			  &gamma_x0_dp);

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
		design_recv_muls0(1.0, BLAS_TRANS, design, &dp2_j, 0.0, &dst);
	}

	/* Part I: X0' * [ Diag(p) - p p' ] * X0
	 */

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
		design_recv_muls0(alpha, BLAS_TRANS, design, &x0_dp2_k,
				  1.0, &dst);
	}

	/* for (i = 0; i < dim; i++) { assert(matrix_item(y, i, i) * alpha > -1e-5); } */

	ssize_t dyn_off = design_recv_dyn_index(design);
	ssize_t dyn_dim = design_recv_dyn_dim(design);
	struct matrix y_1 = matrix_slice(y, 0, dyn_off, dim, dyn_dim);
	struct matrix y1_ = matrix_slice(y, dyn_off, 0, dyn_dim, dim);
	struct matrix y11 = matrix_slice(y, dyn_off, dyn_off, dyn_dim, dyn_dim);

	/*   Part II: dX' * [ Diag(p) - p p' ] * dX
	 */

	matrix_axpy(alpha, &imat->var_dx, &y11);
	/* for (i = 0; i < dim; i++) { assert(matrix_item(y, i, i) * alpha > -1e-5); } */

	/*   Part III:   X0' * [ Diag(p) - p p' ] * dX
	 *             + dX' * [ Diag(p0 - p p' ] * X0
	 */

	struct vector x0_j;
	struct svector e_j;

	vector_init(&x0_j, dim);
	svector_init(&e_j, design_recv_count(design));

	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);

		svector_set_basis(&e_j, jrecv);
		design_recv_muls0(1.0, BLAS_TRANS, design, &e_j, 0.0, &x0_j);

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

	/* for (i = 0; i < dim; i++) { assert(matrix_item(y, i, i) * alpha > -1e-5); } */

}

static void info_init(struct recv_loglik_info *info, const struct recv_model *m)
{
	assert(info);
	assert(m);

	ssize_t dim = recv_model_dim(m);

	matrix_init(&info->imat, dim, dim);
	vector_init(&info->score, dim);
	vector_init(&info->mean, dim);
	info->dev = 0;
	info->nsend = 0;
	info->nrecv = 0;
}

static void info_deinit(struct recv_loglik_info *info)
{
	vector_deinit(&info->mean);
	vector_deinit(&info->score);
	matrix_deinit(&info->imat);
}

static void info_clear(struct recv_loglik_info *info)
{
	matrix_fill(&info->imat, 0);
	vector_fill(&info->score, 0);
	vector_fill(&info->mean, 0);
	info->dev = 0;
}

void sender_init(struct recv_loglik_sender *ll, const struct recv_model *model,
		 ssize_t isend)
{
	assert(ll);
	assert(model);
	assert(0 <= isend
	       && isend < design_send_count(recv_model_design(model)));

	ll->model = (struct recv_model *)model;
	ll->isend = isend;

	ll->n = 0;
	ll->n_last = 0;
	ll->dev_avg = 0.0;
	ll->dev_last = 0.0;
	array_init(&ll->active, sizeof(ssize_t));
	score_init(&ll->score_last, model, isend);
	score_init(&ll->score_avg, model, isend);
	imat_init(&ll->imat_last, model, isend);
	imat_init(&ll->imat_avg, model, isend);
}

void sender_deinit(struct recv_loglik_sender *ll)
{
	assert(ll);

	imat_deinit(&ll->imat_avg);
	imat_deinit(&ll->imat_last);
	score_deinit(&ll->score_avg);
	score_deinit(&ll->score_last);
	array_deinit(&ll->active);
}

void sender_clear(struct recv_loglik_sender *ll)
{
	assert(ll);

	ll->n_last = 0;
	ll->n = 0;
	ll->dev_last = 0.0;
	ll->dev_avg = 0.0;
	array_clear(&ll->active);
	score_clear(&ll->score_last);
	score_clear(&ll->score_avg);
	imat_clear(&ll->imat_last);
	imat_clear(&ll->imat_avg);
}

static void sender_update_active(struct recv_loglik_sender *ll,
				 const ssize_t *active, ssize_t n)
{
	assert(ll);
	assert(active);

	if (n > array_count(&ll->active)) {
		const ssize_t n0 = array_count(&ll->active);
		const ssize_t *begin0 = array_to_ptr(&ll->active);
		const ssize_t *end0 = begin0 + n0;
		const ssize_t *begin1 = active;
		const ssize_t *end1 = begin1 + n;
		const ssize_t *i0, *i1;

		for (i0 = begin0, i1 = begin1; i1 < end1; i1++) {
			if (i0 < end0 && *i0 == *i1) {
				i0++;
			} else {
				assert(i0 == end0 || *i1 < *i0);
				score_insert_active(&ll->score_avg,
						    i1 - begin1);
				imat_insert_active(&ll->imat_avg, i1 - begin1);
			}
		}
		assert(i0 == end0);
		assert(score_active_count(&ll->score_avg) == n);
		assert(imat_active_count(&ll->imat_avg) == n);

		array_clear(&ll->active);
		array_add_range(&ll->active, active, n);
	}
}

void sender_add(struct recv_loglik_sender *ll,
		const struct frame *f, const ssize_t *jrecv, ssize_t n)
{
	ssize_t isend = ll->isend;
	const struct recv_model *model = ll->model;
#ifndef NDEBUG
	ssize_t nreceiver = recv_model_count(model);
#endif

	double ntot = ll->n + n;
	double scale = n / ntot;
	ssize_t i;

	ssize_t *active, nactive;
	recv_model_get_active(model, isend, &active, &nactive);
	sender_update_active(ll, active, nactive);

	ll->dev_last = 0.0;
	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);

		double lp = recv_model_logprob(model, isend, jrecv[i]);
		ll->dev_last += -2 * lp;
	}
	ll->dev_avg += scale * (ll->dev_last / n - ll->dev_avg);

	score_set(&ll->score_last, f, isend, jrecv, n, model);
	score_update(n, &ll->score_last, ll->n, &ll->score_avg);

	imat_set(&ll->imat_last, f, model, isend, &ll->score_last);
	imat_update(n, &ll->imat_last, ll->n, &ll->imat_avg);

	ll->n_last = n;
	ll->n += n;
}

ssize_t sender_count(const struct recv_loglik_sender *sll)
{
	assert(sll);
	return sll->n;
}

ssize_t sender_last_count(const struct recv_loglik_sender *sll)
{
	assert(sll);
	return sll->n_last;
}

double sender_avg_dev(const struct recv_loglik_sender *sll)
{
	assert(sll);
	return sll->dev_avg;
}

double sender_last_dev(const struct recv_loglik_sender *sll)
{
	assert(sll);
	return sll->dev_last;
}

void sender_axpy_avg_mean(double alpha, const struct recv_loglik_sender *sll,
			  struct vector *y)
{
	assert(sll);
	assert(y);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	ssize_t isend = sll->isend;

	score_axpy_mean(alpha, &sll->score_avg, &sll->active, design, isend, y);
}

void sender_axpy_last_mean(double alpha, const struct recv_loglik_sender *sll,
			   struct vector *y)
{
	assert(sll);
	assert(y);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	ssize_t isend = sll->isend;

	score_axpy_mean(sll->n_last * alpha, &sll->score_last, &sll->active,
			design, isend, y);
}

void sender_axpy_avg_score(double alpha, const struct recv_loglik_sender *sll,
			   struct vector *y)
{
	assert(sll);
	assert(y);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	ssize_t isend = sll->isend;

	score_axpy(alpha, &sll->score_avg, sll->n, &sll->active, design, isend,
		   y);
}

void sender_axpy_last_score(double alpha, const struct recv_loglik_sender *sll,
			    struct vector *y)
{
	assert(sll);
	assert(y);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	ssize_t isend = sll->isend;

	score_axpy(sll->n_last * alpha, &sll->score_last, sll->n_last,
		   &sll->active, design, isend, y);
}

void sender_axpy_avg_imat(double alpha, const struct recv_loglik_sender *sll,
			  struct matrix *y)
{
	assert(sll);
	assert(y);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	ssize_t isend = sll->isend;

	imat_axpy(alpha, &sll->imat_avg, &sll->score_avg, &sll->active, design,
		  isend, y);
}

void sender_axpy_last_imat(double alpha, const struct recv_loglik_sender *sll,
			   struct matrix *y)
{
	assert(sll);
	assert(y);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	ssize_t isend = sll->isend;

	imat_axpy(sll->n_last * alpha, &sll->imat_last, &sll->score_last,
		  &sll->active, design, isend, y);
}

/* recv_loglik */

void cohort_init(struct recv_loglik_cohort *cll, const struct recv_model *m,
		 ssize_t c)
{
	(void)c; // unused
	assert(cll);
	assert(m);
	assert(0 <= c && c < recv_model_cohort_count(m));

	info_init(&cll->info, m);
	cll->info_cached = true;
}

void recv_loglik_init(struct recv_loglik *ll, struct recv_model *m)
{
	assert(ll);
	assert(m);

	ll->model = m;

	ssize_t ic, nc = recv_model_cohort_count(m);
	struct recv_loglik_cohort *cohorts = xcalloc(nc, sizeof(*cohorts));
	for (ic = 0; ic < nc; ic++) {
		cohort_init(&cohorts[ic], m, ic);
	}
	ll->cohorts = cohorts;

	ssize_t isend, nsend = recv_model_send_count(m);
	struct recv_loglik_sender *senders = xcalloc(nsend, sizeof(*senders));
	for (isend = 0; isend < nsend; isend++) {
		sender_init(&senders[isend], m, isend);
	}
	ll->senders = senders;

	ll->last = NULL;
}

void cohort_deinit(struct recv_loglik_cohort *cll)
{
	assert(cll);
	info_deinit(&cll->info);
}

void recv_loglik_deinit(struct recv_loglik *ll)
{
	assert(ll);

	struct recv_loglik_sender *senders = ll->senders;
	ssize_t isend, nsend = recv_model_send_count(ll->model);
	for (isend = 0; isend < nsend; isend++) {
		sender_deinit(&senders[isend]);
	}
	free(senders);

	struct recv_loglik_cohort *cohorts = ll->cohorts;
	ssize_t ic, nc = recv_model_cohort_count(ll->model);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cohorts[ic]);
	}
	free(cohorts);
}

void cohort_clear(struct recv_loglik_cohort *cll)
{
	assert(cll);
	info_clear(&cll->info);
	cll->info.dev = 0.0;
	cll->info.nsend = 0;
	cll->info.nrecv = 0;
	cll->info_cached = true;
}

void recv_loglik_clear(struct recv_loglik *ll)
{
	assert(ll);

	struct recv_loglik_cohort *cohorts = ll->cohorts;
	ssize_t ic, nc = recv_model_cohort_count(ll->model);
	for (ic = 0; ic < nc; ic++) {
		cohort_clear(&cohorts[ic]);
	}

	struct recv_loglik_sender *senders = ll->senders;
	ssize_t isend, nsend = recv_model_send_count(ll->model);
	for (isend = 0; isend < nsend; isend++) {
		sender_clear(&senders[isend]);
	}

	ll->last = NULL;
}

void cohort_add(struct recv_loglik_cohort *ll, ssize_t nto)
{
	ll->info.nsend += 1;
	ll->info.nrecv += nto;
	ll->info_cached = false;
}

void recv_loglik_add(struct recv_loglik *ll,
		     const struct frame *f, const struct message *msg)
{
	ssize_t isend = msg->from;
	ssize_t c = recv_model_cohort(ll->model, isend);
	sender_add(&ll->senders[isend], f, msg->to, msg->nto);
	cohort_add(&ll->cohorts[c], msg->nto);

	ll->last = &ll->senders[isend];
}

void recv_loglik_add_all(struct recv_loglik *ll,
			 struct frame *f, const struct messages *msgs)
{
	struct messages_iter it;
	const struct message *msg;
	ssize_t i, n;

	MESSAGES_FOREACH(it, msgs) {
		double t = MESSAGES_TIME(it);
		frame_advance(f, t);

		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			recv_loglik_add(ll, f, msg);
		}

		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(f, msg);
		}
	}

}

ssize_t recv_loglik_count(const struct recv_loglik *ll, ssize_t c)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));
	return ll->cohorts[c].info.nrecv;
}

ssize_t recv_loglik_count_sum(const struct recv_loglik *ll)
{
	assert(ll);

	ssize_t n = 0;
	ssize_t ic, nc = recv_model_cohort_count(ll->model);
	for (ic = 0; ic < nc; ic++) {
		n += recv_loglik_count(ll, ic);
	}

	return n;
}

ssize_t recv_loglik_last_count(const struct recv_loglik *ll)
{
	assert(ll);
	assert(ll->last);
	return sender_last_count(ll->last);
}

static double recv_loglik_avg_dev_nocache(const struct recv_loglik *ll,
					  ssize_t c)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));

	const ptrdiff_t *cohorts = recv_model_cohorts(ll->model);
	size_t isend, nsend = recv_model_count(ll->model);
	ssize_t ntot = 0;
	double dev_avg = 0.0;

	for (isend = 0; isend < nsend; isend++) {
		if (cohorts[isend] != c)
			continue;

		struct recv_loglik_sender *sll = &ll->senders[isend];
		ssize_t n = sender_count(sll);
		if (n > 0) {
			double dev = sender_avg_dev(sll);
			dev_avg += n * (dev - dev_avg) / (n + ntot);
			ntot += n;
		}
	}
	assert(ntot == recv_loglik_count(ll, c));

	return dev_avg;
}

static void recv_loglik_axpy_avg_mean_nocache(double alpha,
					      const struct recv_loglik *ll,
					      ssize_t c, struct vector *y)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));
	assert(y);
	assert(vector_dim(y) == recv_model_dim(ll->model));

	const ptrdiff_t *cohorts = recv_model_cohorts(ll->model);
	size_t isend, nsend = recv_model_count(ll->model);

	struct vector avg_mean, diff;
	vector_init(&avg_mean, recv_model_dim(ll->model));
	vector_init(&diff, recv_model_dim(ll->model));
	ssize_t ntot, n;

	ntot = 0;

	for (isend = 0; isend < nsend; isend++) {
		if (cohorts[isend] != c)
			continue;

		struct recv_loglik_sender *sll = &ll->senders[isend];

		n = sender_count(sll);
		if (n > 0) {
			ntot += n;
			vector_assign_copy(&diff, &avg_mean);
			sender_axpy_avg_mean(-1.0, sll, &diff);
			vector_axpy(-((double)n) / ntot, &diff, &avg_mean);
		}
	}
	assert(ntot == recv_loglik_count(ll, c));

	vector_axpy(alpha, &avg_mean, y);
	vector_deinit(&diff);
	vector_deinit(&avg_mean);
}

static void recv_loglik_axpy_avg_score_nocache(double alpha,
					       const struct recv_loglik *ll,
					       ssize_t c, struct vector *y)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));
	assert(y);
	assert(vector_dim(y) == recv_model_dim(ll->model));

	const ptrdiff_t *cohorts = recv_model_cohorts(ll->model);
	size_t isend, nsend = recv_model_count(ll->model);

	struct vector avg_score, diff;
	vector_init(&avg_score, recv_model_dim(ll->model));
	vector_init(&diff, recv_model_dim(ll->model));
	ssize_t ntot, n;

	ntot = 0;

	for (isend = 0; isend < nsend; isend++) {
		if (cohorts[isend] != c)
			continue;

		struct recv_loglik_sender *sll = &ll->senders[isend];

		n = sender_count(sll);
		if (n > 0) {
			ntot += n;
			vector_assign_copy(&diff, &avg_score);
			sender_axpy_avg_score(-1.0, sll, &diff);
			vector_axpy(-((double)n) / ntot, &diff, &avg_score);
		}
	}
	assert(ntot == recv_loglik_count(ll, c));

	vector_axpy(alpha, &avg_score, y);
	vector_deinit(&diff);
	vector_deinit(&avg_score);
}

static void recv_loglik_axpy_avg_imat_nocache(double alpha,
					      const struct recv_loglik *ll,
					      ssize_t c, struct matrix *y)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));
	assert(y);
	assert(matrix_nrow(y) == recv_model_dim(ll->model));
	assert(matrix_ncol(y) == recv_model_dim(ll->model));

	const ptrdiff_t *cohorts = recv_model_cohorts(ll->model);
	size_t isend, nsend = recv_model_count(ll->model);

	struct matrix avg_imat, new_imat, diff;
	ssize_t dim = recv_model_dim(ll->model);
	matrix_init(&avg_imat, dim, dim);
	matrix_init(&new_imat, dim, dim);
	matrix_init(&diff, dim, dim);
	ssize_t ntot, n;

	ntot = 0;

	for (isend = 0; isend < nsend; isend++) {
		if (cohorts[isend] != c)
			continue;

		struct recv_loglik_sender *sll = &ll->senders[isend];

		n = sender_count(sll);
		if (n > 0) {
			ntot += n;

			matrix_fill(&new_imat, 0.0);
			sender_axpy_avg_imat(1.0, sll, &new_imat);
			matrix_assign_copy(&diff, BLAS_NOTRANS, &avg_imat);
			matrix_axpy(-1.0, &new_imat, &diff);
			//sender_axpy_avg_imat(-1.0, sll, &diff);
			matrix_axpy(-((double)n) / ntot, &diff, &avg_imat);
		}
	}
	assert(ntot == recv_loglik_count(ll, c));

	matrix_axpy(alpha, &avg_imat, y);
	matrix_deinit(&diff);
	matrix_deinit(&new_imat);
	matrix_deinit(&avg_imat);
}

static void cache_info(struct recv_loglik *ll, ssize_t c)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));

	info_clear(&ll->cohorts[c].info);
	ll->cohorts[c].info.dev += recv_loglik_avg_dev_nocache(ll, c);
	recv_loglik_axpy_avg_mean_nocache(1, ll, c, &ll->cohorts[c].info.mean);
	recv_loglik_axpy_avg_score_nocache(1, ll, c,
					   &ll->cohorts[c].info.score);
	recv_loglik_axpy_avg_imat_nocache(1, ll, c, &ll->cohorts[c].info.imat);
}

struct recv_loglik_info *recv_loglik_info(const struct recv_loglik *ll,
					  ssize_t c)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));

	struct recv_loglik *mll = (struct recv_loglik *)ll;	// no const

	if (!mll->cohorts[c].info_cached) {
		cache_info(mll, c);
		mll->cohorts[c].info_cached = true;
	}

	return &mll->cohorts[c].info;
}

double recv_loglik_avg_dev(const struct recv_loglik *ll, ssize_t c)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));

	const struct recv_loglik_info *info = recv_loglik_info(ll, c);
	return info->dev;
}

void recv_loglik_axpy_avg_mean(double alpha, const struct recv_loglik *ll,
			       ssize_t c, struct vector *y)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));
	assert(y);
	assert(vector_dim(y) == recv_model_dim(ll->model));

	const struct recv_loglik_info *info = recv_loglik_info(ll, c);
	vector_axpy(alpha, &info->mean, y);
}

void recv_loglik_axpy_avg_score(double alpha, const struct recv_loglik *ll,
				ssize_t c, struct vector *y)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));
	assert(y);
	assert(vector_dim(y) == recv_model_dim(ll->model));

	const struct recv_loglik_info *info = recv_loglik_info(ll, c);
	vector_axpy(alpha, &info->score, y);
}

void recv_loglik_axpy_avg_imat(double alpha, const struct recv_loglik *ll,
			       ssize_t c, struct matrix *y)
{
	assert(ll);
	assert(0 <= c && c < recv_model_cohort_count(ll->model));
	assert(y);
	assert(matrix_nrow(y) == recv_model_dim(ll->model));
	assert(matrix_ncol(y) == recv_model_dim(ll->model));

	const struct recv_loglik_info *info = recv_loglik_info(ll, c);
	matrix_axpy(alpha, &info->imat, y);
}

double recv_loglik_last_dev(const struct recv_loglik *ll)
{
	assert(ll);
	assert(ll->last);

	return sender_last_dev(ll->last);
}

void recv_loglik_axpy_last_mean(double alpha, const struct recv_loglik *ll,
				struct vector *y)
{
	assert(ll);
	assert(ll->last);
	sender_axpy_last_mean(alpha, ll->last, y);
}

void recv_loglik_axpy_last_score(double alpha, const struct recv_loglik *ll,
				 struct vector *y)
{
	assert(ll);
	assert(ll->last);
	sender_axpy_last_score(alpha, ll->last, y);
}

void recv_loglik_axpy_last_imat(double alpha, const struct recv_loglik *ll,
				struct matrix *y)
{
	assert(ll);
	assert(ll->last);
	sender_axpy_last_imat(alpha, ll->last, y);
}
