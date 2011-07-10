#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "util.h"
#include "recv_loglik.h"

static void recv_sloglik_init(struct recv_sloglik *ll,
			      const struct model *model, ssize_t isend);
static void recv_sloglik_deinit(struct recv_sloglik *ll);

static void recv_sloglik_add(struct recv_sloglik *ll,
			     const struct frame *f, const ssize_t *jrecv,
			     ssize_t n);
static void recv_sloglik_clear(struct recv_sloglik *ll);

static ssize_t recv_sloglik_count(const struct recv_sloglik *sll);
static double recv_sloglik_avg_dev(const struct recv_sloglik *sll);
static void recv_sloglik_axpy_avg_mean(double alpha,
				       const struct recv_sloglik *sll,
				       struct vector *y);
static void recv_sloglik_axpy_avg_score(double alpha,
					const struct recv_sloglik *sll,
					struct vector *y);
static void recv_sloglik_axpy_avg_imat(double alpha,
				       const struct recv_sloglik *sll,
				       struct matrix *y);

static ssize_t recv_sloglik_last_count(const struct recv_sloglik *sll);
static double recv_sloglik_last_dev(const struct recv_sloglik *sll);
static void recv_sloglik_axpy_last_mean(double alpha,
					const struct recv_sloglik *sll,
					struct vector *y);
static void recv_sloglik_axpy_last_score(double alpha,
					 const struct recv_sloglik *sll,
					 struct vector *y);
static void recv_sloglik_axpy_last_imat(double alpha,
					const struct recv_sloglik *sll,
					struct matrix *y);

static void score_init(struct recv_sloglik_score *score,
		       const struct model *model, ssize_t isend)
{
	assert(score);

	const struct recv_model *rm = model_recv_model(model, isend);
	const struct vector *mean0 = recv_model_mean0(rm);
	const struct design *design = model_design(model);
	ssize_t nrecv = design_recv_count(design);
	ssize_t dyn_dim = design_recv_dyn_dim(design);

	score->mean0 = mean0;
	svector_init(&score->nrecv, nrecv);
	vector_init(&score->mean_obs_dx, dyn_dim);
	score->gamma = 1.0;
	vector_init(&score->dp, 0);
	vector_init(&score->mean_dx, dyn_dim);
}

static void imat_init(struct recv_sloglik_imat *imat,
		      const struct model *model, ssize_t isend)
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

static void score_deinit(struct recv_sloglik_score *score)
{
	assert(score);

	vector_deinit(&score->mean_dx);
	vector_deinit(&score->dp);
	vector_deinit(&score->mean_obs_dx);
	svector_deinit(&score->nrecv);
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

static void score_clear(struct recv_sloglik_score *score)
{
	assert(score);
	svector_clear(&score->nrecv);
	vector_fill(&score->mean_obs_dx, 0.0);
	score->gamma = 1.0;
	vector_reinit(&score->dp, 0);
	vector_fill(&score->mean_dx, 0.0);
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

static ssize_t score_active_count(const struct recv_sloglik_score *score)
{
	assert(score);
	return vector_dim(&score->dp);
}

static ssize_t imat_active_count(const struct recv_sloglik_imat *imat)
{
	assert(imat);
	return vector_dim(&imat->gamma_dp);
}

static void score_insert_active(struct recv_sloglik_score *score, ssize_t i)
{
	assert(score);
	assert(0 <= i && i <= score_active_count(score));

	ssize_t n0 = score_active_count(score);
	ssize_t n1 = n0 + 1;
	struct vector dp, src, dst;
	vector_init(&dp, n1);

	src = vector_slice(&score->dp, 0, i);
	dst = vector_slice(&dp, 0, i);
	vector_assign_copy(&dst, &src);

	src = vector_slice(&score->dp, i, n0 - i);
	dst = vector_slice(&dp, i + 1, n0 - i);
	vector_assign_copy(&dst, &src);

	vector_deinit(&score->dp);
	score->dp = dp;

	assert(score_active_count(score) == n1);
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
	matrix_assign_copy(&dst, TRANS_NOTRANS, &src);

	src = matrix_slice_cols(&imat->dx_p, i, n0 - i);
	dst = matrix_slice_cols(&dx_p, i + 1, n0 - i);
	matrix_assign_copy(&dst, TRANS_NOTRANS, &src);

	matrix_deinit(&imat->dx_p);
	imat->dx_p = dx_p;

	/* mean_dx_dp */
	struct matrix mean_dx_dp;
	matrix_init(&mean_dx_dp, matrix_nrow(&imat->mean_dx_dp), n1);

	src = matrix_slice_cols(&imat->mean_dx_dp, 0, i);
	dst = matrix_slice_cols(&mean_dx_dp, 0, i);
	matrix_assign_copy(&dst, TRANS_NOTRANS, &src);

	src = matrix_slice_cols(&imat->mean_dx_dp, i, n0 - i);
	dst = matrix_slice_cols(&mean_dx_dp, i + 1, n0 - i);
	matrix_assign_copy(&dst, TRANS_NOTRANS, &src);

	matrix_deinit(&imat->mean_dx_dp);
	imat->mean_dx_dp = mean_dx_dp;

	/* dp2 */
	struct matrix dp2;
	matrix_init(&dp2, n1, n1);

	dst = matrix_slice(&dp2, 0, 0, i, i);
	src = matrix_slice(&imat->dp2, 0, 0, i, i);
	matrix_assign_copy(&dst, TRANS_NOTRANS, &src);

	dst = matrix_slice(&dp2, i + 1, 0, n0 - i, i);
	src = matrix_slice(&imat->dp2, i, 0, n0 - i, i);
	matrix_assign_copy(&dst, TRANS_NOTRANS, &src);

	dst = matrix_slice(&dp2, 0, i + 1, i, n0 - i);
	src = matrix_slice(&imat->dp2, 0, i, i, n0 - i);
	matrix_assign_copy(&dst, TRANS_NOTRANS, &src);

	dst = matrix_slice(&dp2, i + 1, i + 1, n0 - i, n0 - i);
	src = matrix_slice(&imat->dp2, i, i, n0 - i, n0 - i);
	matrix_assign_copy(&dst, TRANS_NOTRANS, &src);
	matrix_deinit(&imat->dp2);
	imat->dp2 = dp2;
}

static void score_set_obs(struct recv_sloglik_score *score,
			  const struct frame *f,
			  ssize_t isend, const ssize_t *jrecv, const ssize_t n)
{
	ssize_t i;
	svector_clear(&score->nrecv);

	for (i = 0; i < n; i++) {
		*svector_item_ptr(&score->nrecv, jrecv[i]) += 1.0;
	}

	frame_recv_dmuls(1.0 / n, TRANS_TRANS, f, isend,
			 &score->nrecv, 0.0, &score->mean_obs_dx);
}

static void score_set_mean(struct recv_sloglik_score *score,
			   const struct frame *f,
			   const struct recv_model *model)
{
	assert(score);
	assert(model);
	assert(f);

	const ssize_t isend = model->isend;
	const struct array *active = &model->active;
	const ssize_t n = array_count(active);
	ssize_t i;

	score->mean0 = recv_model_mean0(model);

	const double gamma = model->gamma;
	score->gamma = gamma;

	if (vector_dim(&score->dp) != n) {
		vector_reinit(&score->dp, n);
	}
	vector_fill(&score->mean_dx, 0.0);

	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);

		/* mean */
		double p = recv_model_prob(model, jrecv);
		double p0 = recv_model_prob0(model, jrecv);
		double dp = p - gamma * p0;
		vector_set_item(&score->dp, i, dp);

		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		vector_axpy(p, dx, &score->mean_dx);
	}
}

static void score_set(struct recv_sloglik_score *score,
		      const struct frame *f,
		      ssize_t isend, const ssize_t *jrecv, ssize_t n,
		      const struct recv_model *model)
{
	score_set_obs(score, f, isend, jrecv, n);
	score_set_mean(score, f, model);
}

static void imat_set(struct recv_sloglik_imat *imat,
		     const struct frame *f,
		     const struct recv_model *model,
		     const struct recv_sloglik_score *score)
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
		ssize_t jrecv = *(ssize_t *)array_item(active, i);
		double dp = vector_item(&score->dp, i);
		double p = recv_model_prob(model, jrecv);
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

	matrix_assign_copy(&diff, TRANS_NOTRANS, val);
	matrix_sub(&diff, mean);
	matrix_axpy(scale, &diff, mean);
}

static void score_update(const struct array *active,
			 ssize_t n0,
			 const struct recv_sloglik_score *score0,
			 ssize_t n1, struct recv_sloglik_score *score1)
{
	assert(active);
	assert(n0 >= 0);
	assert(score0);
	assert(n1 >= 0);
	assert(score1);
	assert(score_active_count(score0) == array_count(active));
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

static void imat_update(const struct array *active,
			ssize_t n0,
			const struct recv_sloglik_imat *imat0,
			ssize_t n1, struct recv_sloglik_imat *imat1)
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
			   const struct recv_sloglik_score *score,
			   ssize_t nsend,
			   const struct design *design,
			   ssize_t isend, struct vector *y)
{
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);
	struct vector ysub = vector_slice(y, off, dim);

	// (X[0,i])^T n[i]
	design_recv_muls0(alpha / nsend, TRANS_TRANS,
			  design, &score->nrecv, 1.0, y);

	// sum{dx[t,i,j]}
	vector_axpy(alpha, &score->mean_obs_dx, &ysub);
}

static void score_axpy_mean(double alpha,
			    const struct recv_sloglik_score *score,
			    const struct array *active,
			    const struct design *design,
			    ssize_t isend, struct vector *y)
{
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

	design_recv_muls0(alpha, TRANS_TRANS, design, &dp, 1.0, y);
	svector_deinit(&dp);

	const struct vector *mean_dx = &score->mean_dx;
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);
	struct vector ysub = vector_slice(y, off, dim);
	vector_axpy(alpha, mean_dx, &ysub);
}

static void score_axpy(double alpha,
		       const struct recv_sloglik_score *score,
		       ssize_t nsend,
		       const struct array *active,
		       const struct design *design,
		       ssize_t isend, struct vector *y)
{
	score_axpy_obs(alpha, score, nsend, design, isend, y);
	score_axpy_mean(-alpha, score, active, design, isend, y);
}

static void imat_axpy(double alpha,
		      const struct recv_sloglik_imat *imat,
		      const struct recv_sloglik_score *score,
		      const struct array *active,
		      const struct design *design,
		      ssize_t isend, struct matrix *y)
{
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
	design_recv_muls0(1.0, TRANS_TRANS, design, &gamma_dp, 0.0,
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
		design_recv_muls0(1.0, TRANS_TRANS, design, &dp2_j, 0.0,
				  &dst);
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
		design_recv_muls0(alpha, TRANS_TRANS, design, &x0_dp2_k,
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
		design_recv_muls0(1.0, TRANS_TRANS, design, &e_j, 0.0, &x0_j);

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

static void info_init(struct recv_loglik_info *info, const struct model *m)
{
	assert(info);
	assert(m);

	ssize_t dim = model_recv_dim(m);

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

void recv_sloglik_init(struct recv_sloglik *ll, const struct model *model,
		       ssize_t isend)
{
	assert(ll);
	assert(model);
	assert(0 <= isend && isend < design_send_count(model_design(model)));

	ll->model = (struct model *)model;
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

void recv_sloglik_deinit(struct recv_sloglik *ll)
{
	assert(ll);

	imat_deinit(&ll->imat_avg);
	imat_deinit(&ll->imat_last);
	score_deinit(&ll->score_avg);
	score_deinit(&ll->score_last);
	array_deinit(&ll->active);
}

void recv_sloglik_clear(struct recv_sloglik *ll)
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
				score_insert_active(&ll->score_avg,
						    i1 - begin1);
				imat_insert_active(&ll->imat_avg, i1 - begin1);
			}
		}
		assert(i0 == end0);
		assert(score_active_count(&ll->score_avg) == n1);
		assert(imat_active_count(&ll->imat_avg) == n1);

		array_assign_copy(&ll->active, active);
	}
}

void recv_sloglik_add(struct recv_sloglik *ll,
		      const struct frame *f, const ssize_t *jrecv, ssize_t n)
{
	ssize_t isend = ll->isend;
	const struct recv_model *model = model_recv_model(ll->model, isend);
	ssize_t nreceiver = recv_model_count(model);

	double ntot = ll->n + n;
	double scale = n / ntot;
	ssize_t i;

	recv_sloglik_update_active(ll, &model->active);

	ll->dev_last = 0.0;
	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);

		double lp = recv_model_logprob(model, jrecv[i]);
		ll->dev_last += -2 * lp;
	}
	ll->dev_avg += scale * (ll->dev_last / n - ll->dev_avg);

	score_set(&ll->score_last, f, isend, jrecv, n, model);
	score_update(&model->active, n, &ll->score_last, ll->n, &ll->score_avg);

	imat_set(&ll->imat_last, f, model, &ll->score_last);
	imat_update(&model->active, n, &ll->imat_last, ll->n, &ll->imat_avg);

	ll->n_last = n;
	ll->n += n;
}

ssize_t recv_sloglik_count(const struct recv_sloglik *sll)
{
	assert(sll);
	return sll->n;
}

ssize_t recv_sloglik_last_count(const struct recv_sloglik *sll)
{
	assert(sll);
	return sll->n_last;
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

void recv_sloglik_axpy_avg_mean(double alpha, const struct recv_sloglik *sll,
				struct vector *y)
{
	assert(sll);
	assert(y);

	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;

	score_axpy_mean(alpha, &sll->score_avg, &sll->active, design, isend, y);
}

void recv_sloglik_axpy_last_mean(double alpha, const struct recv_sloglik *sll,
				 struct vector *y)
{
	assert(sll);
	assert(y);

	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;

	score_axpy_mean(sll->n_last * alpha, &sll->score_last, &sll->active,
			design, isend, y);
}

void recv_sloglik_axpy_avg_score(double alpha, const struct recv_sloglik *sll,
				 struct vector *y)
{
	assert(sll);
	assert(y);

	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;

	score_axpy(alpha, &sll->score_avg, sll->n, &sll->active, design, isend,
		   y);
}

void recv_sloglik_axpy_last_score(double alpha, const struct recv_sloglik *sll,
				  struct vector *y)
{
	assert(sll);
	assert(y);

	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;

	score_axpy(sll->n_last * alpha, &sll->score_last, sll->n_last,
		   &sll->active, design, isend, y);
}

void recv_sloglik_axpy_avg_imat(double alpha, const struct recv_sloglik *sll,
				struct matrix *y)
{
	assert(sll);
	assert(y);

	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;

	imat_axpy(alpha, &sll->imat_avg, &sll->score_avg, &sll->active, design,
		  isend, y);
}

void recv_sloglik_axpy_last_imat(double alpha, const struct recv_sloglik *sll,
				 struct matrix *y)
{
	assert(sll);
	assert(y);

	const struct model *model = sll->model;
	const struct design *design = model_design(model);
	ssize_t isend = sll->isend;

	imat_axpy(sll->n_last * alpha, &sll->imat_last, &sll->score_last,
		  &sll->active, design, isend, y);
}

/* recv_loglik */

void recv_loglik_init(struct recv_loglik *ll, struct model *m)
{
	assert(ll);
	assert(m);

	const struct design *design = model_design(m);
	ssize_t isend, nsend = design_send_count(design);

	ll->model = model_ref(m);

	array_init(&ll->slogliks, sizeof(struct recv_sloglik));
	array_set_capacity(&ll->slogliks, nsend);
	for (isend = 0; isend < nsend; isend++) {
		struct recv_sloglik *sll = array_add(&ll->slogliks, NULL);
		recv_sloglik_init(sll, m, isend);
	}

	ll->last = NULL;

	info_init(&ll->info, m);
	ll->info_cached = false;
}

void recv_loglik_deinit(struct recv_loglik *ll)
{
	assert(ll);

	info_deinit(&ll->info);
	struct recv_sloglik *sll;
	ARRAY_FOREACH(sll, &ll->slogliks) {
		recv_sloglik_deinit(sll);
	}
	array_deinit(&ll->slogliks);
	model_free(ll->model);
}

void recv_loglik_clear(struct recv_loglik *ll)
{
	assert(ll);

	struct recv_sloglik *sll;

	ARRAY_FOREACH(sll, &ll->slogliks) {
		recv_sloglik_clear(sll);
	}

	ll->last = NULL;
	info_clear(&ll->info);
	ll->info.nsend = 0;
	ll->info.nrecv = 0;
	ll->info_cached = false;
}

struct recv_loglik *recv_loglik_alloc(struct model *m, struct messages *msgs)
{
	struct recv_loglik *ll = xcalloc(1, sizeof(*ll));

	struct frame f;

	frame_init(&f, model_design(m));
	recv_loglik_init(ll, m);
	recv_loglik_add_all(ll, &f, msgs);
	frame_deinit(&f);

	return ll;
}

void recv_loglik_free(struct recv_loglik *ll)
{
	if (ll) {
		recv_loglik_deinit(ll);
		xfree(ll);
	}
}

void recv_loglik_add(struct recv_loglik *ll,
		     const struct frame *f, const struct message *msg)
{
	struct recv_sloglik *sll = array_item(&ll->slogliks, msg->from);
	recv_sloglik_add(sll, f, msg->to, msg->nto);
	ll->info.nsend += 1;
	ll->info.nrecv += msg->nto;
	ll->last = sll;
	ll->info_cached = false;

}

void recv_loglik_add_all(struct recv_loglik *ll,
			 struct frame *f, const struct messages *msgs)
{
	struct messages_iter it;
	const struct message *msg;
	ssize_t i, n;

	MESSAGES_FOREACH(it, msgs) {
		double t = MESSAGES_TIME(it);
		frame_advance_to(f, t);

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

ssize_t recv_loglik_count(const struct recv_loglik *ll)
{
	assert(ll);
	return ll->info.nrecv;
}

ssize_t recv_loglik_last_count(const struct recv_loglik *ll)
{
	assert(ll);
	assert(ll->last);
	return recv_sloglik_last_count(ll->last);
}

static double recv_loglik_avg_dev_nocache(const struct recv_loglik *ll)
{
	assert(ll);
	double dev, dev_avg;
	ssize_t ntot, n;
	struct recv_sloglik *sll;

	ntot = 0;
	dev_avg = 0.0;

	ARRAY_FOREACH(sll, &ll->slogliks) {
		n = recv_sloglik_count(sll);
		if (n > 0) {
			dev = recv_sloglik_avg_dev(sll);
			dev_avg += n * (dev - dev_avg) / (n + ntot);
			ntot += n;
		}
	}
	assert(ntot == recv_loglik_count(ll));

	return dev_avg;
}

static void recv_loglik_axpy_avg_mean_nocache(double alpha,
					      const struct recv_loglik *ll,
					      struct vector *y)
{
	struct vector avg_mean, diff;
	vector_init(&avg_mean, model_recv_dim(ll->model));
	vector_init(&diff, model_recv_dim(ll->model));
	ssize_t ntot, n;
	struct recv_sloglik *sll;

	ntot = 0;

	ARRAY_FOREACH(sll, &ll->slogliks) {
		n = recv_sloglik_count(sll);
		if (n > 0) {
			ntot += n;
			vector_assign_copy(&diff, &avg_mean);
			recv_sloglik_axpy_avg_mean(-1.0, sll, &diff);
			vector_axpy(-((double)n) / ntot, &diff, &avg_mean);
		}
	}
	assert(ntot == recv_loglik_count(ll));

	vector_axpy(alpha, &avg_mean, y);
	vector_deinit(&diff);
	vector_deinit(&avg_mean);
}

static void recv_loglik_axpy_avg_score_nocache(double alpha,
					       const struct recv_loglik *ll,
					       struct vector *y)
{
	struct vector avg_score, diff;
	vector_init(&avg_score, model_recv_dim(ll->model));
	vector_init(&diff, model_recv_dim(ll->model));
	ssize_t ntot, n;
	struct recv_sloglik *sll;

	ntot = 0;

	ARRAY_FOREACH(sll, &ll->slogliks) {
		n = recv_sloglik_count(sll);
		if (n > 0) {
			ntot += n;
			vector_assign_copy(&diff, &avg_score);
			recv_sloglik_axpy_avg_score(-1.0, sll, &diff);
			vector_axpy(-((double)n) / ntot, &diff, &avg_score);
		}
	}
	assert(ntot == recv_loglik_count(ll));

	vector_axpy(alpha, &avg_score, y);
	vector_deinit(&diff);
	vector_deinit(&avg_score);
}

static void recv_loglik_axpy_avg_imat_nocache(double alpha,
					      const struct recv_loglik *ll,
					      struct matrix *y)
{
	struct matrix avg_imat, new_imat, diff;
	ssize_t dim = model_recv_dim(ll->model);
	matrix_init(&avg_imat, dim, dim);
	matrix_init(&new_imat, dim, dim);
	matrix_init(&diff, dim, dim);
	ssize_t ntot, n;
	struct recv_sloglik *sll;

	ntot = 0;

	ARRAY_FOREACH(sll, &ll->slogliks) {
		n = recv_sloglik_count(sll);
		if (n > 0) {
			ntot += n;

			matrix_fill(&new_imat, 0.0);
			recv_sloglik_axpy_avg_imat(1.0, sll, &new_imat);
			matrix_assign_copy(&diff, TRANS_NOTRANS, &avg_imat);
			matrix_axpy(-1.0, &new_imat, &diff);
			//recv_sloglik_axpy_avg_imat(-1.0, sll, &diff);
			matrix_axpy(-((double)n) / ntot, &diff, &avg_imat);
		}
	}
	assert(ntot == recv_loglik_count(ll));

	matrix_axpy(alpha, &avg_imat, y);
	matrix_deinit(&diff);
	matrix_deinit(&new_imat);
	matrix_deinit(&avg_imat);
}

static void cache_info(struct recv_loglik *ll)
{
	info_clear(&ll->info);
	ll->info.dev += recv_loglik_avg_dev_nocache(ll);
	recv_loglik_axpy_avg_mean_nocache(1, ll, &ll->info.mean);
	recv_loglik_axpy_avg_score_nocache(1, ll, &ll->info.score);
	recv_loglik_axpy_avg_imat_nocache(1, ll, &ll->info.imat);
}

struct recv_loglik_info *recv_loglik_info(const struct recv_loglik *ll)
{
	struct recv_loglik *mll = (struct recv_loglik *)ll;	// no const

	if (!mll->info_cached) {
		cache_info(mll);
		mll->info_cached = true;
	}

	return &mll->info;
}

double recv_loglik_avg_dev(const struct recv_loglik *ll)
{
	const struct recv_loglik_info *info = recv_loglik_info(ll);
	return info->dev;
}

void recv_loglik_axpy_avg_mean(double alpha, const struct recv_loglik *ll,
			       struct vector *y)
{
	const struct recv_loglik_info *info = recv_loglik_info(ll);
	vector_axpy(alpha, &info->mean, y);
}

void recv_loglik_axpy_avg_score(double alpha, const struct recv_loglik *ll,
				struct vector *y)
{
	const struct recv_loglik_info *info = recv_loglik_info(ll);
	vector_axpy(alpha, &info->score, y);
}

void recv_loglik_axpy_avg_imat(double alpha, const struct recv_loglik *ll,
			       struct matrix *y)
{
	const struct recv_loglik_info *info = recv_loglik_info(ll);
	matrix_axpy(alpha, &info->imat, y);
}

double recv_loglik_last_dev(const struct recv_loglik *ll)
{
	assert(ll);
	assert(ll->last);

	return recv_sloglik_last_dev(ll->last);
}

void recv_loglik_axpy_last_mean(double alpha, const struct recv_loglik *ll,
				struct vector *y)
{
	assert(ll);
	assert(ll->last);
	recv_sloglik_axpy_last_mean(alpha, ll->last, y);
}

void recv_loglik_axpy_last_score(double alpha, const struct recv_loglik *ll,
				 struct vector *y)
{
	assert(ll);
	assert(ll->last);
	recv_sloglik_axpy_last_score(alpha, ll->last, y);
}

void recv_loglik_axpy_last_imat(double alpha, const struct recv_loglik *ll,
				struct matrix *y)
{
	assert(ll);
	assert(ll->last);
	recv_sloglik_axpy_last_imat(alpha, ll->last, y);
}
