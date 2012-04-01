#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sblas.h"
#include "xalloc.h"
#include "recv_loglik.h"

static void cohort_init(struct recv_loglik_cohort *cll,
			const struct recv_model *m, size_t c);
static void cohort_deinit(struct recv_loglik_cohort *cll);
static void cohort_clear(struct recv_loglik_cohort *cll, size_t dim);

static void sender_init(struct recv_loglik_sender *ll,
			const struct recv_model *model, size_t isend);
static void sender_deinit(struct recv_loglik_sender *ll);

static void sender_add(struct recv_loglik_sender *ll,
		       const struct frame *f, const size_t *jrecv, size_t n);
static void sender_clear(struct recv_loglik_sender *ll, size_t dyn_dim);

static size_t sender_count(const struct recv_loglik_sender *sll);
static double sender_dev(const struct recv_loglik_sender *sll);
static void sender_axpy_mean(double alpha,
			     const struct recv_loglik_sender *sll,
			     double *y);
static void sender_axpy_score(double alpha,
			      const struct recv_loglik_sender *sll,
			      double *y);
static void sender_axpy_imat(double alpha,
			     const struct recv_loglik_sender *sll,
			     struct matrix *y);

static size_t sender_last_count(const struct recv_loglik_sender *sll);
static double sender_last_dev(const struct recv_loglik_sender *sll);
static void sender_axpy_last_mean(double alpha,
				  const struct recv_loglik_sender *sll,
				  double *y);
static void sender_axpy_last_score(double alpha,
				   const struct recv_loglik_sender *sll,
				   double *y);
static void sender_axpy_last_imat(double alpha,
				  const struct recv_loglik_sender *sll,
				  struct matrix *y);

static void score_init(struct recv_loglik_sender_score *score,
		       const struct recv_model *model, size_t isend)
{
	assert(score);

	size_t ic = recv_model_cohort(model, isend);
	const double *mean0 = recv_model_mean0(model, ic);
	const struct design *design = recv_model_design(model);
	size_t dyn_dim = design_dvars_dim(design);

	score->mean0 = mean0;
	score->nrecv = NULL;
	vpattern_init(&score->nrecv_pat);
	score->obs_dx = xcalloc(dyn_dim, sizeof(score->obs_dx[0]));
	score->gamma = 0.0;
	score->dp = NULL;
	score->mean_dx = xcalloc(dyn_dim, sizeof(score->mean_dx[0]));
#ifndef NDEBUG
	score->nz = 0;
#endif
}

static void imat_init(struct recv_loglik_sender_imat *imat,
		      const struct recv_model *model, size_t isend)
{
	assert(imat);

	size_t ic = recv_model_cohort(model, isend);
	const struct matrix *imat0 = recv_model_imat0(model, ic);
	const struct design *design = recv_model_design(model);
	size_t dyn_dim = design_dvars_dim(design);

	imat->imat0 = imat0;
	imat->gamma2 = 0.0;
	imat->gamma_dp = NULL;
	imat->gamma_mean_dx = xcalloc(dyn_dim, sizeof(imat->gamma_mean_dx[0]));
	matrix_init(&imat->dx_p, dyn_dim, 0);
	matrix_init(&imat->mean_dx_dp, dyn_dim, 0);
	matrix_init(&imat->dp2, 0, 0);
	matrix_init(&imat->var_dx, dyn_dim, dyn_dim);
#ifndef NDEBUG
	imat->nz = 0;
#endif
}

static void score_deinit(struct recv_loglik_sender_score *score)
{
	assert(score);

	free(score->mean_dx);
	free(score->dp);
	free(score->obs_dx);
	vpattern_deinit(&score->nrecv_pat);
	free(score->nrecv);
}

static void imat_deinit(struct recv_loglik_sender_imat *imat)
{
	assert(imat);

	matrix_deinit(&imat->var_dx);
	matrix_deinit(&imat->dp2);
	matrix_deinit(&imat->mean_dx_dp);
	matrix_deinit(&imat->dx_p);
	free(imat->gamma_mean_dx);
	free(imat->gamma_dp);
}

static void score_clear(struct recv_loglik_sender_score *score, size_t dyn_dim)
{
	assert(score);
	vpattern_clear(&score->nrecv_pat);
	memset(score->obs_dx, 0, dyn_dim * sizeof(score->obs_dx[0]));
	score->gamma = 0.0;
	score->dp = xrealloc(score->dp, 0);
	memset(score->mean_dx, 0, dyn_dim * sizeof(score->mean_dx[0]));
#ifndef NDEBUG
	score->nz = 0;
#endif
}

static void imat_clear(struct recv_loglik_sender_imat *imat, size_t dyn_dim)
{
	assert(imat);
	imat->gamma2 = 0.0;
	imat->gamma_dp = xrealloc(imat->gamma_dp, 0);
	memset(imat->gamma_mean_dx, 0, dyn_dim * sizeof(imat->gamma_mean_dx[0]));
	matrix_reinit(&imat->dx_p, matrix_nrow(&imat->dx_p), 0);
	matrix_reinit(&imat->mean_dx_dp, matrix_nrow(&imat->mean_dx_dp), 0);
	matrix_reinit(&imat->dp2, 0, 0);
	matrix_fill(&imat->var_dx, 0.0);
#ifndef NDEBUG
	imat->nz = 0;
#endif
}

static void score_insert_active(struct recv_loglik_sender_score *score,
				size_t i, size_t nz0)
{
	assert(score->nz == nz0);

	score->dp = xrealloc(score->dp, (nz0 + 1) * sizeof(score->dp[0]));
	memmove(score->dp + i + 1, score->dp + i, (nz0 - i) *
			sizeof(score->dp[0]));
	score->dp[i] = 0.0;
#ifndef NDEBUG
	score->nz = nz0 + 1;
#endif
}

static void imat_insert_active(struct recv_loglik_sender_imat *imat, size_t i,
		size_t nz0)
{
	assert(imat->nz == nz0);

	imat->gamma_dp = xrealloc(imat->gamma_dp, (nz0 + 1) * sizeof(imat->gamma_dp[0]));
	memmove(imat->gamma_dp + i + 1, imat->gamma_dp + i,
		(nz0 - i) * sizeof(imat->gamma_dp[0]));
	imat->gamma_dp[i] = 0.0;
	matrix_insert_col(&imat->dx_p, i);
	matrix_insert_col(&imat->mean_dx_dp, i);
	matrix_insert_row_col(&imat->dp2, i, i);

#ifndef NDEBUG
	imat->nz = nz0 + 1;
#endif
}

static int size_compar(const void *x1, const void *x2)
{
	size_t y1 = *(size_t *)x1;
	size_t y2 = *(size_t *)x2;

	if (y1 < y2) {
		return -1;
	} else if (y1 > y2) {
		return +1;
	} else {
		return 0;
	}
}

static void score_set_obs(struct recv_loglik_sender_score *score,
			  const struct frame *f,
			  size_t isend, const size_t *jrecv, const size_t n)
{
	size_t i;
	size_t nznew;

	vpattern_clear(&score->nrecv_pat);
	if ((nznew = vpattern_grow(&score->nrecv_pat, n))) {
		score->nrecv = xrealloc(score->nrecv,
					nznew * sizeof(score->nrecv[0]));
	}

	for (i = 0; i < n; i++) {
		score->nrecv[i] = 1.0;
	}
	memcpy(score->nrecv_pat.indx, jrecv, n * sizeof(jrecv[0]));
	qsort(score->nrecv_pat.indx, n, sizeof(jrecv[0]), size_compar);
	score->nrecv_pat.nz = n;

	frame_recv_dmuls(1.0, BLAS_TRANS, f, isend,
			 score->nrecv, &score->nrecv_pat, 0.0, score->obs_dx);
}

static void score_set_mean(struct recv_loglik_sender_score *score,
			   const struct frame *f,
			   const struct recv_model *model, size_t isend, const size_t n)
{
	assert(score);
	assert(model);
	assert(f);

	const struct design *d = recv_model_design(model);
	size_t dyn_dim = design_dvars_dim(d);
	size_t *active, iz, nz;
	size_t ic = recv_model_cohort(model, isend);
	recv_model_get_active(model, isend, &active, &nz);

	score->mean0 = recv_model_mean0(model, ic);

	const double gamma = recv_model_invgrow(model, isend);
	score->gamma = n * gamma;

	assert(score->nz == nz);
	memset(score->mean_dx, 0, dyn_dim * sizeof(score->mean_dx[0]));

	for (iz = 0; iz < nz; iz++) {
		size_t jrecv = active[iz];

		/* mean */
		double p = recv_model_prob(model, isend, jrecv);
		double p0 = recv_model_prob0(model, ic, jrecv);
		double dp = p - gamma * p0;
		score->dp[iz] = dp;

		const double *dx = frame_recv_dx(f, isend, jrecv);
		blas_daxpy(dyn_dim, p, dx, 1, score->mean_dx, 1);
	}
	blas_dscal(nz, n, score->dp, 1);
	blas_dscal(dyn_dim, n, score->mean_dx, 1);
}

static void score_set(struct recv_loglik_sender_score *score,
		      const struct frame *f,
		      size_t isend, const size_t *jrecv, size_t n,
		      const struct recv_model *model)
{
	score_set_obs(score, f, isend, jrecv, n);
	score_set_mean(score, f, model, isend, n);
}

static void imat_set(struct recv_loglik_sender_imat *imat,
		     const struct frame *f,
		     const struct recv_model *model,
		     size_t isend, size_t n,
		     const struct recv_loglik_sender_score *score)
{
	assert(imat);
	assert(model);
	assert(f);

	const struct design *d = recv_model_design(model);
	size_t dyn_dim = design_dvars_dim(d);
	size_t ic = recv_model_cohort(model, isend);
	size_t *active, iz, nz;
	recv_model_get_active(model, isend, &active, &nz);

	assert(score->nz == nz);
	assert(imat->nz == nz);

	imat->imat0 = recv_model_imat0(model, ic);

	const double gamma = recv_model_invgrow(model, isend);

	matrix_fill(&imat->dp2, 0.0);
	struct matrix y;
	double ptot = 0.0;

	/* gamma2 */
	imat->gamma2 = n * gamma * (1 - gamma);

	/* gamma_mean_dx */
	blas_dcopy(dyn_dim, score->mean_dx, 1, imat->gamma_mean_dx, 1);
	blas_dscal(dyn_dim, gamma, imat->gamma_mean_dx, 1);

	/* dp2 */
	matrix_update1(&imat->dp2, -1.0 / n, score->dp, score->dp);

	matrix_init(&y, dyn_dim, nz);

	for (iz = 0; iz < nz; iz++) {
		size_t jrecv = active[iz];
		double dp = score->dp[iz];
		double p = recv_model_prob(model, isend, jrecv);
		const double *dx = frame_recv_dx(f, isend, jrecv);

		/* gamma_dp */
		imat->gamma_dp[iz] = gamma * dp;

		/* mean_dx_dp */
		double *mean_dx_dp_j = matrix_col(&imat->mean_dx_dp, iz);
		blas_dcopy(dyn_dim, score->mean_dx, 1, mean_dx_dp_j, 1);
		blas_dscal(dyn_dim, dp / n, mean_dx_dp_j, 1);

		/* dp2 */
		*matrix_item_ptr(&imat->dp2, iz, iz) += dp;

		/* dx_p */
		double *dx_p_j = matrix_col(&imat->dx_p, iz);
		blas_dcopy(dyn_dim, dx, 1, dx_p_j, 1);
		blas_dscal(dyn_dim, n * p, dx_p_j, 1);

		/* var_dx */
		double *y_i = matrix_col(&y, iz);
		blas_dcopy(dyn_dim, dx, 1, y_i, 1);
		blas_daxpy(dyn_dim, -1.0 / n, score->mean_dx, 1, y_i, 1);
		blas_dscal(dyn_dim, sqrt(p), y_i, 1);

		ptot += p;
	}
	blas_dgemm(BLAS_NOTRANS, BLAS_TRANS, dyn_dim, dyn_dim, nz, n,
		   matrix_to_ptr(&y),
		   matrix_lda(&y), matrix_to_ptr(&y), matrix_lda(&y), 0.0,
		   matrix_to_ptr(&imat->var_dx), matrix_lda(&imat->var_dx));
	/* var_dx */
	matrix_update1(&imat->var_dx, (1.0 - ptot) / n, score->mean_dx,
		       score->mean_dx);

	matrix_deinit(&y);
}

static void nrecv_add(const struct recv_loglik_sender_score *score0,
		      struct recv_loglik_sender_score *score1)
{
	size_t i, n = score0->nrecv_pat.nz;
	size_t nzmax = score1->nrecv_pat.nzmax;

	for (i = 0; i < n; i++) {
		size_t ind = score0->nrecv_pat.indx[i];
		int ins;
		size_t ix = vpattern_search(&score1->nrecv_pat, ind, &ins);
		if (ins) {
			if (nzmax != score1->nrecv_pat.nzmax) {
				nzmax = score1->nrecv_pat.nzmax;
				score1->nrecv = xrealloc(score1->nrecv,
							 nzmax * sizeof(score1->nrecv[0]));
			}
			memmove(score1->nrecv + ix + 1, score1->nrecv + ix,
				(score1->nrecv_pat.nz - 1 - ix) * sizeof(score1->nrecv[0]));
			score1->nrecv[ix] = 0.0;
		}
		score1->nrecv[ix] += score0->nrecv[i];
	}
}

static void score_update(const struct recv_loglik_sender_score *score0,
			 struct recv_loglik_sender_score *score1,
			 size_t dyn_dim, size_t nz)
{
	assert(score0->mean0 == score1->mean0);
	assert(score0->nz == nz);
	assert(score1->nz == nz);

	nrecv_add(score0, score1);
	blas_daxpy(dyn_dim, 1.0, score0->obs_dx, 1, score1->obs_dx, 1);

	score1->gamma += score0->gamma;
	blas_daxpy(nz, 1.0, score0->dp, 1, score1->dp, 1);
	blas_daxpy(dyn_dim, 1.0, score0->mean_dx, 1, score1->mean_dx, 1);
}

static void imat_update(const struct recv_loglik_sender_imat *imat0,
			struct recv_loglik_sender_imat *imat1,
			size_t dyn_dim, size_t nz)
{
	assert(imat0);
	assert(imat1);
	assert(imat0->imat0 == imat1->imat0);
	assert(imat0->nz == nz);
	assert(imat1->nz == nz);

	imat1->gamma2 += imat0->gamma2;

	blas_daxpy(nz, 1.0, imat0->gamma_dp, 1, imat1->gamma_dp, 1);
	blas_daxpy(dyn_dim, 1.0, imat0->gamma_mean_dx, 1, imat1->gamma_mean_dx, 1);

	matrix_axpy(1.0, &imat0->dx_p, &imat1->dx_p);
	matrix_axpy(1.0, &imat0->mean_dx_dp, &imat1->mean_dx_dp);
	matrix_axpy(1.0, &imat0->dp2, &imat1->dp2);
	matrix_axpy(1.0, &imat0->var_dx, &imat1->var_dx);
}

static void score_axpy_obs(double alpha,
			   const struct recv_loglik_sender_score *score,
			   const struct design *d,
			   size_t isend, double *y)
{
	(void)isend;		// unused

	// (X[0,i])^T n[i]
	design_muls0(alpha, BLAS_TRANS, d, score->nrecv, &score->nrecv_pat, 1.0, y);

	size_t off = design_dvars_index(d);
	size_t dim = design_dvars_dim(d);

	// sum{dx[t,i,j]}
	blas_daxpy(dim, alpha, score->obs_dx, 1, y + off, 1);
}

static void score_axpy_mean(double alpha,
			    const struct recv_loglik_sender_score *score,
			    const struct vpattern *active,
			    const struct design *d,
			    size_t isend, double *y)
{
	(void)isend;		// unused

	size_t dim = design_dim(d);
	const double *mean0 = score->mean0;
	double gamma = score->gamma;
	blas_daxpy(dim, alpha * gamma, mean0, 1, y, 1);

	design_muls0(alpha, BLAS_TRANS, d, score->dp, active, 1.0, y);

	const double *mean_dx = score->mean_dx;
	size_t dyn_off = design_dvars_index(d);
	size_t dyn_dim = design_dvars_dim(d);
	blas_daxpy(dyn_dim, alpha, mean_dx, 1, y + dyn_off, 1);

}

static void score_axpy(double alpha,
		       const struct recv_loglik_sender_score *score,
		       const struct vpattern *active,
		       const struct design *design,
		       size_t isend, double *y)
{
	score_axpy_obs(alpha, score, design, isend, y);
	score_axpy_mean(-alpha, score, active, design, isend, y);
}

static void imat_axpy(double alpha,
		      const struct recv_loglik_sender_imat *imat,
		      const struct recv_loglik_sender_score *score,
		      const struct vpattern *active,
		      const struct design *design,
		      size_t isend, struct matrix *y)
{
	(void)isend;		// unused

	const size_t dim = design_dim(design);
	size_t dyn_off = design_dvars_index(design);
	size_t dyn_dim = design_dvars_dim(design);
	const double *mean0 = score->mean0;
	const struct matrix *var0 = imat->imat0;
	const double gamma = score->gamma;
	const double gamma2 = imat->gamma2;

	const size_t n = active->nz;
	size_t i, j, k;

	double *gamma_x0_dp = xmalloc(dim * sizeof(gamma_x0_dp[0]));
	design_muls0(1.0, BLAS_TRANS, design, imat->gamma_dp,
		     active, 0.0, gamma_x0_dp);

	struct matrix x0_dp2;

	matrix_init(&x0_dp2, dim, n);
	for (j = 0; j < n; j++) {
		const double *dp2_j = matrix_item_ptr(&imat->dp2, 0, j);
		double *dst = matrix_col(&x0_dp2, j);
		design_muls0(1.0, BLAS_TRANS, design, dp2_j, active, 0.0, dst);
	}

	/* Part I: X0' * [ Diag(p) - p p' ] * X0
	 */

	matrix_axpy(alpha * gamma, var0, y);
	matrix_update1(y, alpha * gamma2, mean0, mean0);
	matrix_update1(y, -alpha, mean0, gamma_x0_dp);
	matrix_update1(y, -alpha, gamma_x0_dp, mean0);

	double *x0_dp2_k = xmalloc(n * sizeof(x0_dp2_k));
	for (k = 0; k < dim; k++) {
		blas_dcopy(n, matrix_item_ptr(&x0_dp2, k, 0),
			   matrix_lda(&x0_dp2), x0_dp2_k, 1);
		double *dst = matrix_col(y, k);
		design_muls0(alpha, BLAS_TRANS, design, x0_dp2_k, active, 1.0, dst);
	}

	/* for (i = 0; i < dim; i++) { assert(matrix_item(y, i, i) * alpha > -1e-5); } */

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

	double one = 1.0;
	struct vpattern pat_j;
	size_t jrecv = 0;

	double *x0_j = xmalloc(dim * sizeof(x0_j[0]));

	pat_j.nz = 1;
	pat_j.indx = &jrecv;

	for (i = 0; i < n; i++) {
		jrecv = active->indx[i];
		design_muls0(1.0, BLAS_TRANS, design, &one, &pat_j, 0.0, x0_j);

		const double *dx_p_j = matrix_col(&imat->dx_p, i);
		matrix_update1(&y_1, alpha, x0_j, dx_p_j);
		matrix_update1(&y1_, alpha, dx_p_j, x0_j);

		const double *mean_dx_dp_j = matrix_col(&imat->mean_dx_dp, i);
		matrix_update1(&y_1, -alpha, x0_j, mean_dx_dp_j);
		matrix_update1(&y1_, -alpha, mean_dx_dp_j, x0_j);
	}

	matrix_update1(&y_1, -alpha, mean0, imat->gamma_mean_dx);
	matrix_update1(&y1_, -alpha, imat->gamma_mean_dx, mean0);

	free(x0_j);
	free(x0_dp2_k);
	matrix_deinit(&x0_dp2);
	free(gamma_x0_dp);
}

static void info_init(struct recv_loglik_info *info, const struct recv_model *m)
{
	assert(info);
	assert(m);

	size_t dim = recv_model_dim(m);

	matrix_init(&info->imat, dim, dim);
	info->score = xcalloc(dim, sizeof(info->score[0]));
	info->mean = xcalloc(dim, sizeof(info->mean[0]));
	info->dev = 0;
	info->nsend = 0;
	info->nrecv = 0;
}

static void info_deinit(struct recv_loglik_info *info)
{
	free(info->mean);
	free(info->score);
	matrix_deinit(&info->imat);
}

static void info_clear(struct recv_loglik_info *info, size_t dim)
{
	assert(info->imat.nrow == dim);
	assert(info->imat.ncol == dim);

	matrix_fill(&info->imat, 0);
	memset(info->score, 0, dim * sizeof(info->score[0]));
	memset(info->mean, 0, dim * sizeof(info->mean[0]));
	info->dev = 0;
}

void sender_init(struct recv_loglik_sender *ll, const struct recv_model *model,
		 size_t isend)
{
	assert(ll);
	assert(model);

	ll->model = (struct recv_model *)model;
	ll->isend = isend;

	ll->n = 0;
	ll->n_last = 0;
	ll->dev = 0.0;
	ll->dev_last = 0.0;
	vpattern_init(&ll->active);
	score_init(&ll->score_last, model, isend);
	score_init(&ll->score, model, isend);
	imat_init(&ll->imat_last, model, isend);
	imat_init(&ll->imat, model, isend);
}

void sender_deinit(struct recv_loglik_sender *ll)
{
	assert(ll);

	imat_deinit(&ll->imat);
	imat_deinit(&ll->imat_last);
	score_deinit(&ll->score);
	score_deinit(&ll->score_last);
	vpattern_deinit(&ll->active);
}

void sender_clear(struct recv_loglik_sender *ll, size_t dyn_dim)
{
	assert(ll);

	ll->n_last = 0;
	ll->n = 0;
	ll->dev_last = 0.0;
	ll->dev = 0.0;
	vpattern_clear(&ll->active);
	score_clear(&ll->score_last, dyn_dim);
	score_clear(&ll->score, dyn_dim);
	imat_clear(&ll->imat_last, dyn_dim);
	imat_clear(&ll->imat, dyn_dim);
}

static void sender_update_active(struct recv_loglik_sender *ll,
				 const size_t *active, size_t n)
{
	assert(ll);
	assert(active);

	size_t n0 = ll->active.nz;

	if (n > n0) {
		vpattern_grow(&ll->active, n - n0);

		const size_t *begin0 = ll->active.indx;
		const size_t *end0 = begin0 + n0;
		const size_t *begin1 = active;
		const size_t *end1 = begin1 + n;
		const size_t *i0, *i1;

		for (i0 = begin0, i1 = begin1; i1 < end1; i1++) {
			if (i0 < end0 && *i0 == *i1) {
				i0++;
			} else {
				assert(i0 == end0 || *i1 < *i0);
				score_insert_active(&ll->score, i1 - begin1, n0);
				score_insert_active(&ll->score_last, i1 - begin1, n0);
				imat_insert_active(&ll->imat, i1 - begin1, n0);
				imat_insert_active(&ll->imat_last, i1 - begin1, n0);
				n0++;
			}
		}
		assert(i0 == end0);
		assert(n0 == n);

		memcpy(ll->active.indx, active, n * sizeof(active[0]));
		ll->active.nz = n;
	}
}

void sender_add(struct recv_loglik_sender *ll,
		const struct frame *f, const size_t *jrecv, size_t n)
{
	size_t isend = ll->isend;
	const struct recv_model *model = ll->model;
	const struct design *d = recv_model_design(model);
	size_t dyn_dim = design_dvars_dim(d);
	size_t i;

	size_t *active, nactive;
	recv_model_get_active(model, isend, &active, &nactive);
	sender_update_active(ll, active, nactive);

	ll->dev_last = 0.0;
	for (i = 0; i < n; i++) {
		double lp = recv_model_logprob(model, isend, jrecv[i]);
		ll->dev_last += -2 * lp;
	}
	ll->dev += ll->dev_last;

	score_set(&ll->score_last, f, isend, jrecv, n, model);
	score_update(&ll->score_last, &ll->score, dyn_dim, nactive);

	imat_set(&ll->imat_last, f, model, isend, n, &ll->score_last);
	imat_update(&ll->imat_last, &ll->imat, dyn_dim, nactive);

	ll->n_last = n;
	ll->n += n;
}

size_t sender_count(const struct recv_loglik_sender *sll)
{
	assert(sll);
	return sll->n;
}

size_t sender_last_count(const struct recv_loglik_sender *sll)
{
	assert(sll);
	return sll->n_last;
}

double sender_dev(const struct recv_loglik_sender *sll)
{
	assert(sll);
	return sll->dev;
}

double sender_last_dev(const struct recv_loglik_sender *sll)
{
	assert(sll);
	return sll->dev_last;
}

void sender_axpy_mean(double alpha, const struct recv_loglik_sender *sll,
		      double *y)
{
	assert(sll);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	size_t isend = sll->isend;

	score_axpy_mean(alpha, &sll->score, &sll->active, design, isend, y);
}

void sender_axpy_last_mean(double alpha, const struct recv_loglik_sender *sll,
			   double *y)
{
	assert(sll);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	size_t isend = sll->isend;

	score_axpy_mean(alpha, &sll->score_last, &sll->active, design, isend, y);
}

void sender_axpy_score(double alpha, const struct recv_loglik_sender *sll,
		       double *y)
{
	assert(sll);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	size_t isend = sll->isend;

	score_axpy(alpha, &sll->score, &sll->active, design, isend, y);
}

void sender_axpy_last_score(double alpha, const struct recv_loglik_sender *sll,
			    double *y)
{
	assert(sll);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	size_t isend = sll->isend;

	score_axpy(alpha, &sll->score_last, &sll->active, design, isend, y);
}

void sender_axpy_imat(double alpha, const struct recv_loglik_sender *sll,
			  struct matrix *y)
{
	assert(sll);
	assert(y);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	size_t isend = sll->isend;

	imat_axpy(alpha, &sll->imat, &sll->score, &sll->active, design, isend, y);
}

void sender_axpy_last_imat(double alpha, const struct recv_loglik_sender *sll,
			   struct matrix *y)
{
	assert(sll);
	assert(y);

	const struct recv_model *model = sll->model;
	const struct design *design = recv_model_design(model);
	size_t isend = sll->isend;

	imat_axpy(alpha, &sll->imat_last, &sll->score_last, &sll->active, design, isend, y);
}

/* recv_loglik */

void cohort_init(struct recv_loglik_cohort *cll, const struct recv_model *m,
		 size_t c)
{
	(void)c;		// unused
	assert(cll);
	assert(m);
	assert(c < recv_model_cohort_count(m));

	info_init(&cll->info, m);
	cll->info_cached = true;
}

void recv_loglik_init(struct recv_loglik *ll, struct recv_model *m)
{
	assert(ll);
	assert(m);

	ll->model = m;

	size_t ic, nc = recv_model_cohort_count(m);
	struct recv_loglik_cohort *cohorts = xcalloc(nc, sizeof(*cohorts));
	for (ic = 0; ic < nc; ic++) {
		cohort_init(&cohorts[ic], m, ic);
	}
	ll->cohorts = cohorts;

	size_t isend, nsend = recv_model_send_count(m);
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
	size_t isend, nsend = recv_model_send_count(ll->model);
	for (isend = 0; isend < nsend; isend++) {
		sender_deinit(&senders[isend]);
	}
	free(senders);

	struct recv_loglik_cohort *cohorts = ll->cohorts;
	size_t ic, nc = recv_model_cohort_count(ll->model);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cohorts[ic]);
	}
	free(cohorts);
}

void cohort_clear(struct recv_loglik_cohort *cll, size_t dim)
{
	assert(cll);
	info_clear(&cll->info, dim);
	cll->info.dev = 0.0;
	cll->info.nsend = 0;
	cll->info.nrecv = 0;
	cll->info_cached = true;
}

void recv_loglik_clear(struct recv_loglik *ll)
{
	assert(ll);

	const struct recv_model *m = ll->model;
	const struct design *design = recv_model_design(m);
	size_t dyn_dim = design_dvars_dim(design);
	size_t dim = recv_model_dim(m);
	struct recv_loglik_cohort *cohorts = ll->cohorts;
	size_t ic, nc = recv_model_cohort_count(m);

	for (ic = 0; ic < nc; ic++) {
		cohort_clear(&cohorts[ic], dim);
	}

	struct recv_loglik_sender *senders = ll->senders;
	size_t isend, nsend = recv_model_send_count(m);
	for (isend = 0; isend < nsend; isend++) {
		sender_clear(&senders[isend], dyn_dim);
	}

	ll->last = NULL;
}

void cohort_add(struct recv_loglik_cohort *ll, size_t nto)
{
	ll->info.nsend += 1;
	ll->info.nrecv += nto;
	ll->info_cached = false;
}

void recv_loglik_add(struct recv_loglik *ll,
		     const struct frame *f, const struct message *msg)
{
	size_t isend = msg->from;
	size_t c = recv_model_cohort(ll->model, isend);
	sender_add(&ll->senders[isend], f, msg->to, msg->nto);
	cohort_add(&ll->cohorts[c], msg->nto);

	ll->last = &ll->senders[isend];
}

void recv_loglik_add_all(struct recv_loglik *ll,
			 struct frame *f, const struct messages *msgs)
{
	struct messages_iter it;
	const struct message *msg;
	size_t i, n;

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

size_t recv_loglik_count(const struct recv_loglik *ll, size_t c)
{
	assert(ll);
	assert(c < recv_model_cohort_count(ll->model));
	return ll->cohorts[c].info.nrecv;
}

size_t recv_loglik_count_sum(const struct recv_loglik *ll)
{
	assert(ll);

	size_t n = 0;
	size_t ic, nc = recv_model_cohort_count(ll->model);
	for (ic = 0; ic < nc; ic++) {
		n += recv_loglik_count(ll, ic);
	}

	return n;
}

size_t recv_loglik_last_count(const struct recv_loglik *ll)
{
	assert(ll);
	assert(ll->last);
	return sender_last_count(ll->last);
}

static double recv_loglik_avg_dev_nocache(const struct recv_loglik *ll,
					  size_t c)
{
	assert(ll);
	assert(c < recv_model_cohort_count(ll->model));

	const size_t *cohorts = recv_model_cohorts(ll->model);
	size_t isend, nsend = recv_model_count(ll->model);
	size_t ntot = 0;
	double dev_avg = 0.0;

	for (isend = 0; isend < nsend; isend++) {
		if (cohorts[isend] != c)
			continue;

		struct recv_loglik_sender *sll = &ll->senders[isend];
		size_t n = sender_count(sll);
		if (n > 0) {
			double dev = sender_dev(sll) / n;
			dev_avg += n * (dev - dev_avg) / (n + ntot);
			ntot += n;
		}
	}
	assert(ntot == recv_loglik_count(ll, c));

	return dev_avg;
}

static void recv_loglik_axpy_avg_mean_nocache(double alpha,
					      const struct recv_loglik *ll,
					      size_t c, double *y)
{
	assert(ll);
	assert(c < recv_model_cohort_count(ll->model));

	const size_t *cohorts = recv_model_cohorts(ll->model);
	size_t isend, nsend = recv_model_count(ll->model);
	size_t dim = recv_model_dim(ll->model);
	double *avg_mean = xcalloc(dim, sizeof(avg_mean[0]));
	double *diff = xmalloc(dim * sizeof(diff[0]));
	size_t ntot, n;

	ntot = 0;

	for (isend = 0; isend < nsend; isend++) {
		if (cohorts[isend] != c)
			continue;

		struct recv_loglik_sender *sll = &ll->senders[isend];

		n = sender_count(sll);
		if (n > 0) {
			ntot += n;
			blas_dcopy(dim, avg_mean, 1, diff, 1);
			sender_axpy_mean(-1.0 / n, sll, diff);
			blas_daxpy(dim, -((double)n) / ntot, diff, 1, avg_mean, 1);
		}
	}
	assert(ntot == recv_loglik_count(ll, c));

	blas_daxpy(dim, alpha, avg_mean, 1, y, 1);
	free(diff);
	free(avg_mean);
}

static void recv_loglik_axpy_avg_score_nocache(double alpha,
					       const struct recv_loglik *ll,
					       size_t c, double *y)
{
	assert(ll);
	assert(c < recv_model_cohort_count(ll->model));

	const size_t *cohorts = recv_model_cohorts(ll->model);
	size_t isend, nsend = recv_model_count(ll->model);

	size_t dim = recv_model_dim(ll->model);
	double *avg_score = xcalloc(dim, sizeof(avg_score[0]));
	double *diff = xmalloc(dim * sizeof(diff[0]));
	size_t ntot, n;

	ntot = 0;

	for (isend = 0; isend < nsend; isend++) {
		if (cohorts[isend] != c)
			continue;

		struct recv_loglik_sender *sll = &ll->senders[isend];

		n = sender_count(sll);
		if (n > 0) {
			ntot += n;
			blas_dcopy(dim, avg_score, 1, diff, 1);
			sender_axpy_score(-1.0 / (double)n, sll, diff);
			blas_daxpy(dim, -((double)n) / ntot, diff, 1,
				avg_score, 1);
		}
	}
	assert(ntot == recv_loglik_count(ll, c));

	blas_daxpy(dim, alpha, avg_score, 1, y, 1);
	free(diff);
	free(avg_score);
}

static void recv_loglik_axpy_avg_imat_nocache(double alpha,
					      const struct recv_loglik *ll,
					      size_t c, struct matrix *y)
{
	assert(ll);
	assert(c < recv_model_cohort_count(ll->model));
	assert(y);
	assert((size_t)matrix_nrow(y) == recv_model_dim(ll->model));
	assert((size_t)matrix_ncol(y) == recv_model_dim(ll->model));

	const size_t *cohorts = recv_model_cohorts(ll->model);
	size_t isend, nsend = recv_model_count(ll->model);

	struct matrix avg_imat, new_imat, diff;
	size_t dim = recv_model_dim(ll->model);
	matrix_init(&avg_imat, dim, dim);
	matrix_init(&new_imat, dim, dim);
	matrix_init(&diff, dim, dim);
	size_t ntot, n;

	ntot = 0;

	for (isend = 0; isend < nsend; isend++) {
		if (cohorts[isend] != c)
			continue;

		struct recv_loglik_sender *sll = &ll->senders[isend];

		n = sender_count(sll);
		if (n > 0) {
			ntot += n;

			matrix_fill(&new_imat, 0.0);
			sender_axpy_imat(1.0 / (double)n, sll, &new_imat);
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

static void cache_info(struct recv_loglik *ll, size_t c)
{
	assert(ll);
	assert(c < recv_model_cohort_count(ll->model));

	size_t dim = recv_model_dim(ll->model);

	info_clear(&ll->cohorts[c].info, dim);
	ll->cohorts[c].info.dev += recv_loglik_avg_dev_nocache(ll, c);
	recv_loglik_axpy_avg_mean_nocache(1, ll, c, ll->cohorts[c].info.mean);
	recv_loglik_axpy_avg_score_nocache(1, ll, c,
					   ll->cohorts[c].info.score);
	recv_loglik_axpy_avg_imat_nocache(1, ll, c, &ll->cohorts[c].info.imat);
}

struct recv_loglik_info *recv_loglik_info(const struct recv_loglik *ll,
					  size_t c)
{
	assert(ll);
	assert(c < recv_model_cohort_count(ll->model));

	struct recv_loglik *mll = (struct recv_loglik *)ll;	// no const

	if (!mll->cohorts[c].info_cached) {
		cache_info(mll, c);
		mll->cohorts[c].info_cached = true;
	}

	return &mll->cohorts[c].info;
}

double recv_loglik_avg_dev(const struct recv_loglik *ll, size_t c)
{
	assert(ll);
	assert(c < recv_model_cohort_count(ll->model));

	const struct recv_loglik_info *info = recv_loglik_info(ll, c);
	return info->dev;
}

void recv_loglik_axpy_avg_mean(double alpha, const struct recv_loglik *ll,
			       size_t c, double *y)
{
	assert(c < recv_model_cohort_count(ll->model));

	size_t n = recv_model_dim(ll->model);
	const struct recv_loglik_info *info = recv_loglik_info(ll, c);
	blas_daxpy(n, alpha, info->mean, 1, y, 1);
}

void recv_loglik_axpy_avg_score(double alpha, const struct recv_loglik *ll,
				size_t c, double *y)
{
	assert(c < recv_model_cohort_count(ll->model));

	size_t n = recv_model_dim(ll->model);
	const struct recv_loglik_info *info = recv_loglik_info(ll, c);
	blas_daxpy(n, alpha, info->score, 1, y, 1);
}

void recv_loglik_axpy_avg_imat(double alpha, const struct recv_loglik *ll,
			       size_t c, struct matrix *y)
{
	assert(c < recv_model_cohort_count(ll->model));
	assert(y);
	assert((size_t)matrix_nrow(y) == recv_model_dim(ll->model));
	assert((size_t)matrix_ncol(y) == recv_model_dim(ll->model));

	const struct recv_loglik_info *info = recv_loglik_info(ll, c);
	matrix_axpy(alpha, &info->imat, y);
}

double recv_loglik_last_dev(const struct recv_loglik *ll)
{
	assert(ll->last);

	return sender_last_dev(ll->last);
}

void recv_loglik_axpy_last_mean(double alpha, const struct recv_loglik *ll,
				double *y)
{
	assert(ll->last);
	sender_axpy_last_mean(alpha, ll->last, y);
}

void recv_loglik_axpy_last_score(double alpha, const struct recv_loglik *ll,
				 double *y)
{
	assert(ll->last);
	sender_axpy_last_score(alpha, ll->last, y);
}

void recv_loglik_axpy_last_imat(double alpha, const struct recv_loglik *ll,
				struct matrix *y)
{
	assert(ll->last);
	sender_axpy_last_imat(alpha, ll->last, y);
}
