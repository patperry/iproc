#ifndef _LINALG_H
#define _LINALG_H

#include "vector.h"
#include "matrix.h"

enum matrix_uplo {
	UPLO_UPPER,
	UPLO_LOWER
};

enum eig_job {
	EIG_NOVEC,
	EIG_VEC
};

struct symeig {
	ssize_t n;
	enum eig_job job;
	double *work;
	f77int lwork;
	f77int *iwork;
	f77int liwork;
	f77int info;
};


// replace b with a \ b; replace a with its cholesky factor
// returns 0 on success, i > 0 if leading minor of order i is not pd
ssize_t chol_solve(enum matrix_uplo uplo, struct matrix *a, struct matrix *b);



void symeig_init(struct symeig *eig, ssize_t n, enum eig_job job);
void symeig_reinit(struct symeig *eig, ssize_t n, enum eig_job job);
void symeig_deinit(struct symeig *eig);

bool symeig_factor(struct symeig *eig, enum matrix_uplo uplo,
		   struct matrix *a, struct vector *w); // destroys a and w

static inline ssize_t symeig_dim(const struct symeig *eig)
{
	assert(eig);	
	return eig->n;
}

static inline enum eig_job symeig_job(const struct symeig *eig)
{
	assert(eig);
	return eig->job;
}

#endif /* _LINALG_H */
