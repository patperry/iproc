#ifndef _IPROC_LOGSUMEXP_H
#define _IPROC_LOGSUMEXP_H

typedef struct _iproc_logsumexp iproc_logsumexp;

struct _iproc_logsumexp {
	double max;
	double sumexpm1;
};

void iproc_logsumexp_init(iproc_logsumexp * lse);

void iproc_logsumexp_insert(iproc_logsumexp * lse, double value);
double iproc_logsumexp_max(iproc_logsumexp * lse);
double iproc_logsumexp_value(iproc_logsumexp * lse);

#endif /* _IPROC_LOGSUMEXP_H */
