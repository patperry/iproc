#ifndef _LOGSUMEXP_H
#define _LOGSUMEXP_H

struct logsumexp {
	double max;
	double sumexpm1;
};

void logsumexp_init(struct logsumexp *lse);
void logsumexp_insert(struct logsumexp *lse, double value);
double logsumexp_max(const struct logsumexp *lse);
double logsumexp_value(const struct logsumexp *lse);

#endif /* _LOGSUMEXP_H */
