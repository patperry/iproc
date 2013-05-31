#ifndef ENRON_MESSAGES_H
#define ENRON_MESSAGES_H

#include <stddef.h>
#include <stdio.h>

int enron_messages_init(size_t *nsend, size_t *nrecv, double **time,
			size_t **from, size_t ***to, size_t **nto,
			intptr_t **attr, size_t *nmsg, size_t maxrecip);
int enron_messages_init_fread(size_t *nsend, size_t *nrecv, double **time,
			      size_t **from, size_t ***to, size_t **nto,
			      intptr_t **attr, size_t *nmsg, size_t maxrecip,
			      FILE *stream);

#endif /* ENRON_MESSAGES_H */
