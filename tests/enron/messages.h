#ifndef ENRON_MESSAGES_H
#define ENRON_MESSAGES_H

#include <stddef.h>
#include <stdio.h>

int enron_messages_init(size_t maxrecip,
			double **time, size_t **from, size_t ***to, size_t **nto,
			intptr_t **attr, size_t *nmsg);
int enron_messages_init_fread(FILE *stream, size_t maxrecip,
			      double **time, size_t **from, size_t ***to, size_t **nto,
			      intptr_t **attr, size_t *nmsg);

#endif /* ENRON_MESSAGES_H */
