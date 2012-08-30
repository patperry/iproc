#ifndef _ENRON_H
#define _ENRON_H

#include <stdio.h>
#include "messages.h"

#define PARENT_DIR "/Users/patperry/Projects/iproc/"
#define ENRON_EMPLOYEES_FILE PARENT_DIR"tests/enron/employees.json"
#define ENRON_MESSAGES_FILE  PARENT_DIR"tests/enron/messages.json"

enum enron_trait {
	ENRON_TRAIT_L = 0,
	ENRON_TRAIT_T,
	ENRON_TRAIT_J,
	ENRON_TRAIT_F,
	ENRON_TRAIT_LJ,
	ENRON_TRAIT_TJ,
	ENRON_TRAIT_LF,
	ENRON_TRAIT_TF,
	ENRON_TRAIT_JF,
	ENRON_TRAIT_LJF,
	ENRON_TRAIT_TJF
};


#define ENRON_TERMS_MAX 3

int enron_employees_init(size_t *nactorp,
			 double **traitsp, size_t *ntraitp,
			 const char * const **trait_namesp,
			 size_t terms);
int enron_employees_init_fread(size_t *nactorp,
			       double **traitsp, size_t *ntraitp,
			       const char * const **trait_namesp,
			       size_t terms, FILE *stream);

int enron_messages_init(struct messages *messages, size_t maxrecip);
int enron_messages_init_fread(struct messages *messages, size_t maxrecip, FILE *stream);

#endif /* _ENRON_H */
