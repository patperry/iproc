#ifndef _ENRON_H
#define _ENRON_H

#include <stdio.h>
#include "actors.h"
#include "matrix.h"
#include "messages.h"

#define PARENT_DIR "/Users/patperry/Projects/iproc/"
#define ENRON_EMPLOYEES_FILE PARENT_DIR"tests/enron/employees.json"
#define ENRON_MESSAGES_FILE  PARENT_DIR"tests/enron/messages.json"

enum enron_cohort {
	ENRON_COHORT_LJF = 0,
	ENRON_COHORT_LJM,
	ENRON_COHORT_LSF,
	ENRON_COHORT_LSM,
	ENRON_COHORT_TJF,
	ENRON_COHORT_TJM,
	ENRON_COHORT_TSF,
	ENRON_COHORT_TSM,
	ENRON_COHORT_OJF,
	ENRON_COHORT_OJM,
	ENRON_COHORT_OSF,
	ENRON_COHORT_OSM
};

#define ENRON_NCOHORT (ENRON_COHORT_OSM + 1)

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
	ENRON_TRAIT_TJF,	
};

#define ENRON_NTRAIT (ENRON_TRAIT_TJF + 1)

extern bool enron_legal[ENRON_NCOHORT];
extern bool enron_trading[ENRON_NCOHORT];
extern bool enron_other[ENRON_NCOHORT];

extern bool enron_junior[ENRON_NCOHORT];
extern bool enron_senior[ENRON_NCOHORT];

extern bool enron_female[ENRON_NCOHORT];
extern bool enron_male[ENRON_NCOHORT];


bool enron_employees_init(struct actors *employees, struct matrix *traits);
bool enron_employees_init_fread(struct actors *employees, struct matrix *traits, FILE *stream);

bool enron_messages_init(struct messages *messages, ssize_t maxrecip);
bool enron_messages_init_fread(struct messages *messages, ssize_t maxrecip, FILE *stream);

#endif /* _ENRON_H */
