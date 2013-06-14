#ifndef ENRON_ACTORS_H
#define ENRON_ACTORS_H

#include <stddef.h>
#include <stdio.h>


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

#define ENRON_ACTOR_COUNT 156
#define ENRON_TERMS_MAX 3

int enron_actors_init(size_t terms, double **trait_x,
		      const char * const **trait_names, size_t *trait_dim);

int enron_actors_init_fread(FILE *stream, size_t terms, double **trait_x,
			    const char * const **trait_names,
			    size_t *trait_dim);


#endif /* ENRON_ACTORS_H */
