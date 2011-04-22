#ifndef _IPROC_R_ACTORS_H
#define _IPROC_R_ACTORS_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "actors.h"

/* Call once to initialize library */
void Riproc_actors_init(DllInfo * info);

/* External functions */
SEXP Riproc_actors_new(SEXP Rtraits_t);
SEXP Riproc_actors_size(SEXP Ractors);
SEXP Riproc_actors_dim(SEXP Ractors);
SEXP Riproc_actors_traits(SEXP Ractors, SEXP Ractor_ids);

SEXP Riproc_actors_mul(SEXP Ractors, SEXP Rmatrix);
SEXP Riproc_actors_tmul(SEXP Ractors, SEXP Rmatrix);

/* Internal use only */
struct actors *Riproc_to_actors(SEXP Ractors);
SEXP Riproc_from_actors(struct actors *actors);

#endif /* _IPROC_R_ACTORS_H */
