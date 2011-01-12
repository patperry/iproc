#ifndef _IPROC_R_ACTORS_H
#define _IPROC_R_ACTORS_H

#include <Rinternals.h>
#include "actors.h"

/* Call once to initialize library */
void           Riproc_actors_init ();

/* External functions */
SEXP           Riproc_actors_new          (SEXP Rclasses,
                                           SEXP Rclass_vectors_t);
SEXP           Riproc_actors_nclass       (SEXP Ractors);
SEXP           Riproc_actors_size         (SEXP Ractors);
SEXP           Riproc_actors_dim          (SEXP Ractors);
SEXP           Riproc_actors_vector       (SEXP Ractors,
                                           SEXP Ractor_ids);
SEXP           Riproc_actors_class        (SEXP Ractors,
                                           SEXP Ractor_ids);
SEXP           Riproc_actors_class_vector (SEXP Ractors,
                                           SEXP Rclass_ids);

/* Internal use only */
iproc_actors * Riproc_to_actors          (SEXP Ractors);
SEXP           Riproc_from_actors        (iproc_actors *actors);

#endif /* _IPROC_R_ACTORS_H */


