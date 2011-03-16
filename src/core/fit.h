#ifndef _IPROC_FIT_H
#define _IPROC_FIT_H

#include "design.h"
#include "loglik.h"
#include "matrix.h"
#include "messages.h"
#include "model.h"
#include "vector.h"

typedef struct _iproc_fit iproc_fit;

struct _iproc_fit {
    iproc_model *model;
    iproc_messages *messages;
    double penalty;
    iproc_loglik *loglik;
    double        value;
    double        value0;
    iproc_vector *grad;
    iproc_vector *grad0;
    iproc_matrix *inv_hess;
    iproc_vector *search_dir;
    double        step;
};

iproc_fit * iproc_fit_new  (iproc_model    *model0,
                            iproc_messages *messages,
                            double          penalty);
void        iproc_fit_free (iproc_fit      *fit);

void        iproc_fit_step (iproc_fit      *fit);

#endif /* _IPROC_FIT_H */