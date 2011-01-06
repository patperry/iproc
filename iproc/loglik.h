#ifndef _IPROC_LOGLIK_H
#define _IPROC_LOGLIK_H

#include <iproc/array.h>
#include <iproc/model.h>

typedef struct _iproc_loglik      iproc_loglik;
typedef struct _iproc_send_loglik iproc_send_loglik;

struct _iproc_loglik {
    iproc_array *send_loglik;
};

struct _iproc_send_loglik {
    iproc_model  *model;
    int64_t       isend;
    int64_t       nsend;
    int64_t       nrecv;
    iproc_vector *ovarsdiff;
    double        value;
    iproc_vector *evarsdiff;
};

#endif /* _IPROC_LOGLIK_H */
