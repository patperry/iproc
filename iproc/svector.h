#ifndef _IPROC_SVECTOR_H
#define _IPROC_SVECTOR_H

#include <stdint.h>
#include <iproc/array.h>
#include <iproc/refcount.h>
#include <iproc/vector.h>

typedef struct _iproc_svector iproc_svector;

struct _iproc_svector {
    int64_t        dim;
    iproc_array   *index;
    iproc_array   *value;
    iproc_refcount refcount;
};

iproc_svector * iproc_svector_new      (int64_t        dim);
iproc_svector * iproc_svector_new_copy (iproc_svector *svector);
void            iproc_svector_clear    (iproc_svector *svector);
iproc_svector * iproc_svector_ref      (iproc_svector *svector);
void            iproc_svector_unref    (iproc_svector *svector);

int64_t         iproc_svector_dim      (iproc_svector *svector);
double          iproc_svector_get      (iproc_svector *svector,
                                        int64_t        i);
void            iproc_svector_set      (iproc_svector *svector,
                                        int64_t        i,
                                        double         value);
void            iproc_svector_inc      (iproc_svector *svector,
                                        int64_t        i,
                                        double         value);
void            iproc_svector_scale    (iproc_svector *svector,
                                        double         scale);

int64_t         iproc_svector_nnz      (iproc_svector *svector);
int64_t         iproc_svector_nz       (iproc_svector *svector,
                                        int64_t        i);
double          iproc_svector_nz_val   (iproc_svector *svector,
                                        int64_t        i);
double          iproc_svector_nz_inc   (iproc_svector *svector,
                                        int64_t        i,
                                        double         value);
double          iproc_svector_nz_shift (iproc_svector *svector,
                                        double         value);
double          iproc_svector_nz_max   (iproc_svector *svector);

double          iproc_vector_sdot      (iproc_vector  *vector,
                                        iproc_svector *svector);
void            iproc_vector_sacc      (iproc_vector  *dst_vector,
                                        double         scale,
                                        iproc_svector *svector);



#endif _IPROC_SVECTOR_H
