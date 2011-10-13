#ifndef _RECV_BOOT_H
#define _RECV_BOOT_H

/* must #include dSFMT/dSFMT.h before including this file */

#define DSFMT_MEXP 19937
#define DSFMT_DO_NOT_USE_OLD_NAMES
#include "dSFMT/dSFMT.h"

#include "messages.h"
#include "design.h"
#include "matrix.h"
#include "recv_model.h"

struct recv_boot {
	struct frame frame;
	struct recv_model model;
	struct messages messages;
};

void recv_boot_init(struct recv_boot *boot,
		    const struct messages *msgs,
		    const struct design *design,
		    size_t ncohort,
		    const size_t *cohorts,
		    const struct matrix *coefs, dsfmt_t * dsfmt);

void recv_boot_deinit(struct recv_boot *boot);

#endif /* _RECV_BOOT_H */
