#ifndef RECV_BOOT_H
#define RECV_BOOT_H

/* must #include dSFMT/dSFMT.h before including this file */

#define DSFMT_MEXP 19937
#define DSFMT_DO_NOT_USE_OLD_NAMES
#include "dSFMT/dSFMT.h"

#include "messages.h"
#include "design.h"
#include "recv_model.h"

struct recv_boot {
	struct recv_model model;
	struct messages messages;
};

void recv_boot_init(struct recv_boot *boot,
		    struct frame *f,
		    const struct messages *msgs,
		    const struct recv_coefs *coefs, dsfmt_t * dsfmt);
void recv_boot_deinit(struct recv_boot *boot);

#endif /* RECV_BOOT_H */