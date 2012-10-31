#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "sblas.h"
#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "vars.h"
#include "frame.h"


void frame_init(struct frame *f, size_t nsend, size_t nrecv, int has_loops,
		const double *intvls, size_t nintvl)
{
	assert(f);
	assert(intvls || !nintvl);

	f->nsend = nsend;
	f->nrecv = nrecv;
	f->has_loops = has_loops;

	history_init(&f->history, nsend, nrecv, intvls, nintvl);
	design_init(&f->send_design, f, nsend);
	design_init(&f->recv_design, f, nrecv);
	design2_init(&f->dyad_design, f, nsend, nrecv);
}

void frame_deinit(struct frame *f)
{
	assert(f);

	design2_deinit(&f->dyad_design);
	design_deinit(&f->recv_design);
	design_deinit(&f->send_design);
	history_deinit(&f->history);
}
