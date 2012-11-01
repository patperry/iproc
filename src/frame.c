#include "port.h"
#include <assert.h>
#include "frame.h"


void frame_init(struct frame *f, size_t nsend, size_t nrecv, int has_loops,
		const double *intvls, size_t nintvl)
{
	assert(f);
	assert(intvls || !nintvl);

	f->has_loops = has_loops;
	history_init(&f->history, nsend, nrecv, intvls, nintvl);
	design_init(&f->send_design, &f->history, nsend);
	design_init(&f->recv_design, &f->history, nrecv);
	design2_init(&f->dyad_design, &f->history, nsend, nrecv);
}

void frame_deinit(struct frame *f)
{
	assert(f);

	design2_deinit(&f->dyad_design);
	design_deinit(&f->recv_design);
	design_deinit(&f->send_design);
	history_deinit(&f->history);
}
