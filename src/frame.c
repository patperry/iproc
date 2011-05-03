#include "port.h"
#include <assert.h>
#include "frame.h"


double frame_time(const struct frame *f)
{
	assert(f);
	return f->time;
}


bool frame_reserve_dyad_events(struct frame *f, ssize_t nadd)
{
	assert(f);
	assert(nadd >= 0);

	return pqueue_reserve_push(&f->dyad_var_diffs, nadd);
}

bool frame_add_dyad_event(struct frame *f, const struct dyad_var_diff *delta)
{
	assert(f);
	assert(delta);
	assert(!(delta->time <= frame_time(f)));
	
	return pqueue_push(&f->dyad_var_diffs, delta);
}
