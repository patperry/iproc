#include "port.h"
#include <assert.h>
#include <string.h>
#include "frame.h"
#include "vars.h"

static void nsendtot_message_add(void *udata, struct history *h,
				 const struct message *msg)
{
	if (!history_interval_count(h))
		return;

	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = v->design;

	size_t isend = msg->from;
	double nto = (double)(msg->nto);
	double dx_data[1] = { nto };
	size_t dx_index[1] = { 0 };
	size_t nz = 1;
	
	design_update(d, v, isend, dx_data, dx_index, nz);
}


static void nsendtot_message_advance(void *udata, struct history *h,
				     const struct message *msg, size_t intvl)
{
	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = v->design;

	size_t isend = msg->from;
	double nto = (double)(msg->nto);
	double dx_data[2] = { -nto, +nto };
	size_t dx_index[2];
	size_t nz = 2;
	
	size_t ix1 = intvl;
	size_t ix0 = ix1 - 1;
		
	dx_index[0] = ix0;
	dx_index[1] = ix1;
		
	design_update(d, v, isend, dx_data, dx_index, nz);
}


static struct history_callbacks nsendtot_history_callbacks = {
	nsendtot_message_add,
	nsendtot_message_advance,
	NULL
};


static void nsendtot_init(struct tvar *tv, struct history *h, va_list ap)
{
	(void)ap;		// unused;

	size_t n = history_interval_count(h);

	tv->var.rank = 1;
	tv->var.dims[0] = n;
	tv->udata = NULL;
	
	history_add_observer(h, tv, &nsendtot_history_callbacks);
}


static void nsendtot_deinit(struct tvar *tv, struct history *h)
{
	history_remove_observer(h, tv);
}


static struct tvar_type VAR_NSENDTOT_REP = {
	nsendtot_init,
	nsendtot_deinit,
	
};


const struct tvar_type *VAR_NSENDTOT = &VAR_NSENDTOT_REP;
