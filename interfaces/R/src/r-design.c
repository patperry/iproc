#include "port.h"
#include <assert.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "actors.h"
#include "vrecv.h"
#include "vnrecv.h"
#include "r-utils.h"
#include "r-actors.h"
#include "r-design.h"

static SEXP Riproc_design_type_tag;

static R_CallMethodDef callMethods[] = {
	{"Riproc_design_new", (DL_FUNC) & Riproc_design_new, 4},
	{"Riproc_design_dim", (DL_FUNC) & Riproc_design_dim, 1},
	{"Riproc_design_nreceiver", (DL_FUNC) & Riproc_design_nreceiver, 1},
	{"Riproc_design_nsender", (DL_FUNC) & Riproc_design_nsender, 1},
	{"Riproc_design_receivers", (DL_FUNC) & Riproc_design_receivers, 1},
	{"Riproc_design_senders", (DL_FUNC) & Riproc_design_senders, 1},
	{NULL, NULL, 0}
};

struct design_params {
	bool loops;
	bool reffects;
	bool vrecv;
	bool vnrecv;
};

static struct design *design_alloc_params(struct actors *senders,
					  struct actors *receivers,
					  const struct vector *intervals,
					  const struct design_params *params)
{
	assert(senders);
	assert(receivers);
	assert(intervals);
	assert(params);

	struct design *design = design_alloc(senders, receivers, intervals);

	design_set_loops(design, params->loops);
	design_set_recv_effects(design, params->reffects);
	if (params->vrecv) {
		design_add_recv_var(design, VAR_TYPE_RECV);
	}
	if (params->vnrecv) {
		design_add_recv_var(design, VAR_TYPE_NRECV);
	}
	return design;
}

void Riproc_design_init(DllInfo * info)
{
	Riproc_design_type_tag = install("Riproc_design_type_tag");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void Riproc_design_free(SEXP Rdesign)
{
	struct design *design =
	    Riproc_sexp2ptr(Rdesign, TRUE, Riproc_design_type_tag,
			    "design");

	design_free(design);
}

struct design *Riproc_to_design(SEXP Rdesign)
{
	struct design *design =
	    Riproc_sexp2ptr(Rdesign, TRUE, Riproc_design_type_tag,
			    "design");
	return design;
}

SEXP Riproc_from_design(struct design * design)
{
	SEXP Rdesign, class;

	design_ref(design);

	PROTECT(Rdesign =
		R_MakeExternalPtr(design, Riproc_design_type_tag, R_NilValue));
	R_RegisterCFinalizer(Rdesign, Riproc_design_free);

	/* set the class of the result */
	PROTECT(class = allocVector(STRSXP, 1));
	SET_STRING_ELT(class, 0, mkChar("iproc.design"));
	classgets(Rdesign, class);

	UNPROTECT(2);
	return Rdesign;
}

SEXP
Riproc_design_new(SEXP Rsenders,
		  SEXP Rreceivers,
		  SEXP Rreceiver_effects, SEXP Rrecip_intervals)
{
	struct actors *senders = Riproc_to_actors(Rsenders);
	struct actors *receivers = Riproc_to_actors(Rreceivers);
	Rboolean receiver_effects = LOGICAL_VALUE(Rreceiver_effects);
	bool has_loops = true;

	if (receiver_effects == NA_LOGICAL) {
		error("'receiver.effects' value is NA");
	}

	struct vector *intervals = malloc(sizeof(*intervals));
	vector_init(intervals, 0);

	struct design_params params;
	params.loops = has_loops;
	params.reffects = receiver_effects;
	params.vrecv = false;
	params.vnrecv = false;

	struct design *design =
	    design_alloc_params(senders, receivers, intervals, &params);

	SEXP Rdesign;

	PROTECT(Rdesign = Riproc_from_design(design));
	design_free(design);

	UNPROTECT(1);
	return Rdesign;
}

SEXP Riproc_design_dim(SEXP Rdesign)
{
	const struct design *design = Riproc_to_design(Rdesign);
	int dim = (int)design_recv_dim(design);
	return ScalarInteger(dim);
}

SEXP Riproc_design_nsender(SEXP Rdesign)
{
	const struct design *design = Riproc_to_design(Rdesign);
	int nsender = (int)design_send_count(design);
	return ScalarInteger(nsender);
}

SEXP Riproc_design_nreceiver(SEXP Rdesign)
{
	const struct design *design = Riproc_to_design(Rdesign);
	int nreceiver = (int)design_recv_count(design);
	return ScalarInteger(nreceiver);
}

SEXP Riproc_design_senders(SEXP Rdesign)
{
	const struct design *design = Riproc_to_design(Rdesign);
	struct actors *senders = design_senders(design);
	return Riproc_from_actors(senders);
}

SEXP Riproc_design_receivers(SEXP Rdesign)
{
	const struct design *design = Riproc_to_design(Rdesign);
	struct actors *receivers = design_receivers(design);
	return Riproc_from_actors(receivers);
}
