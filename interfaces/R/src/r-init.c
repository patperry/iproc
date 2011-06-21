#include "port.h"
#include <assert.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "r-actors.h"
#include "r-design.h"
#include "r-fit.h"
#include "r-messages.h"
#include "r-utils.h"

void R_init_iproc(DllInfo * info)
{
	Riproc_actors_init(info);
	Riproc_design_init(info);
	Riproc_fit_init(info);
	Riproc_messages_init(info);
	Riproc_utils_init(info);
}

void R_unload_iproc(DllInfo * info)
{

}
