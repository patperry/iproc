
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "r-actors.h"
#include "r-model.h"
#include "r-vars.h"


void
R_init_iproc (DllInfo *info)
{
    Riproc_actors_init(info);
    Riproc_model_init(info);
    Riproc_vars_init(info);
}

void
R_unload_iproc (DllInfo *info)
{
    
}
