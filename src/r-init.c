
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "r-actors.h"


void
R_init_iproc (DllInfo *info)
{
    Riproc_actors_init(info);
}

void
R_unload_iproc (DllInfo *info)
{
    
}
