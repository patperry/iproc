
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "r-actors.h"
#include "r-cursor.h"
#include "r-frame.h"
#include "r-messages.h"
#include "r-model.h"
#include "r-utils.h"



void
R_init_iproc (DllInfo *info)
{
    Riproc_actors_init(info);
    Riproc_cursor_init(info);
    Riproc_frame_init(info);
    Riproc_messages_init(info);
    Riproc_model_init(info);
    Riproc_utils_init(info);
}

void
R_unload_iproc (DllInfo *info)
{
    
}
