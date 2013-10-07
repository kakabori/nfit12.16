#include <tcl.h>
#include <tk.h>
#include <string>
#include <cstring>
#include "tcl_utility.h"
#include "globalVariables.h"
#include "fileTools.h"

using std::string;


/******************************************************************************
This function is used in a thread different from the main thread handling
the Tk event loop. g_flag condition ensures that the last call for 
evalTclCommand_async gets processed before a new one overwrites 
g_commandHanlder state.

I should come up with a more efficient way to handle this issue.
******************************************************************************/
void evalTclCommand(const char *comm)
{
  g_string.assign(comm);
  Tcl_AsyncMark(g_commandHandler);
}


/******************************************************************************
This function is used in the main thread handling the Tk event loop.
Do not use it in a child thread b/c g_flag condition will cause a conflict.
It should be used only in the context of Tcl_AsyncCreate.
******************************************************************************/
int evalTclCommand_async(ClientData clientData, Tcl_Interp *interp, int code)
{
  if (interp == NULL) interp = g_interp;
  // copy content of g_string to a writable c-string so that it can be
  // safely passed to Tcl_EvalEx
  char cmd[g_string.length()+1];
  strcpy(cmd, g_string.c_str());
  Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL);
  return code;
}

