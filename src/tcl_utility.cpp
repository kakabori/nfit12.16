#include <tcl.h>
#include <tk.h>
#include <string>
#include <cstring>
#include "tcl_utility.h"
#include "globalVariables.h"

using std::string;


void evalTclCommand(const char *comm)
{
  g_string.assign(comm);
  Tcl_AsyncMark(g_commandHandler);
}


int evalTclCommand_async(ClientData clientData, Tcl_Interp *interp, int code)
{
  interp = g_interp;
  // copy content of g_string to a writable c-string so that it can be
  // safely passed to Tcl_EvalEx
  char cmd[g_string.length()+1];
  strcpy(cmd, g_string.c_str());
  Tcl_EvalEx(interp, cmd, -1, TCL_EVAL_GLOBAL);
  return code;
}

