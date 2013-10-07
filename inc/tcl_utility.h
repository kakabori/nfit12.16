#ifndef GUARD_TCL_UTILITY_H
#define GUARD_TCL_UTILITY_H

#include <tcl.h>
#include <tk.h>
#include <string>

void evalTclCommand(const char *);
int evalTclCommand_async(ClientData, Tcl_Interp *, int);

#endif
