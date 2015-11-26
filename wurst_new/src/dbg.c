/*
 * Debug hook to bring in debugger
 */
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>

#include "dbg.h"

#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: dbg.c,v 1.1.1.1 2001/09/21 05:19:52 torda Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */


#define need_breaker

pid_t getpid (void);
void
breaker (void)
{
#   ifdef need_breaker
    int kill (pid_t x, int sig);/* with ansi flags set, prototype explicitly */
    int sigignore( int sig);
    sigignore (SIGTRAP);
    kill (getpid(), SIGTRAP);
#   endif /* need_breaker */
}
    
