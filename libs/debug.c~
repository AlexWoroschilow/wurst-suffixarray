
/*
 *  debug.c
 *  debug messages
 *
 *  @author Steve Hoffmann
 *  @email steve@bioinf.uni-leipzig.de
 *  @date 08/26/2007 06:49:02 PM CEST
 *  
 */
 
 #include <stdarg.h>
 #include <stdio.h>
 #include <string.h>
 #include "debug.h"

 int
 debugmsg( char *file, 
           int line, 
           const char *fmt, ...) {

   int ret;
   va_list ap;
   va_start(ap, fmt);
   fprintf(stderr, "[%s] file: %s, line: %d: ", "segemehl", file, line);
   ret = vfprintf(stderr, fmt, ap);
   va_end(ap);

   return ret;
 
 }



