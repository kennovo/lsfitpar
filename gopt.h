/* $Id: gopt.h,v 1.1 2015/05/01 23:17:59 kenno Exp $ */
/* slightly modified version of
 * gopt.h version 8.1 by Tom Vajzovic <tom.viza@gmail.com> PUBLIC DOMAIN 2003-8
 * as downloaded from the author's website http://purposeful.co.uk/software/gopt/
 * the 17th of December, 2009 */

#ifndef GOPT_H_INCLUDED
#define GOPT_H_INCLUDED

#define GOPT_ONCE   0
#define GOPT_REPEAT 1
#define GOPT_NOARG  0
#define GOPT_ARG    2

#define gopt_start(...)  (const void*)( const struct { int k; int f; const char *s; const char*const*l; }[]){ __VA_ARGS__, {0}}
#define gopt_option(k,f,s,l)    { k, f, s, l }
#define gopt_shorts( ... )      (const char*)(const char[]){ __VA_ARGS__, 0 }
#define gopt_longs( ... )       (const char**)(const char*[]){ __VA_ARGS__, NULL }

struct opt_s {
  int key;
  const char *arg;
  const char *arg2;
};
typedef struct opt_s opt_t;

void *gopt_sort( int *argc, const char **argv, const void *opt_specs );
/* returns a pointer for use in the following calls
 * prints to stderr and call exit() on error
 */
size_t gopt( const void *opts, int key );
/* returns the number of times the option was specified
 * which will be 0 or 1 unless GOPT_REPEAT was used
 */
size_t gopt_arg( const void *opts, int key, const char **arg );
/* returns the number of times the option was specified
 * writes a pointer to the option argument from the first (or only) occurance to *arg
 */
void gopt_free( void *opts );
/* releases memory allocated in the corresponding call to gopt_sort()
 * opts can no longer be used 
 */
#endif /* GOPT_H_INCLUDED */
