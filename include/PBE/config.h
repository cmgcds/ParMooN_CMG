#ifndef __CONFIG_H
#define __CONFIG_H
/**
 **
 ** provides basic types and configuration settings
 **
 ** arch-tag: 3c370550-398d-496c-ac9e-b4b7cb79c45f
 **/

/*
 * included header
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include "seconds.h"

/*
 * basic types
 */

typedef enum { 
	FEPC_FALSE, 
	FEPC_TRUE 
} bool_t;

typedef double  fepc_real_t;

/*
 * consistency check
 */

#if !defined(NDEBUG)
#  define ASSERT( cond )   assert( cond )
#else
#  define ASSERT( cond )
#endif

/*
 * enables debugging/output
 */

#define DEBUG  if ( false )
#define PRINT  if ( false )

/*
 * inline function decl.
 */

#if defined(__USE_ISOC99) && ( defined(__GNUC__) || defined(__ICC) || defined(__ECC) )
#define _INLINE_  inline
#else
#define _INLINE_  static
#endif

#endif  /* __CONFIG_H */
