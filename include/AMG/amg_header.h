/****************************************************************************/
/*                                                                                                                                                        */
/* File:          amg_header.c                                                                                                        */
/*                                                                                                                                                        */
/* Purpose:   general header for common things (return values, misc..)                */
/*                                                                                                                                                        */
/* Author:          Peter Bastian                                                                                                         */
/*                          Institut fuer Computeranwendungen III                                                 */
/*                          Universitaet Stuttgart                                                                                */
/*                          Pfaffenwaldring 27                                                                                        */
/*                          70550 Stuttgart                                                                                                */
/*                          email: peter@ica3.uni-stuttgart.de                                                        */
/*                          phone: 0049-(0)711-685-7003                                                                        */
/*                          fax  : 0049-(0)711-685-7000                                                                        */
/*                                                                                                                                                        */
/* History:   28 Jan 1996 Begin                                                                                                */
/*            02 Apr 1996 new memory allocation strategy                                        */
/*            30 Sep 1997 redesign                                                                                        */
/*                                                                                                                                                        */
/* Remarks:                                                                                                                                 */
/*                                                                                                                                                        */
/****************************************************************************/

/* RCS_ID
$Header: /homes/matthies/ARCHIVE_MooNMD/MooNMD/include/AMG/amg_header.h,v 1.1.1.1 2001/01/12 12:55:18 matthies Exp $
*/

/****************************************************************************/
/*                                                                                                                                                        */
/* auto include mechanism and other include files                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

#ifndef __AMGHEADER__
#define __AMGHEADER__

/****************************************************************************/
/*                                                                                                                                                        */
/* defines in the following order                                                                                        */
/*                                                                                                                                                        */
/*                  compile time constants defining static data size (i.e. arrays)        */
/*                  other constants                                                                                                        */
/*                  macros                                                                                                                        */
/*                                                                                                                                                        */
/****************************************************************************/

/* general sizes */
#define AMG_NAME_SIZE                        32  /* for names of objects                                        */

/* general return values*/
#define AMG_OK                                        0        /* operation succeded                                        */
#define AMG_NULL                                NULL/* null pointer                                                        */
#define AMG_FATAL                                9999/* fatal error                                                        */

/* misc macros */
#define AMG_MIN(x,y)            (((x)<(y)) ? (x) : (y))
#define AMG_MAX(x,y)            (((x)>(y)) ? (x) : (y))
/* #define AMG_ABS(x)                    (((x)>=0) ?  (x) : (-x))*/
#define AMG_ABS(i) (((i)<0) ? (-(i)) : (i))

#endif

