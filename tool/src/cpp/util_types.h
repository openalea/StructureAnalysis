/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture 
 *
 *       Copyright 1995-2010 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): C. Nouguier & F. Boudon (frederic.boudon@cirad.fr) nouguier 
 *
 *       $Source$
 *       $Id: util_types.h 9869 2010-11-04 17:31:41Z dbarbeau $
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
 *               
 *  ----------------------------------------------------------------------------
 * 
 *                      GNU General Public Licence
 *           
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */				


				





#ifndef __util_types_h__
#define __util_types_h__
/*! 
	\file util_types.h
    \brief File that contains some utility about type. define the type real_t.
*/

#include "config.h"
#include "tools_namespace.h"


#ifdef _WIN64
#define __WORDSIZE  64
#endif

// Basic types redefinition
#ifdef __GNUC__
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#if defined( __MINGW32__ )
#define __int8_t_defined
#define __uint_t_defined
#endif // __MINGW32__
#endif

#include <climits>
/// type uchar_t (unsigned char)
typedef unsigned char uchar_t;

#ifndef __int8_t_defined

/// type int16_t (short)
typedef short int16_t;
/// int16_t min value
#ifndef INT16_MIN 
#define INT16_MIN SHRT_MIN
#endif
/// int16_t max value
#ifndef INT16_MAX
#define INT16_MAX SHRT_MAX
#endif

#endif

#ifndef __uint_t_defined

/// type uint16_t (unsigned short)
typedef unsigned short uint16_t;

/// type int32_t (long)
typedef int int32_t;

/// type uint_t (unsigned long)
typedef unsigned int uint32_t;

#endif

/// uint16_t min value
#ifndef UINT16_MIN 
#define UINT16_MIN USHRT_MIN
#endif
/// uint16_t max value
#ifndef UINT16_MAX
#define UINT16_MAX USHRT_MAX
#endif
/// int32_t min value
#ifndef INT32_MIN 
#define INT32_MIN LONG_MIN
#endif
/// int32_t max value
#ifndef INT32_MAX
#define INT32_MAX LONG_MAX
#endif
/// uint_t min value
#ifndef UINT32_MIN 
#define UINT32_MIN ULONG_MIN
#endif

/// uint_t max value
#ifndef UINT32_MAX
#define UINT32_MAX ULONG_MAX
#endif

typedef unsigned int      uint_t;
typedef int               int_t;

#include <stddef.h>




/* ----------------------------------------------------------------------- */

// __util_types_h__
#endif
