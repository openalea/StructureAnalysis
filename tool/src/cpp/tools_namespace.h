/* -*-c++-*-
 * ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): C. Nouguier & F. Boudon (frederic.boudon@cirad.fr) boudon
 *
 *       $Source$
 *       $Id: tools_namespace.h 9869 2010-11-04 17:31:41Z dbarbeau $
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
 *
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------
 */


#ifndef __vptools_namespace_h__
#define __vptools_namespace_h__

/*! \file tools_namespace.h
    \brief Definition of tools namespace.

  If the macro \b NAMESPACE is defined, \n
  the \b TOOLS Lib will be compiled without \b NAMESPACE, \n
  Else it will be compiled in a namespace named \b TOOLS by default. \n
  To redefined the namespace name, you can redefined \n
  the macro \b TOOLS_NAMESPACE_NAME.
*/

#include "config.h"


#ifdef NO_NAMESPACE

#ifdef VPTOOLS_NAMESPACE
#undef VPTOOLS_NAMESPACE
#endif

#ifdef VPTOOLS_NAMESPACE_NAME
#undef VPTOOLS_NAMESPACE_NAME
#endif

#else

/// Macro that enable the tools namespace
#define VPTOOLS_NAMESPACE

#endif

#ifdef VPTOOLS_NAMESPACE

#ifndef VPTOOLS_NAMESPACE_NAME

/// Macro that contains the tools namespace name
#define VPTOOLS_NAMESPACE_NAME VPTOOLS
#endif

/// Macro for beginning the tools namespace.
#define VPTOOLS_BEGIN_NAMESPACE namespace VPTOOLS_NAMESPACE_NAME {

/// Macro for ending the tools namespace.
#define VPTOOLS_END_NAMESPACE };

/// Macro for using the tools namespace.
#define VPTOOLS_USING_NAMESPACE using namespace VPTOOLS_NAMESPACE_NAME;

/// Macro for using an object of the tools namespace.
#define VPTOOLS_USING(obj) using VPTOOLS_NAMESPACE_NAME::obj;

/// Macro to use an object from the tools namespace.
#define VPTOOLS(obj) VPTOOLS_NAMESPACE_NAME::obj

#else

#ifdef __GNUC__
#warning namespace VPTOOLS not used
#endif

/// Macro for beginning the tools namespace.
#define VPTOOLS_BEGIN_NAMESPACE  

/// Macro for ending the tools namespace.
#define VPTOOLS_END_NAMESPACE  

/// Macro for using the tools namespace.
#define VPTOOLS_USING_NAMESPACE  

/// Macro for using an object of the tools namespace.
#define VPTOOLS_USING(obj)

/// Macro to use an object from the tools namespace.
#define VPTOOLS(obj) obj

#endif


#if defined(_WIN32)

#if defined(VPTOOLS_NODLL)
#undef VPTOOLS_MAKEDLL
#undef VPTOOLS_DLL
#else
#ifndef VPTOOLS_DLL
#define VPTOOLS_DLL
#endif
#endif

#if defined(VPTOOLS_MAKEDLL)
#ifndef VPTOOLS_DLL
#define VPTOOLS_DLL
#endif
#endif

#ifdef VPTOOLS_DLL

#ifdef VPTOOLS_MAKEDLL             /* create a vptool DLL library */
#define VPTOOLS_API  __declspec(dllexport)
#undef VPTOOLS_FWDEF
#else                                                   /* use a vptool DLL library */
#define VPTOOLS_API  __declspec(dllimport)
#endif

#define VPTOOLS_TEMPLATE_API(T) template class VPTOOLS_API T;
#endif

#else // OS != _WIN32

#undef VPTOOLS_MAKEDLL             /* ignore these for other platforms */
#undef VPTOOLS_DLL

#endif

#ifndef VPTOOLS_API
#define VPTOOLS_API
#define VPTOOLS_TEMPLATE_API(T) 
#endif




#endif // __tools_namespace_h__
