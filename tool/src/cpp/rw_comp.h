/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Ch. Pradal (christophe.pradal@cirad.fr)
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

/*! \file rw_comp.h
    \brief File for Rogue Wave object comparison with the STL.
*/

#ifndef __rw_comp_h__
#define __rw_comp_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */

#include "tools_namespace.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class equality

  \brief Predicate.

*/

template< class T >
struct eqptr {
 private:
  const T* _a;

 public:
  /// Compare the 2 string.
  eqptr( const T* a ) : _a(a) { }
  bool operator()(const T* s) const {
    return *s == *_a;
  }
};

template< class T >
struct lessptr {
 public:
  /// Compare the 2 Objects.
  bool operator()(const T* s, const T* v) const {
    return *s < *v;
  }
};

VPTOOLS_END_NAMESPACE

#endif
// __rw_comp_h
#endif
