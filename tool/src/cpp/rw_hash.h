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
/*! \file rw_hash.h
    \brief File for Rogue Wave Hash with the STL.
*/

#ifndef __rw_hash_h__
#define __rw_hash_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */

#include "tools_namespace.h"
#include "util_hashmap.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class hash

  \brief Hash struct.

*/

template< class T >
class rw_hash
#ifdef WIN32_STL_EXTENSION
	: public STDEXT::hash_compare<T>
#endif
{
private:
  std::size_t (*ptr)(const T&);
public:
  rw_hash() {}
  rw_hash(std::size_t (*f)(const T&)) : ptr(f) {}
  std::size_t operator()(const T& x) const { return ptr(x); }
#ifdef WIN32_STL_EXTENSION
  bool operator()(const T& x,const T& y) const { return x < y; }
#endif
};

VPTOOLS_END_NAMESPACE

#endif
// __rw_hash_h
#endif
