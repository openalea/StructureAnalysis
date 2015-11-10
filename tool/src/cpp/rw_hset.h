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

/*! \file rw_hset.h
    \brief File for Rogue Wave Hash Set with the STL.
*/

#ifndef __rw_hset_h__
#define __rw_hset_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */

#include "tools_namespace.h"

#include "util_hashset.h"
#include "rw_hash.h"
#include "rw_defs.h"




/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class hset

  \brief Hash Set.

*/


template < class T >
struct ValHSet : public tool_hash_set< T, rw_hash<T> >
{
  typedef tool_hash_set< T, rw_hash<T> > Base;
  typedef typename Base::key_type key_type;

  ValHSet( size_t (*h)(const T&), size_t n= RWDEFAULT_CAPACITY ) 
#ifdef WIN32_STL_EXTENSION
 // : Base (n )
#else
  : Base (n, rw_hash<T>(h) )
#endif
    {}

  size_t entries() const { return this->size(); }
  bool contains( const key_type& k ) const
   { return ( this->count(k) != 0 ); }
  bool remove( const key_type& k ) { return this->erase(k); }
};

template < class T >
struct ValHSetIt : public tool_hash_set<T, rw_hash<T> >::iterator
{
  typedef typename tool_hash_set<T, rw_hash<T> >::iterator It;
  typedef ValHSet<T> Container;

  Container* c;
  bool first;

  ValHSetIt( Container& container) :
    It(container.begin()),
    c(&container), first(true)
    {}

  virtual ~ValHSetIt() { c= 0; first= true; }

  bool operator++() { return this->operator()(); }

  bool operator()()
    {
    if ( (*this) == c->end() ) return false;
    if ( first )
      {
      first= false;
      return true;
      }
    else
      return ( It::operator++() != c->end() );
    }

  Container* container() const { return c; }
  void reset() {It::operator=(c->begin()); first= true; }
  T key() const { return It::operator*(); }
};

/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#define RWTValHashSet VPTOOLS(ValHSet)
#define RWTValHashSetIterator VPTOOLS(ValHSetIt)
#define RWTValHashTable VPTOOLS(ValHSet)
#define RWTValHashTableIterator VPTOOLS(ValHSetIt)

/* ----------------------------------------------------------------------- */

#else

#include <rw/tvhset.h>

#endif

// __rw_hset_h
#endif
