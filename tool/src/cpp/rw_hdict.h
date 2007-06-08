/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Ch. Pradal (christophe.pradal@cirad.fr)
 *
 *       $Source$
 *       $Id$
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

/*! \file rw_hdict.h
    \brief File for Rogue Wave Hash Dictionnary with the STL.
*/

#ifndef __rw_hdict_h__
#define __rw_hdict_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */

#include <assert.h>
#include "util_hashmap.h"

#include "tools_namespace.h"
#include "rw_hash.h"
#include "rw_defs.h"

/* ----------------------------------------------------------------------- */

TOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class hdict

  \brief hdict (resp. pslist) : Hash Dictionnary.

*/


template < class Key, class T >
struct ValHDict : public STDEXT::hash_map< Key, T, rw_hash<Key> >
{
  typedef STDEXT::hash_map< Key, T, rw_hash<Key> > Base;
  typedef typename Base::key_type key_type;
#ifdef WIN32_STL_EXTENSION
  typedef typename Base::value_type data_type;
#else
  typedef typename Base::data_type data_type;
#endif
  typedef typename Base::const_iterator const_iterator;


  ValHDict( size_t (*h)(const Key&), size_t sz = RWDEFAULT_CAPACITY) :
    STDEXT::hash_map< Key, T, rw_hash<Key> > (
#ifndef WIN32_STL_EXTENSION
		sz,h
#else
		h
#endif
		) {}
  size_t entries() { return this->size(); }
  bool contains( const key_type& k ) const
   { return ( count(k) != 0 ); }
#ifndef WIN32_STL_EXTENSION
  void insertKeyAndValue( const key_type& k, const data_type& a )
    { operator[](k)= a; }
  bool findValue( const key_type& k, data_type& a ) const
    {
    const_iterator it= find(k);
    if( it == this->end() )
      return false;
    assert(count(k));
    a= (*it).second;
    return true;
    }
#else
  void insertKeyAndValue( const key_type& k, const T& a )
    { operator[](k)= a; }
  bool findValue( const key_type& k, T& a ) const
    {
    const_iterator it= find(k);
    if( it == end() )
      return false;
    assert(count(k));
    a= (*it).second;
    return true;
    }
#endif
  bool remove( const key_type& k ) { return erase(k); }
};

template < class Key, class T >
struct ValHDictIt : public STDEXT::hash_map< Key, T, rw_hash<Key> >::iterator
{
  typedef typename STDEXT::hash_map< Key, T, rw_hash<Key> >::iterator It;
  typedef ValHDict<Key,T> Container;

  Container* c;
  bool first;

  ValHDictIt( Container& container) :
    It(container.begin()),
    c(&container), first(true)
    {}

  virtual ~ValHDictIt() { c= 0; first= true; }
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

  Key key() const { return It::operator*().first; }
  T value() const { return It::operator*().second; }

};

/* ----------------------------------------------------------------------- */

TOOLS_END_NAMESPACE

#define RWTValHashDictionary TOOLS(ValHDict)
#define RWTValHashDictionaryIterator TOOLS(ValHDictIt)

/* ----------------------------------------------------------------------- */

#else

#include <rw/tvhdict.h>

#endif

// __rw_hdict_h
#endif
