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

/*! \file rw_queue.h
    \brief File for Rogue Wave queue with the STL.
*/


#ifndef __rw_queue_h__
#define __rw_queue_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */


#include <queue>
//using namespace std;

#include "tools_namespace.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class RWQueue

  \brief RW queue wrapper.

*/


template< class T, class C >
struct RWQueue : public std::queue<T>
{
#ifndef WIN32_STL_EXTENSION
 typedef typename std::queue<T>::const_reference const_reference;
#endif
  bool isEmpty() const { return this->empty();}
  std::size_t entries() const { return this->size();}
#ifndef WIN32_STL_EXTENSION
  const_reference first () const { return this->front();}
  const_reference last () const { return this->back();}
#else
  const T& first () const { return this->front();}
  const T& last () const { return this->back();}
#endif
  void insert( const T& a ) { this->push(a); }
  T get()
    {
    T a= this->front();
    this->pop();
    return a;
    }

  void clear()
    {
    int s= this->size();
    for( int i= 0; i < s; i++) this->pop();
    assert(this->empty());
    }
};

/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#define RWTQueue VPTOOLS(RWQueue)

/* ----------------------------------------------------------------------- */

#else

#include <rw/tqueue.h>

#endif

// __rw_queue_h
#endif
