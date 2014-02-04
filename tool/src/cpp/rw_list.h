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

/*! \file rw_list.h
    \brief File for Rogue Wave List with the STL.
*/

#ifndef __rw_list_h__
#define __rw_list_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */


#include <list>
//using namespace std;

#include "tools_namespace.h"
#include "rw_macro.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class vlist

  \brief vlist (resp. plist) : Double linked list.

*/

#define CORPUS_IT_LIST(T,OBJ,SUPOBJ) \
 CORPUS_IT(T,OBJ,SUPOBJ) \
 void remove() \
   { \
   c->erase(*this); \
   c->previous(*this); \
   }


 #define CORPUS_LIST(T,OBJ,SUPOBJ) \
 CORPUS(T,OBJ,SUPOBJ) \
 const_reference first() const { return this->front(); } \
 const_reference last() const { const_iterator i = this->end(); i--; return *i; } \
 reference removeFirst() { reference elt= this->front(); this->pop_front(); return elt; } \
 void append( const_reference a ) { this->push_back(a); } \
 inline void insert( const_reference a ) { this->append(a); } \
 iterator insert_after(iterator pos, const_reference x){ return Master::insert(pos,x); } \
 template<class InputIterator> \
 void insert_after(iterator pos, InputIterator f, InputIterator l){ Master::insert(pos,f,l); } \
 iterator previous(iterator pos){ pos--; return pos; } \
 void insertAt( size_t i, const_reference a ) \
        { if(i==0) { prepend(a); return; } \
		iterator it= this->begin(); std::advance(it,i-1); insert_after(it,a);  } \
 void removeIn(size_t b, size_t e) \
   { \
   iterator it1= this->begin(); std::advance(it1,b); \
   iterator it2= this->begin(); std::advance(it2,e+1); \
   this->erase(it1,it2); \
   } \
 void prepend( const_reference a ) { this->push_front(a); } \
 const_reference operator[](size_t n) const \
   { const_iterator it= this->begin(); std::advance(it,n); return *it; } \
 reference operator[](size_t n) \
   { iterator it= this->begin(); std::advance(it,n); return *it; }

 #define CTOR_VAL_LIST(T,OBJ) \
 CTOR(T,OBJ)\
 CTOR_VAL(T,OBJ)

 #define CTOR_PTR_LIST(T,OBJ) CTOR_PTR(T,OBJ)

 DEF_RW(T,ValList,std::list,LIST,VAL)
 DEF_RW(T,PtrList,std::list,LIST,PTR)

/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#define RWTValDlist VPTOOLS(ValList)
#define RWTValDlistIterator VPTOOLS(ValListIt)
#define RWTPtrDlist VPTOOLS(PtrList)
#define RWTPtrDlistIterator VPTOOLS(PtrListIt)

/* ----------------------------------------------------------------------- */

#else

#include <rw/tvdlist.h>
#include <rw/tpdlist.h>

#endif

// __rw_list_h
#endif
