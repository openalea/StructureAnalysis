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
/*! \file rw_vector.h
    \brief File for Rogue Wave vector with the STL.
*/

#ifndef __rw_vector_h__
#define __rw_vector_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */


#include <vector>
//using namespace std;

#include "tools_namespace.h"
#include "rw_macro.h"
#include "rw_comp.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class RWVector

  \brief RW vector wrapper.

*/


 #define CORPUS_IT_VECTOR(T,OBJ,SUPOBJ) \
 CORPUS_IT(T,OBJ,SUPOBJ)

 #define CORPUS_IT_ORD_VECTOR(T,OBJ,SUPOBJ) \
 CORPUS_IT(T,OBJ,SUPOBJ)


 #define CORPUS_VECTOR(T,OBJ,SUPOBJ) \
  CORPUS(T,OBJ,SUPOBJ) \
  void insert( const_reference a ) { this->push_back(a); } \
  void insertAt( size_t i, const_reference a ) \
    { Master::insert(this->begin()+i,a); } \
  size_t length() const { return this->size();} \
  void append( const_reference a ) { this->push_back(a); } \
  void resize( size_t n ) { this->reserve(n); } \
  const_reference operator()(size_t i) const { return this->operator[](i);} \
  reference operator()(size_t i) { return this->operator[](i);} \
  reference removeLast() { reference last= this->back(); this->pop_back(); return last; }

 #define CORPUS_ORD_VECTOR(T,OBJ,SUPOBJ) \
 CORPUS_VECTOR(T,OBJ,SUPOBJ)

 #define CORPUS_IT_SORTED_VECTOR(T,OBJ,SUPOBJ) \
 CORPUS_IT(T,OBJ,SUPOBJ)

 #define CORPUS_SORTED_VECTOR(T,OBJ,SUPOBJ) \
  CORPUS(T,OBJ,SUPOBJ) \
  void insert( const_reference a ) \
    { Master::insert(lower_bound(this->begin(),this->end(),a),a); } \
  size_t occurrencesOf(const_reference a) const \
    { \
    const_iterator lower= lower_bound(this->begin(),this->end(),a); \
    return distance(lower, upper_bound(lower,this->end(),a)); \
    }

 #define CORPUS_IT_SORTED_VECTOR_PTR(T,OBJ,SUPOBJ) \
 CORPUS_IT(T,OBJ,SUPOBJ)

 #define CORPUS_SORTED_VECTOR_PTR(T,OBJ,SUPOBJ) \
  CORPUS(T,OBJ,SUPOBJ) \
  void insert( const_reference a ) \
    { Master::insert(lower_bound(this->begin(),this->end(),a,lessptr<T>()),a); } \
  size_t occurrencesOf(const_reference a) const \
    { \
    const_iterator lower= lower_bound(this->begin(),this->end(),a,lessptr<T>()); \
    return distance(lower, upper_bound(lower,this->end(),a,lessptr<T>())); \
    }

 #define CTOR_VAL_VECTOR(T,OBJ) \
 CTOR(T,OBJ)\
 CTOR_VAL(T,OBJ)

 #define CTOR_PTR_VECTOR(T,OBJ) CTOR_PTR(T,OBJ)

 #define CTOR_VAL_ORD_VECTOR(T,OBJ) \
 ORD_CTOR(T,OBJ)\
 CTOR_VAL(T,OBJ)

 #define CTOR_PTR_ORD_VECTOR(T,OBJ) CTOR_PTR(T,OBJ)

 #define CTOR_VAL_SORTED_VECTOR(T,OBJ) \
 CTOR_VAL_VECTOR(T,OBJ)

 #define CTOR_PTR_SORTED_VECTOR_PTR(T,OBJ) CTOR_PTR(T,OBJ)

 DEF_RW(T,ValVector,std::vector,VECTOR,VAL)
 DEF_RW(T,PtrVector,std::vector,VECTOR,PTR)

 DEF_RW(T,ValOrdVector,std::vector,ORD_VECTOR,VAL)
 DEF_RW(T,PtrOrdVector,std::vector,ORD_VECTOR,PTR)

 DEF_RW(T,ValSortedVector,std::vector,SORTED_VECTOR,VAL)
 DEF_RW(T,PtrSortedVector,std::vector,SORTED_VECTOR_PTR,PTR)

/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#define RWTValVector VPTOOLS(ValVector)
#define RWTValVectorIterator VPTOOLS(ValVectorIt)::iterator

#define RWTValOrderedVector VPTOOLS(ValOrdVector)
#define RWTValOrderedVectorIterator VPTOOLS(ValOrdVectorIt)::iterator

#define RWTPtrVector VPTOOLS(PtrVector)
#define RWTPtrVectorIterator VPTOOLS(PtrVectorIt)::iterator

#define RWTPtrOrderedVector VPTOOLS(PtrOrdVector)
#define RWTPtrOrderedVectorIterator VPTOOLS(PtrOrdVectorIt)::iterator

#define RWTValSortedVector VPTOOLS(ValSortedVector)
#define RWTValSortedVectorIterator VPTOOLS(ValSortedVectorIt)::iterator

#define RWTPtrSortedVector VPTOOLS(PtrSortedVector)
#define RWTPtrSortedVectorIterator VPTOOLS(PtrSortedVectorIt)::iterator

/* ----------------------------------------------------------------------- */

#else

#include <rw/tvvector.h>
#include <rw/tvsrtvec.h>
#include <rw/tpsrtvec.h>

#endif

// __rw_vector_h
#endif
