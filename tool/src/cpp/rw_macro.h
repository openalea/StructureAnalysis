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

/*! \file rw_macro.h
    \brief File for macro to define Rogue Wave object with the STL.
*/

#ifndef __rw_macro_h__
#define __rw_macro_h__

#include "config.h"
#include <assert.h>
#include "rw_defs.h"
#include "rw_comp.h"
#include <algorithm>

#ifdef RWOUT

/* ----------------------------------------------------------------------- */


#define BEGIN_STRUCT_VAL_OBJ(T,OBJ,SUPOBJ) \
template < class T > \
struct OBJ : public SUPOBJ<T> \
{ \
  typedef SUPOBJ<T> Master; \
  typedef typename Master::const_reference const_reference; \
  typedef typename Master::reference reference; \
  typedef typename Master::const_iterator const_iterator; \
  typedef typename Master::iterator iterator;

#define BEGIN_STRUCT_PTR_OBJ(T,OBJ,SUPOBJ) \
template < class T > \
struct OBJ : public SUPOBJ<T*> \
{ \
  typedef SUPOBJ<T*> Master; \
  typedef typename Master::const_reference const_reference; \
  typedef typename Master::reference reference; \
  typedef typename Master::const_iterator const_iterator; \
  typedef typename Master::iterator iterator;

#define CTOR(T,OBJ) \
  OBJ( size_t n ) : Master(n,T()) {}
#define ORD_CTOR(T,OBJ) \
  OBJ( size_t n ) {}

#define CTOR_VAL(T,OBJ) \
  OBJ() : Master() {} \
  OBJ( size_t n, const_reference t ) : Master(n,t) {} \
  ~OBJ() {} \
  const T* data() const { return & (this->front()); }	\
  size_t index( const_reference a) const \
    { \
    const_iterator it= std::find(this->begin(),this->end(),a); \
    return ( it == this->end() ) ? RW_NPOS : std::distance(this->begin(),it); \
    } \
  bool contains( const_reference a ) const \
    { return !( std::find(this->begin(),this->end(),a) == this->end() ); }

#define CTOR_PTR(T,OBJ) \
  OBJ() {} \
  OBJ( size_t n ) {} \
  OBJ( size_t n, const_reference t ) {} \
  ~OBJ() {} \
  T* const* data() const { return &(this->front()); }	\
  void clearAndDestroy() \
    { \
    for(iterator it= this->begin(); it != this->end(); it++) delete(*it); \
    this->clear(); \
    } \
  size_t index( const_reference a) const \
    { \
    const_iterator it= std::find_if(this->begin(),this->end(),eqptr<T>(a)); \
    return ( it == this->end() ) ? RW_NPOS : distance(this->begin(),it); \
    } \
  bool contains( const_reference a ) const \
    { return !( std::find_if(this->begin(),this->end(),eqptr<T>(a)) == this->end() ); }

#define BEGIN_STRUCT_VAL_IT(T,OBJ,SUPOBJ) \
template < class T > \
struct OBJ##It : public SUPOBJ<T>::iterator \
{ \
  typedef typename SUPOBJ<T>::iterator It; \
  typedef OBJ<T> Container; \
  Container* c; \
  typedef typename Container::const_reference const_reference; \
  bool first;

#define BEGIN_STRUCT_PTR_IT(T,OBJ,SUPOBJ) \
template < class T > \
struct OBJ##It : public SUPOBJ<T*>::iterator \
{ \
  typedef typename SUPOBJ<T*>::iterator It; \
  typedef OBJ<T> Container; \
  typedef typename Container::const_reference const_reference; \
  Container* c; \
  bool first;

#define CORPUS(T,OBJ,SUPOBJ) \
  size_t entries() const { return this->size(); } \
  const_reference at(size_t i) const \
    { \
    assert(i<this->size()); \
    const_iterator it= this->begin(); \
    std::advance(it,i); \
    return *it; \
    } \
  reference at(size_t i)\
    { \
    assert(i<this->size()); \
    iterator it= this->begin(); \
    std::advance(it,i); \
    return *it; \
    } \
  void reshape( size_t n ) { Master::resize(n); assert(n==this->size()); } \
  bool isEmpty() const {return this->empty(); }

#define CORPUS_IT(T,OBJ,SUPOBJ) \
  OBJ##It( Container& container) : \
    It(container.begin()), \
    c(&container), first(true) \
    {} \
  virtual ~OBJ##It() { c= 0; first= true; } \
  bool operator()() \
    { \
    if ( (*this) == c->end() ) return false; \
    if( first )  { first= false; return true; } \
    else return ( It::operator++() != c->end() ); \
    } \
  bool operator++() { return this->operator()(); } \
  Container* container() const { return c; } \
  void reset() {It::operator=(c->begin());first= true;} \
  void reset(Container& container) \
    { \
    c= &container; \
    It::operator=(c->begin()); \
    first= true; \
    } \
  const_reference key() const { return It::operator*(); }

#define END_STRUCT };

// TYPE= VAL or PTR
#define DEF_RW(T,OBJ,SUPOBJ,NAME,TYPE) \
 BEGIN_STRUCT_##TYPE##_OBJ(T,OBJ,SUPOBJ) \
 CTOR_##TYPE##_##NAME(T,OBJ) \
 CORPUS_##NAME(T,OBJ,SUPOBJ) \
 END_STRUCT \
 BEGIN_STRUCT_##TYPE##_IT(T,OBJ,SUPOBJ) \
 CORPUS_IT_##NAME(T,OBJ,SUPOBJ) \
 END_STRUCT


/* ----------------------------------------------------------------------- */

#endif
// __rw_macro_h
#endif
