/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
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
 CORPUS_IT(T,OBJ,SUPOBJ)

 #define CORPUS_LIST(T,OBJ,SUPOBJ) \
 CORPUS(T,OBJ,SUPOBJ) \
 const_reference first() const { return this->front(); } \
 void removeFirst() { this->pop_front(); } \
 void append( const_reference a ) { this->push_back( a ); } \
 inline void insert( const_reference a ) { this->append(a); } \

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
