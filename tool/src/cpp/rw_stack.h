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

/*! \file rw_stack.h
    \brief File for Rogue Wave stack with the STL.
*/

#ifndef __rw_stack_h__
#define __rw_stack_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */


#include <stack>
//using namespace std;

#include "tools_namespace.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class RWStack

  \brief RW stack wrapper.

*/
using std::stack;

template< class T, class C >
struct RwStack : public stack<T,C>
{
  bool isEmpty() const { return this->empty();}
  size_t entries() const { return this->size();}
  void clear()
    {
    int s= this->size();
    for( int i= 0; i < s; i++)
      stack<T,C>::pop();
    }
  T pop()
    {
    T t= this->top();
    stack<T,C>::pop();
    return t;
    }
};

/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#define RWTStack VPTOOLS(RwStack)

/* ----------------------------------------------------------------------- */

#else

#include <rw/tstack.h>

#endif

// __rw_stack_h
#endif
