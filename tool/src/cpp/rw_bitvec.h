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

#ifndef __rw_bitvec_h__
#define __rw_bitvec_h__

/*! \file rw_bitvec.h
    \brief File for Rogue Wave Bit Vector with the STL.
*/

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */


#include <vector>

#include "tools_namespace.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class RWBitVec

  \brief RW bit vector.

*/

#ifndef STL_EXTENSION
using std::bit_vector;
#else
typedef std::vector<bool> bit_vector;
#endif

struct BitVec : public bit_vector
{
  BitVec( size_t n, bool t) : bit_vector(n,t) {}
  bool testBit( size_t i ) const { return operator[] (i);}
  void setBit( size_t i ) { operator[] (i)= true; }
};


/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#define RWBitVec VPTOOLS(BitVec)

/* ----------------------------------------------------------------------- */

#else

#include <rw/bitvec.h>

#endif

// __rw_stack_h
#endif
