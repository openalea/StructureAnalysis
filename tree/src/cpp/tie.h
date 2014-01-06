/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2003 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): C. Pradal (christophe.pradal@cirad.fr)
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

#ifndef TREE_TIE_H
#define TREE_TIE_H

/* ----------------------------------------------------------------------- */

#include <utility>
namespace Tree_tie
{

/**
   \class refpair utility.h
   \brief Like a pair, but contains references which can be assigned to.
*/

template< typename T, typename U >
struct refpair
{
   typedef T first_type;
   typedef U second_type;

   /// Construct a pair of references to \c x and \c y.
   refpair( T& x, U& y ) : first(x), second(y) {}
   /// Construct a copy.
   refpair(refpair const& rp) : first(rp.first), second(rp.second) {}

   /// Assign the values of \c p to the references in this pair.
   refpair& operator=( std::pair<T, U> const& p )
     {
     first = p.first;
     second = p.second;
     return *this;
     }

   /// The first member of the pair.
   T& first;
   /// The second member of the pair.
   U& second;
};

/// Creates a refpair.
template< typename T, typename U >
inline refpair<T, U> tie( T& x, U& y ) { return refpair<T, U>(x, y); }}


#endif
// TREE_TIE_H
