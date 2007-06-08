/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
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


#ifndef SB_SELKOW_HEADER
#define SB_SELKOW_HEADER
#include "matching.h"
class Selkow : public Matching
{
public:

//     //Constructor
    Selkow() {};
    Selkow(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    void make(TreeGraph& , TreeGraph& ,NodeCost& ) ;
    //Destructor
    ~Selkow();

    DistanceType match();

    //Distance Between Trees
    virtual DistanceType distanceBetweenTree(int ,int ) ;
    //Distances Betweeen Forest
    virtual DistanceType distanceBetweenForest(int ,int ) ;

    void ForestList(int ,int,Sequence& );
  void TreeList(int ,int ,Sequence& );

    DistanceType getDBF(int ,int ) ;

    DistanceType getInBF(int ,int ) ;

    DistanceType getDeBF(int ,int ) ;

    DistanceType getSuBF(int ,int ) ;

  
};


#endif

