/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/Ordered/matching_O.h,v $
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



#ifndef SB_MATCHING_GOSDE_HEADER
#define SB_MATCHING_GOSDE_HEADER

#include "MS_O_Matching.h"


/**
 *\class Matching_O 
 *\brief Algorithm for evaluating ordered distance 
 *\par Presentation
 * In order to compute distance between two tree graphs, we have implemented
 * the algorithm proposed by Zhang \cite{Zha89} to define a distance between 
 * ordered trees
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Pascal Ferraro
 *\date 2005
 */

class Matching_O
{
public:
  Matching_O(TreeGraph&, TreeGraph&, NodeCost&);
  virtual ~Matching_O();
  virtual DistanceType run();
  virtual DistanceType match(){return run();}
  virtual DistanceType getInsertionCost();
  virtual DistanceType getDeletionCost();
  virtual DistanceType getMatchingCost();
  virtual Sequence* getSequence(); 
  virtual void TreeList(int ,int ,Sequence& );
  MS_O_Matching* getMSM(){return _msm;};
  virtual DistanceType getDBT(int, int);
  virtual DistanceType getInsertCost(int j){return 1;}
  virtual DistanceType getDeleteCost(int i ){return 1;}
  virtual DistanceType getSubstitutionCost(int i,int j){return 0;}
protected:
  MS_O_Matching* _msm;
  
};


#endif

