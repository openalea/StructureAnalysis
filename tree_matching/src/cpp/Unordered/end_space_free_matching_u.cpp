/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
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


#include "end_space_free_matching_u.h"


EndSpaceFreeMatching_U::EndSpaceFreeMatching_U(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
  Matching_U(input,reference,nodeDistance)
{
  _spaceOptimization = 0;
  _sumNbCaseVector = CaseVector(3,0);
  _nbCaseVector = CaseVector(7,0);
  _distanceMatrix = DistanceVectorTable(T1->getNbVertex());
}


DistanceType EndSpaceFreeMatching_U::match()
{
	const int size1 = T1->getNbVertex();
	const int size2 = T2->getNbVertex();
	DistanceType D=MAXDIST;
	_distanceMatrix.resize(size1);
	//cerr << "" << "Already computed :     ";
	if ((!T1->isNull())&&(!T2->isNull()))
	{
		for (int input_vertex=size1-1;input_vertex>=0;input_vertex--)
		{
			_distances.openDistancesVector(input_vertex);
			if (!_spaceOptimization){
			_insertCost.openDistancesVector(input_vertex);
 			_deleteCost.openDistancesVector(input_vertex);
 			_substitutionCost.openDistancesVector(input_vertex);
			}
			_distanceMatrix[input_vertex].resize(size2);
			for (int reference_vertex=size2-1;reference_vertex>=0;reference_vertex--)
			{
				distanceBetweenForest(input_vertex,reference_vertex);
				DistanceType d =distanceBetweenTree(input_vertex,reference_vertex);
				
				(_distanceMatrix[input_vertex])[reference_vertex]=d;
				if ((input_vertex == 0) || (reference_vertex == 0)){
				  if (d<D){
				    D = d;
				    _i_v = input_vertex;
				    _r_v = reference_vertex; 
				  }
				}
			}
 			if (int(100. - 100*input_vertex/size1)%10 == 0)
 				cerr << "\x0d" << "Already computed : "<<int(100. - 100*input_vertex/size1) <<"% " <<" matched ... " << flush;
			for (int i=1;i<=T1->getNbChild(input_vertex);i++)
			{
				_distances.closeDistancesVector(T1->child(input_vertex,i));
				if (!_spaceOptimization){
					_insertCost.closeDistancesVector(T1->child(input_vertex,i));
					_deleteCost.closeDistancesVector(T1->child(input_vertex,i));
					_substitutionCost.closeDistancesVector(T1->child(input_vertex,i));
				}
			}
		}
	}
  else
    {
      D = 0.;
    }
  //cerr<< "\x0d"<<endl;
  return(D);
}

DistanceType EndSpaceFreeMatching_U::match_begin_T2()
{
	const int size1 = T1->getNbVertex();
	const int size2 = T2->getNbVertex();
	DistanceType D=MAXDIST;
	DistanceType d ;
	_distanceMatrix.resize(size1);
	//cerr << "" << "Already computed :     ";
	if ((!T1->isNull())&&(!T2->isNull()))
	{
		for (int input_vertex=size1-1;input_vertex>=0;input_vertex--)
		{
			_distances.openDistancesVector(input_vertex);
			if (!_spaceOptimization){
			_insertCost.openDistancesVector(input_vertex);
 			_deleteCost.openDistancesVector(input_vertex);
 			_substitutionCost.openDistancesVector(input_vertex);
			}
			_distanceMatrix[input_vertex].resize(size2);
			for (int reference_vertex=size2-1;reference_vertex>=0;reference_vertex--)
			{
				distanceBetweenForest(input_vertex,reference_vertex);
				d =distanceBetweenTree(input_vertex,reference_vertex);
				
				(_distanceMatrix[input_vertex])[reference_vertex]=d;
			}
			if (d<D){
			  D = d;
			  _i_v = input_vertex;
			  _r_v = 0; /* r_v = 0 of course */ 
			}
						
// 			if (int(100. - 100*input_vertex/size1)%10 == 0)
// 				cerr << "\x0d" << "Already computed : "<<int(100. - 100*input_vertex/size1) <<"% " <<" matched ... " << flush<<endl;
			for (int i=1;i<=T1->getNbChild(input_vertex);i++)
			{
				_distances.closeDistancesVector(T1->child(input_vertex,i));
				if (!_spaceOptimization){
					_insertCost.closeDistancesVector(T1->child(input_vertex,i));
					_deleteCost.closeDistancesVector(T1->child(input_vertex,i));
					_substitutionCost.closeDistancesVector(T1->child(input_vertex,i));
				}
			}
		}
	}
  else
    {
      D = 0.;
    }
  //cerr<< "\x0d"<<endl;
  return(D);
}



DistanceType EndSpaceFreeMatching_U::match_begin_end_T2()
{
	const int size1 = T1->getNbVertex();
	const int size2 = T2->getNbVertex();
	DistanceType D=MAXDIST;
	DistanceType d ;
	_distanceMatrix.resize(size1);
	//cerr << "" << "Already computed :     ";
	if ((!T1->isNull())&&(!T2->isNull()))
	{
		for (int input_vertex=size1-1;input_vertex>=0;input_vertex--)
		{
			_distances.openDistancesVector(input_vertex);
			if (!_spaceOptimization){
			_insertCost.openDistancesVector(input_vertex);
 			_deleteCost.openDistancesVector(input_vertex);
 			_substitutionCost.openDistancesVector(input_vertex);
			}
			_distanceMatrix[input_vertex].resize(size2);
			for (int reference_vertex=size2-1;reference_vertex>=0;reference_vertex--)
			{
				distanceBetweenForest(input_vertex,reference_vertex);
				d =distanceBetweenTree(input_vertex,reference_vertex);
				
				(_distanceMatrix[input_vertex])[reference_vertex]=d;
				DistanceType deleteT2 = getDBT(EMPTY_TREE,reference_vertex);
				cerr<<"reference_vertex = "<<reference_vertex<<" - delete T2 = "<<deleteT2<<endl;
				cerr<<d-deleteT2<<endl;
				if (d-deleteT2<D){
					D = d-deleteT2;
					_i_v = input_vertex;
					_r_v = reference_vertex; /* r_v = 0 of course */ 
				}
			}
			
// 			if (int(100. - 100*input_vertex/size1)%10 == 0)
// 				cerr << "\x0d" << "Already computed : "<<int(100. - 100*input_vertex/size1) <<"% " <<" matched ... " << flush<<endl;
			for (int i=1;i<=T1->getNbChild(input_vertex);i++)
			{
				_distances.closeDistancesVector(T1->child(input_vertex,i));
				if (!_spaceOptimization){
					_insertCost.closeDistancesVector(T1->child(input_vertex,i));
					_deleteCost.closeDistancesVector(T1->child(input_vertex,i));
					_substitutionCost.closeDistancesVector(T1->child(input_vertex,i));
				}
			}
		}
		D += getDBT(EMPTY_TREE,0);
		cerr<<"D = "<<D<<endl;
	}
  else
    {
      D = 0.;
    }
  //cerr<< "\x0d"<<endl;
  return(D);
}







