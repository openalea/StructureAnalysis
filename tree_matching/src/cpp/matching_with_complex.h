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



#ifndef SB_MATCHING_WITH_COMPLEX_HEADER
#define SB_MATCHING_WITH_COMPLEX_HEADER

#include "definitions.h" 
#include "matchpath.h" 
#include "choicetable.h" 
#include "mdtable.h"
#include "inttable.h"
#include "treegraph.h"
#include "sequence.h" 
#include "nodecost.h"
#include "wnodecost.h"
#include "mnodecost.h"



/**
 *\class Matching With Complex
 *\brief Algorithm for comparing two unordered tree graph using the complex
 *\par Presentation
 * In order to compute a distance between two multiscale tree graph, we have extended 
 * the algorithm proposed by Zhang \cite{Zha93}. The distance is computed as the minimum 
 * cost of sequence of edit operations needed totransform one MTG into another one.
 * In fact we use only two scales of decomposition. We have added two new constraints in 
 * order to respect the MTG structures.
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance), this version take only into account
 * a topological cost.
 *\author Pascal ferraro
 *\date 04/2000
 */

class MatchingWithComplex
{
	
  public :
    /** Default constructor.*/
    MatchingWithComplex() {}

    /** Constructs a MatchingWithComplex using two TreeGraphs defined at the same scale and a NodeCost. */
    MatchingWithComplex(TreeGraph& , TreeGraph& ,NodeCost& ) ;

    /** Constructs a MatchingWithComplex using two TreeGraphs defined at the same scale and a NodeCost. */
    void make(TreeGraph& , TreeGraph& ,NodeCost& ) ;

    /** Destructor. */
    ~MatchingWithComplex();
    
    /** Computes the distances between trees referenced by their roots.
     *  Five different types of distance are defined (see {fer00}) */
    virtual DistanceType distanceBetweenTree(int ,int ) ; 

    /** Computes the distances between forests referenced by their roots. */
    virtual DistanceType distanceBetweenForest(int ,int ) ;
    


    /**\par Operators */

    /** Returns the distance between two trees */
    DistanceType getDBT(int ,int ) const;

    /** Returns the distance between two forests */
    DistanceType getDBF(int ,int ) const;
 
   /** Computes recursively all the distance between subtrees and sub forests
    *  of the original tree graphs  */
    DistanceType  match();

   /** Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between trees which are referenced by their roots.  */
    void getList(int ,int ,Sequence*);

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined an image for the complex of \e v and no image for the complex of \e w.*/
    void ForestList_v_l(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined  the complex of \e v is the image  of the  complex of \e w.*/
    void ForestList_v_w(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined an image for the complex of \e w and no image for the complex of \e v.*/
    void ForestList_l_w(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined no image for the complex of \e v and n the complex of \e w.*/
    void ForestList_l_l(int ,int ,Sequence& );


    int Lat(ChoiceList* L, int vertex);


    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between trees  which are referenced by their roots.*/
    void TreeList(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between trees which are referenced by their roots which
    *  defined an image for the complex of \e v and no image for the complex of \e w.*/
    void TreeList_v_l(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between trees which are referenced by their roots in which
    *   the complex of \e v and the complex of \e w are image on each other .*/
    void TreeList_v_w(int ,int ,Sequence& );
 
    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined an image for the complex of \e w and no image for the complex of \e v.*/
   void TreeList_l_w(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots in which
    *  the complex and \e v and \e w are no image*/
    void TreeList_l_l(int ,int ,Sequence& );

    int operator()(int );

  protected :
    TreeGraph* T1;
    TreeGraph* T2;
    MatchingDistanceTable _distances;  
    MatchingDistanceTable _d_v_w;
    MatchingDistanceTable _d_v_l;
    MatchingDistanceTable _d_l_w;
    MatchingDistanceTable _d_l_l;
    MatchingDistanceTable _restrDistances_l_l;
    MatchingDistanceTable _restrDistances_v_w;
    MatchingDistanceTable _restrDistances;
    ChoiceTable _choices;
    ChoiceTable _choices_v_w;
    ChoiceTable _choices_v_l;
    ChoiceTable _choices_l_w;
    ChoiceTable _choices_l_l;
    NodeCost* ND;
    MatchPath _restrMapp_v_w;
    MatchPath _restrMapp_l_l;
    VertexVector _restrMappList_v_w;
    VertexVector _restrMappList_l_l;
    int M(int,int);
    int M_v_l(int,int);
    int M_l_w(int,int);
    int M_v_w(int,int);
    int M_l_l(int,int);
};


#endif

