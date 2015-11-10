/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *       $Source$
 *       $Id: matchpath.h 3258 2007-06-06 13:18:26Z dufourko $
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


// ----------------------------------------------
// Class:   MatchPath
// ----------------------------------------------

#ifndef SB_MATCH_PATH_HEADER
#define SB_MATCH_PATH_HEADER

#include "definitions.h"
#include "treegraph.h"
#include "heap.h"
#include "mdtable.h"
#include <vector>

typedef int FlowType;

typedef std::vector<FlowType> CapacityVector;
typedef std::vector<DistanceType> CostVector;
typedef std::vector<int> IntVector;
typedef std::vector<int> VertexVector;
typedef std::vector<int> EdgeList;


/* -----------------------------------------------------------------*/

// An evaluator of edge cost for the MatchPath algorithm
class TREEMATCH_API MatchEdgeCost
{
public :

  /** Default constructor. */
  MatchEdgeCost(){}

  /** Destructor. */
  virtual ~MatchEdgeCost() { }

  // should be redefined
  virtual DistanceType edgeCost(int , int ) = 0;
};

// A smart pointer on a MatchEdgeCost
typedef boost::shared_ptr<MatchEdgeCost> MatchEdgeCostPtr;

/* -----------------------------------------------------------------*/

// EdgeCost evaluator that get cost from a table
class TREEMATCH_API TableEdgeCost : public MatchEdgeCost
{
public :

  /** Default constructor. */
  TableEdgeCost(MatchingDistanceTable* mdtable = NULL);

  /** Destructor. */
  virtual ~TableEdgeCost() { }

  // redefined to lookup in the table
  virtual DistanceType edgeCost(int , int );

protected:
	// The distance table stored as a matrix
    MatchingDistanceTable*     _mdtable;  
};


/* -----------------------------------------------------------------*/


class TREEMATCH_API MatchPath
{
  public :
    //Constructor
    MatchPath();
    MatchPath(const NodeList& ,const NodeList& );
    //Destructor
    virtual ~MatchPath();

    // Cree un graphe de flot a partir de deux listes de noeuds,
    // une contenant les noeuds initiaux et l'autre les noeuds de reference
    void make(const NodeList& ,const NodeList& );


    void link(int, MatchingDistanceTable*);
    // Recherche le plus court chemin au sens du cout ameliore de
    // Edmons et Karp
    bool findPath(VertexVector&,EdgeList& );

    // Resout le probleme de flot de cout minimum avec l'algorithme de
    // Busacker et Gowen ameliore par Edmons et Karp
    DistanceType minCostFlow(VertexVector&); //Au lieu de FlowType ...
    // Functions used to get the edges' cost
    DistanceType length(int,int,int);
    bool saturated(int);
    bool empty(int);
    bool reverse(int );
    bool direct(int );
    void initCost();
    // Useful function to simulate the flowgraph structure
    int next_edge(int,int);
    int next_vertex(int,int);
    int nbOut(int);
    int who(int );
    int capacity(int);

	// call the evaluator. Still can be redefined.
    virtual DistanceType edgeCost(int a, int b)
	{ assert(_edgecostEvaluator);
	  return _edgecostEvaluator->edgeCost(a,b); }
	  
 

  protected :
    MatchEdgeCostPtr   _edgecostEvaluator;
    NodeList*          _inputList;     //Liste des noeuds de la foret initiale
    NodeList*          _referenceList; //Liste des noeuds de la foret de reference
    CapacityVector     flow; // Vecteur flot
    CostVector         cost;     // Vecteur cout
    int                nbEdge;          // Nombre d'arcs dans le graphe de flot.
    int                nbVertex;        // Nombre de sommets dans le graphe de flot.
};

#endif



