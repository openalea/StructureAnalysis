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


// ----------------------------------------------
// Class:   MatchPath
// ----------------------------------------------

#ifndef SB_MATCH_PATH_HEADER
#define SB_MATCH_PATH_HEADER

#include "definitions.h"
#include "treegraph.h"
#include "heap.h"
#include "mdtable.h"
#include "inttable.h"
#include <vector>

typedef int FlowType;

typedef std::vector<FlowType> CapacityVector;
typedef std::vector<DistanceType> CostVector;
typedef std::vector<int> IntVector;
typedef std::vector<int> VertexVector;
typedef std::vector<int> EdgeList;

class MatchPath
{
  public :
    //Constructor
    MatchPath(){};
    //Destructor
    ~MatchPath();
    // Cree un graphe de flot a partir de deux listes de noeuds,
    // une contenant les noeuds initiaux et l'autre les noeuds de reference
    void make(NodeList& ,NodeList& );
    void make(int ni, int nj, int *input_list,int *reference_list,int *capacity, double **distances);
    void link(int, DistanceTable&);
    void link(int, DistanceTable&, IntTable&);
    // Recherche le plus court chemin au sens du cout ameliore de
    // Edmons et Karp
    Boolean findPath(VertexVector&,EdgeList& );
    Boolean findPathWithComponents(VertexVector&,EdgeList& );
    // Resout le probleme de flot de cout minimum avec l'algorithme de
    // Busacker et Gowen ameliore par Edmons et Karp
    DistanceType minCostFlow(VertexVector&); //Au lieu de FlowType ...
    DistanceType minCostFlow(int*);
   DistanceType minCostFlowWithComponents(VertexVector&);
    // Functions used to get the edges' cost
    DistanceType length(int,int,int);
    int connected_component(int,int,int);
    Boolean saturated(int);
    Boolean empty(int);
    Boolean reverse(int );
    Boolean direct(int );
    void initCost();
    // Useful function to simulate the flowgraph structure
    int next_edge(int,int);
    int next_vertex(int,int);
    int nbOut(int);
    int who(int );
    int capacity(int);
    DistanceType edgeCost(int , int );
    int edgeComponents(int , int );

  protected :
    DistanceTable*     _treeDistances;  //
    IntTable*          _treeComponents;  //
    NodeList*          _inputList;     //Liste des noeuds de la foret initiale
    NodeList*         _referenceList; //Liste des noeuds de la foret de reference
    CapacityVector     flow; // Vecteur flot
    CostVector         cost;     // Vecteur cout
    IntVector          components;     // Vecteur cout
    IntVector          _capacityList;     // Vecteur cout
    int                _generalCapacity;
    int                nbEdge;          // Nombre d'arcs dans le graphe de flot.
    int                nbVertex;        // Nombre de sommets dans le graphe de flot.
};

#endif



