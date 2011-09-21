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
 *       $Id: matchpath.cpp 3258 2007-06-06 13:18:26Z dufourko $
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


// ---------------------------------------------
// La classe GeneralMatchPath permet d'implementer la
// resolution du probleme de flot maximum et de 
// cout minimum.
// ---------------------------------------------

#include "GeneralMatchPath.h"


/* -----------------------------------------------------------------*/

GeneralMatchPath::GeneralMatchPath():
  _inputList(0), _referenceList(0)  {} 

GeneralMatchPath::GeneralMatchPath(const NodeList& input_list,const NodeList& reference_list, const VertexArray inputSuccessors, const VertexArray referencePredecessors):
  _inputList(0), _referenceList(0), _inputSuccessors(0), _referencePredecessors(0){
  make(input_list,reference_list,inputSuccessors,referencePredecessors);
  int deg_max = I_MAX(input_list.size(),reference_list.size());
  // le degre max correspond aux noeuds vides
  // Une optimissation doit etre possible sur le vecteur flow (a revoir)
  flow.resize(nbEdge);
  cost.resize(nbEdge-input_list.size()-reference_list.size()-2);
}


void GeneralMatchPath::link(int deg_max,MatchingDistanceTable* mdtable)
{
  _edgecostEvaluator = MatchEdgeCostPtr(new TableEdgeCost(mdtable));
  //  _mdtable=mdtable;
  //flow.resize(deg_max*deg_max+3*deg_max);
  //cost.resize(2*deg_max+3);
}




// ---------------------------------------------------------
// On initialise le graphe de flot necessaire a l'algorithme
// d'alignement restreint.
// ---------------------------------------------------------
		
void GeneralMatchPath::make(const NodeList& input_list,const NodeList& reference_list, const VertexArray inputSuccessors, const VertexArray referencePredecessors)
{
  if(_inputList) delete _inputList;
  _inputList= new NodeList(input_list);
  if(_referenceList) delete _referenceList;
  _referenceList=new NodeList(reference_list);

  _inputSuccessors = VertexArray(inputSuccessors);
    _referencePredecessors = VertexArray(referencePredecessors);

  // On recupere le nombre d'arbres des forets initiales et finales
  int ni=_inputList->size();
  int nj=_referenceList->size();

  // Le nombre de vertex est egal a la somme des deux ensembles plus deux noeuds 
  // plus deux noeuds vides plus la source et le puis
  nbVertex=ni+nj+4;
  // le nombre d'arcs depend de la liste des successeurs somme des élements de _inputSuccessors
  nbEdge = 0;
  for (int i = 0; i < _inputSuccessors.size(); i++)
    nbEdge += _inputSuccessors[i].size();
  nbEdge += 2*(ni + nj + 1);

  flow = CapacityVector(flow.size(),0);
  cost = CostVector(cost.size(),0.0);
}


	 
		
// -----------
// Destructeur
// -----------
GeneralMatchPath::~GeneralMatchPath()
{
  delete _inputList;
  delete _referenceList;
}

// ----------------------------------------
// FONTIONS USED TO HANDLE THE FLOW GRAPH
// ----------------------------------------

// -------------------------------------------
// Cette fonction vérifie si un arc est saturé
// ------------------------------------------

bool GeneralMatchPath::saturated(int flow_edge)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  if (flow_edge == ni ) 
    return flow[flow_edge] == nj; // correspond au noeud vide connecte aux refs
  if (flow_edge == nbEdge-1)
    return flow[flow_edge] == ni; // correspond au noeud vide connecte aux input
  
 
  return(flow[flow_edge]==1); // tous les autres arcs sont de capacites 1
}


// --------------------------------------------
// On renvoie la capacite de l'arc no flow_edge
// --------------------------------------------

int GeneralMatchPath::capacity(int flow_edge)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  if (flow_edge == ni ) 
    return  nj; // correspond au noeud vide connecte aux refs
  if (flow_edge == nbEdge-1)
    return  ni; // correspond au noeud vide connecte aux input
  
 
  
  return 1;
}


// ---------------------------------------------
// Verifie si le flot de l'arc flow_edge est nul
// ---------------------------------------------
bool GeneralMatchPath::empty(int flow_edge)
{
  assert(flow_edge!=-1);
  return(flow[flow_edge]==MINDIST);
} 

bool GeneralMatchPath::reverse(int residual_edge)
{
  assert(residual_edge!=-1);
  return((residual_edge%2)==1);
}

bool GeneralMatchPath::direct(int residual_edge)
{
  assert(residual_edge!=-1);
  return((residual_edge%2)==0);
}



// -----------------------------------------
// DIJKSTRA'S SHORTEST PATH ALGORITHM
// On implemente l'algorithme de recherche
// de plus court chemin tel que le presente 
// Tarjan, avec les ameliorations de 
// Edmons et Karp
// -----------------------------------------
bool GeneralMatchPath::findPath(VertexVector& VertexOnThePath,EdgeList& EdgeOnThePath)
{
  // On numerote les sommets
  int source=0;
  int sink=nbVertex-1;
  int current_out_vertex;
  int current_out_edge;
  
  // On utilise un tas d'ordres 2+m/n
  Heap path_heap(2.0+((float) nbEdge)/((float) nbVertex));

  // Le sommet courant est la source
  int current_vertex=source;
  
  vector<bool> heap_index(nbVertex,0);
  
  // Chaque sommet est value
  // On initialise la valeur de tous les sommets avec +l'infini
  CostVector distance(nbVertex,MAXDIST);
  
  // sauf la source qui est valuee avec 0
  distance[0]=MINDIST;
  DistanceType t_dist;
  DistanceType epsilon=1e-8;
  


  int item_pos;
  DistanceType tmp_dist;
  bool path_found=0;
  do
    {
      // Pour tous les sommets relie au sommets courant:
      for (int i=1;i<=nbOut(current_vertex);i++)
	{
	  // On regarde les premier cotes et sommets adjacents du sommet courant
	  // We look at the outgoing edges and vertices of current vertex
	  current_out_vertex=next_vertex(current_vertex,i);
	  current_out_edge=next_edge(current_vertex,i);
	  // On applique la transformation d'Edmonds and Karp pour que l'arc est une valeur non negative
	  t_dist=cost[current_vertex]+length(current_out_edge,current_vertex,current_out_vertex)-cost[current_out_vertex]; 
	  if ((t_dist<epsilon)&&(t_dist>-1.0*epsilon)) t_dist=0.0;
	  assert(t_dist>=0);
	  // We evalute the tentative distance
	  // On evalue la nouvelle distance 
	  tmp_dist=distance[current_vertex]+t_dist;
	  
	  // Si cette nouvelle distance est plus petite que la precedente alors
	  if (tmp_dist<distance[current_out_vertex])
	    {
	      // On met a jour la valeur du noeud,
	      distance[current_out_vertex]=tmp_dist;
	      // et on met le sommet et l'arc dans la liste du chemin. 
	      assert(VertexOnThePath.size() > current_out_vertex);
	      VertexOnThePath[current_out_vertex]=current_vertex;
	      if (EdgeOnThePath.size() <= current_out_vertex)
		cout<<"Probleme acces memoire"<<endl;
	      assert(EdgeOnThePath.size() > current_out_vertex);
	      EdgeOnThePath[current_out_vertex]=current_out_edge;
	      
	      // Si de plus, l'index dans le tas du sommet est faux, alors
	      // ce sommet n'a pas encore etait marque,il ne se 
	      // trouve pas dans le tas, donc on le marque et on l'insere
	      // dans le tas avec sa valeur
	      if (!heap_index[current_out_vertex])
		{
		  heap_index[current_out_vertex]=true;
		  int heap_pos=path_heap.insertItem(distance[current_out_vertex],current_out_vertex);
		}
	      // Sinon c'est qu'il se trouvait deja dans le tas
	      // On le recherche, on lui attribue sa nouvelle valeur
	      // et on reorganise le tas.
	      else
		{
		  item_pos=path_heap.position(current_out_vertex);
		  path_heap.at(item_pos)->putKey(distance[current_out_vertex]);
		  path_heap.siftUp(*path_heap.at(item_pos),item_pos);
		}
	    }
	  // Si ce dernier noeud etait le puit, on a trouve un chemin de la source au puit
	  if (current_out_vertex==sink) {path_found=true;};
	}
      // On prend le noeud du tas dont la valeur est minimum et on l'enleve du tas.
      current_vertex=path_heap.deleteMin();
    }
  // On recommence jusqu'a ce qu'il n'y ait plus d'element dans le tas.
  while(current_vertex!=-1);
  
  // Si on a trouve un chemin,
  // on reevalue le cout pour joindre chaque sommet a la source.
  if (path_found) 
    { 
      for (int i=0;i!=nbVertex;i++) 
	{
	  
	  if (distance[i] >= MAXDIST) 
	    {
	      cost[i] = MAXDIST;
	    } 
	  else 
	    { 
	      cost[i]=cost[i]+distance[i];
	    }
	}
    }
	
  return(path_found);
}



// -------------------------------------------------
// On resout ici le probleme du flot de cout minimum
// -------------------------------------------------

//Méthode générale pour le calcul du cout minimum de flot maximum
// DistanceType GeneralMatchPath::minCostFlow(int *int_list)
// {
//   int ni=_inputList->size();
//   int nj=_referenceList->size();  
//   VertexVector map_list;
//   map_list.resize(ni+nj+2);
//   minCostFlow(map_list);
//   for (int i = 0; i<ni+nj+2;i++)
//     int_list[i] = map_list[i];
  
// }


// Le graphe de flot utilise est tel quel:
// une source s, ni sommets representant les arbres de la foret initiale,
// nj sommet representant la foret de reference, un puit et si ni!=nj, on 
// rajoute un sommet du cote de inf(ni,nj) representant l'arbre vide.
// Chaque arc entre la source et le input a une capacite de 1 et un cout nul,
// de meme pour les arcs des ref au puits. Entre les input et les ref, l'arc
// a une capacite de 1 et un cout la distance d'un arbre a l'autre.
// la capacite entre la source (ou le puit) et le noeud vide est de |ni-nj|
// et le cout nul, entre le ref (ou les input) et le vide la capcite est de
// 1 et le cout celui de la transformation.


DistanceType GeneralMatchPath::minCostFlow(VertexVector& map_list)
{
  int current_vertex;
  VertexVector PredOnThePath(nbVertex,-1);
  EdgeList EdgeOfThePath(nbVertex,-1);

  int ni=_inputList->size();
  int nj=_referenceList->size();

  int source=0;
  int sink=nbVertex-1;

  int nb_input=ni;
  //if (ni<nj) { nb_input=ni+1;}
  nb_input=ni+1;
  // La valeur du flot initialement est de 0
  DistanceType flow_value=0;
  // Le flot maximum est le max de ni, nj.
  DistanceType flow_max=D_MAX(ni,nj);

  bool path = true ;

  for (int f=1;(f<=flow_max)&&(path);f++)
    {
      //On cherche le plus court chemin avec les poids de EDMONS AND KARP"<<endl;
      path=findPath(PredOnThePath,EdgeOfThePath);
      current_vertex = sink;
      // Si on a trouve un chemin, on cree le graphe residuel avec le flot augmentant
      // on modifie le flot et les arcs ...
      if (path)
	{
	  do
	    {
	      int residual_edge=EdgeOfThePath[current_vertex];
	      int flow_edge=(int) residual_edge/2;
	      int pred=PredOnThePath[current_vertex];
	      flow_value=flow_value+length(residual_edge,pred,current_vertex);
	      if ((PredOnThePath[current_vertex]<=nb_input)&&(PredOnThePath[current_vertex]!=source))
		{
		  map_list[PredOnThePath[current_vertex]]=current_vertex;
		  map_list[current_vertex]=PredOnThePath[current_vertex];
		} 
	      // Si l'arc considere est un arc de renversement alors, on diminue le flot de 
	      // une unite sur cet arc,
	      if (reverse(residual_edge))
		{
		  flow[flow_edge]=flow[flow_edge]-1;

		  assert(flow[flow_edge]>=0);
		}
	      // sinon on l'augmente de une unite
	      // en verifiant toujours que le flot reste compatible
	      else
		{
		  flow[flow_edge]=flow[flow_edge]+1;

		  assert(flow[flow_edge]<=capacity(flow_edge));
		}
	      current_vertex=PredOnThePath[current_vertex];
	    }
	  while(current_vertex!=source);
	}
    }
  return(flow_value);	
}


// Calcule le cout du passage du sommet vertex1 au sommet vertex2 sur l'arc 
// residual_edge qui relie les deux sommets
DistanceType GeneralMatchPath::length(int residual_edge,int vertex1,int vertex2)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  int flow_edge=(int) residual_edge/2;
  int source=0;
  int sink=nbVertex-1;
  
  
  if ((direct(residual_edge))&&(saturated(flow_edge))) 	{return(2*MAXDIST);};
  if ((reverse(residual_edge))&&(empty(flow_edge))) 	{return(2*MAXDIST);};
  
  if ((vertex1 == source)||(vertex2==source)||(vertex1==sink)||(vertex2==sink))
    {
      return 0;
    }
  else
    {
      int treevertex1 = who(vertex1);
      int treevertex2 = who(vertex2);
      int empty=ni+1;
      if (vertex1<=empty)
	return((capacity(flow_edge)-flow[flow_edge])*edgeCost(treevertex1,treevertex2));
      else
	return(-1*flow[flow_edge]*edgeCost(treevertex2,treevertex1));
    }				
}


int GeneralMatchPath::who(int vertex)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  int source=0;
  int sink=nbVertex-1;
  int EMPTY=-1;
  
  if (vertex==-1) return EMPTY;
  assert((vertex!=source)&&(vertex!=sink));
  int empty_input=ni+1;
  int empty_ref=ni+nj+2;
  if ((vertex == empty_input) || (vertex == empty_ref))
    return EMPTY;
  if (vertex < ni) 
    {
      NodeList::iterator begin;
      begin = _inputList->begin();
      //	  for (int i=0;i<vertex-1;i++)
      for (int i=1;i<vertex;i++)
	begin++;
      return(*begin);
    }
  else 
    {
      NodeList::iterator begin;
      begin = _referenceList->begin();
      //	  for (int i=0;i<vertex-ni-1;i++)
      for (int i=ni+1;i<vertex;i++)
	begin++;
      return(*begin);
    }
}

//--------------------------------------------------------------------------
// FUNCTIONS USED TO SIMULATE THE FLOWGRAPH DATA STRUCTURE
//--------------------------------------------------------------------------

// Cette fonction renvoie le nombre de sommets relie au sommet n

int GeneralMatchPath::nbOut(int n)
{

  int ni=_inputList->size();
  int nj=_referenceList->size();
  ni++;
  nj++;
  if (n==0)  
    return ni; // connection de la source aux autres noeuds
  if (n < ni) 
    return _inputSuccessors[n].size();
  if (n == ni+1)
    return nj;
  if (n < ni + nj +1)
    return _referencePredecessors[n].size();
  else
    return 1;
}

int GeneralMatchPath::next_edge(int n,int i)
{
  // n represente l'index d'un noeud, i est l'index de l'arc sortant
  // la fonction retourne l'index de l'arc
  int ni=_inputList->size();
  int nj=_referenceList->size();
  ni++;
  nj++;
  // source
  if (n==0) 
    return(2*(i-1));
  // sink
  if (n==ni+nj+1) 
    return(2*ni+2*ni*nj+2*(i-1)+1);
  if (n<ni)
    {
      if (i==nbOut(n))
	return 2*(n-1)+1;
      else{
	int successor = _inputSuccessors[n-1][i];
	return 2*ni+2*nj*(n-1)+2*(successor-1);
      }
    }
  if (n == ni) // empty node
    if (i==nbOut(n))
      return 2*(n-1)+1;
    else
 	return 2*ni+2*nj*(n-1)+2*(i-1);
  if (n < ni+nj)
    {
      if (i==nbOut(n))
	return (2*(ni+ni*nj+(n-ni-1)));
      else{
	int predecessor = _referencePredecessors[n-ni-2][i];
	return 2*ni+2*nj*(predecessor-1)+2*(n-ni-1)+1;
      }
    }
  else{ //empty node
      if (i==nbOut(n))
	{
	  return(2*(ni+ni*nj+(n-ni-1)));
	}
      else // Pas sur de ça
	{
	  return(2*ni+2*nj*(i-1)+2*(n-ni-1)+1);	
	}
  }
}

int GeneralMatchPath::next_vertex(int n,int i)
{
		
  int ni=_inputList->size();
  int nj=_referenceList->size();
  ni++;
  nj++;
  
  if (n==0)
    return i;
  if (n==ni+nj+1) 
    return ni+i ;
  if (n<ni)
    {      
      if (i==nbOut(n))        
	return 0;      
      else    
	//	return ni+i;        
	return ni+_inputSuccessors[n-1][i]+1;
    }
  if (n==ni){
    if (i==nbOut(n))        
      return 0;      
    else    
      return ni+i;        
  }  
  if (n < ni+nj)
    {
      if (i==nbOut(n))
	return ni+nj+1;
      else
	return _referencePredecessors[n-1][i]+1;;
    }
  else {
    if (i==nbOut(n))        
      return 0;      
    else    
      return i;        
  }  
}
//--------------------------------------------------------------------------















