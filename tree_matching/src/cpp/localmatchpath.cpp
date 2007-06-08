/* -*-c++-*- 
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture 
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr) 
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
// La classe MatchPath permet d'implementer la
// resolution du probleme de flot maximum et de 
// cout minimum.
// ---------------------------------------------

#include "localmatchpath.h"


// -------------------------------------------------
// Renvoie le cout du passage d'un sommet initial au
// sommet de reference et le place dans le tableau
// des distances
// -------------------------------------------------

DistanceType LocalMatchPath::edgeCost(int input_vertex,int reference_vertex)
{
  if (input_vertex==-1) 
    //input_vertex = _treeDistances->getSimulatedSize()-1; 
    return(_max_distance);
  if (reference_vertex==-1) 
    //reference_vertex = _treeDistances->getColumnSize()-1;
    return(_max_distance);
  return(_max_distance-_treeDistances->getDistance(input_vertex,reference_vertex));
}

// -----------------------------------------
// DIJKSTRA'S SHORTEST PATH ALGORITHM
// On implemente l'algorithme de recherche
// de plus court chemin tel que le presente 
// Tarjan, avec les ameliorations de 
// Edmons et Karp
// -----------------------------------------
Boolean LocalMatchPath::findPath(VertexVector& VertexOnThePath,EdgeList& EdgeOnThePath)
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
  
  vector<Boolean> heap_index(nbVertex,0);
  
  // Chaque sommet est value
  // On initialise la valeur de tous les sommets avec +l'infini
  CostVector distance(nbVertex,MAXDIST);
  
 // sauf la source qui est valuee avec 0
  distance[0]=MINDIST;
  DistanceType t_dist;
  DistanceType epsilon=1e-8;
  


  int item_pos;
  DistanceType tmp_dist;
  Boolean path_found=0;
  
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
	      VertexOnThePath[current_out_vertex]=current_vertex;
	      EdgeOnThePath[current_out_vertex]=current_out_edge;
	      
	      // Si de plus, l'index dans le tas du sommet est faux, alors
	      // ce sommet n'a pas encore jamais etait marque,il ne se 
	      // trouve pas dans le tas, donc on le marque et on l'insere
	      // dans le tas avec sa valeur
	      if (!heap_index[current_out_vertex])
		{
		  heap_index[current_out_vertex]=1;
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
	  if (current_out_vertex==sink) {path_found=1;};
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
// On resout ici le probleme du flot de cout maximum
// -------------------------------------------------


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


DistanceType LocalMatchPath::maxCostFlow(VertexVector& map_list)
{
  int current_vertex;
  VertexVector PredOnThePath(nbVertex,-1);
  EdgeList EdgeOfThePath(nbVertex,-1);
  int ni=_inputList->size();
  int nj=_referenceList->size();
  map_list.resize(ni+nj+3,-1);	

	
  int source=0;
  int sink=nbVertex-1;

  int nb_input=ni;
  if (ni<nj) 
    nb_input=ni+1;
// La valeur du flot initialement est de 0
  DistanceType flow_value=0;
// Le flot maximum est le max de ni, nj.
  DistanceType flow_max=D_MAX(ni,nj);
  
  _max_distance = 0.;
  NodeList::iterator begin_input,begin_reference;
  begin_input = _inputList->begin();
  
  for (int i = 0;i<ni;i++)
   {
     begin_reference = _referenceList->begin();
     for ( int j = 0;j<nj;j++)
       {
	 if (_max_distance<_treeDistances->getDistance(*begin_input,*begin_reference))
	   _max_distance = _treeDistances->getDistance(*begin_input,*begin_reference);
	 begin_reference++;
       }
     begin_input++;
   }
  //_max_distance =_max_distance+1.;

  Boolean path=1;
  for (int f=1;(f<=flow_max)&&(path);f++)
  {
// On cherche le plus court chemin avec les poids de EDMONS AND KARP
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
	flow_value=flow_value+reverseLength(residual_edge,pred,current_vertex);
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
DistanceType LocalMatchPath::length(int residual_edge,int vertex1,int vertex2)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  int flow_edge=(int) residual_edge/2;
  int source=0;
  int sink=nbVertex-1;
  
	
  if ((direct(residual_edge))&&(saturated(flow_edge))) 	{return(2*MAXDIST);};
  if ((reverse(residual_edge))&&(empty(flow_edge))) 	{return(2*MAXDIST);};
  
  if ((vertex1==source)||(vertex2==source)||(vertex1==sink)||(vertex2==sink))
    {
      return(0);
    }
  else
    {
      int treevertex1=who(vertex1);
      int treevertex2=who(vertex2);
      int empty=ni+1;
      if (ni<nj)
	{
	  if (vertex1<=empty)
	    {
	      return((capacity(flow_edge)-flow[flow_edge])*edgeCost(treevertex1,treevertex2));
	    }
	  else
	    {
	      return(-1*flow[flow_edge]*edgeCost(treevertex2,treevertex1));
	    }
	}
      else
	{
	  if (vertex1<empty)
	    {
	      return((capacity(flow_edge)-flow[flow_edge])*edgeCost(treevertex1,treevertex2));
	    }
	  else 
	    {
	      return(-1*flow[flow_edge]*edgeCost(treevertex2,treevertex1));
	    }
	}
    }				
}

// Calcule le cout du passage du sommet vertex1 au sommet vertex2 sur l'arc 
// residual_edge qui relie les deux sommets
DistanceType LocalMatchPath::reverseLength(int residual_edge,int vertex1,int vertex2)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  int flow_edge=(int) residual_edge/2;
  int source=0;
  int sink=nbVertex-1;
  
	
  if ((direct(residual_edge))&&(saturated(flow_edge))) 	{return(2*MAXDIST);};
  if ((reverse(residual_edge))&&(empty(flow_edge))) 	{return(2*MAXDIST);};
  
  if ((vertex1==source)||(vertex2==source)||(vertex1==sink)||(vertex2==sink))
    {
      return(0);
    }
  else
    {
      int treevertex1=who(vertex1);
      int treevertex2=who(vertex2);
      int empty=ni+1;
      if (ni<nj)
	{
	  if (vertex1<=empty)
	    {
	      return((capacity(flow_edge)-flow[flow_edge])*(_max_distance-edgeCost(treevertex1,treevertex2)));
	    }
	  else
	    {
	      return(-1*flow[flow_edge]*(_max_distance-edgeCost(treevertex2,treevertex1)));
	    }
	}
      else
	{
	  if (vertex1<empty)
	    {
	      return((capacity(flow_edge)-flow[flow_edge])*(_max_distance-edgeCost(treevertex1,treevertex2)));
	    }
	  else 
	    {
	      return(-1*flow[flow_edge]*(_max_distance-edgeCost(treevertex2,treevertex1)));
	    }
	}
    }				
}
