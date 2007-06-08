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


// ---------------------------------------------
// La classe MatchPath permet d'implementer la
// resolution du probleme de flot maximum et de 
// cout minimum.
// ---------------------------------------------

#include "matchpath.h"

void MatchPath::link(int deg_max,DistanceTable& tree_distances)
{
  _treeDistances=&tree_distances;
  flow.resize(deg_max*deg_max+3*deg_max);
  cost.resize(2*deg_max+3);
}



void MatchPath::link(int deg_max,DistanceTable& tree_distances,IntTable& tree_components)
{
  _treeDistances=&tree_distances;
  _treeComponents=&tree_components;
  flow.resize(deg_max*deg_max+3*deg_max);
  cost.resize(2*deg_max+3);
  components.resize(2*deg_max+3);

}

  // ---------------------------------------------------------
  // On initialise le graphe de flot necessaire a l'algorithme
  // d'alignement restreint.
  // ---------------------------------------------------------
		
void MatchPath::make(NodeList& input_list,NodeList& reference_list)
{
  
  _inputList=&input_list;
  _referenceList=&reference_list;
  _generalCapacity =0;
  // On recupere le nombre d'arbres des forets initiales et finales
  int ni=_inputList->size();
  int nj=_referenceList->size();
  
  
  if (ni!=nj) 
    {
  // Les sommets du graphe de flot sont ni + nj + 3:
  // Une source et un puit, ni sommets representant les arbres initiaux,
  // nj sommets representant les arbres finaux plus un noeud representant
  // l'arbre vide.
      nbVertex=ni+nj+3;
      if (ni<nj) 
	{
  // Si ni<nj, le noeud representant l'arbre vide est du cote des noeuds 
  // initiaux donc le nombre d'arc est:
  //    ni entre la source et les init,
  // +  ni*nj entre les init et les ref,
  // +  nj entre les ref et le puits,
  // +  1 entre la source et le vide,
  // +  nj entre le vide et les ref. 
  // d'ou nbEdge = ni+nj*ni+nj+nj+1 = ni + ni*nj + 2*nj +1 !!!!!

	  nbEdge=ni+ni*nj+2*nj;
	}
      else
	{
	  nbEdge=2*ni+ni*nj+nj;
	}
      
    }
  else
    {
      nbVertex=ni+nj+2;
      nbEdge=(ni+ni*nj+nj);
    }

  // On initialise le flot et le cout a 0

  flow = CapacityVector(flow.size(),0);
  cost = CostVector(cost.size(),0.0);
  components= IntVector(components.size(),(int) 0);		
}


  // ---------------------------------------------------------
  // On initialise le graphe de flot necessaire a l'algorithme
  // d'alignement restreint Spécial pour Fred
  // ---------------------------------------------------------
		
void MatchPath::make(int ni, int nj, int *input_list,int *reference_list,int *capacity, double **distances)
{
  int s1,s2;

  _treeDistances = new DistanceTable(2*ni,nj,ni);

    //_treeDistances->resize(2*ni,nj,ni);

  _generalCapacity = 1;
  _inputList=new NodeList();
  _referenceList=new NodeList();
  _capacityList.resize(nj);

  for (s1=0;s1<ni;s1++) 
    {
      _inputList->push_back(input_list[s1]); 
      _treeDistances->openDistanceVector(s1);
    }
  for (s2=0;s2<nj;s2++) 
    {
      _referenceList->push_back(reference_list[s2]); 
      _capacityList[s2] = capacity[s2];
   }
  
  int deg_max = I_MAX (ni,nj);
  flow.resize(deg_max*deg_max+3*deg_max);
  cost.resize(2*deg_max+3);

  for (s1=0;s1<ni;s1++) 
    {
      for (s2=0;s2<nj;s2++)
	_treeDistances->putDistance(distances[s1][s2],s1,s2);
    }
 
  nbEdge = ni + ni*nj + nj;
  nbVertex = ni+nj+2;


  // On initialise le flot et le cout a 0

  flow = CapacityVector(flow.size(),0);
  cost = CostVector(cost.size(),0.0);
  components= IntVector(components.size(),(int) 0);		
}
	 
		
// -----------
// Destructeur
// -----------
MatchPath::~MatchPath()
{
  if (_generalCapacity==1)
    {
      int ni=_inputList->size();
      for (int i =0; i<ni; i++)
      _treeDistances->closeDistanceVector(i);
    }
     
}

// ----------------------------------------
// FONTIONS USED TO HANDLE THE FLOW GRAPH
// ----------------------------------------

// -------------------------------------------
// Cette fonction vérifie si un arc est saturé
// ------------------------------------------

Boolean MatchPath::saturated(int flow_edge)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  if (_generalCapacity==0)
    {       
      if ((ni<nj)&&(flow_edge==ni)) { return(flow[flow_edge]==nj-ni); } 
      
      if ((ni>nj)&&(flow_edge==(2*ni+ni*nj))) { return(flow[flow_edge]==ni-nj); } 
      
      return(flow[flow_edge]==1);
    }
  else
    {
      if (flow_edge>=ni+ni*nj)
	return(flow[flow_edge]==_capacityList[flow_edge-ni+ni*nj]);
      else
	return(flow[flow_edge]==1);
    }
	
}


// --------------------------------------------
// On renvoie la capacite de l'arc no flow_edge
// --------------------------------------------

int MatchPath::capacity(int flow_edge)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  if (_generalCapacity==0)
    {
      if ((ni<nj)&&(flow_edge==ni)) { return(nj-ni); } 
      
      if ((ni>nj)&&(flow_edge==(2*ni+ni*nj))) { return(ni-nj); } 
      
      return(1);
    }
  else
    {
      if (flow_edge>=ni+ni*nj)
	return(_capacityList[flow_edge-ni+ni*nj]);
      else
	return(1);
    }
}


// ---------------------------------------------
// Verifie si le flot de l'arc flow_edge est nul
// ---------------------------------------------
Boolean MatchPath::empty(int flow_edge)
{
  assert(flow_edge!=-1);
  return(flow[flow_edge]==MINDIST);
} 

Boolean MatchPath::reverse(int residual_edge)
{
  assert(residual_edge!=-1);
  return((residual_edge%2)==1);
}

Boolean MatchPath::direct(int residual_edge)
{
  assert(residual_edge!=-1);
  return((residual_edge%2)==0);
}

// -------------------------------------------------
// Renvoie le cout du passage d'un sommet initial au
// sommet de reference et le place dans le tableau
// des distances
// -------------------------------------------------

DistanceType MatchPath::edgeCost(int input_vertex,int reference_vertex)
{
  if (input_vertex==-1) { input_vertex = _treeDistances->getSimulatedSize()-1; }
  if (reference_vertex==-1) { reference_vertex = _treeDistances->getColumnSize()-1; }
  return(_treeDistances->getDistance(input_vertex,reference_vertex));
}

int MatchPath::edgeComponents(int input_vertex,int reference_vertex)
{
  if (input_vertex==-1) { input_vertex = _treeComponents->getSimulatedSize()-1; }
  if (reference_vertex==-1) { reference_vertex = _treeComponents->getColumnSize()-1; }
  return(_treeComponents->getInt(input_vertex,reference_vertex));
}

// -----------------------------------------
// DIJKSTRA'S SHORTEST PATH ALGORITHM
// On implemente l'algorithme de recherche
// de plus court chemin tel que le presente 
// Tarjan, avec les ameliorations de 
// Edmons et Karp
// -----------------------------------------
Boolean MatchPath::findPath(VertexVector& VertexOnThePath,EdgeList& EdgeOnThePath)
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

// -----------------------------------------
// DIJKSTRA'S SHORTEST PATH ALGORITHM
// On implemente l'algorithme de recherche
// de plus court chemin tel que le presente 
// Tarjan, avec les ameliorations de 
// Edmons et Karp et avec les composantes connexes
// -----------------------------------------
Boolean MatchPath::findPathWithComponents(VertexVector& VertexOnThePath,EdgeList& EdgeOnThePath)
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
  IntVector compo(nbVertex,MAXINT);
  
  // sauf la source qui est valuee avec 0
  distance[0]=MINDIST;
  compo[0]=0;

  DistanceType t_dist;
  int t_comp;

  DistanceType epsilon=1e-8;
  


  int item_pos;
  DistanceType tmp_dist;
  int tmp_comp;
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
	  t_comp=components[current_vertex]+connected_component(current_out_edge,current_vertex,current_out_vertex)-components[current_out_vertex]; 
	  if ((t_dist<epsilon)&&(t_dist>-1.0*epsilon)) t_dist=0.0;
	  assert(t_dist>=0);
	  // We evalute the tentative distance
	  // On evalue la nouvelle distance 
	  tmp_dist=distance[current_vertex]+t_dist;
	  tmp_comp=compo[current_vertex]+t_comp;
	  
	  // Si cette nouvelle distance est plus petite que la precedente alors
	  if ((tmp_dist<distance[current_out_vertex])||((tmp_dist==distance[current_out_vertex])&&(tmp_comp<compo[current_out_vertex])))
	    {
	      // On met a jour la valeur du noeud,
	      distance[current_out_vertex]=tmp_dist;
	      compo[current_out_vertex]=tmp_comp;
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
	      components[i]= MAXINT;
	    } 
	  else 
	    { 
	      cost[i]=cost[i]+distance[i];
	      components[i]=components[i]+compo[i];
	    }
	}
    }
	
  return(path_found);
}



// -------------------------------------------------
// On resout ici le probleme du flot de cout minimum
// -------------------------------------------------

//Méthode générale pour le calcul du cout minimum de flot maximum
// DistanceType MatchPath::minCostFlow(int *int_list)
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


DistanceType MatchPath::minCostFlow(VertexVector& map_list)
{
  int current_vertex;
  VertexVector PredOnThePath(nbVertex,-1);
  EdgeList EdgeOfThePath(nbVertex,-1);

  int ni=_inputList->size();
  int nj=_referenceList->size();
	
  int source=0;
  int sink=nbVertex-1;

  int nb_input=ni;
  if (ni<nj) { nb_input=ni+1;};
// La valeur du flot initialement est de 0
  DistanceType flow_value=0;
// Le flot maximum est le max de ni, nj.
  DistanceType flow_max=D_MAX(ni,nj);

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

// -------------------------------------------------
// On resout ici le probleme du flot de cout minimum
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


DistanceType MatchPath::minCostFlowWithComponents(VertexVector& map_list)
{
  int current_vertex;
  VertexVector PredOnThePath(nbVertex,-1);
  EdgeList EdgeOfThePath(nbVertex,-1);

  int ni=_inputList->size();
  int nj=_referenceList->size();
	
  int source=0;
  int sink=nbVertex-1;

  int nb_input=ni;
  if (ni<nj) { nb_input=ni+1;};
// La valeur du flot initialement est de 0
  DistanceType flow_value=0;
// Le flot maximum est le max de ni, nj.
  DistanceType flow_max=D_MAX(ni,nj);

  Boolean path=1;
   	
  for (int f=1;(f<=flow_max)&&(path);f++)
  {
// On cherche le plus court chemin avec les poids de EDMONS AND KARP
    path=findPathWithComponents(PredOnThePath,EdgeOfThePath);
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
DistanceType MatchPath::length(int residual_edge,int vertex1,int vertex2)
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

int MatchPath::connected_component(int residual_edge,int vertex1,int vertex2)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  int flow_edge=(int) residual_edge/2;
  int source=0;
  int sink=nbVertex-1;
  
	
  if ((direct(residual_edge))&&(saturated(flow_edge))) 	{return(MAXINT);};
  if ((reverse(residual_edge))&&(empty(flow_edge))) 	{return(MAXINT);};
  
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
	      return((capacity(flow_edge)-flow[flow_edge])*edgeComponents(treevertex1,treevertex2));
	    }
	  else
	    {
	      return(-1*flow[flow_edge]*edgeComponents(treevertex2,treevertex1));
	    }
	}
      else
	{
	  if (vertex1<empty)
	    {
	      return((capacity(flow_edge)-flow[flow_edge])*edgeComponents(treevertex1,treevertex2));
	    }
	  else 
	    {
	      return(-1*flow[flow_edge]*edgeComponents(treevertex2,treevertex1));
	    }
	}
    }				
}

int MatchPath::who(int vertex)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  int source=0;
  int sink=nbVertex-1;
  int EMPTY=-1;
  
  if (vertex==-1) return EMPTY;
  assert((vertex!=source)&&(vertex!=sink));
  if (ni==nj)
    {
      if (vertex<=ni) 
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
  else
    { 
      int empty=ni+1;
      if (vertex==empty) 
	return(EMPTY);
      if (vertex<empty) 
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
// 	  for (int i=0;i<vertex-ni-2-1;i++)
	  for (int i=empty+1;i<vertex;i++)
	    begin++;
	  return(*begin);
	}
    }
	
}

//--------------------------------------------------------------------------
// FUNCTIONS USED TO SIMULATE THE FLOWGRAPH DATA STRUCTURE
//--------------------------------------------------------------------------

// Cette fonction renvoie le nombre de sommets relie au sommet n

int MatchPath::nbOut(int n)
{

  int ni=_inputList->size();
  int nj=_referenceList->size();
  if (_generalCapacity==0)
    {
      if (ni<nj) 
	ni++;
      if (ni>nj)
	nj++;
    }
  if (n==0) 
    return(ni);
  else
    {
      if (n<=ni) 
	return(nj+1);
      else
	{
	  if (n<=ni+nj) 
	    return(ni+1);
	  else
	    return(nj);
	}
      //  if (n==ni+nj+1) {return(nj);}
    }
 
}

int MatchPath::next_edge(int n,int i)
{
  int ni=_inputList->size();
  int nj=_referenceList->size();
  if (_generalCapacity==0)
    {
      if (ni<nj) ni++;
      if (ni>nj) nj++;
    }

  if (n==0) {return(2*(i-1));}
  if (n==ni+nj+1) {return(2*ni+2*ni*nj+2*(i-1)+1);}
  if (n<=ni)
    {
      if (i==nbOut(n))
	{
	  return(2*(n-1)+1);
	}
      else
	{
	  return(2*ni+2*nj*(n-1)+2*(i-1));
	}
    }
  else
    {
      if (i==nbOut(n))
	{
	  return(2*(ni+ni*nj+(n-ni-1)));
	}
      else
	{
	  return(2*ni+2*nj*(i-1)+2*(n-ni-1)+1);	
	}
    }
}

int MatchPath::next_vertex(int n,int i)
{
		
  int ni=_inputList->size();
  int nj=_referenceList->size();
  if (_generalCapacity==0)
    {
      if (ni<nj) ni++;
      if (ni>nj) nj++;
    }
  
  if (n==0) return(i);
  if (n==ni+nj+1) return(ni+i);
  if (n<=ni)
    {      
      if (i==nbOut(n))        
	{       
	  return(0);      
	}       
      else    
	{       
	  return(ni+i);        
	}
    }
  else
    {
      if (i==nbOut(n))
	{
	  return(ni+nj+1);
	}
      else
	{
	  return(i);
	}
    }
  
}
//--------------------------------------------------------------------------















