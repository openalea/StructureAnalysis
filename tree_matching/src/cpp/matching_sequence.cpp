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


#include "matching_sequence.h"


  // -------------
  // Constructeur
  // -------------
MatchingSequence::MatchingSequence(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance, int s)
{
  T1=&input;
  T2=&reference;
  ND=&nodeDistance;
  _scale = s;
  _distances.make(*T1,*T2,nodeDistance);
  _insertCost.make(*T1,*T2,nodeDistance);
  _deleteCost.make(*T1,*T2,nodeDistance);
  _substitutionCost.make(*T1,*T2,nodeDistance);
  _sumNbCaseVector = CaseVector(3,0);
  _nbCaseVector = CaseVector(7,0);
  _distanceMatrix = DistanceVectorTable(T1->getNbVertex());
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  // constante qui va permettre de calculer l'alignement restreint
  _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
}
// -------------
// Destructeur
  // -------------
MatchingSequence::~MatchingSequence()
{
  int size1 = T1->getNbVertex();
  int size2 = T2->getNbVertex();
//   for (int iv=0;iv<=size1-1;iv++)
//     {
//       delete (DistanceType*) _distanceMatrix[iv];
//     }
//   delete (DistanceType**) _distanceMatrix;
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux sequences S1[input_vertex] et S2[reference_vertex]
// ----------------------------------------------------------------------------------
DistanceType MatchingSequence::distanceBetweenTree(int input_vertex,int reference_vertex)
{
 
  DistanceType insert,del,sub;
  DistanceType min,MIN=2*MAXDIST;
  int MTC=0;
  int i;

// Case 1
  int input_child=T1->child(input_vertex,1);
  min = getDBT(input_child,reference_vertex) + getDCost(input_vertex);
  if (min<MIN) { MIN=min; MTC=1; }

// Case 2  
  int reference_child=T2->child(reference_vertex,1);
  min = getDBT(input_vertex,reference_child) + getICost(reference_vertex);
  if (min<MIN) { MIN=min; MTC=2; }

//Case3
  ChoiceTable c = ChoiceTable();
  min = getDBT(input_child,reference_child) + getCCost(input_vertex,reference_vertex,c);

  if (min<MIN) { MIN=min; MTC=3;  }

  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
    //-------------------------------
  //We maintain the matching lists
  // On maintient les listes d'alignement
  //-------------------------------
  _choices.createList(input_vertex,reference_vertex);
  _choices.putFirst(input_vertex,reference_vertex,MTC);

    switch (MTC)
      {
      case 1 :{
	_choices.putFirst(input_vertex,reference_vertex,input_child);
	_choices.putLast(input_vertex,reference_vertex,-1);
      }break;
      case 2 :{
	_choices.putFirst(input_vertex,reference_vertex,reference_child);
	_choices.putLast(input_vertex,reference_vertex,M(input_vertex,reference_child));
      }break;
      case 3 :{
	_choices.putFirst(input_vertex,reference_vertex,-1);
	_choices.putLast(input_vertex,reference_vertex,reference_vertex);
      }break;
      default :   assert(0);break;
      }
    _choices.putFirst(input_vertex,reference_vertex,MTC);
    _insertCost.putDBT(input_vertex,reference_vertex,insert);
    _deleteCost.putDBT(input_vertex,reference_vertex,del);
    _substitutionCost.putDBT(input_vertex,reference_vertex,sub);
    _distances.putDBT(input_vertex,reference_vertex,MIN);
  //  cout<<endl<<"Distance Between "<<input_vertex<<" and "<<reference_vertex<<" = "<<MIN<<endl;
    //    delete (Sequence*) s;
  return(MIN);
}

void MatchingSequence::getList(int input_vertex, int reference_vertex, Sequence* sequence)
{
  TreeList(input_vertex,reference_vertex,*sequence);
}

int MatchingSequence::Lat(ChoiceList* L, int vertex)
{
  ChoiceList::iterator begin;
  begin = L->begin();
  for (int i=0;i<vertex;i++)
    begin++;
  return(*begin);
}



void MatchingSequence::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      switch(tree_choice)
        {
        case 1:
          {
            TreeList(Lat(L,1),reference_vertex,sequence);
          }
          break;
        case 2:
          {
            TreeList(input_vertex,Lat(L,1),sequence);
          }
          break;
        case 3: {
          sequence.append(input_vertex,reference_vertex,_distances.getCCost(input_vertex,reference_vertex));
          int i_node=T1->child(input_vertex,1);
          int r_node=T2->child(reference_vertex,1);
	  if (i_node != -1)
	    TreeList(i_node,r_node,sequence);
        }break;
        default : break;
        }
    }
}


DistanceType MatchingSequence::match()
{
  const int size1 = T1->getNbVertex();
  const int size2 = T2->getNbVertex();
  DistanceType D=0;
  _distanceMatrix.resize(size1);
  if ( _scale==2)
    cerr << "\x0d" << "Already computed : 0% matched ...                                   " << flush;
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      for (int input_vertex=size1-1;input_vertex>=0;input_vertex--)
        {
          _distances.openDistancesVector(input_vertex);
          _insertCost.openDistancesVector(input_vertex);
          _deleteCost.openDistancesVector(input_vertex);
          _substitutionCost.openDistancesVector(input_vertex);
    	  _distanceMatrix[input_vertex].resize(size2);
          for (int reference_vertex=size2-1;reference_vertex>=0;reference_vertex--)
            {
	      DistanceType d = distanceBetweenTree(input_vertex,reference_vertex);
	      //(_distanceMatrix[input_vertex])[reference_vertex]=d;
            }
	  if ((int(100. - 100*input_vertex/size1)%5 == 0)&&(_scale==2))
			cerr << "\x0d" << "Already computed : "<<int(100. - 100*input_vertex/size1) <<"% " <<" matched ...                                   " << flush;
          for (int i=1;i<=T1->getNbChild(input_vertex);i++)
            {
              _distances.closeDistancesVector(T1->child(input_vertex,i));
              _insertCost.closeDistancesVector(T1->child(input_vertex,i));
              _deleteCost.closeDistancesVector(T1->child(input_vertex,i));
              _substitutionCost.closeDistancesVector(T1->child(input_vertex,i));
            }
        }
      if ( _scale==2)
	cerr<<endl;
		D=getDBT(0,0);
		//cerr<<"Distance between = "<<D<<endl;
    }
  else
    {
      if (T1->isNull())
        {
          if (!T2->isNull()) {D=_distances.referenceTreeFromEmpty(0);}
        }
      else
        {
          D=_distances.inputTreeToEmpty(0);
        }
    }
  return(D);
}

DistanceType  MatchingSequence::getInsertCost()
{
  return(getInBT(0,0));
}

DistanceType  MatchingSequence::getDeleteCost()
{
  return(getDeBT(0,0));
}

DistanceType  MatchingSequence::getSubstitutionCost()
{
  return(getSuBT(0,0));
}

// --------------------------------------------
// Renvoie les distances entre arbres et forets
// --------------------------------------------
DistanceType MatchingSequence::getDistanceMatrix(int iv,int rv) const
{
  return(_distanceMatrix[iv][rv]);
}

DistanceType MatchingSequence::getDBT(int input_vertex,int reference_vertex) const
{
  return(_distances.getDBT(input_vertex,reference_vertex));
}

DistanceType MatchingSequence::getInBT(int input_vertex,int reference_vertex) const
{
  return(_insertCost.getDBT(input_vertex,reference_vertex));
}

DistanceType MatchingSequence::getDeBT(int input_vertex,int reference_vertex) const
{
  return(_deleteCost.getDBT(input_vertex,reference_vertex));
}


DistanceType MatchingSequence::getSuBT(int input_vertex,int reference_vertex) const
{
  return(_substitutionCost.getDBT(input_vertex,reference_vertex));
}


// renvoie le dernier element de la liste de la case node du tableau maintenant les listes d'alignement
int MatchingSequence::M(int i_node,int r_node)
{
  return(_choices.getList(i_node,r_node)->back());
}

DistanceType MatchingSequence::getDCost(int input_vertex) const
{
  if (_scale == 1) 
    return _distances.getDCost(input_vertex);
  else
    {
      int root1 = (T1->getNode(input_vertex))->getVertex();
      TreeGraph* tree1 = new TreeGraph(*(T1->getMTG()),root1,COMPO,ANY);
      TreeGraph* tree2 = new TreeGraph();
      NodeCost* MCF=new NodeCost(TOPOLOGIC);
      MatchingSequence* M=new MatchingSequence(*tree1,*tree2,*MCF,_scale-1);
      //effective matching : the result is stored in D
      DistanceType D=M->match();
      Sequence* s=new Sequence();
      M->TreeList(0,0,*s);
      delete (Sequence*) s;
      delete (NodeCost*) MCF;
      delete (MatchingSequence*) M;
      return D;
    }
}

DistanceType MatchingSequence::getICost(int reference_vertex) const
{
  if (_scale == 1) 
    return _distances.getICost(reference_vertex);
  else
    {
      VId root2 = (T2->getNode(reference_vertex))->getVertex();
      TreeGraph* tree1 = new TreeGraph();
      TreeGraph* tree2 = new TreeGraph(*(T2->getMTG()),root2,COMPO,ANY);
      NodeCost* MCF=new NodeCost(TOPOLOGIC);
      MatchingSequence* M=new MatchingSequence(*tree1,*tree2,*MCF,_scale-1);
      //effective matching : the result is stored in D
      DistanceType D=M->match();
      Sequence* s=new Sequence();
      M->TreeList(0,0,*s);
      delete (Sequence*) s;
      delete (NodeCost*) MCF;
      delete (MatchingSequence*) M;
      return D;
    }
}


DistanceType MatchingSequence::getCCost(int input_vertex,int reference_vertex,ChoiceTable c) const
{
  if (_scale == 1) 
    return _distances.getCCost(input_vertex,reference_vertex);
  else
    {
      VId root1 = (T1->getNode(input_vertex))->getVertex();
      TreeGraph* tree1 = new TreeGraph(*(T1->getMTG()),root1,COMPO,ANY);
      VId root2 = (T2->getNode(reference_vertex))->getVertex();
      TreeGraph* tree2 = new TreeGraph(*(T2->getMTG()),root2,COMPO,ANY);
      NodeCost* MCF=new NodeCost(TOPOLOGIC);
      MatchingSequence* M=new MatchingSequence(*tree1,*tree2,*MCF,_scale-1);
      //effective matching : the result is stored in D
      DistanceType D=M->match();
      c = M->getChoiceTable();
      delete (NodeCost*) MCF;
      delete (MatchingSequence*) M;
      return D;
    }
}

ChoiceTable MatchingSequence::getChoiceTable() const
{
  return _choices;
}





