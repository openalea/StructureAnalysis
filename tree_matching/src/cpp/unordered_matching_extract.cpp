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


#include "unordered_treematch.h"
#include "treenode.h"
#include <math.h>
#include "Tools/timer.h"

#ifdef SYSTEM_IS__Linux
// AML2
// Apparament, il y a un probleme avec la fonction truncate,
// qui pourtant est bien defini dans mathcall.h inclu par math.h
// Je ne sais pas faire mieux que de la definir "a la main".
extern "C" double   trunc(double);
#endif


//------------------------------------------------------------------------------------------------------
// DISPLAY THE RESULTS OF THE MATCHING
//------------------------------------------------------------------------------------------------------
ostream& TreeMatch_U::viewAllMatching(ostream& out_fich) const
{
  char TAB='\t';
  out_fich<<TAB;
  int length=_roots->entries();
  int i,j,k;
  int nb_matching=length*(length-1)/2;
  DistanceType max_distance=MINDIST;
  DistanceType min_distance=MAXDIST;
  DistanceType average_distance=0.0;

  float ins_rate=0.0;
  float sub_rate=0.0;
  float del_rate=0.0;
  float mat_rate=0.0;


  out_fich<<" MATCHING OF "<<length<<" TREES"<<endl;
  out_fich<<endl<<" DISTANCE MATRIX "<<endl;
  out_fich<<TAB;
  for (k=0;k<2*length;k++) { out_fich <<"--------"; };out_fich<<endl;
  out_fich<<TAB;
  for (i=0; i<length;i++)
    {

      out_fich<<"|"<<TAB<<"T"<<i<<TAB;
    }
  out_fich<<"|"<<endl;

  for (i=0; i<length;i++)
    {
      for (k=0;k<=2*length;k++) { out_fich <<"--------"; };out_fich<<endl;
      out_fich<<"|"<<"  T"<<i<<TAB;
      for(j=0;j<length;j++)
        {
          DistanceType D=0;
          if (i!=j) D=0.01*((float) trunc(100.0*getDistance(i,j)));
          out_fich<<"|"<<TAB<<D<<TAB;
        }
      out_fich<<"|"<<endl;
    }
  for (k=0;k<=2*length;k++) { out_fich <<"--------"; };out_fich<<endl;

  out_fich<<endl<<" CONSERVATION RATE MATRIX "<<endl;
  out_fich<<TAB;
  for (k=0;k<length;k++) { out_fich <<"--------"; };out_fich<<endl;
  for (i=0; i<length;i++)
    {
      out_fich<<TAB<<"|"<<"T"<<i;
    }
  out_fich<<TAB<<"|"<<endl;

  for (i=0; i<length;i++)
    {
      for (k=0;k<=length;k++) { out_fich <<"--------"; };out_fich<<endl;
      out_fich<<"|"<<"T"<<i<<"";
      for(int j=0;j<length;j++)
        {

          if (i!=j)
            {
              float T=0;
              int nb_inp_vertex=_trees[i]->getNbVertex();
              int s_size=getSequence(i,j)->getSize();
              T=0.01*trunc(100.0*((float) s_size)/((float) nb_inp_vertex));
              out_fich<<TAB<<"|"<<T;
            }
          else
            {
              out_fich<<TAB<<"|   *";
            };
        }
      out_fich<<TAB<<"|"<<endl;
    }
  for (k=0;k<=length;k++) { out_fich <<"--------"; };out_fich<<endl;



  for (i=0; i<length-1;i++)
    {
      int tree_size1=_trees[i]->getNbVertex();
      for(int j=i+1;j<length;j++)
        {
          Sequence* mat_seq=getSequence(i,j);
          int seq_size=mat_seq->getSize();
          int sub_number=0;
          int mat_number=0;
          mat_seq->reset();
          do
            {
              if (mat_seq->getCurrent()->getCost()==0.0)
                {
                  mat_number++;
                }
              else
                {
                  sub_number++;
                };
            } while(mat_seq->next());
          mat_seq->putNbMat(mat_number);
          mat_seq->putNbSub(sub_number);
          int tree_size2=_trees[j]->getNbVertex();
          mat_seq->putNbDel(tree_size1-seq_size);
          mat_seq->putNbIns(tree_size2-seq_size);
          int nb_edit_op=tree_size1+tree_size2-seq_size;
          del_rate = del_rate + ((float) tree_size1-seq_size)/((float) nb_edit_op);
          ins_rate = ins_rate + ((float) tree_size2-seq_size)/((float) nb_edit_op);
          sub_rate = sub_rate + ((float) sub_number)/((float) nb_edit_op);
          mat_rate = mat_rate + ((float) mat_number)/((float) nb_edit_op);
          DistanceType distance=getDistance(i,j);
          max_distance=I_MAX(max_distance,distance);
          min_distance=I_MIN(min_distance,distance);
          average_distance=average_distance+distance;
        }
    }
  average_distance=average_distance/((float) nb_matching);
  del_rate = del_rate/((float) nb_matching);
  sub_rate = sub_rate/((float) nb_matching);
  ins_rate = ins_rate/((float) nb_matching);
  mat_rate = mat_rate/((float) nb_matching);
  out_fich<<endl;
  out_fich<<"MAXIMUM DISTANCE         : "<<max_distance<<endl;
  out_fich<<"MINIMUM DISTANCE         : "<<min_distance<<endl;
  out_fich<<"AVERAGE DISTANCE         : "<<average_distance<<endl;
  out_fich<<endl;
  out_fich<<"AVERAGE INSERTION RATE   : "<<ins_rate<<endl;
  out_fich<<"AVERAGE DELETION  RATE   : "<<del_rate<<endl;
  out_fich<<"AVERAGE MATCHING RATE    : "<<mat_rate<<endl;
  out_fich<<"AVERAGE SUBSTITUTION RATE : "<<sub_rate<<endl;
  out_fich<<endl;

  if (_fileName)
    {
      fstream outfile(_fileName,ios::out);
      outfile<<out_fich<<endl;
      outfile.close();
    }
  return(out_fich);
}

//------------------------------------------------------------------------------------------------------
// GET THE MATRIX OF THE MATCHING
//------------------------------------------------------------------------------------------------------

Distance_matrix* TreeMatch_U::getMatrix()
{
  int nb_tree ;

  nb_tree = getNbTree();

  int nb_del;
  int nb_ins;
  double del_cost;
  double ins_cost;
  Distance_matrix* dmatrix = new Distance_matrix(nb_tree, -1, -1, "ARBORESCENCE");

  for (int inp_tree=0;inp_tree<nb_tree-1;inp_tree++)
    {
      for (int ref_tree=inp_tree+1;ref_tree<nb_tree;ref_tree++)
        {
	  Sequence* s=new Sequence();
	  s=getSequence(inp_tree,ref_tree);
	  nb_del = s->getNbDel();
	  del_cost =s->getDelCost();
	  nb_ins = s->getNbIns();
	  ins_cost = s->getInsCost();
	  
	  dmatrix->update(inp_tree+1,
			  ref_tree+1,
			  getDistance(inp_tree,ref_tree),
			  s->getSize()*2+s->getNbDel() + s->getNbIns(),
			  del_cost,nb_del,
			  ins_cost,nb_ins,
			  s->getNbMat(),getDistance(inp_tree,ref_tree) -
			  ins_cost-del_cost, s->getNbSub());
	  
	  // On remplit l'autre moitie de la matrice triangulaire
	  dmatrix->update(ref_tree+1,
			  inp_tree+1,
			  getDistance(inp_tree,ref_tree),
			  s->getSize()*2+s->getNbDel() + s->getNbIns(),
                          ins_cost , s->getNbIns(),
			  del_cost, s->getNbDel(),
			  s->getNbMat(),getDistance(inp_tree,ref_tree) -
			  ins_cost-del_cost, s->getNbSub()); 
        }
    }
  return(dmatrix);
}






//------------------------------------------------------------------------------------------------------
// DISPLAY THE RESULTS OF THE MATCHING
//------------------------------------------------------------------------------------------------------
ostream& TreeMatch_U::introduce(ostream& out_info) const
{
  int nb_tree=_trees.size();
  int m_nb_vertex=0;
  int m_degree=0;
  for (int i=0;i<nb_tree;i++)
    {
      m_nb_vertex=m_nb_vertex+_trees[i]->getNbVertex();
      m_degree=m_degree+_trees[i]->getDegree();
    }
  out_info<<"NB TREE = "<<nb_tree<<" AVERAGE NB VERTEX = "<<SINLT(m_nb_vertex/nb_tree)<<" AVERAGE DEGREE = "<<SINLT(m_degree/nb_tree)<<endl;
  return out_info;
};

//------------------------------------------------------------------------------------------------------
//  RETURN A MATCHING LIST
//------------------------------------------------------------------------------------------------------
SLArray* TreeMatch_U::getList(int i_tree,int r_tree)
{
  SLArray* matched_vtx=new SLArray;
  SLArray* imp_vtx=new SLArray;
  SLArray* ref_vtx=new SLArray;
  Sequence* list_vtx = getSequence(i_tree,r_tree);
  list_vtx->reset();
  do
    {
      *imp_vtx += AMObj(AMObjType::VTX, list_vtx->getCurrent()->getIV());
      *ref_vtx += AMObj(AMObjType::VTX, list_vtx->getCurrent()->getRV());
    } while(list_vtx->next());
  *matched_vtx += AMObj(AMObjType::ARRAY, imp_vtx);
  *matched_vtx += AMObj(AMObjType::ARRAY, ref_vtx);
  return(matched_vtx);
}

//------------------------------------------------------------------------------------------------------
// RETURN A DISTANCE
//------------------------------------------------------------------------------------------------------
DistanceType TreeMatch_U::getDist(int i_tree,int r_tree)
{
  return(getDistance(i_tree,r_tree));
}

//------------------------------------------------------------------------------------------------------
// RETURN TIME COMPUTATION
//------------------------------------------------------------------------------------------------------
DistanceType TreeMatch_U::getTime(int inp_tree,int ref_tree) const
{
  if (inp_tree==ref_tree) return(0.0);

  if (inp_tree>ref_tree)
    {
      return(_time[ref_tree][inp_tree-ref_tree-1]);
    }
  else
    {
      return(_time[inp_tree][ref_tree-inp_tree-1]);
    }
}


//------------------------------------------------------------------------------------------------------
// RETURN A NORMALIZED DISTANCE
//------------------------------------------------------------------------------------------------------
DistanceType TreeMatch_U::getNormalizedDistance(int i_tree,int r_tree)
{

  int nb_del;
  int nb_ins;
  double del_cost;
  double ins_cost;
  Sequence* s=new Sequence();
  s=getSequence(i_tree,r_tree);
  nb_del = s->getNbDel();
  del_cost =s->getDelCost();
  nb_ins = s->getNbIns();
  ins_cost = s->getInsCost();
  DistanceType d = getDist(i_tree,r_tree);
  return(d / (double(s->getSize()*2+nb_del+nb_ins)));
}

//------------------------------------------------------------------------------------------------------
// RETURN A NORMALIZED DISTANCE
//------------------------------------------------------------------------------------------------------
DistanceType TreeMatch_U::getSelfSimilarNormalizedDistance(int i_tree,int r_tree)
{

  DistanceType d = getDist(i_tree,r_tree);

  int size1 = _trees[0]->getNbDesc(i_tree);
  int size2 = _trees[1]->getNbDesc(r_tree);

  return d/(double(size1+size2));
}

//------------------------------------------------------------------------------------------------------
// RETURN A SelfSimilar DISTANCE
//------------------------------------------------------------------------------------------------------
DistanceType TreeMatch_U::getSelfSimilarDistance(int i_tree,int r_tree) const
{

	int reference = _trees[0]->getNumber(r_tree);
	int input = _trees[0]->getNumber(i_tree);
	DistanceType d = getDistance(input,reference);
    return d;
}


//------------------------------------------------------------------------------------------------------
// RETURN A PARACLADIAL COEFFICIENT
//------------------------------------------------------------------------------------------------------

DistanceType TreeMatch_U::getC(ostream& out_info,int i_tree,int r_tree) const
{
	int father = _trees[0]->getNumber(r_tree);
	int input = _trees[0]->getNumber(i_tree);
	int nb_child = _trees[0]->getNbChild(father);
	if ((father!=-1)&&(input!=-1)){
		DistanceType min = getDistance(input,father);
		int child = _trees[0]->getChildLess(father);
		while (child!=-1)
		{
			//int i = 0;
			//while ((i<nb_child)&&(!_trees[0]->childIsInAxis(father,i)))
			//	i++;
			double d = getDistance(input,child);
			if (min > d)
				min = d;
			father = child;
			child = _trees[0]->getChildLess(father);
		}
		return min;
	}
	else
		return -1.;
}

//------------------------------------------------------------------------------------------------------
// DISPLAY THE RESULTS OF THE MATCHING
//------------------------------------------------------------------------------------------------------
ostream& TreeMatch_U::viewOneMatching(ostream& out_info,int imp_tree,int ref_tree) const
{
  char TAB='\t';
  int i=0,j=0,k=0;
  TreeGraph* Tree1=_trees[imp_tree];
  TreeGraph* Tree2=_trees[ref_tree];
  int nb_imp_vertex=Tree1->getNbVertex();
  int nb_ref_vertex=Tree2->getNbVertex();
  Sequence* s= getSequence(imp_tree,ref_tree);
  VertexList imp_vtx;
  VertexList ref_vtx;
  int nb_del=0;
  int nb_ins=0;
  int nb_mat=0;
  int nb_sub=0;
  int nb_ope=0;
  DistanceType sub_cost=0;
  DistanceType ins_cost=0;
  DistanceType del_cost=0;
  DistanceVector imp_rate_by_class(_mtg->classNb(),0);
  DistanceVector ref_rate_by_class(_mtg->classNb(),0);
  DistanceVector class_imp_nb_vertex(_mtg->classNb(),0);
  DistanceVector class_ref_nb_vertex(_mtg->classNb(),0);
  DistanceVector imp_rate_by_order(1+_maxOrder,0);
  DistanceVector ref_rate_by_order(1+_maxOrder,0);
  DistanceVector order_imp_nb_vertex(1+_maxOrder,0);
  DistanceVector order_ref_nb_vertex(1+_maxOrder,0);


  s->reset();
  do
    {
      VId imp_current=s->getCurrent()->getIV();
      VId ref_current=s->getCurrent()->getRV();
      DistanceType cost=s->getCurrent()->getCost();
      if (cost==0.0) { nb_mat++; } else { nb_sub++; };
      sub_cost=sub_cost+cost;
      imp_vtx.push_back(imp_current);
      ref_vtx.push_back(ref_current);
    } while(s->next());

  del_cost = s->getDelCost();
  ins_cost = s->getInsCost();

  nb_del=nb_imp_vertex-s->getSize();
  nb_ins=nb_ref_vertex-s->getSize();
  nb_ope=nb_del+nb_ins+nb_mat+nb_sub;


  // input_tree
  for (i=0;i<Tree1->getNbVertex();i++)
    {
      VId vertex=Tree1->getNode(i)->getVertex();
      VClass v_class=_mtg->vclass(vertex);
      int order=Tree1->getNode(i)->getOrder();
      int contains=0;
      VertexList::iterator begin = imp_vtx.begin();
      VertexList::iterator end = imp_vtx.end();
      while((begin!=end)||(*begin!=vertex))
        {
          begin++;
          if (*begin==vertex) contains=1;
        }
      if (contains)
        {
          imp_rate_by_class[v_class] = imp_rate_by_class[v_class]+1.0;
          imp_rate_by_order[order]   = imp_rate_by_order[order]+1.0;
        }
      class_imp_nb_vertex[v_class] = class_imp_nb_vertex[v_class]+1.0;
      order_imp_nb_vertex[order]   = order_imp_nb_vertex[order]+1.0;

    }

  // reference_tree
  for (j=0;j<Tree2->getNbVertex();j++)
    {
      VId vertex=Tree2->getNode(j)->getVertex();
      VClass v_class=_mtg->vclass(vertex);
      int order=Tree2->getNode(j)->getOrder();
      int contains=0;
      VertexList::iterator begin = ref_vtx.begin();
      VertexList::iterator end = ref_vtx.end();
      while((begin!=end)||(*begin!=vertex))
        {
          begin++;
          if (*begin==vertex) contains=1;
        }
      if (contains)
        {
          ref_rate_by_order[order]   = ref_rate_by_order[order]+1.0;
          ref_rate_by_class[v_class] = ref_rate_by_class[v_class]+1.0;
        };
      class_ref_nb_vertex[v_class] = class_ref_nb_vertex[v_class]+1.0;
      order_ref_nb_vertex[order]   = order_ref_nb_vertex[order]+1.0;

    }

  // printing of the results
  out_info<<" MATCHING PLANT "<<imp_tree<<" WITH PLANT "<<ref_tree<<endl<<endl;
  out_info<<"DISTANCE : "<<getDistance(imp_tree,ref_tree)<<endl;
  out_info<<"EDIT OPERATIONS "<<endl;
  out_info<<TAB<<TAB;
  for (k=0;k<3;k++) { out_info <<"----------------"; };out_info<<endl;
  out_info<<TAB<<TAB<<"|"<<" nb operations |"<<TAB
          <<"rate"<<TAB<<"|"<<TAB
          <<"value"<<TAB<<"|"<<endl;
  for (k=0;k<4;k++) { out_info <<"----------------"; };out_info<<endl;
  out_info<<"| "<<"DELETION"<<TAB
          <<"|"<<TAB<<nb_del<<TAB
          <<"|"<<TAB<<0.01*trunc(100.0*((float) nb_del/nb_ope))<<TAB
          <<"|"<<TAB<<del_cost<<TAB<<"|"<<endl;
  out_info<<"| "<<"INSERTION"<<TAB
          <<"|"<<TAB<<nb_ins<<TAB
          <<"|"<<TAB<<0.01*trunc(100.0*((float) nb_ins/nb_ope))<<TAB
          <<"|"<<TAB<<ins_cost<<TAB<<"|"<<endl;
  out_info<<"| "<<"MATCHING"<<TAB
          <<"|"<<TAB<<nb_mat<<TAB
          <<"|"<<TAB<<0.01*trunc(100.0*((float) nb_mat/nb_ope))<<TAB
          <<"|"<<TAB<<"0"<<TAB<<"|"<<endl;
  out_info<<"| "<<"SUBSTITUTION"<<TAB
          <<"|"<<TAB<<nb_sub<<TAB
          <<"|"<<TAB<<0.01*trunc(100.0*((float) nb_sub/nb_ope))<<TAB
          <<"|"<<TAB<<sub_cost<<TAB<<"|"<<endl;
  for (k=0;k<4;k++) { out_info <<"----------------"; };out_info<<endl;


  out_info<<endl<<"IMPUT TREE : CLASS CONSERVATION RATE "<<endl;
  for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
  out_info<<"|   CLASS"<<TAB<<"|"<<TAB
          <<"rate"<<TAB<<"|"<<endl;
  for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
  for (i=0;i<_mtg->classNb();i++)
    {
      if (class_imp_nb_vertex[i] != 0.0 )
        {
          out_info<<"|"<<TAB<<i<<TAB
                  <<"|"<<TAB<<(int) 100.0*imp_rate_by_class[i]/class_imp_nb_vertex[i]<<TAB<<"|"<<endl;
          for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
        };
    };


  out_info<<endl<<"REFERENCE TREE : CLASS CONSERVATION RATE "<<endl;
  for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
  out_info<<"|   CLASS"<<TAB<<"|"<<TAB<<"rate"<<TAB<<"|"<<endl;
  for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
  for (i=0;i<_mtg->classNb();i++)
    {
      if (class_ref_nb_vertex[i] != 0.0 )
        {
          out_info<<"|"<<TAB<<i<<TAB
                  <<"|"<<TAB<<(int) 100.0*ref_rate_by_class[i]/class_ref_nb_vertex[i]<<TAB<<"|"<<endl;
          for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
        };
    };

  out_info<<endl<<"IMPUT TREE : ORDER CONSERVATION RATE "<<endl;
  for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
  out_info<<"|   ORDER"<<TAB<<"|"<<TAB
          <<"rate"<<TAB<<"|"<<endl;
  for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
  for (i=0;i<=_maxOrder;i++)
    {
      if (order_imp_nb_vertex[i] != 0.0 )
        {
          out_info<<"|"<<TAB<<i<<TAB
                  <<"|"<<TAB<<(int) 100.0*imp_rate_by_order[i]/order_imp_nb_vertex[i]<<TAB<<"|"<<endl;
          for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
        };
    };


  out_info<<endl<<"REFERENCE TREE : ORDER CONSERVATION RATE "<<endl;
  for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
  out_info<<"|   ORDER"<<TAB<<"|"<<TAB<<"rate"<<TAB<<"|"<<endl;
  for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
  for (i=0;i<=_maxOrder;i++)
    {
      if (order_ref_nb_vertex[i] != 0.0 )
        {
          out_info<<"|"<<TAB<<i<<TAB
                  <<"|"<<TAB<<(int) 100.0*ref_rate_by_order[i]/order_ref_nb_vertex[i]<<TAB<<"|"<<endl;
          for (k=0;k<2;k++) { out_info <<"----------------"; };out_info<<endl;
        };
    };
  return(out_info);
}
//------------------------------------------------------------------------------------------------------
// DISPLAY THE DISTANCE BETWEEN PLANT
//------------------------------------------------------------------------------------------------------
DistanceType TreeMatch_U::viewDistanceMatching(ostream& out_info,int inp_tree,int ref_tree) const
{
  if (inp_tree == ref_tree)
    {
      //out_info<<0.<<endl;
      return ( 0.);
    }
  else
    {
	  //cout<<inp_tree<<" - "<<ref_tree<<" = "<<getDistance(inp_tree,ref_tree)<<endl;
      return(getDistance(inp_tree,ref_tree));
    }
}

DistanceType TreeMatch_U::viewNormalizedDistance(ostream& out_info,int i_tree,int r_tree)
{
  if (i_tree == r_tree)
    {
      return ( 0.);
    }
  else
    {
      int nb_del;
      int nb_ins;
      double del_cost;
      double ins_cost;
      Sequence* s=new Sequence();
      s=getSequence(i_tree,r_tree);
      nb_del = s->getNbDel();
      del_cost =s->getDelCost();
      nb_ins = s->getNbIns();
      ins_cost = s->getInsCost();
      DistanceType d = getDist(i_tree,r_tree);
      //out_info<<d<<endl;
      return(d / (double(s->getSize()*2+nb_del+nb_ins)));
    }
}

int TreeMatch_U::viewVertex(int i_tree)
{
  TreeNode* tree_node=_trees[0]->getNode(i_tree);
  int vertex = tree_node->getVertex();
  //  delete (tree_node);
    //out_info<<vertex<<endl;
  return(vertex);
}

void TreeMatch_U::putDistance(DistanceType distance,int inp_tree,int ref_tree)
{
  if (inp_tree>ref_tree)
    {
      _distances[ref_tree][inp_tree-ref_tree-1]=distance;
    }
  else
    {
      _distances[inp_tree][ref_tree-inp_tree-1]=distance;
    }
}

void TreeMatch_U::putTime(DistanceType time,int inp_tree,int ref_tree)
{
  if (inp_tree>ref_tree)
    {
      _time[ref_tree][inp_tree-ref_tree-1]=time;
    }
  else
    {
      _time[inp_tree][ref_tree-inp_tree-1]=time;
    }
}



void TreeMatch_U::putDistanceMatrix(MatchingDistanceTable matrix)
{
  int nbTree = _trees[0]->getNbVertex();
  _distances.resize(nbTree);
  _sequences.resize(nbTree);
  for (int tree=0;tree<nbTree;tree++)
    {
      //cout<<tree<<endl;
       _sequences[tree]=new SequenceVector(nbTree-tree,(Sequence*) NULL);
       //cout<<tree<<endl;
      _distances[tree].resize(nbTree-tree,0.0);

    }
  for (int i_vertex = 0;i_vertex<nbTree;i_vertex++)
    {
      for (int r_vertex = i_vertex+1;r_vertex<nbTree;r_vertex++)
        {
          if (i_vertex>r_vertex)
            {
              _distances[r_vertex][i_vertex-r_vertex-1]=matrix.getDBT(r_vertex,i_vertex);
            }
          else
            {
              _distances[i_vertex][r_vertex-i_vertex-1]=matrix.getDBT(i_vertex,r_vertex);
              cout<<i_vertex<<" - "<<r_vertex<<" D = "<<matrix.getDBT(i_vertex,r_vertex)<<" - "<<getDistance(i_vertex,r_vertex)<<endl;
            }
        }
    }
}

DistanceType TreeMatch_U::getDistance(int inp_tree,int ref_tree) const
{
  if (inp_tree==ref_tree) return(0.0);

  if (inp_tree>ref_tree)
    {
      return(_distances[ref_tree][inp_tree-ref_tree-1]);
    }
  else
    {
      return(_distances[inp_tree][ref_tree-inp_tree-1]);
    }
}

void TreeMatch_U::putSequence(Sequence* sequence,int inp_tree,int ref_tree)
{
  if (inp_tree>ref_tree)
    {
      (*_sequences[ref_tree])[inp_tree-ref_tree-1]=sequence;
    }
  else
    {
      (*_sequences[inp_tree])[ref_tree-inp_tree-1]=sequence;
    }
}

void TreeMatch_U::putSequenceTab(Sequence*** sequence_tab,int inp_tree,int ref_tree)
{
  if (inp_tree>ref_tree)
    {
      (*_sequences_tab[ref_tree])[inp_tree-ref_tree-1]=sequence_tab;
    }
  else
    {
      (*_sequences_tab[inp_tree])[ref_tree-inp_tree-1]=sequence_tab;
    }
}
Sequence* TreeMatch_U::getSequence(int inp_tree,int ref_tree) const
{
  if (inp_tree>ref_tree)
    {
      return((*_sequences[ref_tree])[inp_tree-ref_tree-1]);
    }
  else
    {
      return((*_sequences[inp_tree])[ref_tree-inp_tree-1]);
    }
}
Sequence*** TreeMatch_U::getSequenceTab(int inp_tree,int ref_tree) const
{
  if (inp_tree>ref_tree)
    {
      return((*_sequences_tab[ref_tree])[inp_tree-ref_tree-1]);
    }
  else
    {
      return((*_sequences_tab[inp_tree])[ref_tree-inp_tree-1]);
    }
}

TreeGraph* TreeMatch_U::getTree(int tree) const
{
  return(_trees[tree]);
}

int  TreeMatch_U::getNbTree() const
{
  return(_nbTree);
}

MatchingType TreeMatch_U::getMatchingType() const
{
  return(_matchingType);
}



int  TreeMatch_U::getMaxOrder() const
{
  return(_maxOrder);
}

int  TreeMatch_U::getNbVertexMax() const
{
  return(_maxNbVertex);
}



bool TreeMatch_U::mtg_write( const char *path ) const
{
  return _trees[0]->mtg_write(path);
}

