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


#include "treematchext.h"
#include <math.h>
#include <fstream>
using namespace std;
	 
TreeMatchExtract::TreeMatchExtract(	MTG& mtg,
			TreeMatch& treematch)
{
	_mtg=&mtg;
	_treematch=&treematch;
}


void TreeMatchExtract::statistics()
{
  int nb_tree     = _treematch->getNbTree();
  int nb_matching = ((nb_tree-1)*nb_tree)/2;
  int inp_tree    = 0;
  int ref_tree    = 0;
  int v=0;

// Order Rates
  
  _delVertexRate.reshape(_treematch->getMaxOrder()+1);
  _delVertexRate=(Rate) 0;
  _insVertexRate.reshape(_treematch->getMaxOrder()+1);
  _insVertexRate=(Rate) 0;
  _subVertexRate.reshape(_treematch->getMaxOrder()+1);
  _subVertexRate=(Rate) 0;
  _matVertexRate.reshape(_treematch->getMaxOrder()+1);
  _matVertexRate=(Rate) 0;
 
  for (inp_tree=0; inp_tree<nb_tree-1;inp_tree++)
  {
    for(ref_tree=inp_tree+1;ref_tree<nb_tree;ref_tree++)
    {
      TreeGraph* T1=_treematch->getTree(inp_tree);
      TreeGraph* T2=_treematch->getTree(ref_tree);
      RWTValVector<int>* del_vertex_vector= new RWTValVector<int>(_treematch->getMaxOrder()+1,0);
      RWTValVector<int>* ins_vertex_vector= new RWTValVector<int>(_treematch->getMaxOrder()+1,0);
      RWTValVector<int>* sub_vertex_vector= new RWTValVector<int>(_treematch->getMaxOrder()+1,0);
      RWTValVector<int>* mat_vertex_vector= new RWTValVector<int>(_treematch->getMaxOrder()+1,0);
      RWTValVector<int>* nb_vertex_by_order = new RWTValVector<int>(_treematch->getMaxOrder()+1,0);
 

      // int nb_imp_vertex=0;
      // int nb_ref_vertex=0;	
	//nb_imp_vertex=T1->getNbVertex();
	//nb_ref_vertex=T2->getNbVertex();

      Sequence* s=_treematch->getSequence(inp_tree,ref_tree);
	  
      // input tree
      for (v=0;v<T1->getNbVertex();v++)
      {
	VId vertex=T1->getNode(v)->getVertex();
	int order=T1->getNode(v)->getOrder();

	  switch(editOp(s,vertex,0))
	  {
	           case DEL : (*del_vertex_vector)[order]++;break;
		   case INS : (*ins_vertex_vector)[order]++;break;
		   case MAT : (*mat_vertex_vector)[order]++;break;
		   case SUB : (*sub_vertex_vector)[order]++;break;
	  }
	  (*nb_vertex_by_order)[order]++;
     }

	
      // reference tree
      for (v=0;v<T2->getNbVertex();v++)
      {
	VId vertex=T2->getNode(v)->getVertex();
	int order=T2->getNode(v)->getOrder();

	  switch(editOp(s,vertex,1))
	  {
	           case DEL : (*del_vertex_vector)[order]++;break;
		   case INS : (*ins_vertex_vector)[order]++;break;
		   case MAT : (*mat_vertex_vector)[order]++;break;
		   case SUB : (*sub_vertex_vector)[order]++;break;
	  }
	  (*nb_vertex_by_order)[order]++;
     }
         
      for (int order=0;order<=_treematch->getMaxOrder();order++)
      {
	_delVertexRate[order] = _delVertexRate[order]+((float) (*del_vertex_vector)[order])/((float) (*nb_vertex_by_order)[order]);
        _insVertexRate[order] = _insVertexRate[order]+((float) (*ins_vertex_vector)[order])/((float) (*nb_vertex_by_order)[order]);
	_matVertexRate[order] = _matVertexRate[order]+((float) (*mat_vertex_vector)[order])/((float) (*nb_vertex_by_order)[order]);
	_subVertexRate[order] = _subVertexRate[order]+((float) (*sub_vertex_vector)[order])/((float) (*nb_vertex_by_order)[order]);
      };
	
      delete (RWTValVector<int>*) del_vertex_vector ;
      delete (RWTValVector<int>*) ins_vertex_vector ;
      delete (RWTValVector<int>*) sub_vertex_vector ;
      delete (RWTValVector<int>*) mat_vertex_vector ;
      delete (RWTValVector<int>*) nb_vertex_by_order ;
    };
  };

 for (int order=0;order<=_treematch->getMaxOrder();order++)
 {
   _delVertexRate[order] = _delVertexRate[order]/nb_matching;
   _insVertexRate[order] = _insVertexRate[order]/nb_matching;
   _matVertexRate[order] = _matVertexRate[order]/nb_matching;
   _subVertexRate[order] = _subVertexRate[order]/nb_matching;
 };
	
}

void TreeMatchExtract::plot_statistics(const char* prefix_file_name)
{
        statistics();
	RWCString name(prefix_file_name);
	int nb_tree=_treematch->getNbTree();
	int i=0,k=0;
	
	// writing of the rate file
 
	RWCString rate_dat_name=name+"_RATE.dat";
	fstream rate_dat_file(rate_dat_name.data(),ios::out);
	int order_max=_treematch->getMaxOrder()+1;
	for (k=0;k<order_max;k++)
	{
	  rate_dat_file<<k+0.0<<" "<<_delVertexRate[k]<<" "
	               <<k+0.2<<" "<<_insVertexRate[k]<<" "
		       <<k+0.4<<" "<<_matVertexRate[k]<<" "
		       <<k+0.6<<" "<<_subVertexRate[k]<<endl;
	};
	rate_dat_file.close();
		
	int tree=0;
	RWCString torder_dat_name=name+"_TORDER.dat";
	fstream torder_dat_file(torder_dat_name.data(),ios::out);
	for (tree=0;tree<nb_tree;tree++) 
	{ 
	  torder_dat_file<<tree<<" "<<_treematch->getTree(tree)->getOrder()<<endl; 
	};
	torder_dat_file.close();
	RWCString tvertex_dat_name=name+"_TVERTEX.dat";
	fstream tvertex_dat_file(tvertex_dat_name.data(),ios::out);
	for (tree=0;tree<nb_tree;tree++) 
	{ 
	  tvertex_dat_file<<tree<<" "<<_treematch->getTree(tree)->getNbVertex()<<endl; 
	};
	tvertex_dat_file.close();

	int itree;
	int rtree;
	RWCString tdist_dat_name1=name+"_TDIST1.dat";
	fstream tdist_dat_file1(tdist_dat_name1.data(),ios::out);
	DistanceType distance_max=1e-5;
	DistanceVector average_distance(nb_tree+1,0.0);
	for (itree=0;itree<nb_tree;itree++) 
	{ 
	  for (rtree=0;rtree<nb_tree;rtree++) 
	  { 
	    average_distance[itree]= average_distance[itree]+(_treematch->getDistance(itree,rtree))/((float) nb_tree);
	    tdist_dat_file1<<itree<<" "<<rtree<<" "<<_treematch->getDistance(itree,rtree)<<endl; 
	    distance_max=D_MAX(distance_max,_treematch->getDistance(itree,rtree));
	  }
	};
	tdist_dat_file1.close();

		
	RWCString tdist_dat_name2=name+"_TDIST2.dat";
	fstream tdist_dat_file2(tdist_dat_name2.data(),ios::out);
	RWTValVector<int> nb_couple(100,0);
	DistanceType dist_max=1.01*distance_max;
	for (itree=0;itree<nb_tree-1;itree++) 
	{ 
	  for (rtree=itree+1;rtree<nb_tree;rtree++) 
	  { 
	    int index = LINGT(100.0*_treematch->getDistance(itree,rtree)/dist_max);
	    nb_couple[index]++;
	  }
	};
	for (i=0;i<100;i++)
	{
	  tdist_dat_file2<<(((float) i*dist_max)/100.0)<<" "<<nb_couple[i]<<endl; 
	}
	tdist_dat_file2.close();


	RWCString tdist_dat_name3=name+"_TDIST3.dat";
	fstream tdist_dat_file3(tdist_dat_name3.data(),ios::out);
	average_distance[nb_tree]=MAXDIST;
	_newOrder.insert(nb_tree);
	for (itree=0;itree<nb_tree;itree++) 
	{  
	  int nb_elem=_newOrder.entries();
	  AmlBoolean inserted=FALSE;
	  for (int n=0;n<nb_elem;n++)
	  {
	     if ((!inserted)&&(average_distance[itree]<average_distance[_newOrder.at(n)]))
	     {
	       _newOrder.insertAt(n,itree);
	       inserted=TRUE;
	     } 
	    	    
	   }
	};
	
	for (itree=0;itree<nb_tree;itree++) 
	{
	  tdist_dat_file3<<average_distance[_newOrder.at(itree)]<<endl; 
	}
	tdist_dat_file3.close();


	fstream out_file;
	RWCString file_name;
	  
	for (int file_type = 1 ; file_type <= 2 ; file_type++) 
	{
	  switch (file_type) 
	  {
	    case PLOT : 
	    {
	      file_name=name+".plot";
	      out_file.open(file_name.data(), ios::out ) ; break;
	    };

	    case PRINT : 
	    {
	      file_name=name+".print";
	      out_file.open(file_name.data(), ios::out ) ; break;
	    };
	    default : break;
	  }
		  
	  if (file_type == PRINT )
	  {
	    out_file << "set terminal postscript" << endl;
	    out_file << "set output \"" << name.data() <<".eps"<< "\"\n\n";
	  }
       
	  out_file << "set noborder\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
	  out_file << "\nset title \"NUMBER OF VERTEX PER TREE \"";
	  out_file << "\nset xtics 0,1" << endl;
	  out_file << "\nset ytics " << endl;
	  out_file << "set grid" << endl;
	  out_file << "\nset xlabel \" TREE \"" << endl;
	  out_file << "\nset ylabel \" NB VERTEX \"" << endl;
	  out_file << "plot [0:" << nb_tree << "] [0:"
	           << _treematch->getNbVertexMax()+1 << "] \"" << tvertex_dat_name.data()
		   << "\" using 1:2 title "<<"\" NB VERTEX "
		   << "\" with impulses";
	  out_file << "\npause -1 \"<Return>: continue.\"" << endl;

	  out_file << "set noborder\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
	  out_file << "\nset title \" MAXIMUM ORDER PER TREE \"";
	  out_file << "\nset xtics 0,1" << endl;
	  out_file << "\nset ytics " << endl;
	  out_file << "set grid" << endl;
	  out_file << "\nset xlabel \" TREE \"" << endl;
	  out_file << "\nset ylabel \" ORDER MAX \"" << endl;
	  out_file << "plot [0:" << nb_tree << "] [0:"
		   << _treematch->getMaxOrder()+1 << "] \"" << torder_dat_name.data()
		   << "\" using 1:2 title "<<"\" ORDER "
		   << "\" with impulses";
	  out_file << "\npause -1 \"<Return>: continue.\"" << endl;
		  

	  out_file << "set border\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
	  if (nb_tree<10){ out_file << "\nset xtics 0,1 "<<endl; } else { out_file << "\nset xtics "<<endl; }
	  if (nb_tree<10){ out_file << "\nset ytics 0,1 "<<endl; } else { out_file << "\nset ytics "<<endl; }
	  out_file << "\nset xlabel \" TREE \"" << endl;
	  out_file << "\nset ylabel \" TREE \"" << endl;
	  out_file << "\nset title \" DISTANCE DISTRIBUTION 1\""<<endl;
	  out_file << "set dummy u,v"<<endl;
	  out_file << "set parametric"<<endl;  
	  out_file << "\nsplot [:] [:] "; 
	  out_file << "[0:" << nb_tree-1<<"] "; 
	  out_file << "[0:" << nb_tree-1<<"] "; 
	  int scale = (int) (0.5*powf(10.0 ,(int) log10(distance_max)));
	  DistanceType d_max=scale*(((int) distance_max/scale) + 1);
	  out_file << "["<<distance_max/3.0<<":"<<d_max<<"] "; 
	  out_file << "\""<<tdist_dat_name1.data()<<"\" using 1:2:3 title \"DISTANCE BETWEEN TREES\" with impulses"<<endl; 
	  out_file << "\npause -1 \"<Return>: continue.\"" << endl;
		
	  out_file << "set border\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
	  out_file << "\nset ytics "<<endl;
	  if (nb_tree<10){ out_file << "\nset xtics 0,1 "<<endl; } else { out_file << "\nset xtics "<<endl; }
	  out_file << "\nset xlabel \" DISTANCE \"" << endl;
	  out_file << "\nset ylabel \" TREE \"" << endl;
	  out_file << "set grid" << endl;
	  out_file << "\nset title \" AVERAGE DISTANCE DISTRIBUTION \""<<endl;
	  for (int nb=0;nb<nb_tree;nb++)
	  {
	    out_file<<"set label \""<<_newOrder.at(nb)+1<<"\" at "<<nb<<","<<average_distance[_newOrder.at(nb)]<<endl;
	  }
	  out_file << "plot [0:" << nb_tree << "] [0:"
	           <<  dist_max<< "] \"" << tdist_dat_name3.data()
		   << "\" using 1:2 title "<<"\"TREE "
		   << "\" with impulses"; 
	  out_file << "\npause -1 \"<Return>: continue.\"" << endl;
 
	  out_file << "set border\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
	  out_file << "\nset xtics "<<endl;
	  if (nb_tree<10){ out_file << "\nset ytics 0,1 "<<endl; } else { out_file << "\nset ytics "<<endl; }
	  out_file << "\nset xlabel \" DISTANCE \"" << endl;
	  out_file << "\nset ylabel \" NB TREES \"" << endl;
	  out_file << "set grid" << endl;
	  out_file << "\nset title \" DISTANCE DISTRIBUTION 3\""<<endl;
	  out_file << "plot [0:" << dist_max << "] [0:"
	           <<  nb_tree << "] \"" << tdist_dat_name2.data()
		   << "\" using 1:2 title "<<"\" NB TREES "
		   << "\" with impulses"; 
	  out_file << "\npause -1 \"<Return>: continue.\"" << endl;

	  out_file << "set noborder\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
	  out_file << "\nset title \"EDIT OPERATION RATES PER ORDER \"";
	  out_file << "\nset xtics 0,1" << endl;
	  out_file << "\nset ytics " << endl;
	  out_file << "set grid" << endl;
	  out_file << "\nset xlabel \" ORDER \"" << endl;
	  out_file << "\nset ylabel \" RATE \"" << endl;
	  out_file << "plot [0:"<<_treematch->getMaxOrder()+1<<"] [0:1]"
	           << "\"" <<rate_dat_name.data()<< "\" using 1:2 title "<<"\" DELETION RATE     "<< "\" with impulses,\\"<<endl
		   << "\"" <<rate_dat_name.data()<< "\" using 3:4 title "<<"\" INSERTION RATE    "<< "\" with impulses,\\"<<endl
		   << "\"" <<rate_dat_name.data()<< "\" using 5:6 title "<<"\" MATCHING RATE     "<< "\" with impulses,\\"<<endl
		   << "\"" <<rate_dat_name.data()<< "\" using 7:8 title "<<"\" SUBSTITUTION RATE "<< "\" with impulses"<<endl;
	  out_file << "\npause -1 \"<Return>: continue.\"" << endl;

	  /*
	     out_file << "set noborder\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
	     out_file << "\nset title \"CLASS\"";
	     out_file << "\nset xtics 0,1" << endl;
	     out_file << "\nset ytics " << endl;
	     out_file << "set grid" << endl;
	     out_file << "\nset xlabel \" CLASS \"" << endl;
	     out_file << "\nset ylabel \" RATE \"" << endl;
	     out_file << "plot [0:" << _mtg->classNb()<< "] [0:"
	              << 100 << "] \"" << class_dat_name.data()
		      << "\" using 1:2 title "<<"\" RATE "
		      << "\" with impulses";
	 */
	  out_file << "\npause 0 \"End.\"" << endl;
	  out_file.close();
	};
	
};


void TreeMatchExtract::plot_statistic(const char* prefix_file_name,int inp_tree,int ref_tree)
{
  TreeGraph* T1 = _treematch->getTree(inp_tree);
  TreeGraph* T2 = _treematch->getTree(ref_tree);
  int max_order = _treematch->getMaxOrder()+1;
  RWTValVector<int>* nb_del_vertex = new RWTValVector<int>(max_order);
  RWTValVector<int>* nb_ins_vertex = new RWTValVector<int>(max_order);
  RWTValVector<int>* nb_mat_vertex = new RWTValVector<int>(max_order);
  RWTValVector<int>* nb_sub_vertex = new RWTValVector<int>(max_order);
  RWTValVector<int>* nb_vertex_by_order = new RWTValVector<int>(max_order);
  RateVector* del_vertex_rate = new RateVector(2*max_order,0);
  RateVector* ins_vertex_rate = new RateVector(2*max_order,0);
  RateVector* sub_vertex_rate = new RateVector(2*max_order,0);
  RateVector* mat_vertex_rate = new RateVector(2*max_order,0);
  
  Sequence* s=_treematch->getSequence(inp_tree,ref_tree);
 
  int v=0;
  int order=0;

  // input tree
  (*nb_del_vertex)=0;
  (*nb_ins_vertex)=0;
  (*nb_mat_vertex)=0;
  (*nb_sub_vertex)=0;
  (*nb_vertex_by_order)=0;

  
  for (v=0;v<T1->getNbVertex();v++)
  {
    VId vertex=T1->getNode(v)->getVertex();
    int order=T1->getNode(v)->getOrder();
      switch(editOp(s,vertex,0))
      {
	           case DEL : (*nb_del_vertex)[order]++;break;
		   case INS : (*nb_ins_vertex)[order]++;break;
		   case MAT : (*nb_mat_vertex)[order]++;break;
		   case SUB : (*nb_sub_vertex)[order]++;break;
      }
      (*nb_vertex_by_order)[order]++;
    
  }
  for (order=0;order<max_order;order++)
  {
    (*del_vertex_rate)[order] = ((float) (*nb_del_vertex)[order])/((float) (*nb_vertex_by_order)[order]);
    (*ins_vertex_rate)[order] = ((float) (*nb_ins_vertex)[order])/((float) (*nb_vertex_by_order)[order]);
    (*mat_vertex_rate)[order] = ((float) (*nb_mat_vertex)[order])/((float) (*nb_vertex_by_order)[order]);
    (*sub_vertex_rate)[order] = ((float) (*nb_sub_vertex)[order])/((float) (*nb_vertex_by_order)[order]);
  };

  // reference tree
  (*nb_del_vertex)=0;
  (*nb_ins_vertex)=0;
  (*nb_mat_vertex)=0;
  (*nb_sub_vertex)=0;
  (*nb_vertex_by_order)=0;

  for (v=0;v<T2->getNbVertex();v++)
  {
    VId vertex=T2->getNode(v)->getVertex();
    int order=T2->getNode(v)->getOrder();
    switch(editOp(s,vertex,1))
      {
      case DEL : (*nb_del_vertex)[order]++;break;
      case INS : (*nb_ins_vertex)[order]++;break;
      case MAT : (*nb_mat_vertex)[order]++;break;
      case SUB : (*nb_sub_vertex)[order]++;break;
      }
    (*nb_vertex_by_order)[order]++;
    
  }
         
  for (order=0;order<max_order;order++)
  {
	(*del_vertex_rate)[order+max_order] = ((float) (*nb_del_vertex)[order])/((float) (*nb_vertex_by_order)[order]);
        (*ins_vertex_rate)[order+max_order] = ((float) (*nb_ins_vertex)[order])/((float) (*nb_vertex_by_order)[order]);
	(*mat_vertex_rate)[order+max_order] = ((float) (*nb_mat_vertex)[order])/((float) (*nb_vertex_by_order)[order]);
	(*sub_vertex_rate)[order+max_order] = ((float) (*nb_sub_vertex)[order])/((float) (*nb_vertex_by_order)[order]);
  };

  RWCString name(prefix_file_name);
  
  RWCString inp_edop_dat_name=name+"_INP_EDOP_RATE.dat";
  fstream inp_edop_dat_file(inp_edop_dat_name.data(),ios::out);
  for (order=0;order<max_order;order++)
  {
    inp_edop_dat_file<<order+0.0<<" "<<(*del_vertex_rate)[order]<<" "
                     <<order+0.2<<" "<<(*ins_vertex_rate)[order]<<" " 
		     <<order+0.4<<" "<<(*mat_vertex_rate)[order]<<" "  
		     <<order+0.6<<" "<<(*sub_vertex_rate)[order]<<endl; 
  };
  inp_edop_dat_file.close();

  RWCString ref_edop_dat_name=name+"_REF_EDOP_RATE.dat";
  fstream ref_edop_dat_file(ref_edop_dat_name.data(),ios::out);
  for (order=0;order<max_order;order++)
  {
    ref_edop_dat_file<<order+0.0<<" "<<(*del_vertex_rate)[order+max_order]<<" "
                     <<order+0.2<<" "<<(*ins_vertex_rate)[order+max_order]<<" " 
		     <<order+0.4<<" "<<(*mat_vertex_rate)[order+max_order]<<" "  
		     <<order+0.6<<" "<<(*sub_vertex_rate)[order+max_order]<<endl; 
  };
  ref_edop_dat_file.close();
	
  delete (RWTValVector<int>*) nb_del_vertex ;
  delete (RWTValVector<int>*) nb_ins_vertex ;
  delete (RWTValVector<int>*) nb_mat_vertex ;
  delete (RWTValVector<int>*) nb_sub_vertex ;
  delete (RWTValVector<int>*) nb_vertex_by_order ;
  delete (RateVector*) del_vertex_rate ;
  delete (RateVector*) ins_vertex_rate ;
  delete (RateVector*) sub_vertex_rate ;
  delete (RateVector*) mat_vertex_rate ;


  // Making of the GNUPLOT files
  fstream out_file;
  RWCString file_name;

  for (int file_type = 1; file_type <= 2; file_type++) 
  {
    switch (file_type) 
    {
      case PLOT : 
      {
	file_name=name+".plot";
	out_file.open(file_name.data(), ios::out ) ; break;
      };
      case PRINT : 
      {
	file_name=name+".print";
	out_file.open(file_name.data(), ios::out ) ; break;
      };
      default : break;
    }

    if (file_type == PRINT )
    {
      out_file << "set terminal postscript" << endl;
      out_file << "set output \"" << name.data() <<".eps"<< "\"\n\n";
    }
	
    out_file << "set noborder\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
    out_file << "\nset xtics 0,1" << endl;
    out_file << "set ytics " << endl;
    out_file << "set grid" << endl;
    out_file << "\nset title \"" <<" MATCHING "<<inp_tree<<" WITH "<<ref_tree<<" : EDIT OPERATIONS \"" << endl;
    out_file << "\nset xlabel \"ORDER : INPUT TREE\"";
    out_file << "\nset ylabel \" RATE \"" << endl;
    out_file << "\nset grid" << endl;
    out_file << "plot [0:"<<max_order<<"] [0:1]"
	     << "\"" <<inp_edop_dat_name.data()<< "\" using 1:2 title "<<"\" DELETION RATE     "<< "\" with impulses,\\"<<endl
	     << "\"" <<inp_edop_dat_name.data()<< "\" using 3:4 title "<<"\" INSERTION RATE    "<< "\" with impulses,\\"<<endl
	     << "\"" <<inp_edop_dat_name.data()<< "\" using 5:6 title "<<"\" MATCHING RATE     "<< "\" with impulses,\\"<<endl
	     << "\"" <<inp_edop_dat_name.data()<< "\" using 7:8 title "<<"\" SUBSTITUTION RATE "<< "\" with impulses"<<endl;
    out_file << "\npause -1 \"<Return>: continue.\"" << endl;

    out_file << "set noborder\n" << "set tics out\n" << "set noparametric\n" << "set nolabel\n";
    out_file << "\nset xtics 0,1" << endl;
    out_file << "set ytics " << endl;
    out_file << "set grid" << endl;
    out_file << "\nset title \"" <<" MATCHING "<<inp_tree<<" WITH "<<ref_tree<<" : EDIT OPERATIONS \"" << endl;
    out_file << "\nset xlabel \"ORDER : REFERENCE TREE\"";
    out_file << "\nset ylabel \" RATE \"" << endl;
    out_file << "plot [0:"<<max_order<<"] [0:1]"
	     << "\"" <<ref_edop_dat_name.data()<< "\" using 1:2 title "<<"\" DELETION RATE     "<< "\" with impulses,\\"<<endl
	     << "\"" <<ref_edop_dat_name.data()<< "\" using 3:4 title "<<"\" INSERTION RATE    "<< "\" with impulses,\\"<<endl
	     << "\"" <<ref_edop_dat_name.data()<< "\" using 5:6 title "<<"\" MATCHING RATE     "<< "\" with impulses,\\"<<endl
	     << "\"" <<ref_edop_dat_name.data()<< "\" using 7:8 title "<<"\" SUBSTITUTION RATE "<< "\" with impulses"<<endl;
    out_file << "\npause -1 \"<Return>: continue.\"" << endl;

    out_file << "\npause 0 \"End.\"" << endl;
    out_file.close();
  };

};


EditOperation TreeMatchExtract::editOp(Sequence* s, VId vertex, int side)
{
  s->reset();
  AmlBoolean contains=FALSE;
  DistanceType cost=0.0; 
  do
  {
    VId current;
    if (side)
    {
      current=s->getCurrent()->getRV();
    }
    else
    {
      current=s->getCurrent()->getIV();
    };
    
    if (current==vertex) 
    {
      cost=s->getCurrent()->getCost();
      contains=1;
    }
  }while((s->next())&&(contains==FALSE));
  
  if (contains==TRUE)
  {
    if (cost==0.0)
    {
      return(MAT);
    }
    else
    {
      return(SUB);
    }
  }
  else
  {
    if (side)
    {
      return(INS);
    }
    else
    {
      return(DEL);
    };
  };
}




