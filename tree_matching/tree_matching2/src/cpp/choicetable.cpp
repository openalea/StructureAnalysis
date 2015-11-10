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
 *       $Id: choicetable.cpp 3258 2007-06-06 13:18:26Z dufourko $
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


#include "choicetable.h"
#include <assert.h>

ChoiceTable::ChoiceTable(int i_size,int r_size)
{
  _i_size=i_size;
  _r_size=r_size;
  _data.resize(_i_size);
  for (int i=0;i<_i_size;i++)
    {
      // Attention, la liste n'est pas initialisée à 0.
      _data[i].resize(_r_size);
    }
}

ChoiceTable::~ChoiceTable()
{
}

void ChoiceTable::resize(int i_size,int r_size)
{
  _i_size=i_size;
  _r_size=r_size;
  _data.resize(_i_size);
  for (int i=0;i<_i_size;i++)
    {
      // Attention, la liste n'est pas initialisée à 0.
      (_data[i]).resize(_r_size);
    }
}

void ChoiceTable::print(){
  std::cerr<<"Choice Table Size =>"<<_data.size()<<std::endl;
  for (int i=0;i<_data.size();i++)
    {
      std::cerr<<i<<"->"<<_data[i].size()<<std::endl;
    }
}

void ChoiceTable::putLast(int i_node,int r_node,int elem)
{
  (_data[i_node])[r_node].push_back(elem);
}

int ChoiceTable::getFirst(int i_node,int r_node) const 
{
  return((_data[i_node])[r_node].front());
}

void ChoiceTable::putFirst(int i_node,int r_node,int choice)
{
  (_data[i_node])[r_node].push_front(choice);
}

void ChoiceTable::createList(int i_node ,int r_node)
{
  (_data[i_node])[r_node]= ChoiceList();
}

void ChoiceTable::destroyList(int i_node ,int r_node)
{
  //if ((_data[i_node])[r_node]) delete (ChoiceList*) (*_data[i_node])[r_node];
}

ChoiceList* ChoiceTable::getList(int i_node,int r_node) 
{
  //std::cout<<"i_node : "<<i_node<<" - r_node : "<<r_node<<std::endl;
  assert((i_node < _i_size)&&(r_node < _r_size));
  return(&((_data[i_node])[r_node]));
}





void ChoiceTable::getList(int input_vertex, int reference_vertex, TreeGraphPtr T1, TreeGraphPtr T2, MatchRecordList* sequence){
  TreeList(input_vertex,reference_vertex,T1,T2,*sequence);
}

int ChoiceTable::Lat(ChoiceList* L, int vertex){
  assert(L->size() > vertex); // check for wrong access in the list.
  ChoiceList::iterator begin = L->begin();
  for (int i=0 ; i<vertex ; ++i) ++begin;
  return (*begin); 
}



void ChoiceTable::TreeList(int input_vertex, int reference_vertex, TreeGraphPtr T1, TreeGraphPtr T2, MatchRecordList& sequence){
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      ChoiceList* L=getList(input_vertex,reference_vertex);
      int tree_choice=L->front();
      switch(tree_choice)
	{
	case 1:
	  {
	    TreeList(Lat(L,1),reference_vertex,T1,T2,sequence);
	  }
	  break;
	case 2:
	  {
	    TreeList(input_vertex,Lat(L,1),T1,T2,sequence);
	  }
	  break;
	case 3: 
	  {
	    sequence.push_back(MatchRecord(input_vertex,reference_vertex));
	    ForestList(input_vertex,reference_vertex,T1,T2,sequence);
	  }
	  break;
	case 4: 
	  {
	    int size = Lat(L,1);
	    int nbChoices = L->size();
	    /*
	      Fred temptative of fix:
	      The actual list of id of matched elements seems to be the n last value of the choice list.
	    */
	    int input_vertex;
	    for (int i=0;i<size;i++){
	      input_vertex = Lat(L,nbChoices-size+i);
	      sequence.push_back(MatchRecord(input_vertex,reference_vertex));
	    }
	    // question : is input_vertex to continue the fist or the last of the list. I choose the first
	    input_vertex = Lat(L,nbChoices-size); 
	    ForestList(input_vertex,reference_vertex,T1,T2,sequence);
	  }
	  break;
	case 5: 
	  {
	    int size = Lat(L,1);
	    int nbChoices = L->size();
	    /*
	      Fred temptative of fix:
	      The actual list of id of matched elements seems to be the n last value of the choice list.
	    */
	    int reference_vertex;
	    for (int i=0;i<size;i++){
	      reference_vertex = Lat(L,nbChoices-size+i);
	      sequence.push_back(MatchRecord(input_vertex,reference_vertex));
	    }
	    // question : is reference_vertex to continue the fist or the last of the list. I choose the first
	    reference_vertex = Lat(L,nbChoices-size); 
	    ForestList(input_vertex,reference_vertex,T1,T2,sequence);
	  }
	  break;
	default : break;
	}
    }
}

void ChoiceTable::ForestList(int input_vertex,int reference_vertex,TreeGraphPtr T1, TreeGraphPtr T2, MatchRecordList& sequence)
{
  ChoiceList* L=getList(input_vertex,reference_vertex);
  int forest_choice=Lat(L,2);
  
  switch(forest_choice)
    {
    case 1: ForestList(Lat(L,3),reference_vertex,T1,T2,sequence);break;
    case 2: ForestList(input_vertex,Lat(L,3),T1,T2,sequence);break;
    case 3: 
      {
	for (int i=0;i<T1->getNbChild(input_vertex);i++)
	  {
	    int i_node=T1->child(input_vertex,i);
	    int r_node=Lat(L,3+i);
	    if (r_node!=-1) TreeList(i_node,r_node,T1,T2,sequence);
	  }
      }
      break;
    default : break;
    }
}

#include <iostream>

template<class A>
void bin_write(ofstream& stream, const A& value){ stream.write((const char*)&value,sizeof(A)); }

void ChoiceTable::dump(const std::string& fname) const
{

  std::ofstream stream(fname.c_str(),std::ios::out | std::ios::binary);
  if (stream) {
    bin_write(stream,1.0f); // version
    bin_write(stream,_i_size); // sizes
    bin_write(stream,_r_size);
    for (ChoiceListArray::const_iterator itcr = _data.begin(); itcr != _data.end(); ++itcr){
      for (ChoiceListVector::const_iterator itcc = itcr->begin(); itcc != itcr->end(); ++itcc){
	// size of a record
	bin_write(stream,itcc->size());
	for (ChoiceList::const_iterator itc = itcc->begin(); itc != itcc->end(); ++itc){
	  // values of a record
	  bin_write(stream,*itc);
	}
      }
    }

  }
}

template<class A>
void bin_read(ifstream& stream, A& value){ stream.read((char*)&value,sizeof(A)); }

ChoiceTable ChoiceTable::load(const std::string& fname)
{
  std::ifstream stream(fname.c_str(),std::ios::in | std::ios::binary);
  if (stream) {
    float version;
    bin_read(stream,version); // version
    if (version != 1.0) { printf("cannot read. Wrong version nb.\n"); return ChoiceTable(0,0); }
    int i_size, r_size;
    bin_read(stream,i_size); // sizes
    bin_read(stream,r_size);
    ChoiceTable result(i_size,r_size);
    for (ChoiceListArray::iterator itcr = result._data.begin(); itcr != result._data.end(); ++itcr){
      for (ChoiceListVector::iterator itcc = itcr->begin(); itcc != itcr->end(); ++itcc){
	size_t sizerecord;
	// read size of a record
	bin_read(stream,sizerecord);
	itcc->resize(sizerecord);
	for (ChoiceList::iterator itc = itcc->begin(); itc != itcc->end(); ++itc){
	  // read values of a record
	  bin_read(stream,*itc);
	}
      }
    }
    return result;
  }
  return ChoiceTable(0,0);
}
