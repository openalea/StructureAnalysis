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
      _data[i].resize(_r_size);
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



