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


#include "emchoicetable.h"

ExtendedMatchingChoiceTable::ExtendedMatchingChoiceTable(int i_size,int r_size)
{
  _i_size=i_size;
  _r_size=r_size;
  _matchingChoice.resize(_i_size);
  _ipfrpfMatchingChoice.resize(_i_size);
  _imfrmfMatchingChoice.resize(_i_size);
  _imfrfMatchingChoice.resize(_i_size);
  _ifrmfMatchingChoice.resize(_i_size);
  for (int i=0;i<_i_size;i++)
    {
      _matchingChoice[i]       = new ChoiceListVector(_r_size);
      _ipfrpfMatchingChoice[i] = new ChoiceListVector(_r_size);
      _imfrmfMatchingChoice[i] = new ChoiceListVector(_r_size);
      _imfrfMatchingChoice[i]  = new ChoiceListVector(_r_size);
      _ifrmfMatchingChoice[i]  = new ChoiceListVector(_r_size);
    }
}

ExtendedMatchingChoiceTable::~ExtendedMatchingChoiceTable()
{
  for (int i=0;i<_i_size;i++)
    {
      for (int j=0;j<_r_size;j++)
	{
	  if ((*_matchingChoice[i])[j])       delete (ChoiceList*) (*_matchingChoice[i])[j];
	  if ((*_ipfrpfMatchingChoice[i])[j]) delete (ChoiceList*) (*_ipfrpfMatchingChoice[i])[j];
	  if ((*_imfrmfMatchingChoice[i])[j]) delete (ChoiceList*) (*_imfrmfMatchingChoice[i])[j];
	  if ((*_imfrfMatchingChoice[i])[j])  delete (ChoiceList*) (*_imfrfMatchingChoice[i])[j];
	  if ((*_ifrmfMatchingChoice[i])[j])  delete (ChoiceList*) (*_ifrmfMatchingChoice[i])[j];
	}
      delete (ChoiceListVector*) _matchingChoice[i];
      delete (ChoiceListVector*) _ipfrpfMatchingChoice[i];
      delete (ChoiceListVector*) _imfrmfMatchingChoice[i];
      delete (ChoiceListVector*) _imfrfMatchingChoice[i];
      delete (ChoiceListVector*) _ifrmfMatchingChoice[i];
    }
}

void ExtendedMatchingChoiceTable::resize(int i_size,int r_size)
{
  _i_size=i_size;
  _r_size=r_size;
  _matchingChoice.resize(_i_size);
  _ipfrpfMatchingChoice.resize(_i_size);
  _imfrmfMatchingChoice.resize(_i_size);
  _imfrfMatchingChoice.resize(_i_size);
  _ifrmfMatchingChoice.resize(_i_size);
  for (int i=0;i<_i_size;i++)
    {
      _matchingChoice[i]       = new ChoiceListVector(_r_size);
      _ipfrpfMatchingChoice[i] = new ChoiceListVector(_r_size);
      _imfrmfMatchingChoice[i] = new ChoiceListVector(_r_size);
      _imfrfMatchingChoice[i]  = new ChoiceListVector(_r_size);
      _ifrmfMatchingChoice[i]  = new ChoiceListVector(_r_size);
    }	
}

void ExtendedMatchingChoiceTable::putLast(int i_node,int r_node,int elem)
{
  (*_matchingChoice[i_node])[r_node]->push_back(elem);
}

int ExtendedMatchingChoiceTable::getFirst(int i_node,int r_node) const 
{
  return((*_matchingChoice[i_node])[r_node]->front());
}

void ExtendedMatchingChoiceTable::putFirst(int i_node,int r_node,int choice)
{
  ChoiceList::iterator first = ((*_matchingChoice[i_node])[r_node])->begin();
  (*_matchingChoice[i_node])[r_node]->insert(first,choice);
  //(*_matchingChoice[i_node])[r_node]->prepend(choice);
}

void ExtendedMatchingChoiceTable::putIMFRMFLast(int i_node,int r_node,int elem)
{
  (*_imfrmfMatchingChoice[i_node])[r_node]->push_back(elem);
}

void ExtendedMatchingChoiceTable::putIPFRPFLast(int i_node,int r_node,int elem)
{
  (*_ipfrpfMatchingChoice[i_node])[r_node]->push_back(elem);
}

void ExtendedMatchingChoiceTable::putIMFRFLast(int i_node,int r_node,int elem)
{
  (*_imfrfMatchingChoice[i_node])[r_node]->push_back(elem);
}

void ExtendedMatchingChoiceTable::putIFRMFLast(int i_node,int r_node,int elem)
{
  (*_ifrmfMatchingChoice[i_node])[r_node]->push_back(elem);
}

void ExtendedMatchingChoiceTable::createList(int i_node ,int r_node)
{
  (*_matchingChoice[i_node])[r_node]       = new ChoiceList();
  (*_imfrmfMatchingChoice[i_node])[r_node] = new ChoiceList();
  (*_ipfrpfMatchingChoice[i_node])[r_node] = new ChoiceList();
  (*_imfrfMatchingChoice[i_node])[r_node]  = new ChoiceList();
  (*_ifrmfMatchingChoice[i_node])[r_node]  = new ChoiceList();
}

void ExtendedMatchingChoiceTable::destoyList(int i_node ,int r_node)
{
  if ((*_matchingChoice[i_node])[r_node])       delete (ChoiceList*) (*_matchingChoice[i_node])[r_node];
  if ((*_ipfrpfMatchingChoice[i_node])[r_node]) delete (ChoiceList*) (*_ipfrpfMatchingChoice[i_node])[r_node];
  if ((*_imfrmfMatchingChoice[i_node])[r_node]) delete (ChoiceList*) (*_imfrmfMatchingChoice[i_node])[r_node];
  if ((*_imfrfMatchingChoice[i_node])[r_node])  delete (ChoiceList*) (*_imfrfMatchingChoice[i_node])[r_node];
  if ((*_ifrmfMatchingChoice[i_node])[r_node])  delete (ChoiceList*) (*_ifrmfMatchingChoice[i_node])[r_node];
}

ChoiceList* ExtendedMatchingChoiceTable::getList(int i_node,int r_node) const
{
  return((*_matchingChoice[i_node])[r_node]);
}

ChoiceList* ExtendedMatchingChoiceTable::getIMFRMFList(int i_node,int r_node) const
{
  return((*_imfrmfMatchingChoice[i_node])[r_node]);
}

ChoiceList* ExtendedMatchingChoiceTable::getIPFRPFList(int i_node,int r_node) const
{
  return((*_ipfrpfMatchingChoice[i_node])[r_node]);
}

ChoiceList* ExtendedMatchingChoiceTable::getIMFRFList(int i_node,int r_node) const
{
  return((*_imfrfMatchingChoice[i_node])[r_node]);
}

ChoiceList* ExtendedMatchingChoiceTable::getIFRMFList(int i_node,int r_node) const
{
  return((*_ifrmfMatchingChoice[i_node])[r_node]);
}

