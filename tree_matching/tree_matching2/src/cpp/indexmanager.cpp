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
 *       $Id: indexmanager.cpp 3258 2007-06-06 13:18:26Z dufourko $
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


#include "indexmanager.h"

const int NO_INDEX = -1;

IndexManager::IndexManager(int r_size,int s_size)
{
  resize(r_size,s_size);
}

void IndexManager::resize(int r_size,int s_size)
{
  _realSize      = r_size;
  _simulatedSize = s_size;
  _index.resize(_simulatedSize);
  for (int i = 0 ; i<_simulatedSize; i++)
    _index[i] =NO_INDEX;
  for (int j=0;j<_realSize;j++) 
  { 
    _freeIndices.push_back(j);
  }
}

IndexManager::~IndexManager()
{

}

void IndexManager::open(int s_index)  
{
  assert((s_index>=0)&&(s_index<_simulatedSize));
   assert(_freeIndices.size());
  _index[s_index]=_freeIndices.front();
  _freeIndices.pop_front();
}

void IndexManager::close(int s_index) 
{
  assert((s_index>=0)&&(s_index<_simulatedSize));
  assert(_index[s_index] != NO_INDEX);
  _freeIndices.push_back(_index[s_index]); 
}

int IndexManager::getIndex(int s_index) const 
{ 
  assert((s_index>=0)&&(s_index<_simulatedSize)); 
  return(_index[s_index]); 
}

