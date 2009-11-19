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
 *       $Id: distancetable.cpp 3258 2007-06-06 13:18:26Z dufourko $
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


#include "distancetable.h"

// ----------------
// Constructeur
// ----------------

DistanceTable::DistanceTable(int r_size,int c_size,int s_size)
{
  _rowSize    = r_size;
  _columnSize = c_size;
  _simulatedSize  = s_size;
  _indices.resize(_rowSize , _simulatedSize);
  _distances.resize(_rowSize);
  for (int i=0;i<_rowSize;i++) 
  { 
    _distances[i].resize(_columnSize); 
  }
}

void DistanceTable::resize(int r_size,int c_size,int s_size)
{
  _rowSize    = r_size;
  _columnSize = c_size;
  _simulatedSize  = s_size;
  _indices.resize(_rowSize , _simulatedSize);
  _distances.resize(_rowSize);
  for (int i=0;i<_rowSize;i++) 
  { 
    _distances[i].resize(_columnSize); 
  }
}

DistanceTable::~DistanceTable()
{
}

// renvoie la valeur de la colonne reference_vertex et de la ligne
// correspondant a l'index du sommet initial, gere par le manager 
// d'index dans le tableau des distances.
DistanceType DistanceTable::getDistance(int input_vertex ,int reference_vertex ) const 
{
  assert((input_vertex>=0));
  assert((input_vertex<_simulatedSize));
  int index = _indices.getIndex(input_vertex);
  return(_distances[index][reference_vertex]);
}

// insere la valeur de la distance entre les deux sommets dans la
// colonne reference_vertex et de la ligne correspondant a l'index
// du sommet initial, gere par le manager d'index du tableau de distance.
void DistanceTable::putDistance(DistanceType distance ,int input_vertex ,int reference_vertex )
{
  assert((input_vertex>=0));
  assert((input_vertex<_simulatedSize));
  int index = _indices.getIndex(input_vertex);
  _distances[index][reference_vertex]=distance;  
}


void DistanceTable::openDistanceVector(int input_vertex)
{
  assert((input_vertex>=0));
  assert((input_vertex<_simulatedSize));
  //cout<<"distancetable s_size "<<_simulatedSize<<endl;
  _indices.open(input_vertex); 
}

void DistanceTable::closeDistanceVector(int input_vertex)
{
  assert((input_vertex>=0));
  assert((input_vertex<_simulatedSize));
  _indices.close(input_vertex);
}













