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


#include "relation.h"
     
Relation::Relation()
{
  _IV = 0;
  _RV = 0;
  _cost = DIST_UNDEF;
  _next = (Relation*) NULL;
}


Relation::Relation(int iv, int rv,DistanceType cost,Relation* relation)
{
  _IV=iv;
  _RV=rv;
  _cost=cost;
  _next= relation;	
}

Relation::~Relation()
{
}

int 	Relation::getIV() const 
{
  return(_IV);
}
int 	Relation::getRV() const 
{
  return(_RV);
}
void 	Relation::putIV(int vertex)
{
  _IV=vertex;
}
void 	Relation::putRV(int vertex)
{
  _RV=vertex;
}

Relation* Relation::getNext() const
{
  return(_next);
}

void Relation::putNext(Relation* new_next)
{
  _next=new_next;
}

int 	Relation::operator==(const Relation relation) 
{
  return((_IV==relation.getIV())&&(_RV==relation.getRV())&&(_cost==relation.getCost()));
}



