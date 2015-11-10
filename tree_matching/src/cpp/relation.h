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



#ifndef SB_RELATION_HEADER
#define SB_RELATION_HEADER

#include <iostream>
#include"definitions.h"

/**
 *\class Relation
 *\brief Retains the mapped vertices
 *\author Pascal ferraro
 *\date 1999
 */

class Relation 
{
  friend class Sequence;
  public :
    Relation();
  Relation(int ,int,DistanceType,Relation* );
 
    ~Relation();
    Relation* getNext() const;
    void putNext(Relation* ) ;
    int getIV() const ;
    int getRV() const ;
    DistanceType getCost() const { return(_cost); }
    void putIV(int );
    void putRV(int);
    void putCost(DistanceType cost) { _cost=cost;}
    int operator==(const Relation relation) ;
    
		
  private :
    int _IV;
    int _RV;
    DistanceType _cost;
    Relation* _next;
};

#endif

