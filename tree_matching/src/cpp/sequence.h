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



#ifndef SB_SEQUENCE_HEADER
#define SB_SEQUENCE_HEADER

#include "relation.h"

/**
 *\class Sequence
 *\brief Matching List between two trees referenced by their roots.
 *\author Pascal ferraro
 *\date 1999
 */

class Sequence 
{
  public :
    Sequence();
    ~Sequence();
    void append(int ,int ,DistanceType ) ;
 

    void append(Relation* ) ;
    void clear();
    void link(Sequence* );
    void add(Sequence* );
    Relation* getFirst() const ;
    Relation* getLast() const ;
    int getSize() const ;
    void operator=(const Sequence sequence);
    void reset();
    int next();
    void print();
    Relation* getCurrent() const ;
    int isUsed() const { return(_used);} ;
    void keep() { _used=1; };
    void free() { _used=0; };
    
    int getNbIns() const { return(_nbIns); };
    void putNbIns(int nb_ins) { _nbIns = nb_ins; };
    int getNbDel() const { return(_nbDel); };
    void putNbDel(int nb_del) { _nbDel = nb_del; };
    int getNbMat() const { return(_nbMat); };
    void putNbMat(int nb_mat) { _nbMat = nb_mat; };
    int getNbSub() const { return(_nbSub); };
    void putNbSub(int nb_sub) { _nbSub = nb_sub; };


    DistanceType getInsCost() const { return(_insertCost); };
    void putInsCost(DistanceType insertCost) { _insertCost = insertCost; };
    DistanceType getDelCost() const { return(_deleteCost); };
    void putDelCost(DistanceType deleteCost) { _deleteCost = deleteCost; };
    DistanceType getSubCost() const { return(_substitutionCost); };
    void putSubCost(DistanceType substitutionCost) { _substitutionCost = substitutionCost; };
    
    
    private :
    int _size;
    Relation* _first;
    Relation* _last;
    Relation* _current;
    int _used;
    int _nbIns;
    int _nbDel;
    int _nbMat;
    int _nbSub;
    DistanceType _insertCost;
    DistanceType _deleteCost;
    DistanceType _substitutionCost;
};

#endif

