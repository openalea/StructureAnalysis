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
 *       $Id: indexmanager.h 3258 2007-06-06 13:18:26Z dufourko $
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


#ifndef SB_INDEX_MANAGER_HEADER
#define SB_INDEX_MANAGER_HEADER

#include <iostream>
#include "definitions.h"
#include <vector>
#include <list>

typedef std::vector<int> IndexVector;
typedef std::list<int> IndexList;

/**
 *\class IndexManager
 *\brief
 *\author Pascal ferraro
 *\date 1999
 */

class TREEMATCH_API IndexManager
{

 public :
 IndexManager(){};
 IndexManager(int ,int );
 ~IndexManager();
 void resize(int ,int );
 void open(int );
 void close(int );
 int getIndex(int) const ;
 int getRealSize() const { return(_realSize); };
 int getSimulatedSize() const { return(_simulatedSize); };
 protected :
 IndexVector _index;
 IndexList _freeIndices;
 int _realSize;
 int _simulatedSize;

};

#endif



