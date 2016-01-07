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
 *       $Id: distancetable.h 3258 2007-06-06 13:18:26Z dufourko $
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


// -------------------------------
// Classe: DistanceTable
// Cette classe gere un tableau de
// de distance de r ligne et c
// colonnes
// -------------------------------

#ifndef SB_DISTANCE_TABLE_HEADER
#define SB_DISTANCE_TABLE_HEADER

#include "indexmanager.h"

typedef std::vector<DistanceType> DistanceVector;
typedef std::vector<DistanceVector> DistanceVectorTable;


class TREEMATCH_API DistanceTable
{

 public:
  DistanceTable(){};
  DistanceTable(int ,int ,int );
  void resize(int ,int ,int );
  ~DistanceTable();
  DistanceType getDistance(int , int ) const ;
  void putDistance(DistanceType , int , int );
  void openDistanceVector(int );
  void closeDistanceVector(int );
  int getColumnSize() const { return(_columnSize); };
  int getRowSize() const { return(_rowSize); };
  int getSimulatedSize() const { return(_simulatedSize); };

 private:
  IndexManager _indices;
  DistanceVectorTable _distances;
  int _columnSize;
  int _rowSize;
  int _simulatedSize;
};

#endif



