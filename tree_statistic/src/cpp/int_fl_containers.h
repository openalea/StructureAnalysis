/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *       many thanks to: P. Barbier de Reuille
 *
 *       $Source
 *       $Id: int_fl_containers.h 3186 2007-05-25 15:10:30Z dufourko $
 *
 *       Forum for V-Plants developers:
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
 *       MERCHANTABILITY or FITNESS for A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */

#ifndef INT_FL_CONTAINTERS_H
#define INT_FL_CONTAINTERS_H

/****************************************************************
 *
 *  Purpose :
 *     provides containers for the storage of multidimensional
 *  mixed data (i.e. vectors of integers/floats)
 */

#include <iostream>
#include <vector>
#include <assert.h>
#include <cstdlib>

namespace Stat_trees {

/***************************************************************/
// TOOLS_USING_NAMESPACE
/***************************************************************/

// class foo;

/****************************************************************
 *
 *  Type definitions
 */

// typedef bar;

/****************************************************************
 *
 *  Class definitions :
 */

class Int_fl_container
{ // a vector of int and double

   friend std::ostream& operator<<(std::ostream& os,
                                   const  Int_fl_container& ifc);

protected :

   int _nb_integral;
   int _nb_float;
   typedef union
   {
      int  i;
      double f;
   } int_or_float;
   std::vector<int_or_float> data;

public :

   Int_fl_container();
   Int_fl_container(int inb_integral, int inb_float);
   Int_fl_container(const Int_fl_container& ifc);
   ~Int_fl_container();
   Int_fl_container& operator=(const Int_fl_container& ifc);

   void reset(int inb_integral= 1, int inb_float= 1);

   const int& Int(int index) const;
   int& Int(int index);

   const double& Double(int index) const;
   double& Double(int index);

   int nb_int() const;
   int nb_float() const;
};

class One_int_container : public Int_fl_container
{
   friend std::ostream& operator<<(std::ostream& os,
                                   const One_int_container& ifc);

public :

   One_int_container();
   One_int_container(const One_int_container& ifc);
   ~One_int_container();

   One_int_container& operator=(const One_int_container& ifc);

   const int& Int() const;
   int& Int();

   int nb_int() const;
   int nb_float() const;
};

}; // end namespace
#endif
