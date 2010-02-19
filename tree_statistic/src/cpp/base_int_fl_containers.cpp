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
 *       $Source$
 *       $Id: base_int_fl_containers.cpp 3186 2007-05-25 15:10:30Z dufourko $
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

#include "base_int_fl_containers.h"

using namespace std;
using namespace Stat_trees;

Base_Int_fl_container<1, 0>::Base_Int_fl_container()
{}

Base_Int_fl_container<1, 0>&
Base_Int_fl_container<1, 0>::operator=(int i)
{
   data= i;
   return *this;
}

Base_Int_fl_container<1, 0>&
Base_Int_fl_container<1, 0>::operator=(const Base_Int_fl_container<1, 0>& bifc)
{
   data= bifc.data;
   return *this;
}

const int& Base_Int_fl_container<1, 0>::Int(int i) const
{
  assert(i < 1);
  return (data);
}

int& Base_Int_fl_container<1, 0>::Int(int i)
{
  assert(i < 1);
  return (data);
}

const double& Base_Int_fl_container<1, 0>::Double(int i) const
{
   assert(i < 0);
   return (const double&)(data);
}

double& Base_Int_fl_container<1, 0>::Double(int i)
{
   assert(i < 0);
   return (double&)(data);
}

int Base_Int_fl_container<1, 0>::nb_int() const
{ return 1; }

int Base_Int_fl_container<1, 0>::nb_float() const
{ return 0; }

std::ostream& Stat_trees::operator<<(std::ostream& os,
                                     const Base_Int_fl_container<1, 0>& bifc)
{
   os << bifc.Int(0);
   return os;
}

Base_Int_fl_container<0, 1>::Base_Int_fl_container()
{}

Base_Int_fl_container<0, 1>& Base_Int_fl_container<0, 1>::operator=(double d)
{
   data= d;
   return *this;
}

Base_Int_fl_container<0, 1>& Base_Int_fl_container<0, 1>::operator=(const Base_Int_fl_container<0, 1>& bifc)
{
   data= bifc.data;
   return *this;
}

const double& Base_Int_fl_container<0, 1>::Double(int i) const
{
  assert(i < 1);
  return (data);
}

double& Base_Int_fl_container<0, 1>::Double(int i)
{
  assert(i < 1);
  return (data);
}

const int& Base_Int_fl_container<0, 1>::Int(int i) const
{
   assert(i < 0);
   return (const int&)(data);
}

int& Base_Int_fl_container<0, 1>::Int(int i)
{
   assert(i < 0);
   return (int&)(data);
}

int Base_Int_fl_container<0, 1>::nb_int() const
{ return 0; }

int Base_Int_fl_container<0, 1>::nb_float() const
{ return 1; }

std::ostream& Stat_trees::operator<<(std::ostream& os,
                                     const Base_Int_fl_container<0, 1>& bifc)
{
   os << bifc.Double(0);
   return os;
}

Base_Int_fl_container<0, 0>::Base_Int_fl_container()
 :  _nb_integral(0)
 , _nb_float(0)
{}

int Base_Int_fl_container<0, 0>::nb_int() const
{ return 0; }

int Base_Int_fl_container<0, 0>::nb_float() const
{ return 0; }

std::ostream& Stat_trees::operator<<(std::ostream& os,
                                     const Base_Int_fl_container<0, 0>& bifc)
{
   // int i;

   os << "()" ;

   return os;
}

string Type_Error_Exception::get_string_d2i()
{ return "Warning : used an int instead of a double"; }

string Type_Error_Exception::get_string_i2d()
{ return "Warning : used a double instead of a int"; }

IntDataSet::IntDataSet() { type= INT_VALUE; }
IntDataSet::IntDataSet(int size) : v(size) { type= INT_VALUE; }
const int& IntDataSet::Int(int i) const { return v[i]; }
int& IntDataSet::Int(int i) { return v[i]; }

int& DataSet::Int(int i)
{ throw Type_Error_Exception(); }

const int& DataSet::Int(int i) const
{ throw Type_Error_Exception(); }

DoubleDataSet::DoubleDataSet() { type= REAL_VALUE; }
DoubleDataSet::DoubleDataSet(int size) : v(size) { type= REAL_VALUE; }
const double& DoubleDataSet::Double(int i) const { return v[i]; }
double& DoubleDataSet::Double(int i) { return v[i]; }

const double& DataSet::Double(int i) const
{ throw Type_Error_Exception(); }

double& DataSet::Double(int i)
{ throw Type_Error_Exception(); }

