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
 *       $Id: int_fl_containers.cpp 3186 2007-05-25 15:10:30Z dufourko $
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

#include "int_fl_containers.h"
using namespace Stat_trees;

Int_fl_container::Int_fl_container()
 : data()
{
   _nb_integral= 0;
   _nb_float= 0;
}

Int_fl_container::Int_fl_container(int inb_integral, int inb_float)
 : _nb_integral(inb_integral)
 , _nb_float(inb_float)

{
   int index;
   int_or_float val;

   for(index= 0; index < _nb_integral; index++)
   {
       val.i= 0;
       data.push_back(val);
   }

   for(index= _nb_integral; index < _nb_integral+_nb_float; index++)
   {
       val.f= 0;
       data.push_back(val);
   }
}

Int_fl_container::Int_fl_container(const Int_fl_container& ifc)
 : _nb_integral(ifc._nb_integral)
 , _nb_float(ifc._nb_float)
{ data= ifc.data; }

Int_fl_container::~Int_fl_container()
{}

Int_fl_container& Int_fl_container::operator=(const Int_fl_container& ifc)
{
   int index;

   if ((_nb_integral != ifc._nb_integral) ||
      (_nb_float != ifc._nb_float))
      reset(ifc._nb_integral, ifc._nb_float);

   for(index= 0; index < _nb_integral; index++)
      Int(index)= ifc.Int(index);
   for(index= 0; index < _nb_float; index++)
      Double(index)= ifc.Double(index);
   return *this;
}

void Int_fl_container::reset(int inb_integral, int inb_float)
{
   int index;
   int_or_float val;

   for(index= 0; index < _nb_integral+_nb_float; index++)
       data.pop_back();

   assert(data.size()==0);

   _nb_integral= inb_integral;
   _nb_float= inb_float;

   for(index= 0; index < _nb_integral; index++)
   {
       val.i= 0;
       data.push_back(val);
   }

   for(index= _nb_integral; index < _nb_integral+_nb_float; index++)
   {
       val.f= 0;
       data.push_back(val);
   }

   assert(data.size() == (unsigned int)(_nb_integral+_nb_float));
}

const int& Int_fl_container::Int(int index) const
{
   assert(index < _nb_integral);
   return (data[index].i);
}

int& Int_fl_container::Int(int index)
{
   assert(index < _nb_integral);
   return (data[index].i);
}

const double& Int_fl_container::Double(int index) const
{
   assert(index < _nb_float);
   return (data[index+_nb_integral].f);
}

double& Int_fl_container::Double(int index)
{
   assert(index < _nb_float);
   return (data[index+_nb_integral].f);
}

int Int_fl_container::nb_int() const
{ return _nb_integral; }

int Int_fl_container::nb_float() const
{ return _nb_float; }

std::ostream& Stat_trees::operator<<(std::ostream& os,
                                     const Int_fl_container& ifc)
{
   int index;

   if ((ifc._nb_integral+ifc._nb_float)!= 1)
       os << "(" ;
   for(index= 0; index < ifc._nb_integral-1; index++)
   {
       os << ifc.Int(index);
       if ((ifc._nb_integral+ifc._nb_float)!= 1);
          os<< ", ";
   }
   if (ifc._nb_integral > 0)
   {
       os << ifc.Int(ifc._nb_integral-1);
       if (ifc._nb_float > 0)
          os << ", ";
   }

   for(index= 0; index < ifc._nb_float-1; index++)
   {
       os << ifc.Double(index);
       if ((ifc._nb_integral+ifc._nb_float)!= 1);
          os << ", ";
   }
   if (ifc._nb_float > 0)
     os << ifc.Double(ifc._nb_float-1);
   if ((ifc._nb_integral+ifc._nb_float)!= 1)
     os << ")" ;

   return os;
}

// One_int_container::

One_int_container::One_int_container()
 : Int_fl_container(1, 0)
{}

One_int_container::One_int_container(const One_int_container& ifc)
{ data= ifc.data; }

One_int_container::~One_int_container()
{}

One_int_container& One_int_container::operator=(const One_int_container& ifc)
{
   Int()= ifc.Int();
   return *this;
}

const int& One_int_container::Int() const
{ return (data[0].i); }

int& One_int_container::Int()
{ return (data[0].i); }

int One_int_container::nb_int() const { return 1; }
int One_int_container::nb_float() const { return 0; }

std::ostream& Stat_trees::operator<<(std::ostream& os,
                                     const One_int_container& ifc)
{
   os << ifc.Int();
   return os;
}
