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
 *       $Id: base_int_fl_containers.hpp 3186 2007-05-25 15:10:30Z dufourko $
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

#ifndef BASE_INT_FL_CONTAINTERS_TCC
#define BASE_INT_FL_CONTAINTERS_TCC

// using namespace std;

template<int inb_integral, int inb_float>
Base_Int_fl_container<inb_integral, inb_float>::Base_Int_fl_container()
 : _nb_integral(inb_integral)
 , _nb_float(inb_float)
{
  pair.push_back(new IntDataSet(_nb_integral));
  pair.push_back(new DoubleDataSet(_nb_float));
}

template<int inb_integral, int inb_float>
Base_Int_fl_container<inb_integral, inb_float>&
Base_Int_fl_container<inb_integral, inb_float>::operator=(const Base_Int_fl_container<inb_integral, inb_float>& bifc)
{
   pair= bifc.pair;
   return *this;
}

template<int inb_integral, int inb_float>
const int& Base_Int_fl_container<inb_integral, inb_float>::Int(int i) const
{
  assert(i < _nb_integral);
  return (pair[0]->Int(i));
}

template<int inb_integral, int inb_float>
int& Base_Int_fl_container<inb_integral, inb_float>::Int(int i)
{
  assert(i < _nb_integral);
  return (pair[0]->Int(i));
}

template<int inb_integral, int inb_float>
const double& Base_Int_fl_container<inb_integral, inb_float>::Double(int i) const
{
  assert(i < _nb_float);
  return pair[1]->Double(i);
}

template<int inb_integral, int inb_float>
double& Base_Int_fl_container<inb_integral, inb_float>::Double(int i)
{
  assert(i < _nb_float);
  return pair[1]->Double(i);
}

template<int inb_integral, int inb_float>
int Base_Int_fl_container<inb_integral, inb_float>::nb_int() const
{ return _nb_integral; }

template<int inb_integral, int inb_float>
int Base_Int_fl_container<inb_integral, inb_float>::nb_float() const
{ return _nb_float; }

template<int inb_integral, int inb_float>
std::ostream& operator<<(std::ostream& os,
                         const  Base_Int_fl_container<inb_integral, inb_float>& bifc)
{
   int i;

   os << "(" ;
   for(i= 0; i < bifc._nb_integral; i++)
     os << bifc.Int(i) << ", ";
   for(i= 0; i < bifc._nb_float-1; i++)
     os << bifc.Double(i) << ", ";
   if (bifc._nb_float > 0)
     os << bifc.Double(bifc._nb_float-1);
   os << ")" ;

   return os;
}

template<int inb_integral>
Base_Int_fl_container<inb_integral, 0>::Base_Int_fl_container()
 : _nb_integral(inb_integral)
 , _nb_float(0)
{ data= new int[_nb_integral]; }

template<int inb_integral>
Base_Int_fl_container<inb_integral, 0>&
Base_Int_fl_container<inb_integral, 0>::operator=(const Base_Int_fl_container<inb_integral, 0>& bifc)
{
  int i;

  data= new int[_nb_integral];
  for(i= 0; i < _nb_integral; i++)
    data[i]= bifc.data[i];
  return *this;
}

template<int inb_integral>
const int& Base_Int_fl_container<inb_integral, 0>::Int(int i) const
{
  assert(i < _nb_integral);
  return (data[i]);
}

template<int inb_integral>
int& Base_Int_fl_container<inb_integral, 0>::Int(int i)
{
  assert(i < _nb_integral);
  return (data[i]);
}

template<int inb_integral>
const double& Base_Int_fl_container<inb_integral, 0>::Double(int i) const
{
  assert(i < _nb_float);
#ifdef __GNUC__
  #warning no return statement
#endif
#ifndef __clang__
  return 0;
#endif
}

template<int inb_integral>
double& Base_Int_fl_container<inb_integral, 0>::Double(int i)
{
  assert(i < _nb_float);
#ifdef __GNUC__
  #warning no return statement
#endif
#ifndef __clang__
  return 0;
#endif
}

template<int inb_integral>
int Base_Int_fl_container<inb_integral, 0>::nb_int() const
{ return _nb_integral; }

template<int inb_integral>
int Base_Int_fl_container<inb_integral, 0>::nb_float() const
{ return 0; }

template<int inb_integral>
std::ostream& operator<<(std::ostream& os,
                         const  Base_Int_fl_container<inb_integral, 0>& bifc)
{
   int i;

   os << "(" ;
   for(i= 0; i < bifc._nb_integral-1; i++)
        os << bifc.Int(i) << ", ";
   os << bifc.Int(bifc._nb_integral-1);
   os << ")" ;

   return os;
}

template<int inb_float>
Base_Int_fl_container<0, inb_float>::Base_Int_fl_container()
 : _nb_integral(0)
 , _nb_float(inb_float)
{ data= new double[_nb_float]; }

template<int inb_float>
Base_Int_fl_container<0, inb_float>&
Base_Int_fl_container<0, inb_float>::operator=(const Base_Int_fl_container<0, inb_float>& bifc)
{
  int i;

  data= new double[_nb_float];
  for(i= 0; i < _nb_float; i++)
    data[i]= bifc.data[i];
  return *this;
}

template<int inb_float>
const double& Base_Int_fl_container<0, inb_float>::Double(int i) const
{
  assert(i < _nb_float);
  return (data[i]);
}

template<int inb_float>
double& Base_Int_fl_container<0, inb_float>::Double(int i)
{
  assert(i < _nb_float);
  return (data[i]);
  // #warning no return statement
  // return 0;
}

template<int inb_float>
const int& Base_Int_fl_container<0, inb_float>::Int(int i) const
{
  assert(i < 0);
#ifdef __GNUC__
  #warning no return statement
#endif
#ifndef __clang__
  return 0;
#endif
}

template<int inb_float>
int& Base_Int_fl_container<0, inb_float>::Int(int i)
{
  assert(i < 0);
#ifdef __GNUC__
  #warning no return statement
#endif  
#ifndef __clang__
  return 0;
#endif
}

template<int inb_float>
int Base_Int_fl_container<0, inb_float>::nb_int() const
{ return 0; }

template<int inb_float>
int Base_Int_fl_container<0, inb_float>::nb_float() const
{ return _nb_float; }

template<int inb_float>
std::ostream& operator<<(std::ostream& os,
                         const  Base_Int_fl_container<0, inb_float>& bifc)
{
   int i;

   os << "(" ;
   for(i= 0; i < bifc._nb_float-1; i++)
        os << bifc.Double(i) << ", ";
   os << bifc.Double(bifc._nb_float-1);
   os << ")" ;

   return os;
}

#endif
