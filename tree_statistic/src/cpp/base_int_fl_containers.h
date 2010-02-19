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
 *       $Id: base_int_fl_containers.h 3193 2007-05-29 10:03:19Z dufourko $
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

#ifndef BASE_INT_FL_CONTAINTERS_H
#define BASE_INT_FL_CONTAINTERS_H

/****************************************************************
 *
 *  Purpose :
 *     provides containers for the storage of multidimensional
 *  mixed data (i.e. vectors of integers/floats)
 */

#ifdef __MINGW32__
#define DO_NOT_DLLIMPORT_RCOBJECT
#endif

#include "tool/rcobject.h"

// Useful for RefCountObject (<=> DECLARE_PTR(DataSet))
#include <iostream>
#include <vector>
#include <assert.h>
#include <cstdlib>
#include <string>

/***************************************************************/
VPTOOLS_USING_NAMESPACE
/***************************************************************/

namespace Stat_trees
{

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

class Type_Error_Exception
{ // exceptions raised when using floats instead of integers
  // within DataSet objects and vice-versa

public:

  std::string get_string_d2i();
  std::string get_string_i2d();

// use : catch (Type_Error_Exception e)
//             cout e.getstring_d2i() << endl;
};

class DataSet : public RefCountObject
{ // a container of int or double (or something else)

public :

  enum Int_or_fl { INT_VALUE, REAL_VALUE };

  Int_or_fl type;

  virtual const int& Int(int i) const;
  virtual int& Int(int i);
  virtual const double& Double(int i) const;
  virtual double& Double(int i);

  // virtual someting_else Something_else(int i);
};

class IntDataSet : public DataSet
{ // a container of int

protected :

  std::vector<int> v;

public :

  IntDataSet();
  IntDataSet(int size);
  virtual const int& Int(int i) const;
  virtual int& Int(int i);
};

class DoubleDataSet : public DataSet
{ // a container of double

protected :

  std::vector<double> v;

public :

  DoubleDataSet();
  DoubleDataSet(int size);
  virtual const double& Double(int i) const;
  virtual double& Double(int i);
};

DECLARE_PTR(DataSet);

template<int inb_integral, int inb_float>
class Base_Int_fl_container
{ // a pair of vectors of int and double

   template<int i, int f>
   friend std::ostream& operator<<(std::ostream& os,
                                   const Base_Int_fl_container<i, f>& bifc);

protected :

   int _nb_integral;
   int _nb_float;
   std::vector<DataSetPtr> pair;

public :

   Base_Int_fl_container();
   Base_Int_fl_container<inb_integral, inb_float>&
   operator=(const Base_Int_fl_container<inb_integral, inb_float>& bifc);

   const int& Int(int i) const;
   int& Int(int i);

   const double& Double(int i) const;
   double& Double(int i);

   int nb_int() const;
   int nb_float() const;

};

template<int inb_integral>
class Base_Int_fl_container<inb_integral, 0>
{ // a vector of int

   template<int i>
   friend std::ostream& operator<<(std::ostream& os,
                                   const  Base_Int_fl_container<i, 0>& bifc);
protected :

  int _nb_integral;
  int _nb_float;
  int *data;

public :

  Base_Int_fl_container();
  Base_Int_fl_container<inb_integral, 0>&
  operator=(const Base_Int_fl_container<inb_integral, 0>& bifc);

  const int& Int(int i) const;
  int& Int(int i);

  const double& Double(int i) const;
  double& Double(int i);

  int nb_int() const;
  int nb_float() const;

};

template<>
class Base_Int_fl_container<1, 0>
{ // one int

  friend std::ostream& operator<<(std::ostream& os,
                                  const  Base_Int_fl_container<1, 0>& bifc);

protected :

  int data;

public :

  Base_Int_fl_container();
  Base_Int_fl_container<1, 0>& operator=(int i);
  Base_Int_fl_container<1, 0>& operator=(const Base_Int_fl_container<1, 0>& bifc);

  // int& operator=(const Base_Int_fl_container<1, 0>&)
  // { return (data); }

  const int& Int(int i) const;
  int& Int(int i);

  const double& Double(int i) const;
  double& Double(int i);

  int nb_int() const;
  int nb_float() const;

};


template<int inb_float>
class Base_Int_fl_container<0, inb_float>
{ // a vector of double

   template<int f>
   friend std::ostream& operator<<(std::ostream& os,
                                   const  Base_Int_fl_container<0, f>& bifc);

protected :

   int _nb_integral;
   int _nb_float;
   double *data;

public :

   Base_Int_fl_container();
   Base_Int_fl_container<0, inb_float>&
   operator=(const Base_Int_fl_container<0, inb_float>& bifc);

   const double& Double(int i) const;
   double& Double(int i);

   const int& Int(int i) const;
   int& Int(int i);

   int nb_int() const;
   int nb_float() const;

};

template<>
class Base_Int_fl_container<0, 1>
{ // one double

   friend std::ostream& operator<<(std::ostream& os,
                                   const  Base_Int_fl_container<0, 1>& bifc);

protected :

   double data;

public :

   Base_Int_fl_container();
   Base_Int_fl_container<0, 1>& operator=(double d);
   Base_Int_fl_container<0, 1>& operator=(const Base_Int_fl_container<0, 1>& bifc);

   const double& Double(int i) const;
   double& Double(int i);

   const int& Int(int i) const;
   int& Int(int i);

   int nb_int() const;
   int nb_float() const;

};

template<>
class Base_Int_fl_container<0, 0>
{ // an empty container

   friend std::ostream& operator<<(std::ostream& os,
                                   const  Base_Int_fl_container<0, 0>& bifc);

protected :

   int _nb_integral;
   int _nb_float;

public :

   Base_Int_fl_container();

   int nb_int() const;
   int nb_float() const;

};

#include "base_int_fl_containers.hpp"
}; // end namespace

#endif
