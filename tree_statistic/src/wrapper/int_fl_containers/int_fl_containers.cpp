
// Includes ====================================================================
#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "tree_statistic/int_fl_containers.h"

#include <boost/python.hpp>
#include <boost/python/make_constructor.hpp>
#include <sstream>

// Using =======================================================================
using namespace boost::python;
using namespace Stat_trees;

// Declarations ================================================================
namespace  {

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Int_fl_container_reset_overloads_0_2, reset, 0, 2)

void SetInt_wrapper(Int_fl_container& i,  int index, int v)
{ i.Int(index)= v; }

void SetDouble_wrapper(Int_fl_container& i,  int index, double v)
{ i.Double(index)= v; }

Int_fl_container* wrapper_init1(object tree_value)
{
   Int_fl_container *res;
   boost::python::list values, types;
   int nb_integral, nb_float, index, type,
       int_index, float_index; // length,
   object current_object;
   ostringstream error_message;
   // bool status= true;

   values= extract<boost::python::list>(tree_value.attr("Values")());
   types= extract<boost::python::list>(tree_value.attr("Types")());
   nb_integral= extract<int>(tree_value.attr("NbInt")());
   nb_float= extract<int>(tree_value.attr("NbFloat")());
   int_index= 0;
   float_index= 0;

   res= new Int_fl_container(nb_integral, nb_float);
   for(index=0; index < nb_integral+nb_float; index++)
   {
      type= extract<int>(types[index]);
      if (type==REAL_VALUE)
         res->Double(float_index++)= extract<double>(values[index]);
      else
      {
         if ((type==INT_VALUE) || (type==STATE) || (type==NB_INTERNODE))
            res->Int(int_index++)= extract<int>(values[index]);
         else
         {
            PyErr_SetString(PyExc_IndexError, "unknown variable type");
            throw_error_already_set();
         }
      }
   }
   return res;
}

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(int_fl_containers)
{
    class_< Int_fl_container >("Int_fl_container", init< const Int_fl_container& >())
        .def(init< >())
        .def(init< int, int >())
        .def("__init__", make_constructor(wrapper_init1))
        .def("Reset", &Int_fl_container::reset, Int_fl_container_reset_overloads_0_2())
        .def("Int", (const int& (Int_fl_container::*)(int) const)&Int_fl_container::Int, return_value_policy< copy_const_reference >())
        .def("Double", (const double& (Int_fl_container::*)(int) const)&Int_fl_container::Double, return_value_policy< copy_const_reference >())
        .def("SetInt", &SetInt_wrapper)
        .def("SetDouble", &SetDouble_wrapper)
        .def("NbInt", &Int_fl_container::nb_int)
        .def("NbFloat", &Int_fl_container::nb_float)
        .def(self_ns::str(self))
    ;

}

