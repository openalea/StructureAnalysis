/*------------------------------------------------------------------------------
 *                                                                              
 *        VPlants.Stat_Tool : VPlants Statistics module
 *                                                                              
 *        Copyright 2006-2007 INRIA - CIRAD - INRA                      
 *                                                                              
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>         
 *                                                                              
 *        Distributed under the GPL 2.0 License.                               
 *        See accompanying file LICENSE.txt or copy at                          
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *                                                                              
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr                    
 *       
 *        $Id: scenecontainer_wrap.cpp 559 2007-05-25 12:25:30Z dufourko $
 *                                                                       
 *-----------------------------------------------------------------------------*/

#include "export_distribution.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"


#include <boost/python.hpp>
// definition of boost::python::len
#include <boost/python/detail/api_placeholder.hpp>
// definition of boost::python::make_constructor
#include <boost/python/make_constructor.hpp>

using namespace boost::python;


// Wrappers

// Histogram constructor 
Distribution_data* distribution_data_wrapper_init(list int_list)
{
   int *pelement= NULL;
   object o;
   ostringstream error_message;
   bool status= true;
   Distribution_data *histo= NULL;
   // string s;

   int nb_element= boost::python::len(int_list);

   if (not nb_element)
   {
     status= false;
     error_message << "At least one observation is required to initialize Histogram"; // << endl;
     PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
     throw_error_already_set();
   }

   pelement= new int[nb_element];
   for (int i= 0; i < nb_element; i++)
     {
       try
         {
	   extract<int> x(int_list[i]);
	   if (x.check())
	     pelement[i]= x();
	   else
	     status=false;
         }
       catch (...)
         {
	   status= false;
         }
       if (!status)
	 error_message << "Incorrect type for element " << i
		       << " of argument list: expecting an INT" << endl;
     }
   if (!status)
     {
       delete [] pelement;
       PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
       throw_error_already_set();
     }

   // create new distribution_data object
   histo= new Distribution_data(nb_element, pelement);
   delete [] pelement;
   if (not histo)
     {
       error_message << "Could not initialize Histogram from argument"; // << endl;
       PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
       throw_error_already_set();
     }
   
   return histo;
}



// Overloads





// Boost.Python Wrapper export function
void class_distribution()
{

  
  // _Parametric Model
  class_< Parametric_model >("_Parametric_model",  init< optional< int, int, int, int, double, double > >())
    .def(init< int, int, int, double, double, optional< double > >())
    .def(init< const Histogram& >())
    .def(init< const Distribution& >())
    .def(init< const Parametric& >())
    //.def(init< const Distribution&, const Histogram& >())
    //.def(init< const Parametric&,  const Histogram& >())
    .def(init< const Parametric_model&, optional< bool > > ())
    .def("extract_data", &Parametric_model::extract_data, return_value_policy< manage_new_object >())
    .def(self_ns::str(self))
    ;


  // _Distribution_data
  class_<Distribution_data, bases<Histogram> >("_Distribution_data", init< const Histogram& >())
    .def(init< optional< int > >())
    .def(init< const Distribution& >())
    .def(init< const Distribution_data& >())
    .def(init< const Histogram&, char, int >())
    .def("__init__", make_constructor(distribution_data_wrapper_init))
    ;


  //standalone functions
  def("_histogram_ascii_read", histogram_ascii_read, return_value_policy< manage_new_object >() );
}
