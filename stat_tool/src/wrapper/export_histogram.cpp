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

#include "export_histogram.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"

#include <boost/python.hpp>
using namespace boost::python;


// Boost.Python Wrapper export function
void class_histogram()
{

  // _Histogram
  class_<Histogram>("_Histogram", init<const Histogram&>());

  // _Distribution_data
  class_<Distribution_data, bases<Histogram> >("_Distribution_data", init< const Histogram& >())
    .def(init< optional< int > >())
    .def(init< const Distribution& >())
    .def(init< const Distribution_data& >())
    .def(init< int, int* >())
    .def(init< int, const Histogram** >())
    .def(init< const Histogram&, char, int >());


}

