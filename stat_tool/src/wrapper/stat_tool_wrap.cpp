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



/* WRAPPER Boost.python for stat_tool class */

#include "export_base.h"
#include "export_distribution.h"
#include "export_mixture.h"

#include <boost/python.hpp>
using namespace boost::python;


// Define python module "_stat_tool"
BOOST_PYTHON_MODULE(_stat_tool)
{
  class_base();
  class_distribution();
  class_mixture();
}
