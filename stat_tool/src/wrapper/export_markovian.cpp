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
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_base.cpp 5432 2008-08-25 16:51:24Z jbdurand $
 *
 *-----------------------------------------------------------------------------*/



#include "export_markovian.h"
#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/markovian.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;
using namespace stat_tool;



void class_markovian()
{
    enum_<stat_tool::wrap_util::UniqueInt<8, 1> >("RestorationAlgorithm")
      .value("NO_COMPUTATION", 0)
      .value("FORWARD", FORWARD)
      .value("FORWARD_BACKWARD", FORWARD_BACKWARD)
      .value("VITERBI", VITERBI)
      .value("GENERALIZED_VITERBI", GENERALIZED_VITERBI)
      .value("FORWARD_BACKWARD_SAMPLING", FORWARD_BACKWARD_SAMPLING)
      .value("GIBBS_SAMPLING", GIBBS_SAMPLING)
      .value("FORWARD_DYNAMIC_PROGRAMMING", FORWARD_DYNAMIC_PROGRAMMING)
      .export_values()
    ;

}
