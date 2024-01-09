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
 *        $Id: stat_tool_wrap.cpp 18041 2015-04-23 07:17:30Z guedon $
 *
 *-----------------------------------------------------------------------------*/



/* WRAPPER Boost.python for stat_tool class */
#include "export_plotable.h"

#include "export_base.h"
#include "export_curves.h"
#include "export_distribution.h"
#include "export_frequency_distribution.h"
#include "export_vectors.h"
#include "export_mixture.h"
#include "export_convolution.h"
#include "export_compound.h"
#include "export_distance_matrix.h"
#include "export_markovian.h"
#include "export_regression.h"
#include "export_categorical_process.h"
#include "export_chain.h"
#include "export_chain_reestimation.h"

#include "wrapper_util.h"

#include <boost/python.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 103400
#include <boost/python/docstring_options.hpp>
#endif

using namespace boost::python;
using namespace stat_tool::wrap_util;



// Define python module "_stat_tool"
BOOST_PYTHON_MODULE(_stat_tool)
{
  //show_user_defined : true
  //show_signatures : false
#if BOOST_VERSION >= 103400
  docstring_options doc_options(true, false);
#endif

  StatErrorClass = object(handle<>(PyErr_NewException("_stat_tool.StatError",NULL,NULL)));
  scope().attr("StatError") = StatErrorClass;

  class_constant();
  class_stat_error();
  class_stat_interface();

  class_regression_kernel();
  class_regression();

  class_vectors();
  class_vectordistance();

  class_distance_matrix();
  class_cluster();
  class_dendrogram();

  class_distribution();
  class_discrete_parametric();
  class_discrete_parametric_model();

  // depends on discrete_parametric
  class_forward();

  class_frequency_distribution();
  class_distribution_data();

  class_mixture();
  class_mixture_data();

  class_multivariate_mixture();
  class_multivariate_mixture_data();

  class_convolution();
  class_convolution_data();

  class_compound();
  class_compound_data();

  class_markovian();

  class_plotable();

  class_curves();

  class_non_parametric_process();
  class_chain();
  class_chain_reestimation();
}

