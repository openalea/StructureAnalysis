/*------------------------------------------------------------------------------
 *
 *        VPlants.Sequence_analysis : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id:  $
 *
 *-----------------------------------------------------------------------------*/



/* WRAPPER Boost.python for sequences class */
#include <cstdlib>
#include "stat_tool/stat_tools.h"

#include "export_base.h"
#include "export_function.h"
// #include "export_tops.h"
#include "export_sequences.h"
#include "export_correlation.h"
#include "export_markovian_sequences.h"
#include "export_nonhomogeneous_markov.h"
#include "export_categorical_sequence_process.h"
#include "export_renewal.h"
#include "export_time_events.h"
#include "export_semi_markov.h"
#include "export_hidden_semi_markov.h"
#include "export_variable_order_markov.h"
#include "export_hidden_variable_order_markov.h"

#include <boost/python.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 103400
#include <boost/python/docstring_options.hpp>
#endif

using namespace boost::python;


// Define python module "_sequence_analysis"
BOOST_PYTHON_MODULE(_sequence_analysis)
{
#if BOOST_VERSION >= 103400
  docstring_options doc_options(true, false);
#endif

  class_constant_sequence();

  class_function();

  class_sequences();
  class_sequence_characteristics();

  class_markovian_sequences();
  class_self_transition();

  class_correlation();

//  class_tops();
//  class_top_parameters();

  class_nonhomogeneous_markov();
  class_nonhomogeneous_markov_data();

  class_categorical_sequence_process();

  class_time_events();

  class_renewal();
  class_renewal_data();
  class_renewal_iterator();

  class_semi_markov();
  class_semi_markov_data();
  class_semi_markov_iterator();
  class_hidden_semi_markov();

  class_variable_order_markov();
  class_variable_order_markov_data();
  class_variable_order_markov_iterator();

  class_hidden_variable_order_markov();
  def( "srand", srand );
}

