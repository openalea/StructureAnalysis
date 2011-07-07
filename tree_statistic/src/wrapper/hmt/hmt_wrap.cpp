/*------------------------------------------------------------------------------
 *
 *        VPlants.Tree-Statistic : VPlants Tree-Statistic module
 *        HiddenMarkovTrees
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: hmt_wrap.cpp 9099 2010-06-08 09:03:00Z pradal $
 *
 *-----------------------------------------------------------------------------*/



/* WRAPPER Boost.python for hmt classes */
#include "export_hmt.h"

#include "stat_tool/wrapper_util.h"

#include <boost/python.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 103400
#include <boost/python/docstring_options.hpp>
#endif

using namespace boost::python;
using namespace stat_tool::wrap_util;



// Define python module "_stat_tool"
BOOST_PYTHON_MODULE(_hmt)
{
  //show_user_defined : true
  //show_signatures : false
#if BOOST_VERSION >= 103400
  docstring_options doc_options(true, false);
#endif

  StatErrorClass = object(handle<>(PyErr_NewException("_stat_tool.StatError",NULL,NULL)));
  scope().attr("StatError") = StatErrorClass;

  class_errors();
  enum_entropy_algorithms();

  class_hmt();
  class_hmiot();
  class_hmt_data();
}

