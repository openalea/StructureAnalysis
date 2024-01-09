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



#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequence_analysis/sequences.h"
#include "sequence_analysis/sequence_label.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"


using namespace boost::python;
using namespace boost;
using namespace stat_tool;
using namespace sequence_analysis;



class CategoricalSequenceProcessWrap {

public:

/*  static Distribution*
  get_first_occurrence(const CategoricalSequenceProcess &input, int index)
  {
    HEADER(Distribution);
    ret = input.first_occurrence(index);
    FOOTER;
  }

    static Distribution*
  get_recurrence_time(const CategoricalSequenceProcess &input, int index)
  {
    HEADER(Distribution);
    ret = input.get_recurrence_time(index);
    FOOTER;
  }

  static Distribution*
  get_nb_run(const CategoricalSequenceProcess &input, int index)
  {
    HEADER(Distribution);
    ret = input.get_nb_run(index);
    FOOTER;
  }

  static Distribution*
  get_nb_occurrence(const CategoricalSequenceProcess &input, int index)
  {
    HEADER(Distribution);
    ret = input.get_nb_occurrence(index);
    FOOTER;
  }

  static DiscreteParametric*
  get_sojourn_time(const CategoricalSequenceProcess &input, int index)
  {
    HEADER(DiscreteParametric);
    ret = input.get_sojourn_time(index);
    FOOTER;
  } */

};

// Boost declaration
//class CategoricalSequenceProcess : public CategoricalProcess {  // processus d'observation categoriel

void class_categorical_sequence_process() {

  class_<CategoricalSequenceProcess, bases<CategoricalProcess> > ("_CategoricalSequenceProcess", "CategoricalSequenceProcess")
    .def(init <int, int, int>())

/*    DEF_RETURN_VALUE_NO_ARGS("get_length", &CategoricalSequenceProcess::get_length,  "get length")
    DEF_RETURN_VALUE_NO_ARGS("get_index_value", &CategoricalSequenceProcess::get_index_value, "get index value")


    .def("get_absorption", &CategoricalSequenceProcess::get_absorption, args("value"), "get absorption")
    .def("get_leave", &CategoricalSequenceProcess::get_leave, args("value"),  "get leave")
    .def("get_no_occurrence", &CategoricalSequenceProcess::get_no_occurrence, args("value"),  "get no occurrence")

    DEF_RETURN_VALUE("get_nb_occurrence", &CategoricalSequenceProcessWrap::get_nb_occurrence, args("value"), "get nb absorption")
    DEF_RETURN_VALUE("get_nb_run", &CategoricalSequenceProcessWrap::get_nb_run, args("value"),  "get nb run")
    DEF_RETURN_VALUE("get_sojourn_time", &CategoricalSequenceProcessWrap::get_sojourn_time, args("value"), "get sojourn time")
    DEF_RETURN_VALUE("get_recurrence_time", &CategoricalSequenceProcessWrap::get_recurrence_time, args("value"), "get reccurrence time")
    DEF_RETURN_VALUE("get_first_occurrence", &CategoricalSequenceProcessWrap::get_first_occurrence, args("value"), "get firs occurrence") */

    ;

/*
    CategoricalSequenceProcess(int inb_state , DiscreteParametric **occupancy);
    CategoricalSequenceProcess(const CategoricalProcess &process);
    CategoricalSequenceProcess(const CategoricalSequenceProcess &process , char manip = 'c' , int param = true);
    //

*/

}
