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
#include "stat_tool/markovian.h"
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


class NonParametricSequenceProcessWrap {

public:

  static Distribution*
  get_first_occurrence(const Nonparametric_sequence_process& input, int index)
  {
    HEADER(Distribution);
    ret = input.get_first_occurrence(index);
    FOOTER;
  }

    static Distribution*
  get_recurrence_time(const Nonparametric_sequence_process& input, int index)
  {
    HEADER(Distribution);
    ret = input.get_recurrence_time(index);
    FOOTER;
  }

  static Distribution*
  get_nb_run(const Nonparametric_sequence_process& input, int index)
  {
    HEADER(Distribution);
    ret = input.get_nb_run(index);
    FOOTER;
  }

  static Distribution*
  get_nb_occurrence(const Nonparametric_sequence_process& input, int index)
  {
    HEADER(Distribution);
    ret = input.get_nb_occurrence(index);
    FOOTER;
  }

  static Parametric*
  get_sojourn_time(const Nonparametric_sequence_process& input, int index)
  {
    HEADER(Parametric);
    ret = input.get_sojourn_time(index);
    FOOTER;
  }

};

// Boost declaration
//class Nonparametric_sequence_process : public Nonparametric_process {  // processus d'observation non-parametrique

void class_nonparametric_sequence_process() {

  class_<Nonparametric_sequence_process, bases<Nonparametric_process> > ("_Nonparametric_sequence_process", "Nonparametric_sequence_process")
    .def(init <int, int, int>())

    DEF_RETURN_VALUE_NO_ARGS("get_length", &Nonparametric_sequence_process::get_length,  "get length")
    DEF_RETURN_VALUE_NO_ARGS("get_index_value", &Nonparametric_sequence_process::get_index_value, "get index value")


    .def("get_absorption", &Nonparametric_sequence_process::get_absorption, args("value"), "get absorption")
    .def("get_leave", &Nonparametric_sequence_process::get_leave, args("value"),  "get leave")
    .def("get_no_occurrence", &Nonparametric_sequence_process::get_no_occurrence, args("value"),  "get no occurrence")

    DEF_RETURN_VALUE("get_nb_occurrence", &NonParametricSequenceProcessWrap::get_nb_occurrence, args("value"), "get nb absorption")
    DEF_RETURN_VALUE("get_nb_run", &NonParametricSequenceProcessWrap::get_nb_run, args("value"),  "get nb run")
    DEF_RETURN_VALUE("get_sojourn_time", &NonParametricSequenceProcessWrap::get_sojourn_time, args("value"), "get sojourn time")
    DEF_RETURN_VALUE("get_recurrence_time", &NonParametricSequenceProcessWrap::get_recurrence_time, args("value"), "get reccurrence time")
    DEF_RETURN_VALUE("get_first_occurrence", &NonParametricSequenceProcessWrap::get_first_occurrence, args("value"), "get firs occurrence")

    ;

/*
    Nonparametric_sequence_process(int inb_state , Parametric **occupancy);
    Nonparametric_sequence_process(const Nonparametric_process &process);
    Nonparametric_sequence_process(const Nonparametric_sequence_process &process , char manip = 'c' , int param = true);
    //

*/

}
