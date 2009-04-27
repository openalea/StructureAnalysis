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
 *        $Id: export_tops.cpp 6169 2009-04-01 16:42:59Z cokelaer $
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

using namespace boost::python;
using namespace boost;
using namespace stat_tool;

class NonParametricSequenceProcessWrap {

public:

    static Distribution* get_first_occurrence(const Nonparametric_sequence_process& input, int index)
    {
        Format_error error;
        //strcpy(error,"error to be filled");
        Distribution* ret;
        ret = input.get_first_occurrence(index);
        if (!ret)
            stat_tool::wrap_util::throw_error(error);
        return ret;
    }
    static Distribution* get_recurrence_time(const Nonparametric_sequence_process& input, int index)
    {
        Format_error error;
        Distribution* ret;
        ret = input.get_recurrence_time(index);
        if (!ret)
            stat_tool::wrap_util::throw_error(error);
        return ret;
    }
    static Distribution* get_nb_run(const Nonparametric_sequence_process& input, int index)
    {
        Format_error error;
        Distribution* ret;
        ret = input.get_nb_run(index);
        if (!ret)
            stat_tool::wrap_util::throw_error(error);
        return ret;
    }
    static Distribution* get_nb_occurrence(const Nonparametric_sequence_process& input, int index)
    {
        Format_error error;
        Distribution* ret;
        ret = input.get_nb_occurrence(index);
        if (!ret)
            stat_tool::wrap_util::throw_error(error);
        return ret;
    }
    static Parametric* get_sojourn_time(const Nonparametric_sequence_process& input, int index)
    {
        Format_error error;
        Parametric* ret;
        ret = input.get_sojourn_time(index);
        if (!ret)
            stat_tool::wrap_util::throw_error(error);
        return ret;
    }



/*
 * Distribution** get_first_occurrence() const { return first_occurrence; }
    Distribution** get_recurrence_time() const { return recurrence_time; }
    Parametric** get_sojourn_time() const { return sojourn_time; }
    Distribution** get_nb_run() const { return nb_run; }
    Distribution** get_nb_occurrence() const { return nb_occurrence; }
*/
};

// Boost declaration
//class Nonparametric_sequence_process : public Nonparametric_process {  // processus d'observation non-parametrique

void class_nonparametric_sequence_process() {

	class_<Nonparametric_sequence_process, bases<Nonparametric_process> > ("_Nonparametric_sequence_process", "Nonparametric_sequence_process")
	.def(init <int, int, int>())

    .def("get_length", &Nonparametric_sequence_process::get_length,
        return_value_policy<manage_new_object> (),
        "get length")
    .def("get_index_value", &Nonparametric_sequence_process::get_index_value,
        return_value_policy<manage_new_object> (),
        "get index value")
    .def("get_absorption", &Nonparametric_sequence_process::get_absorption,
        python::args("value"),
        "get absorption")
    .def("get_leave", &Nonparametric_sequence_process::get_leave,
        python::args("value"),
        "get leave")
    .def("get_no_occurrence", &Nonparametric_sequence_process::get_no_occurrence,
        python::args("value"),
        "get no occurrence")

    .def("get_nb_occurrence", &NonParametricSequenceProcessWrap::get_nb_occurrence,
        return_value_policy<manage_new_object> (),
        python::args("value"),
        "get nb absorption")
    .def("get_nb_run", &NonParametricSequenceProcessWrap::get_nb_run,
        return_value_policy<manage_new_object> (),
        python::args("value"),
        "get nb run")
    .def("get_sojourn_time", &NonParametricSequenceProcessWrap::get_sojourn_time,
        return_value_policy<manage_new_object> (),
        python::args("value"),
        "get sojourn time")
    .def("get_recurrence_time", &NonParametricSequenceProcessWrap::get_recurrence_time,
        return_value_policy<manage_new_object> (),
        python::args("value"),
        "get reccurrence time")
    .def("get_first_occurrence", &NonParametricSequenceProcessWrap::get_first_occurrence,
        return_value_policy<manage_new_object> (),
        python::args("value"),
        "get firs occurrence")

;

/*

    Nonparametric_sequence_process(int inb_state , Parametric **occupancy);
    Nonparametric_sequence_process(const Nonparametric_process &process);
    Nonparametric_sequence_process(const Nonparametric_sequence_process &process , char manip = 'c' , int param = true);
   

    //
    Distribution** get_first_occurrence() const { return first_occurrence; }
    Distribution** get_recurrence_time() const { return recurrence_time; }
    Parametric** get_sojourn_time() const { return sojourn_time; }
    Distribution** get_nb_run() const { return nb_run; }
    Distribution** get_nb_occurrence() const { return nb_occurrence; }
*/    


}



