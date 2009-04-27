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
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/tops.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;
using namespace stat_tool;

////////////////////// Export tops ////////////////////////////////////////

class TopParametersWrap {

public:

	static boost::shared_ptr<Top_parameters> top_parameters_from_file(
			char* filename, int max_position) {
		Format_error error;
		Top_parameters *top_parameters = NULL;

		top_parameters = top_parameters_ascii_read(error, filename,
				max_position);

		if (!top_parameters) {
			stat_tool::wrap_util::throw_error(error);
		}

		return boost::shared_ptr<Top_parameters>(top_parameters);
	}

	static Parametric_model* extract(const Top_parameters& top, int position) {
		Format_error error;
		Parametric_model* ret;

		ret = top.extract(error, position);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Tops* simulation1(const Top_parameters& top, int nb_top,
			const Distribution &nb_trial, const Distribution &nb_axillary) {
		Format_error error;
		Tops* ret;

		ret = top.simulation(error, nb_top, nb_trial, nb_axillary);

		if (!ret)
		stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Tops* simulation2(const Top_parameters& top, int nb_top,
				int nb_trial, int nb_axillary) {
			Format_error error;
			Tops* ret;

			ret = top.simulation(error, nb_top, nb_trial, nb_axillary);

			if (!ret)
			stat_tool::wrap_util::throw_error(error);

			return ret;
		}

};

// Boost declaration

void class_top_parameters() {

	//todo : constructors

	class_<Top_parameters, bases<STAT_interface> > ("_Top_parameters",
			"Top parameters")
	.def("__init__", make_constructor(TopParametersWrap::top_parameters_from_file))
	.def(init<double, double, double, int>())
	.def(self_ns::str(self)) //__str__
    // difference between add property  and readonly ?
	.add_property("get_probability", &Top_parameters::get_probability)
	.def_readonly("get_axillary_probability",&Top_parameters::get_axillary_probability)
	.def_readonly("get_rhythm_ratio", &Top_parameters::get_rhythm_ratio)
	.def_readonly("get_max_position", &Top_parameters::get_max_position)
	.def("get_tops", &Top_parameters::get_tops,
			return_value_policy<manage_new_object> (),
			"returns tops")
	.def("extract",	&TopParametersWrap::extract,
			return_value_policy<manage_new_object> (),
			python::args("position"))
	.def("simulation1", &TopParametersWrap::simulation1,
			return_value_policy<manage_new_object> (),
			python::args("position"), "simulation type1")
	.def("simulation2",	&TopParametersWrap::simulation2,
			return_value_policy<manage_new_object> (),
			python::args("position"),
			"simulation type2")

	;

}

class TopsWrap {

public:

	static boost::shared_ptr<Tops> tops_from_file(char* filename,
			bool old_format) {
		Format_error error;
		Tops *tops = NULL;

		tops = tops_ascii_read(error, filename, old_format);

		if (!tops) {
			stat_tool::wrap_util::throw_error(error);
		}

		return boost::shared_ptr<Tops>(tops);
	}

	static boost::shared_ptr<Tops> tops_from_lists(
			const boost::python::list& identifiers,
			const boost::python::list& nb_position,
			bool init_flag) {
		Format_error error;
		Tops * ret = NULL;
		int *ids;
		int *pos;

		int nb_id = len(identifiers);
		int nb_pos = len(nb_position);

		ids = new int[nb_id];
		pos = new int[nb_pos];

		for (int i = 0; i < nb_id; i++) {
			ids[i] = boost::python::extract<int>(identifiers[i]);
		}
		for (int i = 0; i < nb_pos; i++) {
			pos[i] = boost::python::extract<int>(nb_position[i]);
		}

		ret = new Tops(nb_id, ids, pos, init_flag);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return boost::shared_ptr<Tops>(ret);

	}




	static Tops* shift(const Tops& top, int nb_internode) {
		Format_error error;
		Tops* ret;

		ret = top.shift(error, nb_internode);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Distribution_data* extract(const Tops& top, int position) {
		Format_error error;
		Distribution_data* ret;

		ret = top.extract(error, position);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Tops* select_individual(const Tops& top,
			const boost::python::list& identifiers, bool keep) {
		Format_error error;
		Tops * ret = NULL;

		int nb_id = len(identifiers);
		stat_tool::wrap_util::auto_ptr_array<int> ids(new int[nb_id]);

		for (int i = 0; i < nb_id; i++) {
			ids[i] = boost::python::extract<int>(identifiers[i]);
		}
		ret = top.select_individual(error, nb_id, ids.get(), keep);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

	static Tops* reverse(const Tops& top) {
		Format_error error;
		Tops* ret;

		ret = top.reverse(error);

		if (!ret)
			stat_tool::wrap_util::throw_error(error);

		return ret;
	}

};

// Boost declaration

void class_tops() {

	//TODO constrcuctors from sequences (seg fault),
	// what about other constructors ? from lists?
	// Tops::Tops(int nb_top , int *iidentifier , int *nb_position , bool init_flag)
	// Tops::Tops(const Tops &tops , int inb_sequence , int *index):Sequences(tops , inb_sequence , index)

	class_<Tops, bases<Sequences> > ("_Tops", "Tops")

	.def("__init__", make_constructor(TopsWrap::tops_from_file))
    .def("__init__", make_constructor(TopsWrap::tops_from_lists))
	.def(init<Sequences> ())
	.def(self_ns::str(self)) //__str__
	.def_readonly("get_max_position", &Tops::get_max_position)
	.def("get_nb_internode", &Tops::get_nb_internode,
			return_value_policy<manage_new_object> (),
			"returns histogram of nb internode")
	.def("get_top_parameters", &Tops::get_top_parameters,
			return_value_policy<manage_new_object> (),
			"returns top parameters")
	.def("get_axillary_nb_internode", &Tops::get_axillary_nb_internode,
			return_value_policy<manage_new_object> (),
			python::args("position"),
			"returns histogram of axillary nb internode")
	.def("extract",	&TopsWrap::extract,
			return_value_policy<manage_new_object> (),
			python::args("position"))
	.def("shift", &TopsWrap::shift,
			return_value_policy<manage_new_object> (),
			python::args("nb_internode"))
	.def("select_individual", &TopsWrap::select_individual,
			return_value_policy<manage_new_object> (),
			python::args("identifiers", "keep"))
	.def("reverse", &TopsWrap::reverse,
			return_value_policy<manage_new_object> ())
	.def("build_nb_internode_histogram", &Tops::build_nb_internode_histogram)

	;
}
