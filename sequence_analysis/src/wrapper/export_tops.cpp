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

#include "boost_python_aliases.h"
using namespace boost::python;
using namespace boost;
//using namespace stat_tool;

////////////////////// Export tops ////////////////////////////////////////

class TopParametersWrap
{

public:

  static boost::shared_ptr<Top_parameters>
  top_parameters_from_file(char* filename, int max_position)
  {
    Format_error error;
    Top_parameters *top_parameters = NULL;

    top_parameters = top_parameters_ascii_read(error, filename, max_position);

    if (!top_parameters)
        sequence_analysis::wrap_util::throw_error(error);

    return boost::shared_ptr<Top_parameters>(top_parameters);
  }

  static Parametric_model*
  extract(const Top_parameters& top, int position)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, extract, Parametric_model, position);
    }

  static Tops*
  simulation_dists(const Top_parameters& top, int nb_top,
      const Distribution &nb_trial, const Distribution &nb_axillary)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, simulation, Tops,
          nb_top, nb_trial, nb_axillary);
  }

  static Tops*
  simulation(const Top_parameters& top, int nb_top, int nb_trial,
      int nb_axillary)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, simulation, Tops,
        nb_top, nb_trial, nb_axillary);
  }

};

// Boost declaration

void class_top_parameters() {

  //todo : constructors

  class_<Top_parameters, bases<STAT_interface> >
  ("_Top_parameters", "Top parameters")
    .def("__init__", make_constructor(TopParametersWrap::top_parameters_from_file))
    .def(init< optional <double, double, double, int> >())
    .def(init< const Top_parameters& ,optional <bool> >())

    .def(self_ns::str(self)) //__str__

    .add_property("probability", &Top_parameters::get_probability)
    .add_property("axillary_probability", &Top_parameters::get_axillary_probability)
    .add_property("rhythm_ratio", &Top_parameters::get_rhythm_ratio)
    .add_property("max_position", &Top_parameters::get_max_position)

    DEF_RETURN_VALUE("get_axillary_nb_internode", &Top_parameters::get_axillary_nb_internode,args("position"), "returns axillary nb internode distribution")
    DEF_RETURN_VALUE("extract", &TopParametersWrap::extract, args("position"), "extract parametric model")
    DEF_RETURN_VALUE("simulation_dists", &TopParametersWrap::simulation_dists,args("position"), "simulation type1")
    DEF_RETURN_VALUE("simulation", &TopParametersWrap::simulation,args("position"), "simulation ")

    DEF_RETURN_VALUE_NO_ARGS("get_tops", &Top_parameters::get_tops,     "returns tops")
    ;


  //herits from Stat interface the following functions. Do we need to expose them ?
  /*
    std::ostream& line_write(std::ostream &os) const;
    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,  const char *title = 0) const;

    void axillary_nb_internode_computation(int imax_position);
  */

}

class TopsWrap
{

public:

  static boost::shared_ptr<Tops>
  tops_from_file(char* filename, bool old_format)
  {
    Format_error error;
    Tops *tops = NULL;
    tops = tops_ascii_read(error, filename, old_format);
    if (!tops)
      sequence_analysis::wrap_util::throw_error(error);
    return boost::shared_ptr<Tops>(tops);
  }

  static boost::shared_ptr<Tops>
  tops_from_lists(const boost::python::list& identifiers,
      const boost::python::list& nb_position, bool init_flag)
  {
    Format_error error;
    Tops * ret = NULL;
    int *ids;
    int *pos;

    int nb_id = len(identifiers);
    int nb_pos = len(nb_position);

    ids = new int[nb_id];
    pos = new int[nb_pos];

    for (int i = 0; i < nb_id; i++)
      {
        ids[i] = boost::python::extract<int>(identifiers[i]);
      }
    for (int i = 0; i < nb_pos; i++)
      {
        pos[i] = boost::python::extract<int>(nb_position[i]);
      }

    ret = new Tops(nb_id, ids, pos, init_flag);

    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);

    return boost::shared_ptr<Tops>(ret);

  }

  static Tops*
  shift(const Tops& top, int nb_internode)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, shift, Tops, nb_internode);
  }

  static Distribution_data*
  extract(const Tops& top, int position)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, extract, Distribution_data, position);
  }

  static Tops*
  select_individual(const Tops& top, const boost::python::list& identifiers,
      bool keep)
  {
    CREATE_ARRAY(identifiers, int);
    SIMPLE_METHOD_TEMPLATE_1(top, select_individual, Tops, size, data.get(),
        keep);
  }

  static Tops*
  reverse(const Tops& top)
  {
    SIMPLE_METHOD_TEMPLATE_0(top, reverse, Tops);
  }

  static Top_parameters*
  estimation(const Tops& top, int imin_position, int imax_position,
      int neighborhood, bool equal_probability)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, estimation, Top_parameters, imin_position,
        imax_position, neighborhood, equal_probability);
  }

  static Top_parameters*
  estimation2(const Tops& top, int neighborhood, bool equal_probability = false)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, estimation, Top_parameters, neighborhood,
        equal_probability);
  }

};

// Boost declaration

void class_tops() {

	//TODO constrcuctors from sequences (seg fault),

	class_<Tops, bases<Sequences> > ("_Tops", "Tops")
	.def(init<Sequences> ())
	.def("__init__", make_constructor(TopsWrap::tops_from_file))
    .def("__init__", make_constructor(TopsWrap::tops_from_lists))

	.def(self_ns::str(self)) //__str__

	.add_property("max_position", &Tops::get_max_position, "returns max position attribute")

	.def("build_nb_internode_histogram", &Tops::build_nb_internode_histogram)

	DEF_RETURN_VALUE("get_axillary_nb_internode", &Tops::get_axillary_nb_internode, args("position"), "returns histogram of axillary nb internode")
	DEF_RETURN_VALUE("select_individual", TopsWrap::select_individual,args("identifiers", "keep"), "select individual")
	DEF_RETURN_VALUE("extract",TopsWrap::extract,args("position"), "extract method")
    DEF_RETURN_VALUE("shift", TopsWrap::shift, args("nb_internode"), "shift method")
    DEF_RETURN_VALUE("estimation", &TopsWrap::estimation, args("min_position","max_position", "neighborhood", "equal_probability"), "estimation method 1")
    DEF_RETURN_VALUE("estimation2", &TopsWrap::estimation2, args("neighborhood", "equal_probability"), "estimation method 2")

    DEF_RETURN_VALUE_NO_ARGS("reverse", TopsWrap::reverse,"reverse method")
	DEF_RETURN_VALUE_NO_ARGS("get_nb_internode", &Tops::get_nb_internode, "returns histogram of nb internode")
	DEF_RETURN_VALUE_NO_ARGS("get_top_parameters", &Tops::get_top_parameters, "returns top parameters")
	;

//herits from Stat interface the following functions. Do we need to expose them ?
/*
  Tops(const Tops &tops , int inb_sequence , int *index);
  Tops(const Tops &tops , bool model_flag = true , bool reverse_flag = false);
  Tops(int nb_sample , const Tops **ptops);

  std::ostream& line_write(std::ostream &os) const;
  std::ostream& ascii_data_write(std::ostream &os , char format = 'c' , bool exhaustive = false) const;
  bool ascii_data_write(Format_error &error , const char *path , char format = 'c' , bool exhaustive = false) const;
  std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
  bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;
  bool spreadsheet_write(Format_error &error , const char *path) const;
  bool plot_write(Format_error &error , const char *prefix, const char *title = 0) const;

  void build_nb_internode_histogram();
*/

}
