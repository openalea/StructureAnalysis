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

#define WRAP TopParametersWrap
class WRAP
{

public:

  static boost::shared_ptr<TopParameters>
  top_parameters_from_file(char *filename, int max_position)
  {
    StatError error;
    TopParameters *top_parameters = NULL;

    top_parameters = top_parameters_ascii_read(error, filename, max_position);
    if (!top_parameters)
        sequence_analysis::wrap_util::throw_error(error);

    return boost::shared_ptr<TopParameters>(top_parameters);
  }

  static DiscreteParametricModel*
  extract(const TopParameters &top, int position)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, extract, DiscreteParametricModel, position);
  }

  static Tops*
  simulation_dists(const TopParameters &top, int nb_top,
      const Distribution &nb_trial, const Distribution &nb_axillary)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, simulation, Tops, nb_top, nb_trial,
        nb_axillary);
  }

  static Tops*
  simulate(const TopParameters &top, int nb_top, int nb_trial,
      int nb_axillary)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, simulation, Tops, nb_top, nb_trial,
        nb_axillary);
  }

  static MultiPlotSet*
  get_plotable(const TopParameters &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret)
      ERROR;
    return ret;
  }

  static Distribution*
  get_axillary_nb_internode(const TopParameters &input, int position)
  {
    Distribution *res;

    if (position <=0 || position > input.get_max_position())
    {
      PyErr_SetString(PyExc_IndexError,
            "position must be positive and less or equal to max_position ");
      boost::python::throw_error_already_set();
    }


    res = new Distribution(*input.get_axillary_nb_internode(position));
    return res;
  }

};

// Boost declaration

void class_top_parameters() {

  //todo : constructors

  class_<TopParameters, bases<StatInterface> >
  ("_TopParameters", "Top parameters")
    .def("__init__", make_constructor(TopParametersWrap::top_parameters_from_file))
    .def(init< optional <double, double, double, int> >())
    .def(init< const TopParameters& ,optional <bool> >())

    .def(self_ns::str(self)) //__str__

    .add_property("probability", &TopParameters::get_probability)
    .add_property("axillary_probability", &TopParameters::get_axillary_probability)
    .add_property("rhythm_ratio", &TopParameters::get_rhythm_ratio)
    .add_property("max_position", &TopParameters::get_max_position)

    DEF_RETURN_VALUE("get_axillary_nb_internode", WRAP::get_axillary_nb_internode, args("position"),
        "returns axillary nb internode distribution")
    DEF_RETURN_VALUE("extract", WRAP::extract, args("position"), "extract parametric model")
    DEF_RETURN_VALUE("simulation_dists", WRAP::simulation_dists,args("position"), "simulation type1")
    DEF_RETURN_VALUE("simulate", WRAP::simulate,args("position"), "simulate ")
    
    DEF_RETURN_VALUE_NO_ARGS("get_tops", &TopParameters::get_tops,     "returns tops")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "Return a plotable")

    ;


  //herits from Stat interface the following functions. Do we need to expose them ?
  /*
    std::ostream& line_write(std::ostream &os) const;
    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,  const char *title = 0) const;

    void axillary_nb_internode_computation(int imax_position);
  */

}
#undef WRAP

#define WRAP TopsWrap
class WRAP
{

public:

  static boost::shared_ptr<Tops>
  tops_from_file(char *filename, bool old_format)
  {
    StatError error;
    Tops *tops = NULL;
    tops = tops_ascii_read(error, filename, old_format);
    if (!tops)
      sequence_analysis::wrap_util::throw_error(error);
    return boost::shared_ptr<Tops>(tops);
  }

  static boost::shared_ptr<Tops>
  tops_from_lists(const boost::python::list& input_identifiers,
      const boost::python::list& input_position, bool init_flag)
  {
    HEADER(Tops);

    CREATE_ARRAY(input_identifiers, int, ids);
    CREATE_ARRAY(input_position, int, pos);

    ret = new Tops(ids_size, ids.get(), pos.get(), init_flag);

    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);

    return boost::shared_ptr<Tops>(ret);

  }

  static Tops*
  shift(const Tops &top, int nb_internode)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, shift, Tops, nb_internode);
  }

  static DiscreteDistributionData*
  extract(const Tops &top, int position)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, extract, DiscreteDistributionData, position);
  }

  static Tops*
  select_individual(const Tops &top, const boost::python::list& identifiers,
      bool keep)
  {
    CREATE_ARRAY(identifiers, int, data);
    SIMPLE_METHOD_TEMPLATE_1(top, select_individual, Tops, data_size, data.get(),
        keep);
  }

  static Tops*
  reverse(const Tops &top)
  {
    SIMPLE_METHOD_TEMPLATE_0(top, reverse, Tops);
  }

  static TopParameters*
  estimation(const Tops &top, int imin_position, int imax_position,
      int neighborhood, bool equal_probability)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, estimation, TopParameters, imin_position,
        imax_position, neighborhood, equal_probability);
  }

  static TopParameters*
  estimation2(const Tops &top, int neighborhood, bool equal_probability = false)
  {
    SIMPLE_METHOD_TEMPLATE_1(top, estimation, TopParameters, neighborhood,
        equal_probability);
  }

  static Tops*
  merge(const Tops &input, const boost::python::list input_timev)
   {
     HEADER(Tops);
     CREATE_ARRAY(input_timev, const Tops*, timev)
     ret = new Tops(timev_size, timev.get());
     FOOTER;
   }
  
  static FrequencyDistribution*
  get_axillary_nb_internode(const Tops &input, int position)
  {
    FrequencyDistribution *res;
    if (position <=0 || position > input.get_max_position())
    {
      PyErr_SetString(PyExc_IndexError,
            "position must be positive and less or equal to max_position ");
      boost::python::throw_error_already_set();
    }

    res = new FrequencyDistribution(*input.get_axillary_nb_internode(position));
    return res;
  }



};

// Boost declaration

void class_tops() {

  //TODO constrcuctors from sequences (seg fault),

  class_<Tops, bases<Sequences> > ("_Tops", "Tops")
    .def(init<Sequences> ())
    .def("__init__", make_constructor(WRAP::tops_from_file))
    .def("__init__", make_constructor(WRAP::tops_from_lists))

    .def(self_ns::str(self)) //__str__

    .add_property("max_position", &Tops::get_max_position, "returns max position attribute")

    .def("build_nb_internode_frequency_distribution", &Tops::build_nb_internode_frequency_distribution)

    DEF_RETURN_VALUE("get_axillary_nb_internode", WRAP::get_axillary_nb_internode, args("position"), "returns histogram of axillary nb internode")
    DEF_RETURN_VALUE("select_individual", WRAP::select_individual,args("identifiers", "keep"), "select individual")
    DEF_RETURN_VALUE("extract", WRAP::extract,args("position"), "extract method")
    DEF_RETURN_VALUE("shift", WRAP::shift, args("nb_internode"), "shift method")
    DEF_RETURN_VALUE("estimation", WRAP::estimation, args("min_position","max_position", "neighborhood", "equal_probability"), "estimation method 1")
    DEF_RETURN_VALUE("estimation2", WRAP::estimation2, args("neighborhood", "equal_probability"), "estimation method 2")

    DEF_RETURN_VALUE_NO_ARGS("reverse", WRAP::reverse,"reverse method")
    DEF_RETURN_VALUE_NO_ARGS("get_nb_internode", &Tops::get_nb_internode, "returns histogram of nb internode")
    DEF_RETURN_VALUE_NO_ARGS("get_top_parameters", &Tops::get_top_parameters, "returns top parameters")
    DEF_RETURN_VALUE("merge", WRAP::merge, args("list of top instances"), "returns top parameters")

    ;

//herits from Stat interface the following functions. Do we need to expose them ?
/*
  Tops(const Tops &tops , int inb_sequence , int *index);
  Tops(const Tops &tops , bool model_flag = true , bool reverse_flag = false);

  std::ostream& line_write(std::ostream &os) const;
  std::ostream& ascii_data_write(std::ostream &os , char format = 'c' , bool exhaustive = false) const;
  bool ascii_data_write(StatError &error , const char *path , char format = 'c' , bool exhaustive = false) const;
  std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
  bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
  bool spreadsheet_write(StatError &error , const char *path) const;
  bool plot_write(StatError &error , const char *prefix, const char *title = 0) const;

  void build_nb_internode_frequency_distribution();
*/

}
#undef WRAP
