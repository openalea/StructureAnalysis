/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Thomas cokelaer <Thomas.Cokelaer@cirad.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id$
 *
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/compound.h"
#include "stat_tool/distribution.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/shared_ptr.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;


class CompoundWrap
{

  public:

  WRAP_METHOD1(Compound, simulation, Compound_data, int);
  WRAP_METHOD0(Compound, extract_data, Compound_data);
  WRAP_METHOD_FILE_ASCII_WRITE( Compound);

  static boost::shared_ptr<Compound>
  compound_from_file(char* filename)
  {
    Format_error error;
    Compound *compound = NULL;
    compound = compound_ascii_read(error, filename);
    if (!compound)
      stat_tool::wrap_util::throw_error(error);
    return boost::shared_ptr<Compound>(compound);
  }
  static boost::shared_ptr<Compound>
  compound_from_dists(boost::python::list& dists)
  {
    double cumul_threshold = COMPOUND_THRESHOLD;

    return boost::shared_ptr<Compound>(compound_from_dists_and_threshold(dists,
        cumul_threshold));
  }

  static boost::shared_ptr<Compound>
  compound_from_dists_and_threshold(boost::python::list& dists,
      double cumul_threshold)
  {
    Format_error error;
    Compound *compound = NULL;

    // Test list length
    if (boost::python::len(dists) != 2)
      {
        stat_tool::wrap_util::throw_error("Input lists must contains 2 objects");
      }

    stat_tool::wrap_util::auto_ptr_array<const Parametric*> dist(
        new const Parametric*[2]);

    int i = 0;
    for (i = 0; i < 2; i++)
      {
        boost::python::extract<Parametric*> get_param(dists[i]);
        if (get_param.check())
          {
            dist[i] = new Parametric(*get_param());
          }
        else
          {
            dist[i] = new Parametric(*boost::python::extract<Distribution*>(
                dists[i])());
          }
      }

    compound = new Compound(*dist[0], *dist[1]);

    if (!compound)
      stat_tool::wrap_util::throw_error(error);

    for (i = 0; i < 2; i++)
      {
        delete dist[i];
      }

    return boost::shared_ptr<Compound>(compound);
  }

  static Parametric_model*
  extract_compound(const Compound& compound)
  {
    Parametric_model* ret;
    Compound_data* compound_data = NULL;
    compound_data = compound.get_compound_data();
    ret = new Parametric_model(*((Distribution*) (&compound)),
        (Histogram*) compound_data);
    return ret;
  }

  static Parametric_model*
  extract_sum(const Compound& compound)
  {
    Parametric_model* ret;
    Compound_data* compound_data = NULL;
    compound_data = compound.get_compound_data();
    ret = new Parametric_model(*(compound.get_sum_distribution()),
        (compound_data ? compound_data->get_sum_histogram() : NULL));
    return ret;
  }

  static Parametric_model*
  extract_elementary(const Compound& compound)
  {
    Parametric_model* ret;
    Compound_data* compound_data = NULL;
    compound_data = compound.get_compound_data();

    ret = new Parametric_model(*(compound.get_distribution()),
        (compound_data ? compound_data->get_histogram() : NULL));
    return ret;
  }

  static MultiPlotSet*
  survival_get_plotable(const Compound& p)
  {
    Format_error error;
    MultiPlotSet* ret = p.survival_get_plotable(error);
    if (!ret)
      ERROR;
    return ret;
  }


};


#define WRAP CompoundWrap
void class_compound()
{
    class_< Compound, bases<STAT_interface, Distribution> >  ("_Compound", "Compound" )

    //"constructor from 2 distribution and an optional cumul threshold")
    // .def(init<Parametric, Parametric, optional<double> >())
    .def("__init__", make_constructor(WRAP::compound_from_file),	"Build from a filename")
    .def("__init__", make_constructor(WRAP::compound_from_dists),"Build from list of objects")
    .def("__init__", make_constructor(WRAP::compound_from_dists_and_threshold),"Build from list of objects")

    .def(self_ns::str(self)) // __str__

    DEF_RETURN_VALUE("simulate", WRAP::simulation,ARGS("nb_element"), "Simulate nb_element elements")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data,"Return the data")
    DEF_RETURN_VALUE_NO_ARGS("extract_compound", WRAP::extract_compound, "Return the compound distribution")
    DEF_RETURN_VALUE_NO_ARGS("extract_sum", WRAP::extract_sum,	"Return the sum distribution")
    DEF_RETURN_VALUE("extract_elementary", WRAP::extract_elementary,	ARGS("index"),	"Return the elementary distribution")
    DEF_RETURN_VALUE_NO_ARGS("file_ascii_write", WRAP::file_ascii_write, "Save Compound into a file")
    DEF_RETURN_VALUE_NO_ARGS("survival_get_plotable", WRAP::survival_get_plotable, "Return a survival plotable")
    ;

    /*
     * Compound(const Parametric &sum_dist , const Parametric &dist , char type);
      Compound(const Compound &compound , bool data_flag = true)
      void computation(int min_nb_value = 1 ,
	                     double cumul_threshold = COMPOUND_THRESHOLD ,
	                     bool sum_flag = true , bool dist_flag = true);
    // done in compound_data so
     * 	*/
}
#undef WRAP


class CompoundDataWrap
{

public:

  WRAP_METHOD1(Compound_data, extract, Distribution_data, char);

  static Distribution_data*
  extract_sum_distribution(const Compound_data& input)
  {
    Distribution_data* ret;
    ret = new Distribution_data(*(input.get_sum_histogram()),
        input.get_compound()->get_sum_distribution());
    return ret;
  }

  static Distribution_data*
  extract_distribution(const Compound_data& input)
  {
    Distribution_data* ret;
    ret = new Distribution_data(*(input.get_histogram()),
        input.get_compound()->get_distribution());
    return ret;
  }
};


#define WRAP CompoundDataWrap
void class_compound_data()
{
  class_< Compound_data, bases< STAT_interface, Histogram > > ("_CompoundData", "Compound data")
  DEF_RETURN_VALUE_NO_ARGS("extract", WRAP::extract, "Return the data")
  DEF_RETURN_VALUE_NO_ARGS("extract_sum", WRAP::extract_sum_distribution, "Return the sum distribution")
  DEF_RETURN_VALUE("extract_elementary", WRAP::extract_distribution,ARGS("index"), "Return the elementary distribution")
    ;


  /*
  Compound_data();
    Compound_data(const Histogram &histo , const Compound &icompound);
    Compound_data(const Compound &icompound);
    Compound_data(const Compound_data &compound_histo , bool model_flag = true)    :Histogram(compound_histo) { copy(compound_histo , model_flag); }

    std::ostream& line_write(std::ostream &os) const;
    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,  bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,  const char *title = 0) const;

    Compound* get_compound() const { return compound; }
*/


}
#undef WRAP

