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
 *-----------------------------------------------------------------------------*/
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

using namespace boost::python;



class CompoundWrap
{

public:

	static boost::shared_ptr<Compound> compound_from_file(char* filename)
	  {
	    Format_error error;
	    Compound *compound = NULL;
	    compound = compound_ascii_read(error, filename);
	    if(!compound)
	    {
	      stat_tool::wrap_util::throw_error(error);
	    }
	    return boost::shared_ptr<Compound>(compound);
	  }

	static boost::shared_ptr<Compound> compound_two_distributions(
			const Parametric &sum_dist,
			const Parametric &dist
			)
	  {
		  Format_error error;
		  Compound *cmpnd = NULL;

		  cmpnd = new Compound(sum_dist, dist, COMPOUND_THRESHOLD);

		  if(!cmpnd)
		      stat_tool::wrap_util::throw_error(error);

		    return boost::shared_ptr<Compound>(cmpnd);
	  }

	static boost::shared_ptr<Compound> compound_two_distributions_and_threshold(
			const Parametric &sum_dist,
			const Parametric &dist,
			double threshold
			)
	  {
		  Format_error error;
		  Compound *cmpnd = NULL;

		  cmpnd = new Compound(sum_dist,dist, threshold);

		  if(!cmpnd)
		      stat_tool::wrap_util::throw_error(error);

		    return boost::shared_ptr<Compound>(cmpnd);
	  }

	static Compound_data* simulation(const Compound& compound, int nb_element)
	  {
	    Format_error error;
	    Compound_data* ret = NULL;

	    ret = compound.simulation(error, nb_element);
	    if(!ret) stat_tool::wrap_util::throw_error(error);

	    return ret;
	  }

	static Compound_data* extract_data(const Compound& compound)
	  {
	    Format_error error;
	    Compound_data* ret = NULL;

	    ret = compound.extract_data(error);
	    if(!ret) stat_tool::wrap_util::throw_error(error);

	    return ret;
	  }

	static Parametric_model* extract_compound(const Compound& compound)
	  {
	    Parametric_model* ret;
	    Compound_data* compound_histo = NULL;

	    compound_histo = compound.get_compound_data();

	    ret = new Parametric_model(compound,
						(compound_histo ? compound_histo->get_compound() : NULL));

	    return ret;
	  }

	static Parametric_model* extract_sum_distribution(const Compound& compound)
	{
	    Parametric_model* ret;
	    Compound_data* compound_data = NULL;

	    compound_data = compound.get_compound_data();
	    ret = new Parametric_model(*(compound.get_sum_distribution()),
	              (compound_data ? compound_data->get_sum_histogram() : NULL));
	    return ret;
	}


	static Parametric_model* extract_distribution(const Compound& compound)
	{
	    Parametric_model* ret;
	    Compound_data* compound_hist = NULL;

	    compound_hist = compound.get_compound_data();
	    ret = new Parametric_model(*(compound.get_distribution()),
	                   (compound_hist ? compound_hist->get_histogram() : NULL));
	    return ret;
	}

};



void class_compound()
{
	class_< Compound, bases<STAT_interface, Distribution> >
    ("_Compound", "Compound")
    .def("__init__", make_constructor(CompoundWrap::compound_from_file),
    "Build from a filename")
    .def("__init__", make_constructor(CompoundWrap::compound_two_distributions),
    "Build from two distributions")
    .def("__init__", make_constructor(CompoundWrap::compound_two_distributions_and_threshold),
    "Build from two distributions and a threshold")
    .def("simulate", CompoundWrap::simulation, return_value_policy< manage_new_object >(),
     boost::python::arg("nb_element"), "Simulate nb_element elements")
     .def("extract_data", CompoundWrap::extract_data,
     return_value_policy< manage_new_object >(), "Return the data")
     .def("extract_compound", CompoundWrap::extract_compound,
     return_value_policy< manage_new_object >(), "Return the compound distribution")
     .def("extract_sum", CompoundWrap::extract_sum_distribution,
     return_value_policy< manage_new_object >(), "Return the sum distribution")
     .def("extract_elementary", CompoundWrap::extract_distribution,
     return_value_policy< manage_new_object >(), "Return the elementary distribution")
	;
}


void class_compound_data()
{
  class_< Compound_data, bases< STAT_interface, Histogram > >
    ("_CompoundData", "Compound data")
    ;
}


