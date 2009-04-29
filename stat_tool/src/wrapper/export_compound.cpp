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

#include "boost_python_aliases.h"

using namespace boost::python;


class CompoundWrap
{

public:

	WRAP_METHOD1(Compound, simulation, Compound_data, int);
	WRAP_METHOD0(Compound, extract_data, Compound_data);
	WRAP_METHOD_FILE_ASCII_WRITE(Compound);

	static boost::shared_ptr<Compound> compound_from_file(char* filename)
	{
		Format_error error;
		Compound *compound = NULL;
		compound = compound_ascii_read(error, filename);
		if(!compound)
		    stat_tool::wrap_util::throw_error(error);
	    return boost::shared_ptr<Compound>(compound);
	}

	static boost::shared_ptr<Compound> compound_two_distributions(
			const Parametric &sum_dist,	const Parametric &dist)
	{
		Format_error error;
		Compound *cmpnd = NULL;
		cmpnd = new Compound(sum_dist, dist, COMPOUND_THRESHOLD);
		if(!cmpnd)
			stat_tool::wrap_util::throw_error(error);
		return boost::shared_ptr<Compound>(cmpnd);
	}

	static boost::shared_ptr<Compound> compound_two_distributions_and_threshold(
			const Parametric &sum_dist,	const Parametric &dist,	double threshold)
	{
		Format_error error;
		Compound *cmpnd = NULL;
		cmpnd = new Compound(sum_dist, dist, threshold);
		if(!cmpnd)
		    stat_tool::wrap_util::throw_error(error);
		return boost::shared_ptr<Compound>(cmpnd);
	}


	static Parametric_model* extract_compound(const Compound& compound)
	{
	    Parametric_model* ret;
	    Compound_data* compound_hist = NULL;
	    compound_hist = compound.get_compound_data();
	    ret = new Parametric_model(compound,
				(compound_hist ? compound_hist->get_compound() : NULL));
	    return ret;
	}

	static Parametric_model* extract_sum_distribution(const Compound& compound)
	{
	    Parametric_model* ret;
	    Compound_data* compound_hist = NULL;
	    compound_hist = compound.get_compound_data();
	    ret = new Parametric_model(*(compound.get_sum_distribution()),
	            (compound_hist ? compound_hist->get_sum_histogram() : NULL));
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


#define WRAP CompoundWrap
void class_compound()
{
	class_< Compound, bases<STAT_interface, Distribution> >
    ("_Compound", "Compound" )
    DEF_INIT_MAKE_CONSTRUCTOR(WRAP::compound_from_file,
    		"Build from a filename")
    DEF_INIT_MAKE_CONSTRUCTOR(WRAP::compound_two_distributions,
    		"Build from two distributions")
    DEF_INIT_MAKE_CONSTRUCTOR(WRAP::compound_two_distributions_and_threshold,
    		"Build from two distributions and a threshold")

    DEF_STR()

    DEF_RETURN_VALUE("simulate", WRAP::simulation,ARGS("nb_element"),
    		"Simulate nb_element elements")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data,"Return the data")
    DEF_RETURN_VALUE_NO_ARGS("extract_compound", WRAP::extract_compound,
    		"Return the compound distribution")
    DEF_RETURN_VALUE_NO_ARGS("extract_sum", WRAP::extract_sum_distribution,
    		"Return the sum distribution")
    DEF_RETURN_VALUE("extract_elementary", WRAP::extract_distribution,
    		ARGS("index"),
    		"Return the elementary distribution")
	DEF_RETURN_VALUE_NO_ARGS("file_ascii_write", WRAP::file_ascii_write,
			"Save Compound into a file")
	;
}
#undef WRAP


class CompoundDataWrap
{

public:

	WRAP_METHOD1(Compound_data, extract, Distribution_data, char);

	static Distribution_data* extract_sum_distribution(const Compound_data& input)
	{
	    Distribution_data* ret;
	    ret = new Distribution_data(*(input.get_sum_histogram()),
	    		input.get_compound()->get_sum_distribution());
	    return ret;
	}


	static Distribution_data* extract_distribution(const Compound_data& input)
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
  class_< Compound_data, bases< STAT_interface, Histogram > >
	 ("_CompoundData", "Compound data")
     DEF_RETURN_VALUE_NO_ARGS("extract", WRAP::extract, "Return the data")
     DEF_RETURN_VALUE_NO_ARGS("extract_sum", WRAP::extract_sum_distribution, "Return the sum distribution")
     DEF_RETURN_VALUE("extract_elementary", WRAP::extract_distribution,ARGS("index"),"Return the elementary distribution")
    ;
}
#undef WRAP

