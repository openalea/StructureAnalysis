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

#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/compound.h"

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

	static boost::shared_ptr<Compound> compound_from_sum_and_elementary(
			boost::python::list& sum_dist,
	        boost::python::list& dist)
	  {
	   /* Format_error error;
	    Compound *compound = NULL;
	    int nb_component1, nb_component2 = 0;

	    nb_component1 = boost::python::len(sum_dist);
	    nb_component2 = boost::python::len(dist);


	    stat_tool::wrap_util::auto_ptr_array<const Parametric *>
	      sum_dist(new const Parametric*[nb_component1]);

	      stat_tool::wrap_util::auto_ptr_array<const Parametric *>
	           dist(new const Parametric*[nb_component2]);

	    compound = Compound(sum_dist, dist, COMPOUND_THRESHOLD);

	    if(!compound)
	      stat_tool::wrap_util::throw_error(error);


	    return boost::shared_ptr<Compound>(compound);
	  */
	  }

};



void class_compound()
{
	class_< Compound, bases<STAT_interface, Distribution> >
    ("_Compound", "Compound")
    .def("__init__", make_constructor(CompoundWrap::compound_from_file),
    "Build from a filename")
    .def("__init__", make_constructor(CompoundWrap::compound_from_sum_and_elementary),
    "Build from a sum distribution and an elementary distribution")
    ;
}


void class_compound_data()
{
  class_< Compound_data, bases< STAT_interface, Histogram > >
    ("_CompoundData", "Compound data")
    ;
}


