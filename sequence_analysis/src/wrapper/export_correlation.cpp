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
#include "stat_tool/vectors.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/sequence_label.h"
#include "tool/config.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;
using namespace stat_tool;

class CorrelationWrap {

public:

/*	static Correlation* merge(int nb_correl, const boost::list::python& input_correlation)
	{
		Format_error error;
		Correlation **correlation = NULL;

        int nb_correlation = boost::python::len(input_correlation);
        correlation = new Correlation*[nb_correlation];

        for (i=0; i<nb_correlation; i++)
        {

            boost::python::list one_correlation = extract<boost::python::list> (input_correlation[i]);
            int nb = boost::python::len(one_correlation);
            
            correlation[i] = new Correlation[nb];
        }

		correlation = merge(error, nb_correl, old_format);

		return correlation;
	}

*/

};

// Boost declaration

void class_correlation() {

	class_<Correlation, bases<STAT_interface> >
	("_Correlation", "Correlation")
	.def(init<int, int, int, int>())
	.def(init<int, int, bool, int>())
	.def(init<Correlation>())
	// Python Operators
	.def(self_ns::str(self)) //__str__
	.def("get_type", &Correlation::get_type)
    .def("get_variable_type", &Correlation::get_variable_type,
        python::args("index"))
    .def("get_variable1", &Correlation::get_variable1,
       python::args("index"))
   .def("get_variable2", &Correlation::get_variable2,
        python::args("index"))
    .def("get_white_noise", &Correlation::get_white_noise,
        python::args("lag"))
	;

    //todo
	/*
	    virtual ~Correlation();
	    Correlation* merge(Format_error &error , int nb_correl , const Correlation **icorrel) const;
	    std::ostream& line_write(std::ostream &os) const;
	    bool white_noise_correlation(Format_error &error , int nb_point , double *filter ,
	                                 int residual = true);
	    bool white_noise_correlation(Format_error &error , const Distribution &dist);
	    bool white_noise_correlation(Format_error &error , int order);
    */
}



