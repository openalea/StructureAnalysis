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
 *        $Id: export_curves.cpp 6168 2009-04-01 16:42:29Z cokelaer $
 *
 *-----------------------------------------------------------------------------*/


#include <stdio.h>
#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost::python;
using namespace stat_tool;


class CurvesWrap
{


};


void class_curves()
{
  class_< Curves >("_Curves", "Curves")
  .def(init<int, int , optional<bool , bool> >())
  .def(init<const Curves &, optional<char, int> >())

  .def(init<Distribution> ())
  .def(init<FrequencyDistribution> ())

  .def_readonly("nb_curve", &Curves::nb_curve)
  .def_readonly("length", &Curves::nb_curve)
  .def_readonly("offset", &Curves::nb_curve)
  ;



	/*
	    int *frequency;         // effectifs correspondant a chaque abscisse
	    double **point;         // points des courbes

	    std::ostream& ascii_print(std::ostream &os , bool file_flag = false ,
	                              const Curves *curves = 0) const;
	    std::ostream& spreadsheet_print(std::ostream &os , const Curves *curves = 0) const;
	    int plot_length_computation() const;
	    bool plot_print(const char *path , int ilength = I_DEFAULT ,
	                    const Curves *curves_0 = 0 , const Curves *curves_1 = 0) const;
	    bool plot_print_standard_residual(const char *path , double *standard_residual = 0) const;  // sequence_analysis
	    void plotable_print(int index , SinglePlot &plot) const;

	    double mean_computation(int index) const;
	    double total_square_sum_computation(int index , double mean) const;

	    int max_frequency_computation() const;
	    int nb_element_computation() const;



	     * */




}




