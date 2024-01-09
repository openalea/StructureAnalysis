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
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_plotable.cpp 18038 2015-04-23 07:16:50Z guedon $
 *
 *-----------------------------------------------------------------------------*/


#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "export_plotable.h"
#include "stat_tool/plotable.h"

using namespace std;
using namespace boost::python;
using namespace stat_tool;


void class_plotable()
{
  class_ < PlotPoint >("PlotPoint", init<float, float>())
    .def_readwrite("x", &PlotPoint::first )
    .def_readwrite("y", &PlotPoint::second )
    .def_readwrite("min", &PlotPoint::first )
    .def_readwrite("max", &PlotPoint::second );


  class_< SinglePlot >("SinglePlot")

    .def("add_point", (void (SinglePlot::*)(const PlotPoint&)) &SinglePlot::add_point)
    .def("add_point", (void (SinglePlot::*)(float, float)) &SinglePlot::add_point)
    .def("add_text", (void (SinglePlot::*)(float, float, const string&)) &SinglePlot::add_text)
    .def_readwrite("legend", &SinglePlot::legend)
    .def_readwrite("style", &SinglePlot::style)
    .def_readwrite("color", &SinglePlot::color)

    // surely the label needs to be refactored: it clashes with the data Point
    // since it is only used by regression.cpp by now, we'll keep it as it is
    // for the time being.
    .def_readwrite("label", &SinglePlot::label)
    .def("get_label_size", &SinglePlot::get_size, "returns number of labels")
    .def("get_label_x", &SinglePlot::get_x, "returns x position of label(i)")
    .def("get_label_y", &SinglePlot::get_y, "returns y position of label(i)")
    .def("get_label_text", &SinglePlot::get_label, "returns label text(i)")


    .def("__iter__", range(&SinglePlot::begin, &SinglePlot::end))
    .def("__len__", &SinglePlot::size)
    ;


  class_< MultiPlot >("MultiPlot", init<int>())

    .def_readwrite("title", &MultiPlot::title)
    .def_readwrite("xtics", &MultiPlot::xtics)
    .def_readwrite("ytics", &MultiPlot::ytics)

    .def_readwrite("xlabel", &MultiPlot::xlabel)
    .def_readwrite("ylabel", &MultiPlot::ylabel)

    .def_readwrite("xrange", &MultiPlot::xrange)
    .def_readwrite("yrange", &MultiPlot::yrange)

    .def_readwrite("grid", &MultiPlot::grid)
    .def_readwrite("group", &MultiPlot::group)

    .def("resize", &MultiPlot::resize)
    .def("__len__", &MultiPlot::size)
    .def("__getitem__", &MultiPlot::operator[],
	 return_internal_reference<>())
    .def("__iter__", range(&MultiPlot::begin, &MultiPlot::end))
    ;


  class_< MultiPlotSet >("MultiPlotSet", init<int>())
    .def_readwrite("title", &MultiPlotSet::title)
    .def_readwrite("border", &MultiPlotSet::border)
    //.def("variable", (std::vector<int>)&MultiPlotSet::variable)
    .def_readwrite("variable", &MultiPlotSet::variable)
    .def_readwrite("variable_nb_viewpoint", &MultiPlotSet::variable_nb_viewpoint)
    .def_readwrite("viewpoint", &MultiPlotSet::viewpoint)


    .def("__len__", &MultiPlotSet::size)
    .def("__getitem__", &MultiPlotSet::operator[], return_internal_reference<>())
    .def("__iter__", range(&MultiPlotSet::begin, &MultiPlotSet::end))
    ;

class_<std::vector<int> > ("std_vector_int")
  .def(vector_indexing_suite<std::vector<int> >())
;
}



