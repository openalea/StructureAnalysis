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
 *        $Id: export_compound.cpp 6168 2009-04-01 16:42:29Z cokelaer $
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

class CurvesWrap
{


};


void class_curves()
{
	class_< Curves >("_Curves", "Curves", init<>())
    .def(init<int, int , bool , bool>())
	.def(init<Distribution> ())
	.def(init<Histogram> ())
    ;

}




