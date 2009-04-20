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

using namespace boost::python;
using namespace boost;




////////////////////// Export tops ////////////////////////////////////////

class TopParametersWrap
{

public:

  static boost::shared_ptr<Top_parameters> top_parameters_from_file(char* filename, int max_position)
  {
    Format_error error;
    Top_parameters *top_parameters = NULL;

    top_parameters = top_parameters_ascii_read(error, filename, max_position);

/*    if(!top_parameters)
    {
	  stat_tool::wrap_util::throw_error(error);
    }
*/
    return boost::shared_ptr<Top_parameters>(top_parameters);
  }

};



// Boost declaration

void class_top_parameters()
{

  class_< Top_parameters, bases< STAT_interface > >
    ("_Top_parameters", "Top Parameters")
    .def("__init__", make_constructor(TopParametersWrap::top_parameters_from_file))
    ;
}



class TopsWrap
{

public:

  static boost::shared_ptr<Tops> tops_from_file(char* filename, bool old_format)
  {
    Format_error error;
    Tops *tops = NULL;

    tops = tops_ascii_read(error, filename, old_format);

/*    if(!top_parameters)
    {
	  stat_tool::wrap_util::throw_error(error);
    }
*/
    return boost::shared_ptr<Tops>(tops);
  }

};


// Boost declaration

void class_tops()
{

  class_< Tops, bases< Sequences > >
    ("_Tops", "Tops")
    .def("__init__", make_constructor(TopsWrap::tops_from_file))
    ;
}



