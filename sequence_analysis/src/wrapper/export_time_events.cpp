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
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/renewal.h"
#include "sequence_analysis/sequence_label.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;
using namespace stat_tool;

void class_time_events() {


     class_<Time_events, bases<STAT_interface> > ("_Time_events", "Time_events")
         .def(init <int>())
    // Python Operators

         .def("get_nb_element", &Time_events::get_nb_element,"nb elements")
         .def("get_nb_class", &Time_events::get_nb_class,"nb class")
    ;


/*
Time_events(int inb_element , int *itime , int *inb_event){ build(inb_element , itime , inb_event); }
Time_events(int nb_sample , const Time_events **ptimev) { merge(nb_sample , ptimev); }
Time_events(const Time_events &timev) { copy(timev); }

Distribution_data* extract(Format_error &error , int histo_type , int itime = I_DEFAULT) const;

Time_events* time_scaling(Format_error &error , int scaling_coeff) const;
Time_events* time_select(Format_error &error , int min_time ,int max_time) const;
Time_events* nb_event_select(Format_error &error , int min_nb_event ,int max_nb_event) const;

double information_computation() const;

Renewal* estimation(Format_error &error , std::ostream &os , char type ,
                        const Parametric &iinter_event , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int equilibrium_estimator = COMPLETE_LIKELIHOOD ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
Renewal* estimation(Format_error &error , std::ostream &os , char type , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int equilibrium_estimator = COMPLETE_LIKELIHOOD ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;

Histogram* get_htime() const { return htime; }
Histogram* get_hnb_event(int itime) const { return hnb_event[itime]; }
Histogram* get_hmixture() const { return hmixture; }
*/

}



