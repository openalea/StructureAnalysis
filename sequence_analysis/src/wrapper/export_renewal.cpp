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

class RenewalWrap {

public:

	static boost::shared_ptr<Renewal> renewal_from_file(char* filename) {
		Format_error error;
		Renewal *renewal = NULL;
		bool old_format = false;

		renewal = renewal_ascii_read(error, filename, old_format);

		/*    if(!top_parameters)
		 {
		 stat_tool::wrap_util::throw_error(error);
		 }
		 */
		return boost::shared_ptr<Renewal>(renewal);
	}



};

// Boost declaration

void class_renewal() {


	class_<Renewal, bases<STAT_interface> > ("_Renewal", "Renewal")
    .def(init <char, Histogram, Parametric>())
    .def(init <char, Distribution, Parametric>())
	// Python Operators

    .def("get_nb_iterator", &Renewal::get_nb_iterator,"nb iterator")
	;
/*

Renewal(const Renewal_data &irenewal_data , const Parametric &iinter_event);
Renewal(const Renewal &renew , bool data_flag = true){ copy(renew , data_flag); }
void conditional_delete();

Parametric_model* extract(Format_error &error , int dist_type , int itime = I_DEFAULT) const;

void computation(bool inter_event_flag = true , char itype = 'v' ,                     const Distribution *dtime = 0);
double likelihood_computation(const Time_events &timev) const;
Renewal_data* simulation(Format_error &error , char itype ,                             const Histogram &ihtime) const;
Renewal_data* simulation(Format_error &error , char itype ,                             int nb_element , int itime) const;
Renewal_data* simulation(Format_error &error , char itype ,                             int nb_element , const Time_events &itimev) const;

int get_nb_iterator() const { return nb_iterator; }
Renewal_data* get_renewal_data() const { return renewal_data; }
char get_type() const { return type; }
Distribution* get_time() const { return time; }
Parametric* get_inter_event() const { return inter_event; }
Length_bias* get_length_bias() const { return length_bias; }
Backward* get_backward() const { return backward; }
Forward* get_forward() const { return forward; }
Parametric* get_nevent_time(int inb_event) const { return nevent_time[inb_event]; }
Nb_event* get_nb_event(int itime) const { return nb_event[itime]; }
Distribution* get_mixture() const { return mixture; }
Curves* get_index_event() const { return index_event; }
*/
}


void class_renewal_data() {


	class_<Renewal_data, bases<Time_events> > ("_Renewal_data", "Renewal_data")
    .def(init <int, int>())
    .def(init <int, Renewal>())
    .def(init <Time_events, int>())
	// Python Operators
    .def("get_renewal", &Renewal_data::get_renewal, 
        return_value_policy<manage_new_object> (),"get renewal")
    .def("get_type", &Renewal_data::get_type, "get type")



/*
    Renewal_data(int nb_sample , const Renewal_data **itimev);
    Renewal_data(const Renewal_data &timev , bool model_flag = true)    :Time_events(timev) { copy(timev , model_flag); }

    Renewal_data* merge(Format_error &error , int nb_sample , const Renewal_data **itimev) const;
    Distribution_data* extract(Format_error &error , int histo_type , int itime = I_DEFAULT) const;

    Renewal* estimation(Format_error &error , std::ostream &os , const Parametric &iinter_event ,
                        int estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Renewal* estimation(Format_error &error , std::ostream &os , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int mean_computation = COMPUTED ,
                        double weight = D_DEFAULT , int penalty_type = SECOND_DIFFERENCE ,
                        int outside = ZERO) const;

    
   int get_length(int index_seq) const { return length[index_seq]; }
   int get_sequence(int index_seq , int index) const
   { return sequence[index_seq][index]; }
   Histogram* get_inter_event() const { return inter_event; }
   Histogram* get_within() const { return within; }
   Histogram* get_length_bias() const { return length_bias; }
   Histogram* get_backward() const { return backward; }
   Histogram* get_forward() const { return forward; }
   Curves* get_index_event() const { return index_event; }
   };
*/
;

}
