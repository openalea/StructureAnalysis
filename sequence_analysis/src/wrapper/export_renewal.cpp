/*------------------------------------------------------------------------------
 *
 *        VPlants.Sequence_analysis : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id:  $
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

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;
using namespace sequence_analysis;


#define WRAP RenewalWrap
class RenewalWrap {

public:

  static boost::shared_ptr<Renewal>
  constructor_from_file(char* filename)
  {
    Format_error error;
    Renewal *renewal = NULL;
    bool old_format = false;

    renewal = renewal_ascii_read(error, filename, old_format);

    return boost::shared_ptr<Renewal>(renewal);
  }

  static boost::shared_ptr<Renewal>
  constructor_from_inter_event(const Parametric &inter_event, char type, int time)
  {
    Format_error error;
    Renewal *renewal = NULL;

    renewal = renewal_building(error, inter_event, type, time);

    return boost::shared_ptr<Renewal>(renewal);
  }


  static void
  file_ascii_write(const Renewal& d, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
  }


  static Distribution*
  get_time(const Renewal& input)
  {
    Distribution* ret;
    ret = new Distribution(*input.get_time());
    return ret;
}

  static Renewal_data*
  get_renewal_data(const Renewal& input)
  {
    Renewal_data* ret;
    ret = new Renewal_data(*input.get_renewal_data());
    return ret;
  }


  static Renewal_data*
  simulation_histogram(const Renewal& input, char itype, const Histogram & ihtime)
  {
    HEADER(Renewal_data);
    ret = input.simulation(error, itype , ihtime);
    FOOTER;

  }

  static Renewal_data*
  simulation_nb_elements(const Renewal& input, char itype , int nb_element , int itime)
  {
    HEADER(Renewal_data);
    ret = input.simulation(error, itype , nb_element, itime);
    FOOTER;
  }

  static Renewal_data*
  simulation_time_events(const Renewal& input, char itype , int nb_element,
		  const Time_events &timev)
  {
    HEADER(Renewal_data);
    ret = input.simulation(error, itype , nb_element, timev);
    FOOTER;
  }

  static Parametric_model*
  extract(const Renewal& seq, int type, int state)
  {
    HEADER(Parametric_model);
    ret = seq.extract(error, type, state);
    FOOTER;
  }

};

// Boost declaration

void class_renewal() {

  class_<Renewal, bases<STAT_interface> > ("_Renewal", "Renewal")
    //type = 'o' or 'e'
    .def("__init__", make_constructor(RenewalWrap::constructor_from_file))
    .def("__init__", make_constructor(RenewalWrap::constructor_from_inter_event))

    .def(init <char, Histogram, Parametric>())
    .def(init <char, Distribution, Parametric>())

    // Python Operators
    .def(self_ns::str(self)) //__str__

    .add_property("nb_iterator", &Renewal::get_nb_iterator,"nb iterator")
    .add_property("type", &Renewal::get_type,"type")

    .def("file_ascii_write", RenewalWrap::file_ascii_write,"Save vector summary into a file")

    DEF_RETURN_VALUE_NO_ARGS("get_renewal_data", RenewalWrap::get_renewal_data,"returns renewal data")
    DEF_RETURN_VALUE_NO_ARGS("get_time", RenewalWrap::get_time, "returns time")
    DEF_RETURN_VALUE("simulation_histogram", RenewalWrap::simulation_histogram, args("todo"), "simulation")
    .def("extract", WRAP::extract, return_value_policy<manage_new_object> (),  python::args("type", "state"), "Extract distribution data")

    DEF_RETURN_VALUE("simulation_nb_elements", RenewalWrap::simulation_nb_elements, args("todo"), "simulation")
    DEF_RETURN_VALUE("simulation_time_events", RenewalWrap::simulation_time_events, args("todo"), "simulation")
    ;

/*
    friend Renewal* renewal_building(Format_error &error ,
    const Parametric &inter_event ,   char type, int time);

   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool ascii_write(Format_error &error , const char *path ,     bool exhaustive = false) const;
   bool spreadsheet_write(Format_error &error , const char *path) const;
   bool plot_write(Format_error &error , const char *prefix ,  const char *title = 0) const;


	void computation(bool inter_event_flag = true , char itype = 'v' , const Distribution *dtime = 0);
	double likelihood_computation(const Time_events &timev) const;

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
    .def("get_renewal", &Renewal_data::get_renewal, return_value_policy<manage_new_object> (),"get renewal")
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
