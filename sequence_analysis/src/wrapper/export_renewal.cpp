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
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequence_analysis/sequences.h"
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
using namespace stat_tool;
using namespace sequence_analysis;



#define WRAP RenewalWrap
class WRAP {

public:

  static boost::shared_ptr<Renewal>
  constructor_from_file(char *filename)
  {
    StatError error;
    Renewal *renewal = NULL;
    bool old_format = false;

    renewal = renewal_ascii_read(error, filename, old_format);

    return boost::shared_ptr<Renewal>(renewal);
  }

  static boost::shared_ptr<Renewal>
  constructor_from_inter_event(const DiscreteParametric &inter_event, char type, int time)
  {
    StatError error;
    Renewal *renewal = NULL;

    renewal = renewal_building(error, inter_event, type, time);

    return boost::shared_ptr<Renewal>(renewal);
  }


  static void
  file_ascii_write(const Renewal &d, const char *path, bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
  }


  static Distribution*
  get_time(const Renewal &input)
  {
    Distribution* ret;
    ret = new Distribution(*input.get_time());
    return ret;
}

  static RenewalData*
  get_renewal_data(const Renewal &input)
  {
    RenewalData* ret;
    ret = new RenewalData(*input.get_renewal_data());
    return ret;
  }


  static RenewalData*
  simulation_histogram(const Renewal &input, char itype, const FrequencyDistribution &ihtime)
  {
    HEADER(RenewalData);
    ret = input.simulation(error, itype , ihtime);
    FOOTER;

  }

  static RenewalData*
  simulation_nb_elements(const Renewal &input, char itype, int nb_element, int itime)
  {
    HEADER(RenewalData);
    ret = input.simulation(error, itype , nb_element, itime);
    FOOTER;
  }

  static RenewalData*
  simulation_time_events(const Renewal &input, char itype , int nb_element,
		  const TimeEvents &timev)
  {
    HEADER(RenewalData);
    ret = input.simulation(error, itype , nb_element, timev);
    FOOTER;
  }

  static DiscreteParametricModel*
  extract(const Renewal &seq, int type, int state)
  {
    HEADER(DiscreteParametricModel);
    ret = seq.extract(error, type, state);
    FOOTER;
  }

  static MultiPlotSet* get_plotable(const Renewal &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret) ERROR;
    return ret;
  }




};

// Boost declaration

void class_renewal() {

  class_<Renewal, bases<StatInterface> > ("_Renewal", "Renewal", no_init)
    //type = 'o' or 'e'
    .def("__init__", make_constructor(WRAP::constructor_from_file))
    .def("__init__", make_constructor(WRAP::constructor_from_inter_event))

    .def(init <char, FrequencyDistribution, DiscreteParametric>())
    .def(init <char, Distribution, DiscreteParametric>())

    // Python Operators
    .def(self_ns::str(self)) //__str__

    .add_property("nb_iterator", &Renewal::get_nb_iterator,"nb iterator")
    .add_property("type", &Renewal::get_type,"type")

    .def("file_ascii_write", WRAP::file_ascii_write,"Save vector summary into a file")

    DEF_RETURN_VALUE_NO_ARGS("get_renewal_data", WRAP::get_renewal_data,"returns renewal data")
    DEF_RETURN_VALUE_NO_ARGS("get_time", WRAP::get_time, "returns time")
    DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("todo"), "simulation")
    .def("extract", WRAP::extract, return_value_policy<manage_new_object> (),  python::args("type", "state"), "Extract distribution data")

    DEF_RETURN_VALUE("simulation_nb_elements", WRAP::simulation_nb_elements, args("todo"), "simulation")
    DEF_RETURN_VALUE("simulation_time_events", WRAP::simulation_time_events, args("Ordinary or Equilibirium", "size", "timev"), "simulation")

    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "Return a plotable")

    ;

/*

   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool spreadsheet_write(StatError &error , const char *path) const;
   bool plot_write(StatError &error , const char *prefix ,  const char *title = 0) const;


	void computation(bool inter_event_flag = true , char itype = 'v' , const Distribution *dtime = 0);
	double likelihood_computation(const TimeEvents &timev) const;

	DiscreteParametric* get_inter_event() const { return inter_event; }
	LengthBias* get_length_bias() const { return length_bias; }
	Backward* get_backward() const { return backward; }
	Forward* get_forward() const { return forward; }
	DiscreteParametric* get_nevent_time(int inb_event) const { return nevent_time[inb_event]; }
	NbEvent* get_nb_event(int itime) const { return nb_event[itime]; }
	Distribution* get_mixture() const { return mixture; }
	Curves* get_index_event() const { return index_event; }
*/
}

#undef WRAP

#define WRAP RenewalDataWrap
class RenewalDataWrap {

public:
  static DiscreteDistributionData*
  extract(const RenewalData &seq, int histo_type, int itime)
  {
    //default itime = I_DEFAULT
    HEADER(DiscreteDistributionData);
    ret = seq.extract(error, histo_type, itime);
    FOOTER;
  }

  //merge inherited from TimeEvents

  static RenewalData*
  merge(const RenewalData &v, const boost::python::list& vecs)
  {
    StatError error;
    RenewalData * ret = NULL;

    int nb_vec = len(vecs);
    sequence_analysis::wrap_util::auto_ptr_array<const RenewalData*> vects(
        new const RenewalData*[nb_vec]);

    for (int i = 0; i < nb_vec; i++)
      vects[i] = boost::python::extract<RenewalData*> (vecs[i]);

    ret = v.merge(error, nb_vec, vects.get());

    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);

    return ret;
  }

  static Renewal*
  estimation(const RenewalData &input, int estimator, int nb_iter,
      int mean_computation_method, double weight, int penalty_type, int outside)
  {
    HEADER_OS(Renewal);
    ret = input.estimation(error, os, estimator, nb_iter, mean_computation_method,
        weight, penalty_type, outside);

    FOOTER_OS;
  }


   static Renewal*
  estimation_inter_event(const RenewalData &input,
      const DiscreteParametric &input_dist, int estimator, int nb_iter,
      int mean_computation_method, double weight, int penalty_type, int outside)
  {
    HEADER_OS(Renewal);

    ret = input.estimation(error, os, input_dist, estimator, nb_iter,
        mean_computation_method, weight, penalty_type, outside);

    FOOTER_OS;
  }




};

void class_renewal_data() {


  class_<RenewalData, bases<TimeEvents> > ("_RenewalData", "RenewalData", no_init)
    .def(init <int, int>())
    .def(init <int, Renewal>())
    .def(init <TimeEvents, int>())
    .def(init <RenewalData, boost::python::optional<bool> >())

    // Python Operators
    .def("get_renewal", &RenewalData::get_renewal, return_value_policy<manage_new_object> (),"get renewal")
    .def("get_type", &RenewalData::get_type, "get type")
    .def("extract", WRAP::extract, return_value_policy<manage_new_object> (),  python::args("type", "state"), "Extract distribution data")
    DEF_RETURN_VALUE_NO_ARGS("merge", WRAP::merge, "Merge renewal_data. Type Merge? for more information")
    DEF_RETURN_VALUE("estimation", WRAP::estimation, args(""), "estimation")
    DEF_RETURN_VALUE("estimation_inter_event", WRAP::estimation_inter_event, args(""), "estimation")



/*
    RenewalData(int nb_sample , const RenewalData **itimev);


   int get_length(int index_seq) const { return length[index_seq]; }
   int get_sequence(int index_seq , int index) const
   { return sequence[index_seq][index]; }
   FrequencyDistribution* get_inter_event() const { return inter_event; }
   FrequencyDistribution* get_within() const { return within; }
   FrequencyDistribution* get_length_bias() const { return length_bias; }
   FrequencyDistribution* get_backward() const { return backward; }
   FrequencyDistribution* get_forward() const { return forward; }
   Curves* get_index_event() const { return index_event; }
   };
*/
;

}



class RenewalIteratorWrap
{

public:

  static boost::python::list
  simulation(RenewalIterator &input, int nb_sequence = 1, char type = 'v')
  {
    int *sequence;

    input.simulation(nb_sequence, type);
    boost::python::list output_sequence;

    for (int i=0; i < input.get_length(); i++)
    {
        output_sequence.append(input.get_sequence(i));
    }
    return output_sequence;
 }
};

void
class_renewal_iterator()
{

  class_<RenewalIterator >
  ("_RenewalIterator", "RenewalIterator", no_init)

  .def(init<Renewal*>()[with_custodian_and_ward_postcall<1, 2>()])

  .add_property("interval", &RenewalIterator::get_interval)
  .add_property("length", &RenewalIterator::get_length)
  .add_property("counter", &RenewalIterator::get_counter)

  .def("get_sequence", &RenewalIterator::get_sequence, args("index"))  // to be done
  .def("simulation", RenewalIteratorWrap::simulation,  "simulation")
;
}

