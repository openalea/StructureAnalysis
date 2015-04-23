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



class TimeEventsWrap
{

public:



  static boost::shared_ptr<TimeEvents>
  time_events_from_file(char *filename)
  {
    StatError error;
    TimeEvents *time_events = NULL;
    time_events = time_events_ascii_read(error, filename);
    if(!time_events)
    {
      sequence_analysis::wrap_util::throw_error(error);
    }
    return boost::shared_ptr<TimeEvents>(time_events);
  }


  static boost::shared_ptr<TimeEvents>
  time_events_from_histogram(FrequencyDistribution &input, int itime)
  {
    StatError error;
    TimeEvents *time_events = NULL;
    time_events = build_time_events(error, input, itime);
    if(!time_events)
    {
      sequence_analysis::wrap_util::throw_error(error);
    }
    return boost::shared_ptr<TimeEvents>(time_events);
  }


  static DiscreteDistributionData*
  extract(const TimeEvents &input, int histo_type, int itime)
  {

    // to finish !!. This function does not work. core dumped,seg fault!!
    StatError error;
    DiscreteDistributionData* ret;
    ret = new DiscreteDistributionData(*input.extract(error, histo_type, itime));
    if (!ret)
       sequence_analysis::wrap_util::throw_error(error);
    return ret;
  }
  
  static FrequencyDistribution*
  get_mixture(const TimeEvents &input)
  {
    FrequencyDistribution *ret;
    ret = new FrequencyDistribution(*input.get_mixture());
    return ret;
   }

  static FrequencyDistribution*
  get_hnb_event(const TimeEvents &input, int index)
  {
    FrequencyDistribution *ret;
    //todo
    //check value index is in time
    //same results as get_hmixture!!
    ret = new FrequencyDistribution(*input.get_hnb_event(index));
    return ret;
  }

  static FrequencyDistribution*
  get_htime(const TimeEvents &input)
  {
    // check that it is a FrequencyDistribution cast or DiscreteDistributionData
    FrequencyDistribution* ret;
    ret = new FrequencyDistribution(*input.get_htime());
    return ret;
  }

  static TimeEvents*
  time_scaling(const TimeEvents &input, int scaling){
    SIMPLE_METHOD_TEMPLATE_1(input, time_scaling, TimeEvents, scaling);
  }

  static TimeEvents*
  time_select(const TimeEvents &input, int min, int max){
    SIMPLE_METHOD_TEMPLATE_1(input, time_select, TimeEvents, min, max);
  }

  static TimeEvents*
  nb_event_select(const TimeEvents &input, int min, int max){
    SIMPLE_METHOD_TEMPLATE_1(input, nb_event_select, TimeEvents, min, max);
  }

  static void
  file_ascii_write(const TimeEvents &d, const char* path, bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
   }


  static TimeEvents*
  merge(const TimeEvents &input, const boost::python::list input_timev)
  {
    HEADER(TimeEvents);
    CREATE_ARRAY(input_timev, const TimeEvents*, timev)

    ret = new TimeEvents(timev_size, timev.get());

    FOOTER;
  }

  static Renewal*
  estimation_type(const TimeEvents &input, char type, int estimator, int nb_iter,
                 int equilibrium_estimator, int mean_computation_method, double weight,
                 int penalty_type, int outside)
  {
    HEADER_OS(Renewal);
    ret = input.estimation(error, os, type, estimator, nb_iter, equilibrium_estimator,
                           mean_computation_method, weight, penalty_type, outside);

    FOOTER_OS;
  }


  static Renewal*
  estimation_inter_event_type(const TimeEvents &input, char type,
      const DiscreteParametric& input_dist, int estimator, int nb_iter,
      int equilibrium_estimator, int mean_computation_method, double weight,
      int penalty_type, int outside)
  {
    HEADER_OS(Renewal);


    ret = input.estimation(error, os, type, input_dist, estimator, nb_iter,
        equilibrium_estimator, mean_computation_method, weight, penalty_type, outside);

    FOOTER_OS;
  }

  static MultiPlotSet* get_plotable(const TimeEvents &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret) ERROR;
    return ret;
  }




};

void class_time_events() {


  class_<TimeEvents, bases<StatInterface> > ("_TimeEvents", "TimeEvents")
    .def("__init__", make_constructor(TimeEventsWrap::time_events_from_file))
    .def("__init__", make_constructor(TimeEventsWrap::time_events_from_histogram))

    .def(init <int>())
    // Python Operators

    .def(self_ns::str(self)) //__str__

    .add_property("nb_element", &TimeEvents::get_nb_element,"nb elements")
    .add_property("nb_class", &TimeEvents::get_nb_class,"nb class")

    DEF_RETURN_VALUE_NO_ARGS("get_htime", &TimeEventsWrap::get_htime, "returns htime histogram")
    DEF_RETURN_VALUE_NO_ARGS("get_mixture", &TimeEventsWrap::get_mixture, "returns mixture Mixture histogram")
    DEF_RETURN_VALUE("get_hnb_event", &TimeEventsWrap::get_hnb_event, args("itime"),"returns hmixture Mixture histogram")

    DEF_RETURN_VALUE("extract", TimeEventsWrap::extract, args("",""), "See ExtractHistogram")

    .def("file_ascii_write", TimeEventsWrap::file_ascii_write,"Save vector summary into a file")

    DEF_RETURN_VALUE("time_scaling", TimeEventsWrap::time_scaling, args("scaling"),"returns a time-scaled TimeEvents")
    DEF_RETURN_VALUE("time_select", TimeEventsWrap::time_select, args("min index", "max index" ),"returns a time-selectd TimeEvents")
    DEF_RETURN_VALUE("nb_event_select", TimeEventsWrap::nb_event_select, args("nb_event_select"),"returns a nb_event-selected TimeEvents")
    DEF_RETURN_VALUE_NO_ARGS("merge", TimeEventsWrap::merge, "Merge tim events")

    DEF_RETURN_VALUE("estimation_type", TimeEventsWrap::estimation_type, args("tobedone"), "estimation")
    DEF_RETURN_VALUE("estimation_inter_event_type", TimeEventsWrap::estimation_inter_event_type, args("tobedone"), "estimation")

    DEF_RETURN_VALUE_NO_ARGS("get_plotable", TimeEventsWrap::get_plotable, "Return a plotable")

    ;

/*protected
      int *time;              // temps d'observation
      int *nb_event;          // nombre d'evenements
      int *frequency;         // effectif de chacune des classes
                              // {temps, nombre d'evenements}
      std::ostream& ascii_file_write(std::ostream &os , bool exhaustive , char type = 'v') const;
      std::ostream& spreadsheet_write(std::ostream &os , char type = 'v') const;

      void nb_element_computation();
      double min_inter_event_computation() const;

      TimeEvents(int inb_element , int *itime , int *inb_event){ build(inb_element , itime , inb_event); }
    TimeEvents(const TimeEvents &timev) { copy(timev); }



double information_computation() const;


*/

}



