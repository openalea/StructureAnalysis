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
//using namespace sequence_analysis;



class TimeEventsWrap
{

public:



  static boost::shared_ptr<Time_events>
  time_events_from_file(char* filename)
  {
    Format_error error;
    Time_events *time_events = NULL;
    time_events = time_events_ascii_read(error, filename);
    if(!time_events)
    {
      sequence_analysis::wrap_util::throw_error(error);
    }
    return boost::shared_ptr<Time_events>(time_events);
  }


  static boost::shared_ptr<Time_events>
  time_events_from_histogram(const Histogram& input, int itime)
  {
    Format_error error;
    Time_events *time_events = NULL;
    time_events = input.build_time_events(error, itime);
    if(!time_events)
    {
      sequence_analysis::wrap_util::throw_error(error);
    }
    return boost::shared_ptr<Time_events>(time_events);
  }


  static Distribution_data*
  extract(const Time_events& input, int histo_type, int itime)
  {

    // to finish !!. This function does not work. core dumped,seg fault!!
    Format_error error;
    Distribution_data* ret;
    ret = new Distribution_data(*input.extract(error, histo_type, itime));
    if (!ret)
       sequence_analysis::wrap_util::throw_error(error);
    return ret;
  }
  
  static Histogram*
  get_mixture(const Time_events& input)
  {
    Histogram* ret;
    ret = new Histogram(*input.get_mixture());
    return ret;
   }

  static Histogram*
  get_hnb_event(const Time_events& input, int index)
  {
    Histogram* ret;
    //todo
    //check value index is in time
    //same results as get_hmixture!!
    ret = new Histogram(*input.get_hnb_event(index));
    return ret;
  }

  static Histogram*
  get_htime(const Time_events& input)
  {
    // check that it is a Histogram cast or Distribution_data
    Histogram* ret;
    ret = new Histogram(*input.get_htime());
    return ret;
  }

  static Time_events*
  time_scaling(const Time_events& input, int scaling){
    SIMPLE_METHOD_TEMPLATE_1(input, time_scaling, Time_events, scaling);
  }

  static Time_events*
  time_select(const Time_events& input, int min, int max){
    SIMPLE_METHOD_TEMPLATE_1(input, time_select, Time_events, min, max);
  }

  static Time_events*
  nb_event_select(const Time_events& input, int min, int max){
    SIMPLE_METHOD_TEMPLATE_1(input, nb_event_select, Time_events, min, max);
  }

  static void
  file_ascii_write(const Time_events& d, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
   }


  static Time_events*
  merge(const Time_events& input, const boost::python::list input_timev)
  {
    HEADER(Time_events);
    CREATE_ARRAY(input_timev, const Time_events *, timev)

    ret = new Time_events(timev_size, timev.get());

    FOOTER;
  }

  static Renewal*
  estimation_type(const Time_events& input, char type, int estimator, int nb_iter,
                 int equilibrium_estimator, int mean_computation, double weight,
                 int penalty_type, int outside)
  {
    HEADER_OS(Renewal);
    ret = input.estimation(error, os, type, estimator, nb_iter, equilibrium_estimator,
                          mean_computation, weight, penalty_type, outside);

    FOOTER_OS;
  }


  static Renewal*
  estimation_inter_event_type(const Time_events& input, char type,
      const Parametric& input_dist, int estimator, int nb_iter,
      int equilibrium_estimator, int mean_computation, double weight,
      int penalty_type, int outside)
  {
    HEADER_OS(Renewal);


    ret = input.estimation(error, os, type, input_dist, estimator, nb_iter,
        equilibrium_estimator, mean_computation, weight, penalty_type, outside);

    FOOTER_OS;
  }

  static MultiPlotSet* get_plotable(const Time_events& p)
  {
    Format_error error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret) ERROR;
    return ret;
  }




};

void class_time_events() {


  class_<Time_events, bases<STAT_interface> > ("_Time_events", "Time_events")
    .def("__init__", make_constructor(TimeEventsWrap::time_events_from_file))
    .def("__init__", make_constructor(TimeEventsWrap::time_events_from_histogram))

    .def(init <int>())
    // Python Operators

    .def(self_ns::str(self)) //__str__

    .add_property("nb_element", &Time_events::get_nb_element,"nb elements")
    .add_property("nb_class", &Time_events::get_nb_class,"nb class")

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

      Time_events(int inb_element , int *itime , int *inb_event){ build(inb_element , itime , inb_event); }
    Time_events(const Time_events &timev) { copy(timev); }



double information_computation() const;


*/

}



