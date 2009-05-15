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
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(extract_overloads, TimeEventsWrap::extract, 1, 2);



  static Histogram*
  get_hmixture(const Time_events& input)
  {
    Histogram* ret;
    ret = new Histogram(*input.get_hmixture());
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
  merge(Time_events& input_seq, const boost::python::list& seqs)
   {
     int nb_seq = len(seqs);
     sequence_analysis::wrap_util::auto_ptr_array<const Time_events *> sequens(
         new const Time_events*[nb_seq]);

     for (int i = 0; i < nb_seq; i++)
       sequens[i] = boost::python::extract<Time_events *> (seqs[i]);

       //todo does not compile because merge is protected
     //Time_events *ret;
     //ret = input_seq.merge(nb_seq, sequens.get());
     //return ret;
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
    DEF_RETURN_VALUE_NO_ARGS("get_hmixture", &TimeEventsWrap::get_hmixture, "returns hmixture Mixture histogram")
    DEF_RETURN_VALUE("get_hnb_event", &TimeEventsWrap::get_hnb_event, args("itime"),"returns hmixture Mixture histogram")

    .def("extract", (Distribution_data *(*)(const Time_events&, int, int))  TimeEventsWrap::extract, return_value_policy< manage_new_object >(), TimeEventsWrap::extract_overloads(), "todo")
    .def("extract", (Distribution_data *(*)(const Time_events& , int))  TimeEventsWrap::extract, return_value_policy< manage_new_object >(), TimeEventsWrap::extract_overloads(), "todp")

    .def("file_ascii_write", TimeEventsWrap::file_ascii_write,"Save vector summary into a file")

    DEF_RETURN_VALUE("time_scaling", TimeEventsWrap::time_scaling, args("scaling"),"returns a time-scaled TimeEvents")
    DEF_RETURN_VALUE("time_select", TimeEventsWrap::time_select, args("min index", "max index" ),"returns a time-selectd TimeEvents")
    DEF_RETURN_VALUE("nb_event_select", TimeEventsWrap::nb_event_select, args("nb_event_select"),"returns a nb_event-selected TimeEvents")
    DEF_RETURN_VALUE_NO_ARGS("merge", TimeEventsWrap::merge, "Merge sequences")

    ;

/*protected
      int *time;              // temps d'observation
      int *nb_event;          // nombre d'evenements
      int *frequency;         // effectif de chacune des classes
                              // {temps, nombre d'evenements}
      void merge(int nb_sample , const Time_events **ptimev);

      std::ostream& ascii_file_write(std::ostream &os , bool exhaustive , char type = 'v') const;
      std::ostream& spreadsheet_write(std::ostream &os , char type = 'v') const;

      void nb_element_computation();
      double min_inter_event_computation() const;
      /*
/*
Time_events(int inb_element , int *itime , int *inb_event){ build(inb_element , itime , inb_event); }
Time_events(int nb_sample , const Time_events **ptimev) { merge(nb_sample , ptimev); }
Time_events(const Time_events &timev) { copy(timev); }

  Renewal* estimation(Format_error &error , std::ostream &os , char type ,
                       const Parametric &iinter_event , int estimator = LIKELIHOOD ,
                          int nb_iter = I_DEFAULT , int equilibrium_estimator = COMPLETE_LIKELIHOOD ,
                          int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                          int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
  Renewal* estimation(Format_error &error , std::ostream &os , char type , int estimator = LIKELIHOOD ,
                          int nb_iter = I_DEFAULT , int equilibrium_estimator = COMPLETE_LIKELIHOOD ,
                          int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                          int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;



double information_computation() const;


*/

}



