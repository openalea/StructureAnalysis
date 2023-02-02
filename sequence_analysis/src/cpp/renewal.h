/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for V-Plants developers:
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */



#ifndef RENEWAL_H
#define RENEWAL_H


#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"


namespace sequence_analysis {



/****************************************************************
 *
 *  Constants
 */


  const int DEFAULT_TIME = 20;           // default observation period
  const int MAX_TIME = 500;              // maximum observation period
  const int PLOT_NEVENT_TIME = 10;       // maximum number of time to the nth event distributions
                                         // plotted (Gnuplot output)
  const int PLOT_NB_TIME = 5;            // maximum number of distributions of the number of events plotted
                                         // with the mixture of the number of events distributions (Gnuplot output)

  const double RENEWAL_THRESHOLD = 0.99999;  // threshold on the cumulative distribution function for determining
                                             // the upper bound of the support of the inter-event distribution
  const double RB_THRESHOLD = 2000.;     // threshold for using the fast computation of the number of events distribution
                                         // from a binomial inter-event distribution
  const double RNB_THRESHOLD = 2000.;    // threshold for using the fast computation of the number of events distribution
                                         // from a negative binomiale inter-event distribution

  enum renewal_distribution {
    INTER_EVENT ,
    WITHIN_OBSERVATION_PERIOD ,
    LENGTH_BIAS ,
    BACKWARD_RECURRENCE_TIME ,
    FORWARD_RECURRENCE_TIME ,
    NB_EVENT ,
    NB_EVENT_MIXTURE
  };

  const double MIN_NB_EVENT = 0.4;       // minimum mean number of events
  const double MIN_INTER_EVENT = 1.;     // minimum mean time interval between events
  const double RENEWAL_INIT_PROBABILITY = 0.001;  // threshold for probability initialization
  const int RENEWAL_COEFF = 10;          // rounding coefficient for the estimator

  const double MEAN_COEFF = 2.;          // coefficient on the mean for compensating the length bias

  const int RENEWAL_NB_ELEMENT = 1000000;  // maximum sample size for simulation



/****************************************************************
 *
 *  Class definition
 */


  /// \brief Length-biased distribution

  class LengthBias : public stat_tool::DiscreteParametric {

  public :

    LengthBias(int inb_value = 0 , stat_tool::discrete_parametric iident = stat_tool::CATEGORICAL ,
               int iinf_bound = stat_tool::I_DEFAULT , int isup_bound = stat_tool::I_DEFAULT ,
               double iparameter = stat_tool::D_DEFAULT, double iprobability = stat_tool::D_DEFAULT)
    :stat_tool::DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    LengthBias(const DiscreteParametric &inter_event)
    :stat_tool::DiscreteParametric(inter_event) { computation(inter_event); }
    LengthBias(const LengthBias &length_bias)
    :stat_tool::DiscreteParametric((DiscreteParametric&)length_bias) {}

    void computation(const stat_tool::DiscreteParametric&);
  };


  /// \brief Backward recurrence time distribution

  class Backward : public stat_tool::DiscreteParametric {

  public :

    Backward(int inb_value = 0 , stat_tool::discrete_parametric iident = stat_tool::CATEGORICAL ,
             int iinf_bound = stat_tool::I_DEFAULT , int isup_bound = stat_tool::I_DEFAULT ,
             double iparameter = stat_tool::D_DEFAULT, double iprobability = stat_tool::D_DEFAULT)
    :stat_tool::DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Backward(const Backward &dist , int ialloc_nb_value = stat_tool::I_DEFAULT)
    :stat_tool::DiscreteParametric(dist , stat_tool::DISTRIBUTION_COPY , ialloc_nb_value) {}

    void computation(const stat_tool::DiscreteParametric &inter_event , const stat_tool::Distribution &time);
  };


  /// \brief Number of events distribution

  class NbEvent : public stat_tool::DiscreteParametric {

  public :

    stat_tool::process_type type;  ///< renewal process type (ORDINARY/EQUILIBRIUM)
    int time;               ///< observation period

    NbEvent(stat_tool::process_type itype = stat_tool::EQUILIBRIUM , int itime = 0 , int inb_value = 0 ,
            stat_tool::discrete_parametric iident = stat_tool::CATEGORICAL ,
            int iinf_bound = stat_tool::I_DEFAULT , int isup_bound = stat_tool::I_DEFAULT ,
            double iparameter = stat_tool::D_DEFAULT , double iprobability = stat_tool::D_DEFAULT);
    NbEvent(stat_tool::process_type itype , int itime , stat_tool::DiscreteParametric &inter_event);
    NbEvent(const NbEvent &nb_event , int ialloc_nb_value = stat_tool::I_DEFAULT);

    void binomial_computation();
    void negative_binomial_computation();

    void ordinary_computation(stat_tool::DiscreteParametric &inter_event);
    void computation(stat_tool::DiscreteParametric &inter_event);
  };


  class RenewalIterator;
  class TimeEvents;
  class RenewalData;

  /// \brief Renewal process

  class Renewal : public stat_tool::StatInterface {

    friend class RenewalIterator;
    friend class TimeEvents;
    friend class RenewalData;

    friend std::ostream& operator<<(std::ostream &os , const Renewal &renew)
    { return renew.ascii_write(os , renew.renewal_data , false , false); }

  private :

    int nb_iterator;        ///< number of iterators pointing on the Renewal object
    RenewalData *renewal_data;  ///< pointer on a RenewalData object
    stat_tool::process_type type;  ///< renewal process type (ORDINARY/EQUILIBRIUM)
    int nb_event_max;       ///< maximum number of events
    stat_tool::Distribution *time;  ///< observation period distribution
    stat_tool::DiscreteParametric *inter_event;  ///< inter-event distribution
    LengthBias *length_bias;  ///< length-biased distribution
    Backward *backward;     ///< backward recurrence time distribution
    stat_tool::Forward *forward;  ///< forward recurrence time distribution
    stat_tool::DiscreteParametric **nevent_time;  ///< time to the nth event distributions
    NbEvent **nb_event;    ///< number of events distributions for the different observation periods
    stat_tool::Distribution *mixture;  ///< mixture of the number of events distributions
    stat_tool::Curves *index_event;  ///< no-event/event probabilities as a function of time

    void init(int inf_bound , int sup_bound , double parameter , double probability);
    void init(stat_tool::discrete_parametric ident , int inf_bound , int sup_bound ,
              double parameter , double probability);
    void copy(const Renewal &renew , bool data_flag = true);
    void remove();
    void type_init(stat_tool::process_type itype);

    std::ostream& ascii_write(std::ostream &os , const RenewalData *timev ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const RenewalData *timev) const;
    bool plot_write(const char *prefix , const char *title ,
                    const RenewalData *timev) const;
    stat_tool::MultiPlotSet* get_plotable(const RenewalData *timev) const;

    void index_event_computation();

    void expectation_step(const TimeEvents &timev , stat_tool::Reestimation<double> *reestim) const;
    void expectation_step(const TimeEvents &timev , stat_tool::Reestimation<double> *inter_event_reestim ,
                          stat_tool::Reestimation<double> *length_bias_reestim ,
                          stat_tool::censoring_estimator estimator , bool combination = false ,
                          stat_tool::duration_distribution_mean_estimator mean_estimator = stat_tool::COMPUTED) const;

  public :

    Renewal();
    Renewal(stat_tool::process_type itype , const stat_tool::FrequencyDistribution &htime ,
            const stat_tool::DiscreteParametric &iinter_event);
    Renewal(stat_tool::process_type itype , const stat_tool::Distribution &itime ,
            const stat_tool::DiscreteParametric &iinter_event);
    Renewal(const RenewalData &irenewal_data ,
            const stat_tool::DiscreteParametric &iinter_event);
    Renewal(const Renewal &renew , bool data_flag = true)
    { copy(renew , data_flag); }
    ~Renewal();
    void conditional_delete();
    Renewal& operator=(const Renewal &renew);

    stat_tool::DiscreteParametricModel* extract(stat_tool::StatError &error ,
                                                renewal_distribution dist_type ,
                                                int itime = stat_tool::I_DEFAULT) const;

    static Renewal* building(stat_tool::StatError &error , const stat_tool::DiscreteParametric &inter_event ,
                             stat_tool::process_type type = stat_tool::EQUILIBRIUM , int time = DEFAULT_TIME);

    static Renewal* ascii_read(stat_tool::StatError& error , const std::string path ,
                               stat_tool::process_type type = stat_tool::EQUILIBRIUM , int time = DEFAULT_TIME ,
                               double cumul_threshold = RENEWAL_THRESHOLD);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const std::string path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    void computation(bool inter_event_flag = true , stat_tool::process_type itype = stat_tool::DEFAULT_TYPE ,
                     const stat_tool::Distribution *dtime = NULL);

    double likelihood_computation(const TimeEvents &timev) const;

    RenewalData* simulation(stat_tool::StatError &error , stat_tool::process_type itype ,
                            const stat_tool::FrequencyDistribution &ihtime) const;
    RenewalData* simulation(stat_tool::StatError &error , stat_tool::process_type itype ,
                            int nb_element , int itime) const;
    RenewalData* simulation(stat_tool::StatError &error , stat_tool::process_type itype ,
                            int nb_element , const TimeEvents &itimev) const;

    // class member access

    int get_nb_iterator() const { return nb_iterator; }
    RenewalData* get_renewal_data() const { return renewal_data; }
    stat_tool::process_type get_type() const { return type; }
    stat_tool::Distribution* get_time() const { return time; }
    stat_tool::DiscreteParametric* get_inter_event() const { return inter_event; }
    LengthBias* get_length_bias() const { return length_bias; }
    Backward* get_backward() const { return backward; }
    stat_tool::Forward* get_forward() const { return forward; }
    stat_tool::DiscreteParametric* get_nevent_time(int inb_event) const { return nevent_time[inb_event]; }
    NbEvent* get_nb_event(int itime) const { return nb_event[itime]; }
    stat_tool::Distribution* get_mixture() const { return mixture; }
    stat_tool::Curves* get_index_event() const { return index_event; }
  };


  /// \brief Renewal process iterator for asynchronous simulation

  class RenewalIterator {

  private :

    Renewal *renewal;       ///< pointer on a Renewal object
    int interval;           ///< time interval between events
    int counter;            ///< counter
    int length;             ///< sequence length
    int *sequence;          ///< sequence of events

    void copy(const RenewalIterator &iter);

  public :

    RenewalIterator(Renewal *irenewal , int ilength = 1);
    RenewalIterator(const RenewalIterator &iter)
    { copy(iter); }
    ~RenewalIterator();
    RenewalIterator& operator=(const RenewalIterator &iter);

    void simulation(int ilength = 1 , stat_tool::process_type type = stat_tool::DEFAULT_TYPE);

    // class member access

    Renewal* get_renewal() const { return renewal; }
    int get_interval() const { return interval; }
    int get_counter() const { return counter; }
    int get_length() const { return length; }
    int get_sequence(int index) const { return sequence[index]; }
  };


  /// \brief Triplets {observation period, number of events, frequency}

  class TimeEvents : public stat_tool::StatInterface {

    friend class stat_tool::FrequencyDistribution;
    friend class Renewal;

    friend std::ostream& operator<<(std::ostream &os , const TimeEvents &timev)
    { return timev.ascii_write(os , true); }

  protected :

    int nb_element;         ///< sample size
    int nb_class;           ///< number of classes
    int *time;              ///< observation period
    int *nb_event;          ///< number of events
    int *frequency;         ///< frequency of each pair {observation period, number of events}
    stat_tool::FrequencyDistribution *htime;  ///< observation period frequency distribution
    stat_tool::FrequencyDistribution **hnb_event;  ///< number of events frequency distributions for the different observation periods
    stat_tool::FrequencyDistribution *mixture;  ///< mixture of the number of events frequency distributions

    void build_frequency_distribution();
    void build_sample();
    void build(int inb_element , int *itime , int *inb_event);
    void copy(const TimeEvents&);
    void merge(int nb_sample , const TimeEvents **ptimev);
    void remove();

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , stat_tool::process_type type) const;
    std::ostream& ascii_file_write(std::ostream &os , bool exhaustive ,
                                   stat_tool::process_type type = stat_tool::EQUILIBRIUM) const;
    std::ostream& spreadsheet_write(std::ostream &os ,
                                    stat_tool::process_type type = stat_tool::EQUILIBRIUM) const;

    void nb_element_computation();
    double min_inter_event_computation() const;

  public :

    TimeEvents(int inb_class = 0);
    TimeEvents(int inb_element , int *itime , int *inb_event)
    { build(inb_element , itime , inb_event); }
    TimeEvents(stat_tool::FrequencyDistribution &inb_event , int itime);
    TimeEvents(int nb_sample , const TimeEvents **ptimev) { merge(nb_sample , ptimev); }
    TimeEvents(const TimeEvents &timev) { copy(timev); }
    ~TimeEvents();
    TimeEvents& operator=(const TimeEvents &timev);

    TimeEvents* merge(int nb_sample , const std::vector<TimeEvents> itimev) const;
    stat_tool::DiscreteDistributionData* extract(stat_tool::StatError &error ,
                                                 renewal_distribution histo_type ,
                                                 int itime = stat_tool::I_DEFAULT) const;

    TimeEvents* time_scaling(stat_tool::StatError &error , int scaling_coeff) const;
    TimeEvents* time_select(stat_tool::StatError &error , int min_time ,
                            int max_time) const;
    TimeEvents* nb_event_select(stat_tool::StatError &error , int min_nb_event ,
                                int max_nb_event) const;

    static TimeEvents* building(stat_tool::StatError &error , stat_tool::FrequencyDistribution &nb_event , int itime);

    static TimeEvents* ascii_read(stat_tool::StatError &error , const std::string path);
    static TimeEvents* old_ascii_read(stat_tool::StatError &error , const std::string path);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = true) const;
    bool ascii_write(stat_tool::StatError &error , const std::string path , bool exhaustive = true) const;
    bool spreadsheet_write(stat_tool::StatError &error , const std::string path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    double information_computation() const;

    Renewal* estimation(stat_tool::StatError &error , bool display , stat_tool::process_type type ,
                        const stat_tool::DiscreteParametric &iinter_event ,
                        stat_tool::estimation_criterion estimator = stat_tool::LIKELIHOOD ,
                        int nb_iter = stat_tool::I_DEFAULT ,
                        stat_tool::censoring_estimator equilibrium_estimator = stat_tool::COMPLETE_LIKELIHOOD ,
                        stat_tool::duration_distribution_mean_estimator mean_estimator = stat_tool::COMPUTED ,
                        double weight = stat_tool::D_DEFAULT ,
                        stat_tool::penalty_type pen_type = stat_tool::SECOND_DIFFERENCE ,
                        stat_tool::side_effect outside = stat_tool::ZERO) const;
    Renewal* estimation(stat_tool::StatError &error , bool display , stat_tool::process_type type ,
                        stat_tool::estimation_criterion estimator = stat_tool::LIKELIHOOD ,
                        int nb_iter = stat_tool::I_DEFAULT ,
                        stat_tool::censoring_estimator equilibrium_estimator = stat_tool::COMPLETE_LIKELIHOOD ,
                        stat_tool::duration_distribution_mean_estimator mean_estimator = stat_tool::COMPUTED ,
                        double weight = stat_tool::D_DEFAULT ,
                        stat_tool::penalty_type pen_type = stat_tool::SECOND_DIFFERENCE ,
                        stat_tool::side_effect outside = stat_tool::ZERO) const;

    // class member access

    int get_nb_element() const { return nb_element; }
    int get_nb_class() const { return nb_class; }
    stat_tool::FrequencyDistribution* get_htime() const { return htime; }
    stat_tool::FrequencyDistribution* get_hnb_event(int itime) const { return hnb_event[itime]; }
    stat_tool::FrequencyDistribution* get_mixture() const { return mixture; }
  };


  /// \brief Data structure corresponding to a renewal process

  class RenewalData : public TimeEvents {

    friend class Renewal;
    friend class Sequences;

    friend std::ostream& operator<<(std::ostream &os , RenewalData &timev)
    { return timev.ascii_write(os , false); }

  private :

    Renewal *renewal;       ///< pointer on a Renewal object
    stat_tool::process_type type;  ///< renewal process type (ORDINARY/EQUILIBRIUM)
    int *length;            ///< sequence length
    int **sequence;         ///< sequences of events
    stat_tool::FrequencyDistribution *inter_event; ///< inter-event frequency distribution
    stat_tool::FrequencyDistribution *within;  ///< frequency distribution of time intervals between events within the observation period
    stat_tool::FrequencyDistribution *length_bias;  ///< length-biased frequency distribution
    stat_tool::FrequencyDistribution *backward;  ///< backward recurrence time frequency distribution
    stat_tool::FrequencyDistribution *forward;  ///< forward recurrence time frequency distribution
    stat_tool::Curves *index_event;  ///< empirical no-event/event probabilities as a function of time

    void copy(const RenewalData &timev , bool model_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os) const;

    void build_index_event(int offset = 1);

  public :

    RenewalData();
    RenewalData(int nb_element , int itime);
    RenewalData(stat_tool::process_type itype , const Renewal &renew);
    RenewalData(const TimeEvents &timev , stat_tool::process_type itype);
    RenewalData(int nb_sample , const RenewalData **itimev);
    RenewalData(const RenewalData &timev , bool model_flag = true)
    :TimeEvents(timev) { copy(timev , model_flag); }
    ~RenewalData();
    RenewalData& operator=(const RenewalData&);

    RenewalData* merge(stat_tool::StatError &error , int nb_sample , const RenewalData **itimev) const;
    RenewalData* merge(stat_tool::StatError &error , int nb_sample , const std::vector<RenewalData> itimev) const;
    stat_tool::DiscreteDistributionData* extract(stat_tool::StatError &error ,
                                                 renewal_distribution histo_type ,
                                                 int itime = stat_tool::I_DEFAULT) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const std::string path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    Renewal* estimation(stat_tool::StatError &error , bool display ,
                        const stat_tool::DiscreteParametric &iinter_event ,
                        stat_tool::estimation_criterion estimator = stat_tool::LIKELIHOOD ,
                        int nb_iter = stat_tool::I_DEFAULT ,
                        stat_tool::duration_distribution_mean_estimator mean_estimator = stat_tool::COMPUTED ,
                        double weight = stat_tool::D_DEFAULT ,
                        stat_tool::penalty_type pen_type = stat_tool::SECOND_DIFFERENCE ,
                        stat_tool::side_effect outside = stat_tool::ZERO) const;
    Renewal* estimation(stat_tool::StatError &error , bool display ,
                        stat_tool::estimation_criterion estimator = stat_tool::LIKELIHOOD ,
                        int nb_iter = stat_tool::I_DEFAULT ,
                        stat_tool::duration_distribution_mean_estimator mean_estimator = stat_tool::COMPUTED ,
                        double weight = stat_tool::D_DEFAULT ,
                        stat_tool::penalty_type pen_type = stat_tool::SECOND_DIFFERENCE ,
                        stat_tool::side_effect outside = stat_tool::ZERO) const;

    // class member access

    Renewal* get_renewal() const { return renewal; }
    stat_tool::process_type get_type() const { return type; }
    int get_length(int index_seq) const { return length[index_seq]; }
    int get_sequence(int index_seq , int index) const
    { return sequence[index_seq][index]; }
    stat_tool::FrequencyDistribution* get_inter_event() const { return inter_event; }
    stat_tool::FrequencyDistribution* get_within() const { return within; }
    stat_tool::FrequencyDistribution* get_length_bias() const { return length_bias; }
    stat_tool::FrequencyDistribution* get_backward() const { return backward; }
    stat_tool::FrequencyDistribution* get_forward() const { return forward; }
    stat_tool::Curves* get_index_event() const { return index_event; }
  };


};  // namespace sequence_analysis




#endif
