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



#ifndef COMPOUND_H
#define COMPOUND_H


#include "stat_tools.h"
#include "distribution.h"


namespace stat_tool {



/****************************************************************
 *
 *  Constants
 */


  const double COMPOUND_THRESHOLD = 0.99999;  // threshold on the cumulative distribution function
                                              // for determining the upper bound of the support

  const double COMPOUND_INIT_PROBABILITY = 0.001;  // threshold for probability initialization
  const double COMPOUND_LIKELIHOOD_DIFF = 1.e-5;  // threshold for stopping EM iterations
  const int COMPOUND_NB_ITER = 10000;    // maximum number of EM iterations
  const double COMPOUND_DIFFERENCE_WEIGHT = 0.5;  // default penalty weight (1st- or 2nd-order difference cases)
  const double COMPOUND_ENTROPY_WEIGHT = 0.1;  // default penalty weight (entropy case)
  const int COMPOUND_COEFF = 10;         // rounding coefficient for the estimator



/****************************************************************
 *
 *  Class definition
 */


  class CompoundData;

  /// \brief Compound distribution

  class Compound : public StatInterface , public Distribution {

    friend class FrequencyDistribution;
    friend class CompoundData;

    friend std::ostream& operator<<(std::ostream &os , const Compound &compound)
    { return compound.ascii_write(os , compound.compound_data , false , false); }

  private :

    CompoundData *compound_data;  ///< pointer on a CompoundData object
    DiscreteParametric *sum_distribution;  ///< sum distribution
    DiscreteParametric *distribution;  ///< basis distribution

    void copy(const Compound &compound , bool data_flag = true);

    std::ostream& ascii_write(std::ostream &os , const CompoundData *compound_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const CompoundData *compound_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const CompoundData *compound_histo) const;
    MultiPlotSet* get_plotable(const CompoundData *compound_histo) const;

    void computation(DiscreteParametric **power_dist , int min_nb_value ,
                     double cumul_threshold , bool sum_flag , bool dist_flag);
    void expectation_step(const FrequencyDistribution &histo , DiscreteParametric **power_dist ,
                          Reestimation<double> *sum_reestim , Reestimation<double> *reestim) const;

  public :

    Compound();
    Compound(const DiscreteParametric &sum_dist , const DiscreteParametric &dist ,
             double cumul_threshold = COMPOUND_THRESHOLD);
    Compound(const DiscreteParametric &sum_dist , const DiscreteParametric &dist ,
             compound_distribution type);
    Compound(const Compound &compound , bool data_flag = true)
    :Distribution(compound) { copy(compound , data_flag); }
    ~Compound();
    Compound& operator=(const Compound &compound);

    CompoundData* extract_data(StatError &error) const;

    static Compound* ascii_read(StatError &error , const std::string path ,
                                double cumul_threshold = COMPOUND_THRESHOLD);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void computation(int min_nb_value = 1 ,
                     double cumul_threshold = COMPOUND_THRESHOLD ,
                     bool sum_flag = true , bool dist_flag = true);
    CompoundData* simulation(StatError &error , int nb_element) const;

    // class member access

    CompoundData* get_compound_data() const { return compound_data; }
    DiscreteParametric* get_sum_distribution() const { return sum_distribution; }
    DiscreteParametric* get_distribution() const { return distribution; }
  };


  /// \brief Data structure corresponding to a compound distribution

  class CompoundData : public StatInterface , public FrequencyDistribution {
    friend class FrequencyDistribution;
    friend class Compound;

    friend std::ostream& operator<<(std::ostream &os , const CompoundData &compound_histo)
    { return compound_histo.ascii_write(os , false); }

  private :

    Compound *compound;     ///< pointer on a Compound object
    FrequencyDistribution *sum_frequency_distribution;   ///< sum frequency distribution
    FrequencyDistribution *frequency_distribution;   ///< basis frequency distribution

    void copy(const CompoundData &compound_histo , bool model_flag = true);

  public :

    CompoundData();
    CompoundData(const FrequencyDistribution &histo , const Compound &icompound);
    CompoundData(const Compound &icompound);
    CompoundData(const CompoundData &compound_histo , bool model_flag = true)
    :FrequencyDistribution(compound_histo) { copy(compound_histo , model_flag); }
    ~CompoundData();
    CompoundData& operator=(const CompoundData &compound_histo);

    DiscreteDistributionData* extract(StatError &error , compound_distribution type) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    // class member access

    Compound* get_compound() const { return compound; }
    FrequencyDistribution* get_sum_frequency_distribution() const { return sum_frequency_distribution; }
    FrequencyDistribution* get_frequency_distribution() const { return frequency_distribution; }
  };


};  // namespace stat_tool



#endif
