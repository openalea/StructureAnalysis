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



#ifndef CONVOLUTION_H
#define CONVOLUTION_H


#include "stat_tools.h"
#include "distribution.h"


namespace stat_tool {



/****************************************************************
 *
 *  Constants
 */


  const int CONVOLUTION_NB_DISTRIBUTION = 10;  // maximum number of elementary distributions
  const double CONVOLUTION_THRESHOLD = 0.9999;  // threshold on the cumulative distribution function
                                                // for determining the upper bound of the support

  const double CONVOLUTION_INIT_PROBABILITY = 0.001;  // threshold for probability initialization
  const double CONVOLUTION_LIKELIHOOD_DIFF = 1.e-5;  // threshold for stopping EM iterations
  const int CONVOLUTION_NB_ITER = 10000;  // maximum number of EM iterations
  const double CONVOLUTION_DIFFERENCE_WEIGHT = 0.5;  // default penalty weight (1st- or 2nd-order difference cases)
  const double CONVOLUTION_ENTROPY_WEIGHT = 0.1;  // default penalty weight (entropy case)
  const int CONVOLUTION_COEFF = 10;      // rounding coefficient for the estimator



/****************************************************************
 *
 *  Class definition
 */


  class ConvolutionData;

  /// \brief Convolution of discrete distributions

  class Convolution : public StatInterface , public Distribution {

    friend class FrequencyDistribution;
    friend class ConvolutionData;

    friend std::ostream& operator<<(std::ostream &os , const Convolution &convol)
    { return convol.ascii_write(os , convol.convolution_data , false , false); }

  private :

    ConvolutionData *convolution_data;  ///< pointer on a ConvolutionData object
    int nb_distribution;    ///< number of elementary distributions
    DiscreteParametric **distribution;  ///< elementary distributions

    void copy(const Convolution &convol , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const ConvolutionData *convol_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const ConvolutionData *convol_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const ConvolutionData *convol_histo) const;
    MultiPlotSet* get_plotable(const ConvolutionData *convol_histo) const;

    void expectation_step(const FrequencyDistribution &histo , const Distribution **partial_convol ,
                          Reestimation<double> **reestim) const;

  public :

    Convolution();
    Convolution(int nb_dist , const DiscreteParametric **pdist);
    Convolution(int nb_dist , const std::vector<DiscreteParametric> idist);
    Convolution(const DiscreteParametric &known_dist , const DiscreteParametric &unknown_dist);
    Convolution(const Convolution &convol , bool data_flag = true)
    :Distribution(convol) { copy(convol , data_flag); }
    ~Convolution();
    Convolution& operator=(const Convolution &convol);

    DiscreteParametricModel* extract(StatError &error , int index) const;
    ConvolutionData* extract_data(StatError &error) const;

    static Convolution* building(StatError &error , int nb_dist ,
                                 const DiscreteParametric **dist);
    static Convolution* building(StatError &error , int nb_dist ,
                                 const std::vector<DiscreteParametric> dist);
    static Convolution* ascii_read(StatError &error , const std::string path ,
                                   double cumul_threshold = CONVOLUTION_THRESHOLD);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void computation(int min_nb_value = 1 , double cumul_threshold = CONVOLUTION_THRESHOLD ,
                     bool *dist_flag = NULL);
    ConvolutionData* simulation(StatError &error , int nb_element) const;

    // class member access

    ConvolutionData* get_convolution_data() const { return convolution_data; }
    int get_nb_distribution() const { return nb_distribution; }
    DiscreteParametric* get_distribution(int index) const { return distribution[index]; }
  };


  /// \brief Data structure corresponding to a convolution of discrete distributions

  class ConvolutionData : public StatInterface , public FrequencyDistribution {

    friend class FrequencyDistribution;
    friend class Convolution;

    friend std::ostream& operator<<(std::ostream &os , const ConvolutionData &convol_histo)
    { return convol_histo.ascii_write(os , false); }

  private :

    Convolution *convolution;  ///< pointer on a Convolution object
    int nb_distribution;       ///< number of elementary frequency distributions
    FrequencyDistribution **frequency_distribution;  ///< elementary frequency distributions

    void copy(const ConvolutionData &convol_histo , bool model_flag = true);
    void remove();

  public :

    ConvolutionData();
    ConvolutionData(const FrequencyDistribution &histo , int nb_dist);
    ConvolutionData(const Convolution &convol);
    ConvolutionData(const ConvolutionData &convol_histo , bool model_flag = true)
    :FrequencyDistribution(convol_histo) { copy(convol_histo , model_flag); }
    ~ConvolutionData();
    ConvolutionData& operator=(const ConvolutionData &convol_histo);

    DiscreteDistributionData* extract(StatError &error , int index) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    // class member access

    Convolution* get_convolution() const { return convolution; }
    int get_nb_distribution() const { return nb_distribution; }
    FrequencyDistribution* get_frequency_distribution(int index) const { return frequency_distribution[index]; }
  };


};  // namespace stat_tool



#endif
