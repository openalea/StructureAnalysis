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
 *       $Id: mixture.h 14946 2013-09-30 12:29:22Z guedon $
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



#ifndef MIXTURE_H
#define MIXTURE_H


#include "vectors.h"


namespace stat_tool {



/****************************************************************
 *
 *  Constants
 */


  const int MIXTURE_NB_COMPONENT = NB_STATE;  // maximum number of components

  const double MIXTURE_LIKELIHOOD_DIFF = 1.e-6;  // threshold for stopping the EM iterations
  const int MIXTURE_NB_ITER = 500;        // maximum number of EM iterations
  const int POSTERIOR_PROBABILITY_NB_VECTOR = 300;  // maximum number of individuals for the output of
                                                    // the posterior probabilities of the most probable assignments

  const int MIXTURE_NB_VECTOR = 50000;   // maximum sample size for simulation



/****************************************************************
 *
 *  Class definition
 */


  class MixtureData;

  /// \brief Multivariate mixture of distributions

  class Mixture : public StatInterface {

    friend class Vectors;
    friend class MixtureData;

    friend std::ostream& operator<<(std::ostream &os , const Mixture &mixt)
    { return mixt.ascii_write(os , mixt.mixture_data); }

  private :

    MixtureData *mixture_data;  ///< pointer on a MixtureData object
    int nb_component;       ///< number of components
    DiscreteParametric *weight;  ///< weight distribution
//    int explanatory_variable;  categorical explanatory variable for the weights
//    DiscreteParametric **category_weight;  component weights for the different categories
    int nb_output_process;  ///< number of observation processes
    CategoricalProcess **categorical_process;  ///< categorical observation processes
    DiscreteParametricProcess **discrete_parametric_process;  ///< discrete parametric observation processes
    ContinuousParametricProcess **continuous_parametric_process;  ///< continuous parametric observation processes

    Mixture(const DiscreteParametric *iweight , int inb_output_process ,
            CategoricalProcess **categorical_observation ,
            DiscreteParametricProcess **discrete_parametric_observation ,
            ContinuousParametricProcess **continuous_parametric_observation);

    void copy(const Mixture &mixt , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const MixtureData *vec ,
                              bool exhaustive = false , bool file_flag  = false) const;
    std::ostream& spreadsheet_write(std::ostream &os , const MixtureData *vec) const;
    bool plot_write(const char *prefix , const char *title ,
                    const MixtureData *vec) const;
    MultiPlotSet* get_plotable(const MixtureData *vec) const;

    int nb_parameter_computation(double min_probability = 0.) const;
    void individual_assignment(MixtureData &vec , bool assignment = true) const;

  public :

    Mixture();
    Mixture(int inb_component , int inb_output_process , int *nb_value);
    Mixture(int inb_component , double offset , double mean , double standard_deviation ,
            bool common_dispersion);
    Mixture(int inb_component , int ident , double mean , double standard_deviation ,
            bool tied_mean , tying_rule variance_factor);
    Mixture(const Mixture &mixt , bool data_flag = true)
    { copy(mixt , data_flag); }
    ~Mixture();
    Mixture& operator=(const Mixture &mixt);

    DiscreteParametricModel* extract(StatError &error , int variable , int index) const;
    MixtureData* extract_data(StatError &error) const;

    Mixture* thresholding(double min_probability = MIN_PROBABILITY) const;

    static Mixture* ascii_read(StatError &error , const std::string path ,
                               double cumul_threshold = CUMUL_THRESHOLD);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    double classification_likelihood_computation(const MixtureData &vec) const;
    double likelihood_computation(const Vectors &vec , int index = I_DEFAULT) const;

    MixtureData* simulation(StatError &error , int nb_vector) const;

    // class member access

    MixtureData* get_mixture_data() const { return mixture_data; }
    int get_nb_component() const { return nb_component; }
    DiscreteParametric* get_weight() const { return weight; }
    int get_nb_output_process() const { return nb_output_process; }
    CategoricalProcess** get_categorical_process()
    const { return categorical_process; }
    CategoricalProcess* get_categorical_process(int variable)
    const { return categorical_process[variable]; }
    DiscreteParametricProcess** get_discrete_parametric_process() const
    { return discrete_parametric_process; }
    DiscreteParametricProcess* get_discrete_parametric_process(int variable) const
    { return discrete_parametric_process[variable]; }
    ContinuousParametricProcess** get_continuous_parametric_process() const
    { return continuous_parametric_process; }
    ContinuousParametricProcess* get_continuous_parametric_process(int variable) const
    { return continuous_parametric_process[variable]; }
  };


  /// \brief Data structure corresponding to a multivariate mixture of distributions

  class MixtureData : public Vectors {

    friend class Vectors;
    friend class Mixture;

    friend std::ostream& operator<<(std::ostream &os , const MixtureData &vec)
    { return vec.ascii_write(os , false); }

  private :

    Mixture *mixture;       ///< pointer on a Mixture object
//    int explanatory_variable;  categorical explanatory variable for the weights
//    FrequencyDistribution **category_weight;  component weights for the different categories
    FrequencyDistribution ***observation_distribution;  ///< observation frequency distributions
    Histogram ***observation_histogram;  ///< observation histograms
    double likelihood;      ///< log-likelihood of the observed data
    double restoration_likelihood;  ///< log-likelihood of the restored data
    double sample_entropy;  ///< sum of the entropy of the individual assignments
    double *posterior_probability;  ///< posterior probabilities of the most probable assignments
    double *entropy;        ///< entropy of the individual assignments

    void copy(const MixtureData &vec , bool model_flag = true);
    void remove();

    void observation_frequency_distribution_computation(int variable , int nb_component);

  public :

    MixtureData();
    MixtureData(int inb_vector , int inb_variable , variable_nature *itype , bool init_flag = false);
    MixtureData(const Vectors &vec , vector_transformation transform = VECTOR_COPY);
    MixtureData(const MixtureData &vec , bool model_flag = true ,
                vector_transformation transform = VECTOR_COPY)
    :Vectors(vec , transform) { copy(vec , model_flag); }
    ~MixtureData();
    MixtureData& operator=(const MixtureData &vec);

    DiscreteDistributionData* extract(StatError &error , int variable , int index) const;

    std::ostream& ascii_data_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_data_write(StatError &error , const std::string path , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void state_variable_init(variable_nature itype = STATE);

    double classification_information_computation() const;
    double information_computation() const;

    void build_observation_frequency_distribution(int nb_component);
    void build_observation_histogram(int variable , int nb_component , double bin_width = D_DEFAULT);
    void build_observation_histogram(int nb_component);
    bool select_bin_width(StatError &error , int variable , double bin_width ,
                          double imin_value = D_INF);

    // class member access

    Mixture* get_mixture() const { return mixture; }
    FrequencyDistribution*** get_observation_distribution() const
    { return observation_distribution; }
    FrequencyDistribution** get_observation_distribution(int variable) const
    { return observation_distribution[variable]; }
    FrequencyDistribution* get_observation_distribution(int variable , int index) const
    { return observation_distribution[variable][index]; }
    Histogram*** get_observation_histogram() const { return observation_histogram; }
    Histogram** get_observation_histogram(int variable) const
    { return observation_histogram[variable]; }
    Histogram* get_observation_histogram(int variable , int index) const
    { return observation_histogram[variable][index]; }
    double get_likelihood() const { return likelihood; }
    double get_restoration_likelihood() const { return restoration_likelihood; }
    double get_sample_entropy() const { return sample_entropy; }
    double get_posterior_probability(int index) const { return posterior_probability[index]; }
    double get_entropy(int index) const { return entropy[index]; }
  };


};  // namespace stat_tool



#endif
