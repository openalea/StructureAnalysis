/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Identifying patterns in plant architecture and development
 *
 *       Copyright 1995-2018 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for StructureAnalysis developers:
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



#ifndef REESTIMATION_H
#define REESTIMATION_H

#include <iostream>

#include "config.h"

namespace stat_tool {



/****************************************************************
 *
 *  Constants
 */


  const double CUMUL_THRESHOLD = 0.999;  // threshold on the cumulative distribution function
                                         // to determine an upper bound of the support

  const double BISECTION_RATIO_THRESHOLD = 1.e-8;  // threshold for stopping the interval bisection method
  const int BISECTION_NB_ITER = 100;     // maximum number of iterations for the interval bisection method

  enum penalty_type {
    FIRST_DIFFERENCE ,
    SECOND_DIFFERENCE ,
    ENTROPY
  };

  enum side_effect {
    ZERO ,
    CONTINUATION
  };



/****************************************************************
 *
 *  Class definition
 */


  class Distribution;
  class DiscreteParametric;
  class ContinuousParametric;

  /// \brief Frequency distribution with integer or real (for EM algorithms) frequencies

  template <typename Type>
  class Reestimation {

    friend std::ostream& operator<<(std::ostream &os , const Reestimation<Type> &histo)
    { return histo.print(os); }

  public :

    int nb_value;           ///< number of values from 0
    int alloc_nb_value;     ///< number of allocated values
    int offset;             ///< number of values of null frequencies from 0
    Type nb_element;        ///< sample size
    Type max;               ///< maximum frequency
    double mean;            ///< mean
    double variance;        ///< variance
    Type *frequency;        ///< frequency for each value

    void init(int inb_value);
    void copy(const Reestimation<Type> &histo);

    Reestimation(int inb_value = 0) { init(inb_value); }
    Reestimation(const Reestimation<Type> &histo);
    Reestimation(int nb_histo , const Reestimation<Type> **histo);
    ~Reestimation();
    Reestimation<Type>& operator=(const Reestimation<Type> &histo);

    std::ostream& ascii_characteristic_print(std::ostream &os , bool shape = false ,
                                             bool comment_flag = false) const;
    std::ostream& ascii_circular_characteristic_print(std::ostream &os ,
                                                      bool comment_flag = false) const;
    std::ostream& print(std::ostream &os) const;

    void nb_value_computation();
    void offset_computation();
    void nb_element_computation();
    void max_computation();
    double mode_computation() const;
    void mean_computation();
    double quantile_computation(double icumul = 0.5) const;
    void variance_computation(bool bias = false);

    double mean_absolute_deviation_computation(double location) const;
    double log_geometric_mean_computation() const;
    double skewness_computation() const;
    double kurtosis_computation() const;

    void mean_direction_computation(double *mean_direction) const;

    double information_computation() const;
    double likelihood_computation(const Distribution &dist) const;

    void distribution_estimation(Distribution *pdist) const;
    void penalized_likelihood_estimation(Distribution *dist , double weight , penalty_type pen_type ,
                                         double *penalty , side_effect outside) const;

    double binomial_estimation(DiscreteParametric *pdist , int min_inf_bound , bool min_inf_bound_flag) const;
    double poisson_estimation(DiscreteParametric *pdist , int min_inf_bound , bool min_inf_bound_flag ,
                              double cumul_threshold) const;
    double negative_binomial_estimation(DiscreteParametric *pdist , int min_inf_bound , bool min_inf_bound_flag ,
                                        double cumul_threshold) const;
    double poisson_geometric_estimation(DiscreteParametric *pdist , int min_inf_bound , bool min_inf_bound_flag ,
                                        double cumul_threshold) const;

    double parametric_estimation(DiscreteParametric *pdist , int min_inf_bound = 0 , bool min_inf_bound_flag = true ,
                                 double cumul_threshold = CUMUL_THRESHOLD , bool poisson_geometric = false) const;
    double type_parametric_estimation(DiscreteParametric *pdist , int min_inf_bound = 0 , bool min_inf_bound_flag = true ,
                                      double cumul_threshold = CUMUL_THRESHOLD , bool poisson_geometric = false) const;

    DiscreteParametric* type_parametric_estimation(int min_inf_bound = 0 , bool min_inf_bound_flag = true ,
                                                   double cumul_threshold = CUMUL_THRESHOLD) const;

    void equilibrium_process_combination(const Reestimation<Type> *length_bias_reestim , double imean);
    void equilibrium_process_estimation(const Reestimation<Type> *length_bias_reestim ,
                                        Distribution *dist , double imean) const;
    void penalized_likelihood_equilibrium_process_estimation(const Reestimation<Type> *length_bias_reestim ,
                                                             Distribution *dist , double imean ,
                                                             double weight , penalty_type pen_type ,
                                                             double *penalty , side_effect outside) const;

    void state_occupancy_estimation(const Reestimation<Type> *final_run ,
                                    Reestimation<double> *occupancy_reestim ,
                                    Type *occupancy_survivor , Type *censored_occupancy_survivor ,
                                    bool characteristic_computation = true);

    void gamma_estimation(ContinuousParametric *dist , int iter) const;
    void zero_inflated_gamma_estimation(ContinuousParametric *dist , int iter) const;
    void inverse_gaussian_estimation(ContinuousParametric *dist) const;
  };


  STAT_TOOL_API double interval_bisection(Reestimation<double> *distribution_reestim ,
                            Reestimation<double> *length_bias_reestim);


};  // namespace stat_tool



#endif
