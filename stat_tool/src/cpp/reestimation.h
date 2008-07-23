/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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



/****************************************************************
 *
 *  Constantes :
 */


const double CUMUL_THRESHOLD = 0.999;  // seuil sur la fonction de repartition
                                       // pour borner une loi

const double BISECTION_RATIO_THRESHOLD = 1.e-8;  // seuil pour stopper la recherche
                                                 // par bissection d'intervalle
const int BISECTION_NB_ITER = 100;     // nombre d'iterations maximum pour la methode
                                       // de bissection d'intervalle



/****************************************************************
 *
 *  Definition de la classe :
 */


class Distribution;
class Parametric;


template <typename Type>
class Reestimation {    // histogramme a frequences entieres ou reelles (estimateur EM)

/*    friend class Parametric;
      friend class Histogram;
      friend class Compound;
      friend class Convolution;
      friend class Renewal;
      friend class Time_events;
      friend class Markovian_sequences; */

    friend std::ostream& operator<<(std::ostream &os , const Reestimation<Type> &histo)
    { return histo.print(os); }

// protected :
public :

    int nb_value;           // nombre de valeurs a partir de 0
    int alloc_nb_value;     // nombre de valeurs alloues
    int offset;             // nombre de valeurs de frequence nulle a partir de 0
    Type nb_element;        // effectif total
    Type max;               // valeur de frequence maximum
    double mean;            // moyenne
    double variance;        // variance
    Type *frequency;        // frequence de chacune des valeurs

    void init(int inb_value);
    void copy(const Reestimation<Type> &histo);

    std::ostream& ascii_characteristic_print(std::ostream &os , bool shape = false ,
                                             bool comment_flag = false) const;
    std::ostream& print(std::ostream &os) const;

    void nb_value_computation();
    void offset_computation();
    void nb_element_computation();
    void max_computation();
    void mean_computation();
    void variance_computation();

    void distribution_estimation(Distribution *pdist) const;
    void penalized_likelihood_estimation(Distribution *dist , double weight , int type ,
                                         double *penalty , int outside) const;

    double binomial_estimation(Parametric *pdist , int min_inf_bound , bool min_inf_bound_flag) const;
    double poisson_estimation(Parametric *pdist , int min_inf_bound , bool min_inf_bound_flag ,
                              double cumul_threshold) const;
    double negative_binomial_estimation(Parametric *pdist , int min_inf_bound , bool min_inf_bound_flag ,
                                        double cumul_threshold) const;

    double parametric_estimation(Parametric *pdist , int min_inf_bound = 0 , bool min_inf_bound_flag = true ,
                                 double cumul_threshold = CUMUL_THRESHOLD) const;
    double type_parametric_estimation(Parametric *pdist , int min_inf_bound = 0 , bool min_inf_bound_flag = true ,
                                      double cumul_threshold = CUMUL_THRESHOLD) const;

    Parametric* type_parametric_estimation(int min_inf_bound = 0 , bool flag = true ,
                                           double cumul_threshold = CUMUL_THRESHOLD) const;

    void equilibrium_process_combination(const Reestimation<Type> *length_bias_reestim , double imean);
    void equilibrium_process_estimation(const Reestimation<Type> *length_bias_reestim ,
                                        Distribution *dist , double imean) const;
    void penalized_likelihood_equilibrium_process_estimation(const Reestimation<Type> *length_bias_reestim ,
                                                             Distribution *dist , double imean ,
                                                             double weight , int type ,
                                                             double *penalty , int outside) const;

    void state_occupancy_estimation(const Reestimation<Type> *final_run ,
                                    Reestimation<double> *occupancy_reestim ,
                                    Type *occupancy_survivor , Type *censored_occupancy_survivor ,
                                    bool characteristic_computation = true);

// public :

    Reestimation(int inb_value = 0) { init(inb_value); }
    Reestimation(const Reestimation<Type> &histo);
    Reestimation(int nb_histo , const Reestimation<Type> **histo);
    ~Reestimation();
    Reestimation<Type>& operator=(const Reestimation<Type> &histo);

    double mean_absolute_deviation_computation() const;
    double skewness_computation() const;
    double kurtosis_computation() const;

    double information_computation() const;
    double likelihood_computation(const Distribution &dist) const;
};



#endif
