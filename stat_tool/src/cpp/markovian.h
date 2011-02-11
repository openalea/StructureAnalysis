/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
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



#ifndef MARKOVIAN_H
#define MARKOVIAN_H



#include "chain_reestimation.h"



/****************************************************************
 *
 *  Constantes :
 */


const int NB_STATE = 100;              // nombre maximum d'etats d'une chaine de Markov
const int ORDER = 8;                   // ordre maximum d'une chaine de Markov
const double MIN_PROBABILITY = 1.e-5;  // probabilite minimum
const double THRESHOLDING_FACTOR = 0.8;  // facteur pour le seuillage des probabilites
const int NB_PARAMETER = 100000;       // nombre maximum de parametres d'une chaine de Markov
const int NB_OUTPUT_PROCESS = 10;      // nombre maximum de processus d'observation
const int NB_OUTPUT = 12;              // nombre maximum d'observations par etat (cas non-parametrique)
const double OBSERVATION_THRESHOLD = 0.999;  // seuil sur la fonction de repartition pour borner
                                             // une loi d'observation discrete parametrique

const double ACCESSIBILITY_THRESHOLD = 1.e-6;  // seuil pour stopper l'algorithme
                                               // de calcul de l'accessibilite des etats
const int ACCESSIBILITY_LENGTH = 100;  // longueur maximum de sequence pour l'algorithme
                                       // de calcul de l'accessibilite des etats

const double NOISE_PROBABILITY = 0.05;  // perturbation des probabilites d'observation
const double MEAN_SHIFT_COEFF = 0.1;   // coefficient pour decaler les lois d'observations
                                       // continus parametriques

const int MIN_NB_ELEMENT = 10;         // taille minimum de l'echantillon construit par arrondi
const int OBSERVATION_COEFF = 10;      // coefficient arrondi estimateur pour les lois
                                       // d'observation parametriques

const int GAUSSIAN_MAX_NB_DECIMAL = 6;  // nombre maximum de decimales pour la simulation
                                        // d'une loi gaussienne
const int DEGREE_DECIMAL_SCALE = 10;   // facteur pour determiner le nombre de decimales
                                       // pour la simulation d'une loi de Von Mises en degrees
const int RADIAN_DECIMAL_SCALE = 1000;  // facteur pour determiner le nombre de decimales
                                        // pour la simulation d'une loi de Von Mises en radians

// const double SELF_TRANSITION = 0.9;    probabilite de rester dans un etat initiale

enum {
  NON_PARAMETRIC ,
  DISCRETE_PARAMETRIC ,
  CONTINUOUS_PARAMETRIC
};

enum {
  FORWARD ,
  FORWARD_BACKWARD ,
  VITERBI ,
//  VITERBI_FORWARD_BACKWARD ,
  GENERALIZED_VITERBI ,
  FORWARD_BACKWARD_SAMPLING ,
  GIBBS_SAMPLING ,
  FORWARD_DYNAMIC_PROGRAMMING
};



/****************************************************************
 *
 *  Definition des classes :
 */


class ChainData;

class Chain {           // chaine de Markov

//    friend class ChainData;
    friend std::ostream& operator<<(std::ostream &os , const Chain &chain)
    { return chain.ascii_print(os , true); }

// protected :
public :

    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int nb_state;           // nombre d'etats
    int nb_row;             // nombre de lignes de la matrice
                            // des probabilites de transition
    bool **accessibility;   // matrice d'accessibilite des etats
    int nb_component;       // nombre de classes
    int *component_nb_state;  // nombre d'etats par classe
    int **component;        // classes
    char *state_type;       // types des etats ('r' : recurrent,
                            // 't' : transitoire, 'a' : absorbant)
    double *initial;        // probabilites initiales
    double *cumul_initial;  // fonction de repartition correspondant
                            // au probabilites initiales
    double **transition;    // matrice des probabilites de transition
    double **cumul_transition;  // fonctions de repartition correspondant aux lignes
                                // de la matrice des probabilites de transition

    void parameter_copy(const Chain&);
    void copy(const Chain&);
    void remove();

    std::ostream& ascii_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;

    void create_cumul();
    void cumul_computation();
    void remove_cumul();
    void log_computation();

    bool** logic_transition_computation() const;
    bool connex_component_research(StatError &error , bool **ilogic_transition = NULL) const;
    void graph_accessibility_computation(bool **ilogic_transition);
    void probability_accessibility_computation();
    void component_computation(bool **ilogic_transition = NULL);

    void thresholding(double min_probability , bool semi_markov = false);

    int nb_parameter_computation(double min_probability = 0.) const;
    double chi2_value_computation(const ChainData &chain_data) const;

// public :

    Chain(char itype = 'o' , int inb_state = 0 , bool init_flag = true);
    Chain(char itype , int inb_state , int inb_row , bool init_flag);
    Chain(const Chain &chain) { copy(chain); }
    ~Chain();
    Chain& operator=(const Chain &chain);

    void init(bool left_right , double self_transition);

    double likelihood_computation(const ChainData &chain_data , bool initial_flag = true) const;
    void chi2_fit(const ChainData &chain_data , Test &test) const;

    // acces membres de la classe

/*    char get_type() const { return type; }
    int get_nb_state() const { return nb_state; }
    int get_nb_row() const { return nb_row; }
    bool get_accessibility(int state1 , int state2) const
    { return accessibility[state1][state2]; }
    int get_nb_component() const { return nb_component; }
    int get_component_nb_state(int icomponent) const
    { return component_nb_state[icomponent]; }
    int get_component(int icomponent , int index) const
    { return component[icomponent][index]; }
    char get_state_type(int state) const { return state_type[state]; }
    double get_initial(int state) const { return initial[state]; }
    double get_transition(int memory , int state) const
    { return transition[memory][state]; } */
};


Chain* chain_parsing(StatError &error , ifstream &in_file , int &line , char type);



class ChainData : public ChainReestimation<int> {  // structure de donnees correspondant a
                                                   // une chaine de Markov

/*    friend class Chain;
    friend class MarkovianSequences;
    friend class VariableOrderMarkov;
    friend class VariableOrderMarkovData;
    friend class SemiMarkovData;
    friend class NonhomogeneousMarkov;
    friend class NonhomogeneousMarkovData; */

// protected :
public :

    int nb_parameter_computation() const;

// public :

    ChainData(char itype , int inb_state , int inb_row , bool init_flag = false)
    :ChainReestimation<int>(itype , inb_state , inb_row , init_flag) {}
    ChainData(const ChainData &chain_data);

    void estimation(Chain &chain) const;

    // acces membres de la classe

/*    char get_type() const { return type; }
    int get_nb_state() const { return nb_state; }
    int get_nb_row() const { return nb_row; }
    int get_initial(int state) const { return initial[state]; }
    int get_transition(int memory , int state) const
    { return transition[memory][state]; } */
};



class NonparametricProcess {  // processus d'observation discret non-parametrique

    friend class MultivariateMixture;

    friend NonparametricProcess* observation_parsing(StatError &error , ifstream &in_file ,
                                                     int &line , int nb_state , bool hidden);

    friend NonparametricProcess** observation_parsing(StatError &error , ifstream &in_file ,
                                                      int &line , int nb_state , int &nb_output_process);
    friend NonparametricProcess** old_observation_parsing(StatError &error , ifstream &in_file ,
                                                          int &line , int nb_state , int &nb_output_process);

protected :

    int nb_state;           // nombre d'etats
    int nb_value;           // nombre de valeurs
    Distribution **observation;  // lois d'observation

    void copy(const NonparametricProcess &process);
    void remove();

    bool test_hidden() const;
    void thresholding(double min_probability);

    void state_permutation(int *perm) const; // permutation des etats

public :

    NonparametricProcess(int inb_state = 0 , int inb_value = 0 , int observation_flag = false);
    NonparametricProcess(int inb_state , int inb_value , double **observation_probability);
    NonparametricProcess(int inb_state , Distribution **pobservation);
    NonparametricProcess(const NonparametricProcess &process)
    { copy(process); }
    ~NonparametricProcess();
    NonparametricProcess& operator=(const NonparametricProcess &process);

    int nb_parameter_computation(double min_probability) const;
    void init();

    // acces membres de la classe

    int get_nb_state() const { return nb_state; }
    int get_nb_value() const { return nb_value; }
    Distribution* get_observation(int state) const
    { return observation[state]; }

    std::ostream& ascii_print(std::ostream &os , FrequencyDistribution **empirical_observation ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_print(std::ostream &os ,
                                    FrequencyDistribution **empirical_observation = NULL) const;
    bool plot_print(const char *prefix , const char *title , int process ,
                    FrequencyDistribution **empirical_observation = NULL) const;
};


NonparametricProcess* observation_parsing(StatError &error , ifstream &in_file ,
                                          int &line , int nb_state , bool hidden);

NonparametricProcess** observation_parsing(StatError &error , ifstream &in_file ,
                                           int &line , int nb_state , int &nb_output_process);
NonparametricProcess** old_observation_parsing(StatError &error , ifstream &in_file ,
                                               int &line , int nb_state , int &nb_output_process);



class DiscreteParametricProcess {  // processus d'observation discret parametrique

/*    friend class MarkovianSequences;
    friend class SemiMarkov;
    friend class SemiMarkovData;
    friend class SemiMarkovIterator;
    friend class HiddenSemiMarkov;
    friend class VariableOrderMarkov;
    friend class VariableOrderMarkovData;
    friend class VariableOrderMarkovIterator;
    friend class HiddenVariableOrderMarkov; */

    friend DiscreteParametricProcess* discrete_observation_parsing(StatError &error , ifstream &in_file ,
                                                                   int &line , int nb_state ,
                                                                   double cumul_threshold);

// private :
public :

    int nb_state;           // nombre d'etats
    int nb_value;           // nombre de valeurs
    DiscreteParametric **observation;  // lois d'observation
    Distribution *weight;   // poids theoriques des lois d'observation
    Distribution *mixture;  // melange de lois d'observation
    Distribution *restoration_weight;  // poids des lois d'observation
                                       // deduits de la restoration
    Distribution *restoration_mixture;  // melange de lois d'observation

    void copy(const DiscreteParametricProcess &process);
    void state_permutation(int *perm) const;  // permutation des etats
    void remove();

    void nb_value_computation();
    Distribution* mixture_computation(Distribution *pweight);

// public :

    DiscreteParametricProcess(int inb_state = 0 , int inb_value = 0);
    DiscreteParametricProcess(int inb_state , DiscreteParametric **pobservation);
    DiscreteParametricProcess(const DiscreteParametricProcess &process)
    { copy(process); }
    ~DiscreteParametricProcess();
    DiscreteParametricProcess& operator=(const DiscreteParametricProcess &process);

    std::ostream& ascii_print(std::ostream &os , FrequencyDistribution **empirical_observation ,
                              FrequencyDistribution *marginal_distribution ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_print(std::ostream &os ,
                                    FrequencyDistribution **empirical_observation = NULL ,
                                    FrequencyDistribution *marginal_distribution = NULL) const;
    bool plot_print(const char *prefix , const char *title , int process ,
                    FrequencyDistribution **empirical_observation = NULL ,
                    FrequencyDistribution *marginal_distribution = NULL) const;
    void plotable_write(MultiPlotSet &plot , int &index , int process ,
                        FrequencyDistribution **empirical_observation = NULL ,
                        FrequencyDistribution *marginal_distribution = NULL) const;

    int nb_parameter_computation() const;
    void init();

    // acces membres de la classe

/*    int get_nb_state() const { return nb_state; }
    int get_nb_value() const { return nb_value; }
    DiscreteParametric* get_observation(int state) const
    { return observation[state]; } */
};


DiscreteParametricProcess* discrete_observation_parsing(StatError &error , ifstream &in_file ,
                                                        int &line , int nb_state ,
                                                        double cumul_threshold = OBSERVATION_THRESHOLD);



class ContinuousParametricProcess {  // processus d'observation continu parametrique

    friend ContinuousParametricProcess* continuous_observation_parsing(StatError &error , ifstream &in_file ,
                                                                       int &line , int nb_state);

// private :
public :

    int nb_state;           // nombre d'etats
    int ident;              // identificateur des lois d'observation
    int unit;               // unite (degre/radian) pour les lois de von Mises
    ContinuousParametric **observation;  // lois d'observation
    Distribution *weight;   // poids theorique des lois d'observation
    Distribution *restoration_weight;  // poids des lois d'observation
                                       // deduits de la restoration

    void copy(const ContinuousParametricProcess &process);
    void remove();

// public :

    ContinuousParametricProcess(int inb_state = 0);
    ContinuousParametricProcess(int inb_state , ContinuousParametric **pobservation);
    ContinuousParametricProcess(const ContinuousParametricProcess &process)
    { copy(process); }
    ~ContinuousParametricProcess();
    ContinuousParametricProcess& operator=(const ContinuousParametricProcess &process);

    std::ostream& ascii_print(std::ostream &os , Histogram **observation_histogram ,
                              FrequencyDistribution **observation_distribution ,
                              Histogram *marginal_histogram ,
                              FrequencyDistribution *marginal_distribution ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_print(std::ostream &os ,
                                    Histogram **observation_histogram = NULL ,
                                    FrequencyDistribution **observation_distribution = NULL ,
                                    Histogram *marginal_histogram = NULL ,
                                    FrequencyDistribution *marginal_distribution = NULL) const;
    bool plot_print(const char *prefix , const char *title ,
                    int process , Histogram **observation_histogram = NULL ,
                    FrequencyDistribution **observation_distribution = NULL ,
                    Histogram *marginal_histogram = NULL ,
                    FrequencyDistribution *marginal_distribution = NULL ,
                    int nb_value = I_DEFAULT , double **empirical_cdf = NULL) const;
    void plotable_write(MultiPlotSet &plot , int &index , int process ,
                        Histogram **observation_histogram = NULL ,
                        FrequencyDistribution **observation_distribution = NULL ,
                        Histogram *marginal_histogram = NULL ,
                        FrequencyDistribution *marginal_distribution = NULL ,
                        int nb_value = I_DEFAULT , double **empirical_cdf = NULL) const;

    int nb_parameter_computation() const;
    void select_unit(int iunit);
    void init(int iident , double min_value , double max_value ,
              double mean , double variance);
};


ContinuousParametricProcess* continuous_observation_parsing(StatError &error , ifstream &in_file ,
                                                            int &line , int nb_state);




# endif
