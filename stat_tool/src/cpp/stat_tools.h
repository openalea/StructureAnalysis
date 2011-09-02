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



#ifndef STAT_TOOLS_H
#define STAT_TOOLS_H



#include <fstream>
#include "reestimation.h"
#include "plotable.h"

using namespace plotable;



/****************************************************************
 *
 *  Macros :
 */


#ifndef MIN
#define MIN(x,y)  ((x) < (y) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x,y)  ((x) > (y) ? (x) : (y))
#endif



/****************************************************************
 *
 *  Constantes :
 */


const int ERROR_LENGTH = 200;

const int I_DEFAULT = -1;              // int par defaut
const double D_DEFAULT = -1.;          // double par defaut
const double D_INF = -1.e37;           // plus petit nombre reel
const double DOUBLE_ERROR = 1.e-6;     // erreur sur une somme de doubles
// const double DOUBLE_ERROR = 5.e-6;     erreur sur une somme de doubles


enum {
  STANDARD_NORMAL ,
  CHI2 ,
  FISHER ,
  STUDENT
};

const int NB_CRITICAL_PROBABILITY = 2;
const double ref_critical_probability[NB_CRITICAL_PROBABILITY] = {0.05 , 0.01};

const int NB_VALUE = 1000;             // nombre de valeurs prises par une v.a.
const int SAMPLE_NB_VALUE = NB_VALUE;  // nombre de valeurs d'un echantillon

enum {
  FLOOR ,
  ROUND ,
  CEIL
};

enum {
  NONPARAMETRIC ,
  BINOMIAL ,
  POISSON ,
  NEGATIVE_BINOMIAL ,
  UNIFORM ,
  MULTINOMIAL                          // ajout par Florence Chaubert
};

enum {
//  GAMMA ,
  GAUSSIAN ,
  VON_MISES
};

enum {
  DEGREE ,
  RADIAN
};

enum {
  SYMBOLIC ,
  ORDINAL ,
  NUMERIC ,
  CIRCULAR
};

enum {
  PEARSON ,
  SPEARMAN ,
  KENDALL ,
  SPEARMAN2
};

enum {
  FIRST_DIFFERENCE ,
  SECOND_DIFFERENCE ,
  ENTROPY
};

enum {
  LIKELIHOOD ,
  PENALIZED_LIKELIHOOD ,
  PARAMETRIC_REGULARIZATION
};

enum {
  ZERO ,
  CONTINUATION
};

enum {
  PARTIAL_LIKELIHOOD ,
  COMPLETE_LIKELIHOOD ,
  KAPLAN_MEIER
};

enum {
  COMPUTED ,
  ESTIMATED ,
  ONE_STEP_LATE
};

enum {
  AIC ,
  AICc ,
  BIC ,
  BICc ,
  ICL ,
  ICLc
};

enum {
  AGGLOMERATIVE ,
  DIVISIVE ,
  ORDERING
};

enum {
  ERROR ,
  WARNING
};

enum {
  INT_VALUE ,                          // observation entiere
  REAL_VALUE ,                         // observation reelle
  STATE ,                              // etat
  OLD_INT_VALUE ,                      // pour compatibilite ascendante
  NB_INTERNODE ,                       // nombre d'entrenoeuds axe porte
  AUXILIARY                            // variable auxiliaire (lissage/moyenne par segment)
};

enum {
  SELF_TRANSITION ,
  OBSERVATION ,
  INTENSITY ,
  FIRST_OCCURRENCE ,
  RECURRENCE_TIME ,
  SOJOURN_TIME ,
  INITIAL_RUN ,
  FINAL_RUN ,
  NB_RUN ,
  NB_OCCURRENCE ,
  COUNTING ,
  LENGTH ,
  SEQUENCE_CUMUL ,
  SEQUENCE_MEAN
};

const int MAX_INF_BOUND = 10000;       // borne inferieure maximum
const int MAX_DIFF_BOUND = 10000;      // difference maximum entre les bornes
                                       // inferieure et superieure
const double MAX_MEAN = 10000.;        // moyenne maximum

const double B_PROBABILITY = 0.8;      // seuil pour utiliser le calcul descendant de la loi binomiale
const double B_THRESHOLD = 1000.;      // seuil pour utiliser le calcul en log de la loi binomiale
const double P_THRESHOLD = 90.;        // seuil pour utiliser le calcul en log de la loi de Poisson
const double NB_THRESHOLD = 500.;      // seuil pour utiliser le calcul en log de la loi binomiale negative

const double SAMPLE_NB_VALUE_COEFF = 5.;  // facteur pour deduire le nombre de valeurs prises
                                          // par une v.a. du nombre de valeurs d'un echantillon
const int INF_BOUND_MARGIN = 5;        // plage de recherche pour la borne inferieure
const int SUP_BOUND_MARGIN = 3;        // plage de recherche pour la borne superieure
const double POISSON_RATIO = 0.7;      // rapport moyenne/variance minimum pour
                                       // estimer une loi de Poisson
const double POISSON_RANGE = 0.1;      // plage de variation pour choisir une loi de
                                       // Poisson par dilatation de l'echelle des temps
const double NB_VALUE_COEFF = 2.;      // facteur pour deduire le nombre de valeurs prises par une v.a.
                                       // du nombre de valeurs prises par une v.a. initiale

const int MIN_RANGE = 10;              // intervalle de definition minimum
                                       // pour appliquer la methode du rejet
const double MAX_SURFACE = 3.;         // surface maximum pour appliquer la methode du rejet
const int DIST_NB_ELEMENT = 1000000;   // taille maximum de l'echantillon pour la simulation

const double GAUSSIAN_TAIL = 5.e-4;    // traine de la loi de Gauss
const int GAUSSIAN_NB_STEP = 1000;     // nombre de pas pour le calcul de la loi de Gauss
const int GAUSSIAN_NB_SUB_STEP = 10;   // nombre de pas pour le calcul de la loi de Gauss
const int VON_MISES_NB_STEP = 3600;    // nombre de pas pour le calcul de la loi de von Mises
const int VON_MISES_NB_SUB_STEP = 10;  // nombre de pas pour le calcul de la loi de von Mises
// const double CONCENTRATION_THRESHOLD = 10.;  seuil sur le parametre de concentration
                                             // pour appliquer l'approximation gaussienne
                                             // pour le calcul de la loi de von Mises

const int CHI2_FREQUENCY = 2;          // effectif theorique minimum pour un
                                       // test d'ajustement du Chi2

const int MARGINAL_DISTRIBUTION_MAX_VALUE = 20000;  // valeur maximum pour la construction
                                                    // de la loi marginale
const int HISTOGRAM_FREQUENCY = 10;    // frequence moyenne pour definir le pas de regroupement
                                       // d'un histogramme
const double SKEWNESS_ROUNDNESS = 1.e-2;  // arrondi sur le coefficient d'asymetrie

const int NB_ERROR = 10;               // nombre maximum d'erreurs prises en compte

const int LINE_NB_CHARACTER = 100;     // nombre de caracteres par ligne pour les sequences

const int ASCII_NB_VALUE = 15;         // nombre de valeurs maximum (sortie ASCII)
const int ASCII_SPACE = 2;             // nombre d'espaces entre 2 colonnes (sortie ASCII)
const double ASCII_ROUNDNESS = 1.e-5;  // arrondi sur la fonction de repartition
                                       // pour borner une loi (sortie ASCII)

const double SPREADSHEET_ROUNDNESS = 1.e-7;  // arrondi sur la fonction de repartition
                                             // pour borner une loi (sortie tableur)

const int DISPLAY_NB_INDIVIDUAL = 50;  // nombre maximum d'individus selectionnes affiches

const int PLOT_NB_DISTRIBUTION = 10;   // nombre maximum de lois affichees (sortie graphique)
const int PLOT_NB_HISTOGRAM = 10;      // nombre maximum d'histogrammes affiches (sortie graphique)
const double PLOT_ROUNDNESS = 1.e-5;   // arrondi sur la fonction de repartition
                                       // pour borner une loi (sortie graphique)
const double PLOT_SHIFT = 0.2;         // decalage entre 2 histogrammes (sortie graphique)
const double PLOT_MAX_SHIFT = 0.5;     // decalage maximum entre le premier et le dernier
                                       // histogramme (sortie graphique)
const int TIC_THRESHOLD = 10;          // nombre de graduations minimum pour echelle
                                       // automatique (sortie graphique)
const double PLOT_MASS_THRESHOLD = 1.e-3;  // valeur minimale pour afficher un 0 apres la derniere
                                           // valeur possible (sortie graphique)
const double YSCALE = 1.4;             // facteur d'echelle axe y (sortie graphique)



/****************************************************************
 *
 *  Definition des classes :
 */


class Test {            // test d'hypothese

/*    friend class Distribution;
    friend class Vectors;
    friend class Regression;
    friend class Chain; */

    friend std::ostream& operator<<(std::ostream &os , const Test &test)
    { return test.ascii_print(os); }

// private :
public :

    int ident;              // identificateur (Normale / Chi2 / F / t)
    bool one_side;          // unilateral / bilateral
    int df1;                // nombre de degres de liberte (Chi2 / F / t)
    int df2;                // nombre de degres de liberte (F)
    double value;           // valeur
    double critical_probability;  // probabilite critique

    void copy(const Test &test);

// public :

    Test(int iident , bool ione_side = true);
    Test(int iident , bool ione_side , int idf1 , int idf2 , double ivalue);
    Test(int iident , bool ione_side , int idf1 , int idf2 , double ivalue ,
         double icritical_probability);
    Test(const Test &test)
    { copy(test); }
    Test& operator=(const Test &test);

    std::ostream& ascii_print(std::ostream &os , bool comment_flag = false ,
                              bool reference_flag = true) const;
    std::ostream& spreadsheet_print(std::ostream &os , bool reference_flag = true) const;

    void standard_normal_critical_probability_computation();
    void standard_normal_value_computation();
    void chi2_critical_probability_computation();
    void chi2_value_computation();
    void F_critical_probability_computation();
    void F_value_computation();
    void t_critical_probability_computation();
    void t_value_computation();

    // acces membres de la classe

/*    int get_ident() const { return ident; }
    bool get_one_side() const { return one_side; }
    int get_df1() const { return df1; }
    int get_df2() const { return df2; }
    double get_value() const { return value; }
    double get_critical_probability() const { return critical_probability; } */
};



class StatError {       // erreurs

    friend std::ostream& operator<<(std::ostream &os , const StatError &error)
    { return error.ascii_write(os , ERROR); }

private :

    int nb_error;           // nombre d'erreurs
    int max_nb_error;       // nombre maximum d'erreurs
    int *line;              // no de ligne
    int *column;            // no de colonne
    char **label;           // messages d'erreur

public :

    StatError(int imax_nb_error = NB_ERROR);
    ~StatError();

    std::ostream& ascii_write(std::ostream &os , int type = ERROR) const;

    void init() { nb_error = 0; }
    void update(const char *ilabel , int iline = 0 , int icolumn = 0);
    void correction_update(const char *ilabel , const char *correction ,
                           int iline = 0 , int icolumn = 0);
    void correction_update(const char *ilabel , int correction ,
                           int iline = 0 , int icolumn = 0);

    // acces membres de la classe

    int get_nb_error() const { return nb_error; }
    int get_max_nb_error() const { return max_nb_error; }
};



class StatInterface {  // classe abstraite pour les interfaces

public :

    virtual std::ostream& line_write(std::ostream &os) const = 0;

    virtual std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const = 0;
    virtual bool ascii_write(StatError &error , const char *path ,
                             bool exhaustive = false) const = 0;
    virtual bool spreadsheet_write(StatError &error , const char *path) const = 0;
    virtual bool plot_write(StatError &error , const char *prefix ,
                            const char *title = NULL) const = 0;

    virtual MultiPlotSet* get_plotable() const { return NULL; };

//    bool binary_write(StatError &error , const char *path) const;
};



class FrequencyDistribution;
class DiscreteParametricModel;

class Distribution {    // loi de probabilite discrete

/*    template <typename Type> friend class Reestimation;  probleme Windows
    friend class FrequencyDistribution;
    friend class DiscreteParametricModel;
    friend class DiscreteDistributionData;
    friend class Convolution;
    friend class CompoundData;
    friend class Vectors;
    friend class VectorDistance;
    friend class Backward;
    friend class Renewal;
    friend class TimeEvents;
    friend class Curves;
    friend class ChainData;
    friend class NonparametricProcess;
    friend class NonparametricSequenceProcess;
    friend class NonhomogeneousMarkov;
    friend class VariableOrderMarkov;
    friend class HiddenVariableOrderMarkov;
    friend class SemiMarkov;
    friend class HiddenSemiMarkov;
    friend class Sequences;
    friend class MarkovianSequences;
    friend class Correlation;
    friend class TopParameters; */

    friend std::ostream& operator<<(std::ostream& , const Distribution&);
    friend bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                           double *scale , int *dist_nb_value , int nb_histo ,
                           const FrequencyDistribution **histo);

// protected :
  public :

    int nb_value;           // nombre de valeurs a partir de 0
    int alloc_nb_value;     // nombre de valeurs allouees
    int offset;             // nombre de valeurs de probabilite nulle a partir de 0
    double max;             // probabilite maximum
    double complement;      // probabilite complementaire
                            // (> 0. dans le cas d'une loi impropre)
    double mean;            // moyenne
    double variance;        // variance
    int nb_parameter;       // nombre de parametres inconnus
    double *mass;           // probabilites de chaque valeur
    double *cumul;          // fonction de repartition

    void max_computation();
    void mean_computation();
    void variance_computation();

    void mass_copy(const Distribution &dist , int inb_value = I_DEFAULT);
    void equal_size_copy(const Distribution &dist);
    void init(int inb_value);
    void copy(const Distribution &dist , int ialloc_nb_value = I_DEFAULT);
    void normalization_copy(const Distribution &dist);

    std::ostream& ascii_characteristic_print(std::ostream &os , bool shape = false ,
                                             bool comment_flag = false) const;
    std::ostream& ascii_print(std::ostream &os , bool comment_flag ,
                              bool cumul_flag , bool nb_value_flag ,
                              const FrequencyDistribution *histo = NULL) const;

    std::ostream& ascii_print(std::ostream &os , int nb_dist , const Distribution **dist ,
                              double *dist_scale , bool comment_flag , bool cumul_flag ,
                              const FrequencyDistribution *histo = NULL) const;

    std::ostream& spreadsheet_characteristic_print(std::ostream &os , bool shape = false) const;

    std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag ,
                                    bool concentration_flag , bool nb_value_flag ,
                                    const FrequencyDistribution *histo = NULL) const;

    std::ostream& spreadsheet_print(std::ostream &os , int nb_dist , const Distribution **dist ,
                                    double *dist_scale , bool cumul_flag ,
                                    const FrequencyDistribution *histo = NULL) const;

    int plot_nb_value_computation(const FrequencyDistribution *histo = NULL) const;
    bool plot_print(const char *path , double *concentration , double scale) const;
    bool plot_print(const char *path , const FrequencyDistribution *histo = NULL) const;
    virtual std::ostream& plot_title_print(std::ostream &os) const
    { return os; }
    bool survival_plot_print(const char *path , double *survivor) const;
    std::ostream& print(std::ostream&) const;

    void plotable_mass_write(SinglePlot &plot , double scale = 1.) const;
    void plotable_cumul_write(SinglePlot &plot) const;
    void plotable_cumul_matching_write(SinglePlot &plot , const Distribution &reference_dist) const;
    void plotable_concentration_write(SinglePlot &plot) const;
    void plotable_survivor_write(SinglePlot &plot) const;

    void convolution(Distribution &dist1 , Distribution &dist2 ,
                     int inb_value = I_DEFAULT);

    void nb_value_computation();
    void offset_computation();
    double concentration_computation() const;

    void cumul_computation();
    double* survivor_function_computation() const;
    double* concentration_function_computation() const;
    void log_computation();

    double survivor_likelihood_computation(const FrequencyDistribution &histo) const;
    double chi2_value_computation(const FrequencyDistribution &histo) const;
    void chi2_degree_of_freedom(const FrequencyDistribution &histo , Test &test) const;

    void penalty_computation(double weight , int type , double *penalty , int outside) const;

// public :

    Distribution(int inb_value = 0);
    Distribution(const Distribution &dist , double scaling_coeff);
    Distribution(const FrequencyDistribution &histo);
    Distribution(const Distribution &dist , char transform = 'c' ,
                 int ialloc_nb_value = I_DEFAULT);
    virtual ~Distribution();
    Distribution& operator=(const Distribution&);
    bool operator==(const Distribution&) const;
    bool operator!=(const Distribution &dist) const { return !(*this == dist); }

    bool plot_write(StatError &error , const char *prefix , int nb_dist ,
                    const Distribution **idist , const char *title) const;
    MultiPlotSet* get_plotable() const;
    MultiPlotSet* get_plotable_distributions(StatError &error , int nb_dist ,
                                             const Distribution **idist) const;

    std::ostream& survival_ascii_write(std::ostream &os) const;
    bool survival_ascii_write(StatError &error , const char *path) const;
    bool survival_spreadsheet_write(StatError &error , const char *path) const;
    bool survival_plot_write(StatError &error , const char *prefix ,
                             const char *title = NULL) const;
    MultiPlotSet* survival_get_plotable(StatError &error) const;

    double mean_absolute_deviation_computation() const;
    double skewness_computation() const;
    double kurtosis_computation() const;
    double information_computation() const;

    double first_difference_norm_computation() const;
    double second_difference_norm_computation() const;

    double likelihood_computation(const Reestimation<int> &histo) const
    { return histo.likelihood_computation(*this); }
    double likelihood_computation(const Reestimation<double> &histo) const
    { return histo.likelihood_computation(*this); }
    void chi2_fit(const FrequencyDistribution &histo , Test &test) const;

    int simulation() const;

    DiscreteParametricModel* truncate(StatError &error , int imax_value) const;

    // acces membres de la classe

/*    int get_offset() const { return offset; }
    int get_nb_value() const { return nb_value; }
    double get_max() const { return max; }
    double get_complement() const { return complement; }
    double get_mean() const { return mean; }
    double get_variance() const { return variance; }
    double& operator[](int index) const { return mass[index]; }
    double get_mass(int index) const { return mass[index]; }
    double get_cumul(int index) const { return cumul[index]; } */
};


bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                double *scale , int *dist_nb_value , int nb_histo ,
                const FrequencyDistribution **histo);



class Forward;

class DiscreteParametric : public Distribution {  // loi de probabilite discrete parametrique

/*    template <typename Type> friend class Reestimation;  probleme Windows
    friend class FrequencyDistribution;
    friend class DiscreteDistributionData;
    friend class Compound;
    friend class Convolution;
    friend class Mixture;
    friend class Backward;
    friend class Forward;
    friend class LengthBias;
    friend class NbEvent;
    friend class Renewal;
    friend class TimeEvents;
    friend class RenewalData;
    friend class NonparametricSequenceProcess;
    friend class DiscreteParametricProcess;
    friend class NonhomogeneousMarkov;
    friend class VariableOrderMarkov;
    friend class SemiMarkov;
    friend class HiddenSemiMarkov;
    friend class MarkovianSequences;
    friend class TopParameters; */

    friend int nb_value_computation(int ident , int inf_bound , int sup_bound ,
                                    double parameter , double probability ,
                                    double cumul_threshold);
    friend DiscreteParametric* discrete_parametric_parsing(StatError &error ,
                                                           std::ifstream &in_file ,
                                                           int &line , int last_ident ,
                                                           double cumul_threshold,
                                                           int min_inf_bound);
    friend std::ostream& operator<<(std::ostream& , const DiscreteParametric&);

// protected :
public :

    int ident;              // identificateur
    int inf_bound;          // borne inferieure
    int sup_bound;          // borne superieure (binomiale / uniforme)
    double parameter;       // parametre (Poisson / binomiale negative)
    double probability;     // probabilite de succes (binomiale / binomiale negative)

    void init(int iinf_bound , int isup_bound , double iparameter , double iprobability);
    void init(int iident , int iinf_bound , int isup_bound , double iparameter , double iprobability);
    void copy(const DiscreteParametric &dist);
    std::ostream& ascii_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    std::ostream& plot_title_print(std::ostream &os) const;

    void nb_parameter_update();

    void binomial_computation(int inb_value , char mode);
    void poisson_computation(int inb_value , double cumul_threshold ,
                             char mode);
    void negative_binomial_computation(int inb_value , double cumul_threshold ,
                                       char mode);
    void uniform_computation();

    double renewal_likelihood_computation(const Forward &forward_dist ,  // sequence_analysis
                                          const FrequencyDistribution &within ,
                                          const FrequencyDistribution &backward ,
                                          const FrequencyDistribution &forward ,
                                          const FrequencyDistribution *no_event) const;
    void expectation_step(const FrequencyDistribution &within ,  // sequence_analysis
                          const FrequencyDistribution &backward ,
                          const FrequencyDistribution &forward ,
                          const FrequencyDistribution *no_event ,
                          Reestimation<double> *inter_event_reestim ,
                          Reestimation<double> *length_bias_reestim , int iter) const;

    void reestimation(const Reestimation<double> *reestim , int nb_estim = 1);

    double state_occupancy_likelihood_computation(const FrequencyDistribution &sojourn_time ,  // sequence_analysis
                                                  const FrequencyDistribution &final_run) const;
    double state_occupancy_likelihood_computation(const Forward &forward ,  // sequence_analysis
                                                  const FrequencyDistribution &sojourn_time ,
                                                  const FrequencyDistribution &final_run ,
                                                  const FrequencyDistribution &initial_run ,
                                                  const FrequencyDistribution &single_run) const;
    void expectation_step(const FrequencyDistribution &sojourn_time ,  // sequence_analysis
                          const FrequencyDistribution &final_run ,
                          Reestimation<double> *occupancy_reestim , int iter) const;
    void expectation_step(const FrequencyDistribution &sojourn_time ,  // sequence_analysis
                          const FrequencyDistribution &final_run ,
                          const FrequencyDistribution &initial_run ,
                          const FrequencyDistribution &single_run ,
                          Reestimation<double> *occupancy_reestim ,
                          Reestimation<double> *length_bias_reestim , int iter ,
                          bool combination = false , int mean_computation_method = COMPUTED) const;

// public :

    DiscreteParametric(int inb_value = 0 , int iident = NONPARAMETRIC ,
                       int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
                       double iparameter = D_DEFAULT , double iprobability = D_DEFAULT);
    DiscreteParametric(int iident , int iinf_bound , int isup_bound , double iparameter ,
                       double iprobability , double cumul_threshold = CUMUL_THRESHOLD);
    DiscreteParametric(const Distribution &dist , int ialloc_nb_value = I_DEFAULT);
    DiscreteParametric(const Distribution &dist , double scaling_coeff);
    DiscreteParametric(const DiscreteParametric &dist , double scaling_coeff);
    DiscreteParametric(const FrequencyDistribution &histo);
    DiscreteParametric(const DiscreteParametric &dist , char transform = 'c' ,
                       int ialloc_nb_value = I_DEFAULT);
    DiscreteParametric& operator=(const DiscreteParametric &dist);

    int nb_parameter_computation();

    double parametric_mean_computation() const;
    double parametric_variance_computation() const;
    double parametric_skewness_computation() const;
    double parametric_kurtosis_computation() const;

    void computation(int min_nb_value = 1 ,
                     double cumul_threshold = CUMUL_THRESHOLD);
    int simulation() const;

    // acces membres de la classe

/*    int get_ident() const { return ident; }
    int get_inf_bound() const { return inf_bound; }
    int get_sup_bound() const { return sup_bound; }
    double get_parameter() const { return parameter; }
    double get_probability() const { return probability; } */
};


int nb_value_computation(int ident , int inf_bound , int sup_bound ,
                         double parameter , double probability ,
                         double cumul_threshold = CUMUL_THRESHOLD);
DiscreteParametric* discrete_parametric_parsing(StatError &error , std::ifstream &in_file ,
                                                int &line , int last_ident = NEGATIVE_BINOMIAL ,
                                                double cumul_threshold = CUMUL_THRESHOLD ,
                                                int min_inf_bound = 0);
std::ostream& operator<<(std::ostream& , const DiscreteParametric&);



class Forward : public DiscreteParametric {  // loi de l'intervalle de temps residuel

/*    friend class Renewal;
    friend class Renewal_data;
    friend class Semi_markov; */

public :

    Forward(int inb_value = 0 , int iident = NONPARAMETRIC ,
            int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
            double iparameter = D_DEFAULT , double iprobability = D_DEFAULT)
    :DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Forward(const DiscreteParametric &dist , int ialloc_nb_value = I_DEFAULT)
    :DiscreteParametric(dist , 'c' , ialloc_nb_value) { computation(dist); }
    Forward(const Forward &forward , int ialloc_nb_value = I_DEFAULT)
    :DiscreteParametric((DiscreteParametric&)forward , 'c' , ialloc_nb_value) {}

    void computation(const DiscreteParametric &dist);
};



class DiscreteDistributionData;
class ContinuousParametric;
class TimeEvents;
class Mixture;
class Convolution;
class Compound;

// class FrequencyDistribution : protected Reestimation<int> {
class FrequencyDistribution : public Reestimation<int> {  // loi discrete empirique

/*    friend class Distribution;
    friend class DiscreteParametric;
    friend class DiscreteParametricModel;
    friend class DiscreteDistributionData;
    friend class Compound;
    friend class CompoundData;
    friend class Convolution;
    friend class ConvolutionData;
    friend class Mixture;
    friend class MixtureData;
    friend class Vectors;
    friend class Regression;
    friend class VectorDistance;
    friend class Renewal;
    friend class TimeEvents;
    friend class RenewalData;
    friend class Curves;
    friend class ChainData;
    friend class NonparametricSequenceProcess;
    friend class Nonparametric_tree_process;
    friend class DiscreteParametricProcess;
    friend class NonhomogeneousMarkov;
    friend class NonhomogeneousMarkovData;
    friend class VariableOrderMarkov;
    friend class VariableOrderMarkovData;
    friend class HiddenVariableOrderMarkov;
    friend class SemiMarkov;
    friend class SemiMarkovData;
    friend class HiddenSemiMarkov;
    friend class Sequences;
    friend class SequenceCharacteristics;
    friend class MarkovianSequences;
    friend class TopParameters;
    friend class Tops; */

    friend DiscreteDistributionData* frequency_distribution_ascii_read(StatError &error ,
                                                                       const char *path);

    friend bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                           double *scale , int *dist_nb_value , int nb_histo ,
                           const FrequencyDistribution **histo);

// protected :
public :

    void shift(const FrequencyDistribution &histo , int shift_param);
    void cluster(const FrequencyDistribution &histo , int step , int mode);

    std::ostream& ascii_print(std::ostream &os , int comment_flag = false , bool cumul_flag = false) const;
    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_characteristic_print(std::ostream &os , bool shape = false) const;
    std::ostream& spreadsheet_circular_characteristic_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag = false ,
                                    bool concentration_flag = false) const;
    bool plot_print(const char *path , int nb_histo = 0 ,
                    const FrequencyDistribution **histo = NULL) const;
    bool plot_print(const char *path , double *cumul , double *concentration ,
                    double shift = 0.) const;
    bool survival_plot_print(const char *path , double *survivor) const;
    std::ostream& plot_title_print(std::ostream &os) const
    { return os; }

    void plotable_frequency_write(SinglePlot &plot) const;
    void plotable_mass_write(SinglePlot &plot) const;
    void plotable_cumul_write(SinglePlot &plot , double *icumul = NULL ,
                              double scale = D_DEFAULT) const;
    void plotable_cumul_matching_write(SinglePlot &plot , int reference_offset ,
                                       int  reference_nb_value , double *reference_cumul ,
                                       double *icumul = NULL) const;
    void plotable_concentration_write(SinglePlot &plot , double *icumul = NULL ,
                                      double scale = D_DEFAULT) const;
    void plotable_survivor_write(SinglePlot &plot) const;

    double* cumul_computation(double scale = D_DEFAULT) const;
    double concentration_computation() const;

    double* survivor_function_computation(double scale = D_DEFAULT) const;
    double* concentration_function_computation(double scale = D_DEFAULT) const;

    Test* kruskal_wallis_test(int nb_histo , const FrequencyDistribution **ihisto) const;

    std::ostream& dissimilarity_ascii_write(std::ostream &os , int nb_histo ,
                                            const FrequencyDistribution **ihisto ,
                                            int type , double **dissimilarity) const;
    bool dissimilarity_ascii_write(StatError &error , const char *path ,
                                   int nb_histo , const FrequencyDistribution **ihisto ,
                                   int type , double **dissimilarity) const;
    bool dissimilarity_spreadsheet_write(StatError &error , const char *path ,
                                         int nb_histo , const FrequencyDistribution **ihisto ,
                                         int type , double **dissimilarity) const;

    void update(const Reestimation<double> *reestim , int inb_element);
    FrequencyDistribution* frequency_scale(int inb_element) const;
    double* rank_computation() const;
    int cumulative_distribution_function_computation(double **cdf) const;
    int min_interval_computation() const;

    DiscreteParametric* parametric_estimation(int ident , int min_inf_bound = 0 , bool flag = true ,
                                              double cumul_threshold = CUMUL_THRESHOLD) const;

    double likelihood_computation(const ContinuousParametric &dist ,
                                  int min_interval = I_DEFAULT) const;

// public :

    FrequencyDistribution(int inb_value = 0)
    :Reestimation<int>(inb_value) {}
    FrequencyDistribution(const Distribution &dist)
    :Reestimation<int>(dist.nb_value) {}
    FrequencyDistribution(int inb_element , int *pelement);
    FrequencyDistribution(int nb_histo , const FrequencyDistribution **histo)
    :Reestimation<int>(nb_histo , (const Reestimation<int>**)histo) {}
    FrequencyDistribution(const FrequencyDistribution &histo , char transform ,
                          int param , int mode = FLOOR);

    bool operator==(const FrequencyDistribution&) const;
    bool operator!=(const FrequencyDistribution &histo) const { return !(*this == histo); }

    DiscreteDistributionData* shift(StatError &error , int shift_param) const;
    DiscreteDistributionData* cluster(StatError &error , int step , int mode = FLOOR) const;
    DiscreteDistributionData* cluster(StatError &error , double ratio , std::ostream &os) const;
    DiscreteDistributionData* cluster(StatError &error , int nb_class , int *ilimit) const;
    DiscreteDistributionData* transcode(StatError &error , int *isymbol) const;
    DiscreteDistributionData* value_select(StatError &error , int min_value ,
                                           int max_value , bool keep = true) const;

    TimeEvents* build_time_events(StatError &error , int itime) const;  // sequence_analysis

    bool ascii_write(StatError &error , const char *path) const;

    bool plot_write(StatError &error , const char *prefix , int nb_histo ,
                    const FrequencyDistribution **ihisto , const char *title) const;

    MultiPlotSet* get_plotable() const;
    MultiPlotSet* get_plotable_frequency_distributions(StatError &error , int nb_histo ,
                                                       const FrequencyDistribution **ihisto) const;

    std::ostream& survival_ascii_write(std::ostream &os) const;
    bool survival_ascii_write(StatError &error , const char *path) const;
    bool survival_spreadsheet_write(StatError &error , const char *path) const;
    bool survival_plot_write(StatError &error , const char *prefix ,
                             const char *title = NULL) const;
    MultiPlotSet* survival_get_plotable(StatError &error) const;

    bool comparison(StatError &error , std::ostream &os , int nb_histo ,
                    const FrequencyDistribution **ihisto , int type , const char *path = NULL ,
                    char format = 'a') const;

    void F_comparison(std::ostream &os , const FrequencyDistribution &histo) const;
    void t_comparison(std::ostream &os , const FrequencyDistribution &histo) const;
    bool wilcoxon_mann_whitney_comparison(StatError &error , std::ostream &os ,
                                          const FrequencyDistribution &ihisto) const;

    DiscreteParametricModel* fit(StatError &error , const DiscreteParametric &idist) const;

    DiscreteParametricModel* parametric_estimation(StatError &error , int ident ,
                                                   int min_inf_bound = 0 , bool flag = true ,
                                                   double cumul_threshold = CUMUL_THRESHOLD) const;
    DiscreteParametricModel* type_parametric_estimation(StatError &error ,
                                                        int min_inf_bound = 0 , bool flag = true ,
                                                        double cumul_threshold = CUMUL_THRESHOLD) const;

    Mixture* mixture_estimation(StatError &error , const Mixture &imixt , bool *estimate ,
                                int min_inf_bound = 0 , bool mixt_flag = true ,
                                bool component_flag = true , double weight_step = 0.1) const;
    Mixture* mixture_estimation(StatError &error , const Mixture &imixt ,
                                int min_inf_bound = 0 , bool mixt_flag = true ,
                                bool component_flag = true , double weight_step = 0.1) const;
    Mixture* mixture_estimation(StatError &error , int nb_component , int *ident ,
                                int min_inf_bound = 0 , bool mixt_flag = true ,
                                bool component_flag = true , double weight_step = 0.1) const;
    Mixture* mixture_estimation(StatError &error , std::ostream &os , int min_nb_component ,
                                int max_nb_component , int *ident , int min_inf_bound = 0 ,
                                bool mixt_flag = true , bool component_flag = true ,
                                int penalty_type = BICc , double weight_step = 0.1) const;

    Convolution* convolution_estimation(StatError &error , std::ostream &os , const DiscreteParametric &known_dist ,
                                        const DiscreteParametric &unknown_dist , int estimator = LIKELIHOOD ,
                                        int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Convolution* convolution_estimation(StatError &error , std::ostream &os , const DiscreteParametric &known_dist ,
                                        int min_inf_bound , int estimator = LIKELIHOOD ,
                                        int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;

    Compound* compound_estimation(StatError &error , std::ostream &os , const DiscreteParametric &sum_dist ,
                                  const DiscreteParametric &dist , char type , int estimator = LIKELIHOOD ,
                                  int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                  int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Compound* compound_estimation(StatError &error , std::ostream &os , const DiscreteParametric &known_dist ,
                                  char type , int min_inf_bound = 0 , int estimator = LIKELIHOOD ,
                                  int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                  int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;

    DiscreteParametricModel* estimation(StatError &error , std::ostream &os ,
                                        const FrequencyDistribution &backward ,
                                        const FrequencyDistribution &forward ,
                                        const FrequencyDistribution *no_event ,
                                        const DiscreteParametric &iinter_event ,
                                        int estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                                        int mean_computation_method = COMPUTED , double weight = D_DEFAULT ,
                                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO ,
                                        double iinter_event_mean = D_DEFAULT) const;
    DiscreteParametricModel* estimation(StatError &error , std::ostream &os ,
                                        const FrequencyDistribution &backward ,
                                        const FrequencyDistribution &forward ,
                                        const FrequencyDistribution *no_event ,
                                        int estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                                        int mean_computation_method = COMPUTED , double weight = D_DEFAULT ,
                                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;

    // acces membres de la classe

/*    int get_nb_value() const { return nb_value; }
    int get_alloc_nb_value() const { return alloc_nb_value; }
    int get_offset() const { return offset; }
    int get_nb_element() const { return nb_element; }
    int get_max() const { return max; }
    double get_mean() const { return mean; }
    double get_variance() const { return variance; }
    int& operator[](int index) const { return frequency[index]; }
    int get_frequency(int index) const { return frequency[index]; } */
};


DiscreteDistributionData* frequency_distribution_ascii_read(StatError &error , const char *path);



class Histogram;

class ContinuousParametric {  // loi de probabilite continue parametrique

    friend ContinuousParametric* continuous_parametric_parsing(StatError &error ,
                                                               std::ifstream &in_file ,
                                                               int &line);
//    friend std::ostream& operator<<(std::ostream &os , const ContinuousParametric &dist)
//    { return dist.ascii_print(os); }

// private :
public :

    int ident;              // identificateur
    double location;        // parametre de moyenne
    double dispersion;      // parametre de dispersion - ecart-type (GAUSSIAN),
                            // concentration (VON_MISES)
    double min_value;       // valeur minimum
    double max_value;       // valeur maximum
    int unit;               // unite (degre/radian) pour la loi de von Mises
    double *cumul;          // fonction de repartition (loi de von Mises)

    void copy(const ContinuousParametric &dist);

    std::ostream& ascii_parameter_print(std::ostream &os) const;
    std::ostream& ascii_characteristic_print(std::ostream &os ,
                                             bool file_flag = false) const;
    std::ostream& ascii_print(std::ostream &os , bool file_flag = false ,
                              bool cumul_flag = false , const Histogram *histo1 = NULL ,
                              const FrequencyDistribution *histo2 = NULL);
    std::ostream& spreadsheet_parameter_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag = false ,
                                    const Histogram *histo1 = NULL ,
                                    const FrequencyDistribution *histo2 = NULL);
    std::ostream& spreadsheet_characteristic_print(std::ostream &os) const;
    std::ostream& plot_title_print(std::ostream &os) const;
    bool plot_print(const char *path , const Histogram *histo1 = NULL ,
                    const FrequencyDistribution *histo2 = NULL);
    bool q_q_plot_print(const char *path , int nb_value ,
                        double **empirical_cdf) const;
    void plotable_write(SinglePlot &plot , const Histogram *histo1 = NULL ,
                        const FrequencyDistribution *histo2 = NULL);
    void q_q_plotable_write(SinglePlot &plot , int nb_value ,
                            double **empirical_cdf) const;

    double** q_q_plot_computation(int nb_value , double **cdf) const;

// public :

    ContinuousParametric(int iident = GAUSSIAN , double ilocation = D_INF ,
                         double idispersion = D_DEFAULT , int iunit = I_DEFAULT);
    ContinuousParametric(const ContinuousParametric &dist)
    { copy(dist); }
    ~ContinuousParametric();
    ContinuousParametric& operator=(const ContinuousParametric&);

    int nb_parameter_computation() const;

    double von_mises_mass_computation(double inf , double sup) const;
    void von_mises_cumul_computation();
    double mass_computation(double inf , double sup) const;

    double likelihood_computation(const FrequencyDistribution &histo , int min_interval) const
    { return histo.likelihood_computation(*this , min_interval); }

    double simulation() const;
};


ContinuousParametric* continuous_parametric_parsing(StatError &error ,
                                                    std::ifstream &in_file ,
                                                    int &line);



class Histogram {       // histogramme

// private :
public :

    int nb_element;         // effectif total
    int nb_category;        // nombre de categories
    double step;            // pas de regroupement
    int max;                // frequence maximum
    int *frequency;         // frequences
    int type;               // type de la variable (INT_VALUE/REAL_VALUE)
    double min_value;       // valeur minimum
    double max_value;       // valeur maximum

    void copy(const Histogram &histo);

    std::ostream& ascii_print(std::ostream &os , bool comment_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    bool plot_print(const char *path) const;
    void plotable_write(SinglePlot &plot) const;

    void max_computation();

// public :

    Histogram(int inb_category = 0 , bool init_flag = true);
    Histogram(const FrequencyDistribution &histo);
    Histogram(const Histogram &histo)
    { copy(histo); }
    ~Histogram();
    Histogram& operator=(const Histogram &histo);

    double* cumul_computation() const;
};



#include "reestimation.cpp"



#endif
