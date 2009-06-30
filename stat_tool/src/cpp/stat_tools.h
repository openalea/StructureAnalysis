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



#ifndef STAT_TOOLS_H
#define STAT_TOOLS_H



#include <fstream>
// #include <rw/collect.h>
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
const double D_INF = -1.e37;           // plus petit nombre flottant
const double DOUBLE_ERROR = 1.e-6;     // erreur sur une somme de doubles

enum {
  STANDARD_NORMAL ,
  CHI2 ,
  FISHER ,
  STUDENT
};

const int NB_CRITICAL_PROBABILITY = 2;
const double ref_critical_probability[NB_CRITICAL_PROBABILITY] = {0.05 , 0.01};

const int NB_VALUE = 1000;             // nombre de valeurs prises par une v.a.
const int SAMPLE_NB_VALUE = NB_VALUE;  // nombre de valeurs d'un histogramme

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

const int MAX_INF_BOUND = 10000;       // borne inferieure maximum
const int MAX_DIFF_BOUND = 10000;      // difference maximum entre les bornes
                                       // inferieure et superieure
const double MAX_MEAN = 10000.;        // moyenne maximum

const double B_PROBABILITY = 0.8;      // seuil pour utiliser le calcul descendant de la loi binomiale
const double B_THRESHOLD = 1000.;      // seuil pour utiliser le calcul en log de la loi binomiale
const double P_THRESHOLD = 90.;        // seuil pour utiliser le calcul en log de la loi de Poisson
const double NB_THRESHOLD = 500.;      // seuil pour utiliser le calcul en log de la loi binomiale negative

const double SAMPLE_NB_VALUE_COEFF = 5.;  // facteur pour deduire le nombre de valeurs prises
                                          // par une v.a. du nombre de valeurs d'un histogramme
const int INF_BOUND_MARGIN = 5;        // plage de recherche pour la borne inferieure
const int SUP_BOUND_MARGIN = 3;        // plage de recherche pour la borne superieure
const double POISSON_RATIO = 0.7;      // rapport moyenne/variance minimum pour      
                                        // estimer une loi de Poisson    
                                        //       
const double POISSON_RANGE = 0.1;      // plage de variation pour choisir une loi de   
                                       // Poisson par dilatation de l'echelle des temps
const double NB_VALUE_COEFF = 2.;      // facteur pour deduire le nombre de valeurs prises par une v.a.
                                       // du nombre de valeurs prises par une v.a. initiale

const int MIN_RANGE = 10;              // intervalle de definition minimum
                                       // pour appliquer la methode du rejet
const double MAX_SURFACE = 3.0;        // surface maximum pour appliquer la methode du rejet
const int DIST_NB_ELEMENT = 1000000;   // taille maximum de l'echantillon pour la simulation

const int CHI2_FREQUENCY = 2;          // effectif theorique minimum pour un
                                       // test d'ajustement du Chi2
const double MIN_T_VALUE = 2.5;        // seuil sur la variable t

const int MARGINAL_MAX_VALUE = 20000;  // valeur maximum pour la construction de la loi marginale
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
 *  Identificateurs des classes pour la persistance :
 */


/* enum {
  STATI_PARAMETRIC_MODEL = 500 ,
  STATI_MIXTURE ,
  STATI_CONVOLUTION ,
  STATI_COMPOUND ,
  STATI_DISTRIBUTION_DATA ,
  STATI_MIXTURE_DATA ,
  STATI_CONVOLUTION_DATA ,
  STATI_COMPOUND_DATA ,

  STATI_VECTORS ,
  STATI_REGRESSION ,

  STATI_VECTOR_DISTANCE ,

  STATI_DISTANCE_MATRIX ,
  STATI_CLUSTERS ,

  STATI_RENEWAL ,
  STATI_TIME_EVENTS ,
  STATI_RENEWAL_DATA ,

  STATI_MARKOV ,
  STATI_SEMI_MARKOV ,
  STATI_HIDDEN_MARKOV ,
  STATI_HIDDEN_SEMI_MARKOV ,
  STATI_SEQUENCES ,
  STATI_MARKOVIAN_SEQUENCES ,
  STATI_MARKOV_DATA ,
  STATI_SEMI_MARKOV_DATA ,
  STATI_CORRELATION ,

  STATI_TOP_PARAMETERS ,
  STATI_TOPS
}; */



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



class Format_error {    // erreurs

    friend std::ostream& operator<<(std::ostream &os , const Format_error &error)
    { return error.ascii_write(os , ERROR); }

private :

    int nb_error;           // nombre d'erreurs
    int max_nb_error;       // nombre maximum d'erreurs
    int *line;              // no de ligne
    int *column;            // no de colonne
    char **label;           // messages d'erreur

public :

    Format_error(int imax_nb_error = NB_ERROR);
    ~Format_error();

    std::ostream& ascii_write(std::ostream &os , int type = ERROR) const;

    void init() { nb_error = 0; }
    void update(const char *ilabel , int iline = 0 , int icolumn = 0);
    void correction_update(const char *ilabel , const char *correction ,
                           int iline  = 0 , int icolumn = 0);
    void correction_update(const char *ilabel , int correction ,
                           int iline  = 0 , int icolumn = 0);

    // acces membres de la classe

    int get_nb_error() const { return nb_error; }
    int get_max_nb_error() const { return max_nb_error; }
};



class STAT_interface {  // classe abstraite pour les interfaces
// class STAT_interface : public RWCollectable {

public :

    virtual std::ostream& line_write(std::ostream &os) const = 0;

    virtual std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const = 0;
//    virtual std::ostream& spreadsheet_write(std::ostream &os) const = 0;

    virtual bool ascii_write(Format_error &error , const char *path ,
                             bool exhaustive = false) const = 0;
    virtual bool spreadsheet_write(Format_error &error , const char *path) const = 0;
    virtual bool plot_write(Format_error &error , const char *prefix ,
                            const char *title = 0) const = 0;

    virtual plotable::MultiPlotSet* get_plotable() const { return NULL; };

//    bool binary_write(Format_error &error , const char *path) const;
};



class Histogram;
class Parametric_model;

class Distribution {    // loi de probabilite discrete

/*    template <typename Type> friend class Reestimation;  probleme Windows
    friend class Histogram;
    friend class Parametric_model;
    friend class Distribution_data;
    friend class Convolution;
    friend class Compound_data;
    friend class Vectors;
    friend class Vector_distance;
    friend class Backward;
    friend class Renewal;
    friend class Time_events;
    friend class Curves;
    friend class Chain_data;
    friend class Nonparametric_process;
    friend class Nonparametric_sequence_process;
    friend class Markov;
    friend class Hidden_markov;
    friend class Variable_order_markov;
    friend class Hidden_variable_order_markov;
    friend class Semi_markov;
    friend class Hidden_semi_markov;
    friend class Sequences;
    friend class Markovian_sequences;
    friend class Correlation;
    friend class Top_parameters; */

    friend std::ostream& operator<<(std::ostream& , const Distribution&);
    friend bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                           double *scale , int *dist_nb_value , int nb_histo ,
                           const Histogram **histo , int *index_dist);

// protected :
  public :

    int nb_value;           // nombre de valeurs a partir de 0
    int alloc_nb_value;     // nombre de valeurs allouees
    int offset;             // nombre de valeurs de probabilite nulle a partir de 0
    double max;             // valeur de probabilite maximum
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
    std::ostream& ascii_print(std::ostream &os , bool comment_flag , bool cumul_flag ,
                              bool nb_value_flag , const Histogram *histo = 0) const;

    std::ostream& ascii_print(std::ostream &os , int nb_dist , const Distribution **dist ,
                              double *dist_scale , bool comment_flag , bool cumul_flag ,
                              const Histogram *histo = 0) const;

    std::ostream& spreadsheet_characteristic_print(std::ostream &os , bool shape = false) const;

    std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag , bool concentration_flag ,
                                    bool nb_value_flag , const Histogram *histo = 0) const;

    std::ostream& spreadsheet_print(std::ostream &os , int nb_dist , const Distribution **dist ,
                                    double *dist_scale , bool cumul_flag ,
                                    const Histogram *histo = 0) const;

    int plot_nb_value_computation(const Histogram *histo = 0) const;
    bool plot_print(const char *path , double *concentration , double scale) const;
    bool plot_print(const char *path , const Histogram *histo = 0) const;
    virtual std::ostream& plot_title_print(std::ostream &os) const
    { return os; }
    bool survival_plot_print(const char *path , double *survivor) const;
    std::ostream& print(std::ostream&) const;

    void plotable_mass_write(SinglePlot &plot , double scale = 1.) const;
    void plotable_cumul_write(SinglePlot &plot) const;
    void plotable_cumul_matching_write(SinglePlot &plot , const Distribution &reference_dist) const;
    void plotable_concentration_write(SinglePlot &plot) const;
    void plotable_survivor_write(SinglePlot &plot) const;

/*    RWspace binaryStoreSize(int ialloc_nb_value = I_DEFAULT) const;
    void restoreGuts(RWvistream &is);
    void restoreGuts(RWFile &file);
    void saveGuts(RWvostream &os , int ialloc_nb_value = I_DEFAULT) const;
    void saveGuts(RWFile &file , int ialloc_nb_value = I_DEFAULT) const; */

    void convolution(Distribution &dist1 , Distribution &dist2 ,
                     int inb_value = I_DEFAULT);

    void nb_value_computation();
    void offset_computation();
    double concentration_computation() const;

    void cumul_computation();
    double* survivor_function_computation() const;
    double* concentration_function_computation() const;
    void log_computation();

    double survivor_likelihood_computation(const Histogram &histo) const;
    double chi2_value_computation(const Histogram &histo) const;
    void chi2_degree_of_freedom(const Histogram &histo , Test &test) const;

    void penalty_computation(double weight , int type , double *penalty , int outside) const;

// public :

    Distribution(int inb_value = 0);
    Distribution(const Distribution &dist , double scaling_coeff);
    Distribution(const Histogram &histo);
    Distribution(const Distribution &dist , char transform = 'c' ,
                 int ialloc_nb_value = I_DEFAULT);
    virtual ~Distribution();
    Distribution& operator=(const Distribution&);
    bool operator==(const Distribution&) const;
    bool operator!=(const Distribution &dist) const { return !(*this == dist); }

    bool plot_write(Format_error &error , const char *prefix , int nb_dist ,
                    const Distribution **idist , const char *title) const;
    MultiPlotSet* get_plotable() const;
    MultiPlotSet* get_plotable_distributions(Format_error &error , int nb_dist ,
                                             const Distribution **idist) const;

    std::ostream& survival_ascii_write(std::ostream &os) const;
    bool survival_ascii_write(Format_error &error , const char *path) const;
    bool survival_spreadsheet_write(Format_error &error , const char *path) const;
    bool survival_plot_write(Format_error &error , const char *prefix ,
                             const char *title = 0) const;
    MultiPlotSet* survival_get_plotable(Format_error &error) const;

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
    void chi2_fit(const Histogram &histo , Test &test) const;
    int simulation() const;

    Parametric_model* truncate(Format_error &error , int imax_value) const;

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
                const Histogram **histo , int *index_dist);



class Forward;

class Parametric : public Distribution {  // loi de probabilite parametrique

/*    template <typename Type> friend class Reestimation;  probleme Windows
    friend class Histogram;
    friend class Distribution_data;
    friend class Compound;
    friend class Convolution;
    friend class Mixture;
    friend class Backward;
    friend class Forward;
    friend class Length_bias;
    friend class Nb_event;
    friend class Renewal;
    friend class Time_events;
    friend class Renewal_data;
    friend class Nonparametric_sequence_process;
    friend class Parametric_process;
    friend class Markov;
    friend class Variable_order_markov;
    friend class Semi_markov;
    friend class Hidden_semi_markov;
    friend class Markovian_sequences;
    friend class Top_parameters; */

    friend int nb_value_computation(int ident , int inf_bound , int sup_bound ,
                                    double parameter , double probability ,
                                    double cumul_threshold);
    friend Parametric* parametric_parsing(Format_error &error , std::ifstream &in_file ,
                                          int &line , int last_ident,
                                          double cumul_threshold,
                                          int min_inf_bound);
    friend std::ostream& operator<<(std::ostream& , const Parametric&);

// protected :
public :

    int ident;              // identificateur
    int inf_bound;          // borne inferieure
    int sup_bound;          // borne superieure (binomiale / uniforme)
    double parameter;       // parametre (Poisson / binomiale negative)
    double probability;     // probabilite de succes (binomiale / binomiale negative)

    void init(int iinf_bound , int isup_bound , double iparameter , double iprobability);
    void init(int iident , int iinf_bound , int isup_bound , double iparameter , double iprobability);
    void copy(const Parametric &dist);
    std::ostream& ascii_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    std::ostream& plot_title_print(std::ostream &os) const;

    void nb_parameter_update();

/*    RWspace binaryStoreSize(int ialloc_nb_value = I_DEFAULT) const;
    void restoreGuts(RWvistream &is);
    void restoreGuts(RWFile &file);
    void saveGuts(RWvostream &os , int ialloc_nb_value = I_DEFAULT) const;
    void saveGuts(RWFile &file , int ialloc_nb_value = I_DEFAULT) const; */

    void binomial_computation(int inb_value , char mode);
    void poisson_computation(int inb_value , double cumul_threshold ,
                             char mode);
    void negative_binomial_computation(int inb_value , double cumul_threshold ,
                                       char mode);
    void uniform_computation();

    double renewal_likelihood_computation(const Forward &forward_dist , const Histogram &within ,  // sequence_analysis
                                          const Histogram &backward , const Histogram &forward ,
                                          const Histogram *no_event) const;
    void expectation_step(const Histogram &within ,const Histogram &backward ,  // sequence_analysis
                          const Histogram &forward , const Histogram *no_event ,
                          Reestimation<double> *inter_event_reestim ,
                          Reestimation<double> *length_bias_reestim , int iter) const;

    void reestimation(const Reestimation<double> *reestim , int nb_estim = 1);

    double state_occupancy_likelihood_computation(const Histogram &sojourn_time ,  // sequence_analysis
                                                  const Histogram &final_run) const;
    double state_occupancy_likelihood_computation(const Forward &forward , const Histogram &sojourn_time ,  // sequence_analysis
                                                  const Histogram &final_run , const Histogram &initial_run ,
                                                  const Histogram &single_run) const;
    void expectation_step(const Histogram &sojourn_time , const Histogram &final_run ,  // sequence_analysis
                          Reestimation<double> *occupancy_reestim , int iter) const;
    void expectation_step(const Histogram &sojourn_time ,const Histogram &final_run ,  // sequence_analysis
                          const Histogram &initial_run , const Histogram &single_run ,
                          Reestimation<double> *occupancy_reestim ,
                          Reestimation<double> *length_bias_reestim , int iter ,
                          bool combination = false , int mean_computation = COMPUTED) const;

// public :

    Parametric(int inb_value = 0 , int iident = NONPARAMETRIC ,
               int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
               double iparameter = D_DEFAULT , double iprobability = D_DEFAULT);
    Parametric(int iident , int iinf_bound , int isup_bound , double iparameter ,
               double iprobability , double cumul_threshold = CUMUL_THRESHOLD);
    Parametric(const Distribution &dist , int ialloc_nb_value = I_DEFAULT);
    Parametric(const Distribution &dist , double scaling_coeff);
    Parametric(const Parametric &dist , double scaling_coeff);
    Parametric(const Histogram &histo);
    Parametric(const Parametric &dist , char transform = 'c' ,
               int ialloc_nb_value = I_DEFAULT);
    Parametric& operator=(const Parametric &dist);

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
Parametric* parametric_parsing(Format_error &error , std::ifstream &in_file ,
                               int &line , int last_ident = NEGATIVE_BINOMIAL ,
                               double cumul_threshold = CUMUL_THRESHOLD ,
                               int min_inf_bound = 0);
std::ostream& operator<<(std::ostream& , const Parametric&);



class Forward : public Parametric {  // loi de l'intervalle de temps residuel

/*    friend class Renewal;
    friend class Renewal_data;
    friend class Semi_markov; */

public :

    Forward(int inb_value = 0 , int iident = NONPARAMETRIC ,
            int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
            double iparameter = D_DEFAULT , double iprobability = D_DEFAULT)
    :Parametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Forward(const Parametric &dist , int ialloc_nb_value = I_DEFAULT)
    :Parametric(dist , 'c' , ialloc_nb_value) { computation(dist); }
    Forward(const Forward &forward , int ialloc_nb_value = I_DEFAULT)
    :Parametric((Parametric&)forward , 'c' , ialloc_nb_value) {}

    void computation(const Parametric &dist);
};



class Distribution_data;
class Time_events;
class Mixture;
class Convolution;
class Compound;

// class Histogram : protected Reestimation<int> {
class Histogram : public Reestimation<int> {  // histogramme

/*    friend class Distribution;
    friend class Parametric;
    friend class Parametric_model;
    friend class Distribution_data;
    friend class Compound;
    friend class Compound_data;
    friend class Convolution;
    friend class Convolution_data;
    friend class Mixture;
    friend class Mixture_data;
    friend class Vectors;
    friend class Regression;
    friend class Vector_distance;
    friend class Renewal;
    friend class Time_events;
    friend class Renewal_data;
    friend class Curves;
    friend class Chain_data;
    friend class Nonparametric_sequence_process;
    friend class Nonparametric_tree_process;
    friend class Parametric_process;
    friend class Markov;
    friend class Markov_data;
    friend class Hidden_markov;
    friend class Variable_order_markov;
    friend class Variable_order_markov_data;
    friend class Hidden_variable_order_markov;
    friend class Semi_markov;
    friend class Semi_markov_data;
    friend class Hidden_semi_markov;
    friend class Sequences;
    friend class Sequence_characteristics;
    friend class Markovian_sequences;
    friend class Top_parameters;
    friend class Tops; */

    friend Distribution_data* histogram_ascii_read(Format_error &error , const char *path);

    friend bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                           double *scale , int *dist_nb_value , int nb_histo ,
                           const Histogram **histo , int *index_dist);

// protected :
public :

    void shift(const Histogram &histo , int shift_param);
    void cluster(const Histogram &histo , int step , int mode);

    std::ostream& ascii_print(std::ostream &os , int comment_flag = false , bool cumul_flag = false) const;
    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_characteristic_print(std::ostream &os , bool shape = false) const;
    std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag = false ,
                                    bool concentration_flag = false) const;
    bool plot_print(const char *path , int nb_histo = 0 ,
                    const Histogram **histo = 0) const;
    bool plot_print(const char *path , double *cumul , double *concentration ,
                    double shift = 0.) const;
    bool survival_plot_print(const char *path , double *survivor) const;
    std::ostream& plot_title_print(std::ostream &os) const
    { return os; }

    void plotable_frequency_write(SinglePlot &plot) const;
    void plotable_mass_write(SinglePlot &plot) const;
    void plotable_cumul_write(SinglePlot &plot , double *icumul = 0 ,
                              double scale = D_DEFAULT) const;
    void plotable_cumul_matching_write(SinglePlot &plot , int reference_offset ,
                                       int  reference_nb_value , double *reference_cumul ,
                                       double *icumul = 0) const;
    void plotable_concentration_write(SinglePlot &plot , double *icumul = 0 ,
                                      double scale = D_DEFAULT) const;
    void plotable_survivor_write(SinglePlot &plot) const;

/*    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    double* cumul_computation(double scale = D_DEFAULT) const;
    double concentration_computation() const;

    double* survivor_function_computation(double scale = D_DEFAULT) const;
    double* concentration_function_computation(double scale = D_DEFAULT) const;

    Test* kruskal_wallis_test(int nb_histo , const Histogram **ihisto) const;

    std::ostream& dissimilarity_ascii_write(std::ostream &os , int nb_histo ,
                                            const Histogram **ihisto ,
                                            int type , double **dissimilarity) const;
    bool dissimilarity_ascii_write(Format_error &error , const char *path ,
                                   int nb_histo , const Histogram **ihisto ,
                                   int type , double **dissimilarity) const;
    bool dissimilarity_spreadsheet_write(Format_error &error , const char *path ,
                                         int nb_histo , const Histogram **ihisto ,
                                         int type , double **dissimilarity) const;

    void update(const Reestimation<double> *reestim , int inb_element);
    Histogram* frequency_scale(int inb_element) const;
    double* rank_computation() const;

    Parametric* parametric_estimation(int ident , int min_inf_bound = 0 , bool flag = true ,
                                      double cumul_threshold = CUMUL_THRESHOLD) const;

// public :

    Histogram(int inb_value = 0)
    :Reestimation<int>(inb_value) {}
    Histogram(const Distribution &dist)
    :Reestimation<int>(dist.nb_value) {}
    Histogram(int inb_element , int *pelement);
    Histogram(int nb_histo , const Histogram **histo)
    :Reestimation<int>(nb_histo , (const Reestimation<int>**)histo) {}
    Histogram(const Histogram &histo , char transform , int param , int mode = FLOOR);

    bool operator==(const Histogram&) const;
    bool operator!=(const Histogram &histo) const { return !(*this == histo); }

    Distribution_data* shift(Format_error &error , int shift_param) const;
    Distribution_data* cluster(Format_error &error , int step , int mode = FLOOR) const;
    Distribution_data* cluster(Format_error &error , double ratio , std::ostream &os) const;
    Distribution_data* cluster(Format_error &error , int nb_class , int *ilimit) const;
    Distribution_data* transcode(Format_error &error , int *isymbol) const;
    Distribution_data* value_select(Format_error &error , int min_value ,
                                    int max_value , bool keep = true) const;

    Time_events* build_time_events(Format_error &error , int itime) const;  // sequence_analysis

    bool ascii_write(Format_error &error , const char *path) const;

    bool plot_write(Format_error &error , const char *prefix , int nb_histo ,
                    const Histogram **ihisto , const char *title) const;

    MultiPlotSet* get_plotable() const;
    MultiPlotSet* get_plotable_histograms(Format_error &error , int nb_histo ,
                                          const Histogram **ihisto) const;

    std::ostream& survival_ascii_write(std::ostream &os) const;
    bool survival_ascii_write(Format_error &error , const char *path) const;
    bool survival_spreadsheet_write(Format_error &error , const char *path) const;
    bool survival_plot_write(Format_error &error , const char *prefix ,
                             const char *title = 0) const;
    MultiPlotSet* survival_get_plotable(Format_error &error) const;

    bool comparison(Format_error &error , std::ostream &os , int nb_histo ,
                    const Histogram **ihisto , int type , const char *path = 0 ,
                    char format = 'a') const;

    void F_comparison(std::ostream &os , const Histogram &histo) const;
    void t_comparison(std::ostream &os , const Histogram &histo) const;
    bool wilcoxon_mann_whitney_comparison(Format_error &error , std::ostream &os ,
                                          const Histogram &ihisto) const;

    Parametric_model* fit(Format_error &error , const Parametric &idist) const;

    Parametric_model* parametric_estimation(Format_error &error , int ident ,
                                            int min_inf_bound = 0 , bool flag = true ,
                                            double cumul_threshold = CUMUL_THRESHOLD) const;
    Parametric_model* type_parametric_estimation(Format_error &error ,
                                                 int min_inf_bound = 0 , bool flag = true ,
                                                 double cumul_threshold = CUMUL_THRESHOLD) const;

    Mixture* mixture_estimation(Format_error &error , const Mixture &imixt , bool *estimate ,
                                int min_inf_bound = 0 , bool mixt_flag = true ,
                                bool component_flag = true , double weight_step = 0.1) const;
    Mixture* mixture_estimation(Format_error &error , const Mixture &imixt ,
                                int min_inf_bound = 0 , bool mixt_flag = true ,
                                bool component_flag = true , double weight_step = 0.1) const;
    Mixture* mixture_estimation(Format_error &error , int nb_component , int *ident ,
                                int min_inf_bound = 0 , bool mixt_flag = true ,
                                bool component_flag = true , double weight_step = 0.1) const;
    Mixture* mixture_estimation(Format_error &error , std::ostream &os , int min_nb_component ,
                                int max_nb_component , int *ident , int min_inf_bound = 0 ,
                                bool mixt_flag = true , bool component_flag = true ,
                                int penalty_type = BICc , double weight_step = 0.1) const;

    Convolution* convolution_estimation(Format_error &error , std::ostream &os , const Parametric &known_dist ,
                                        const Parametric &unknown_dist , int estimator = LIKELIHOOD ,
                                        int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Convolution* convolution_estimation(Format_error &error , std::ostream &os , const Parametric &known_dist ,
                                        int min_inf_bound , int estimator = LIKELIHOOD ,
                                        int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;

    Compound* compound_estimation(Format_error &error , std::ostream &os , const Parametric &sum_dist ,
                                  const Parametric &dist , char type , int estimator = LIKELIHOOD ,
                                  int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                  int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Compound* compound_estimation(Format_error &error , std::ostream &os , const Parametric &known_dist ,
                                  char type , int min_inf_bound = 0 , int estimator = LIKELIHOOD ,
                                  int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                  int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;

    Parametric_model* estimation(Format_error &error , std::ostream &os , const Histogram &backward ,  // sequence_analysis
                                 const Histogram &forward , const Histogram *no_event ,
                                 const Parametric &iinter_event , int estimator = LIKELIHOOD ,
                                 int nb_iter = I_DEFAULT , int mean_computation = COMPUTED ,
                                 double weight = D_DEFAULT , int penalty_type = SECOND_DIFFERENCE ,
                                 int outside = ZERO , double iinter_event_mean = D_DEFAULT) const;
    Parametric_model* estimation(Format_error &error , std::ostream &os , const Histogram &backward ,  // sequence_analysis
                                 const Histogram &forward , const Histogram *no_event ,
                                 int estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                                 int mean_computation = COMPUTED , double weight = D_DEFAULT ,
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


Distribution_data* histogram_ascii_read(Format_error &error , const char *path);

bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                double *scale , int *dist_nb_value , int nb_histo ,
                const Histogram **histo , int *index_dist);



#include "reestimation.cpp"



#endif
