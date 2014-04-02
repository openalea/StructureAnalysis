/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2014 CIRAD/INRA/Inria Virtual Plants
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



#ifndef SEQUENCES_H
#define SEQUENCES_H



/****************************************************************
 *
 *  Constantes :
 */


const int DEFAULT_LENGTH = 20;         // longueur par defaut d'une sequence

const int SEQUENCE_NB_VARIABLE = 30;   // nombre maximum de variables observees

const int PLOT_NB_SEQUENCE = 200;      // nombre maximum de sequences visualisees (sortie Gnuplot)
const int PLOT_LEGEND_NB_SEQUENCE = 15;  // nombre maximum de sequences identifiees (sortie Gnuplot)

enum {
  IMPLICIT_TYPE ,                      // parametre d'index implicite
  TIME ,                               // parametre d'index temps
  TIME_INTERVAL ,                      // parametre d'index intervalle de temps
  POSITION ,                           // parametre d'index position
  POSITION_INTERVAL                    // parametre d'index intervalle inter-position
};

enum {
  OBSERVED_VALUE ,
  OBSERVED_STATE ,
  THEORETICAL_STATE ,
  OBSERVED_OUTPUT ,
  THEORETICAL_OUTPUT
};

enum {
  UNCHANGED ,
  ADD_INITIAL_RUN ,
  REMOVE_INITIAL_RUN
};

enum {
  CTM_BIC ,                            // algorithme Context Tree Maximizing/BIC
  CTM_KT ,                             // algorithme Context Tree Maximizing/Krichevsky-Trofimov
  LOCAL_BIC ,                          // algorithme d'elagage recursif/BIC 
  CONTEXT                              // algorithme Context
};

enum {
  MAXIMUM_LIKELIHOOD ,
  LAPLACE ,
  ADAPTATIVE_LAPLACE ,
  UNIFORM_SUBSET ,
  UNIFORM_CARDINALITY
};

enum {
  CATEGORICAL_CHANGE ,
  POISSON_CHANGE ,
  MULTIVARIATE_POISSON_CHANGE ,
  GEOMETRIC_0_CHANGE ,
  GEOMETRIC_1_CHANGE ,
  MULTIVARIATE_GEOMETRIC_0_CHANGE ,
  ORDINAL_GAUSSIAN_CHANGE ,
  GAUSSIAN_CHANGE ,
  MEAN_CHANGE ,
  VARIANCE_CHANGE ,
  MEAN_VARIANCE_CHANGE ,
  LINEAR_MODEL_CHANGE ,
  INTERCEPT_SLOPE_CHANGE ,
  BAYESIAN_POISSON_CHANGE ,
  BAYESIAN_GAUSSIAN_CHANGE
};

const double PRIOR_VARIANCE_FACTOR = 100.;  // facteur pour deduire la variance
                                            // de la loi a priori gamma
const double PRIOR_SAMPLE_SIZE = 1.;  // hyperparametre de la loi gaussiennne-gamma
const double PRIOR_DEGREES_OF_FREEDOM = 2.;  // hyperparametre de la loi gaussiennne-gamma
const double PRIOR_DISPERSION_FACTOR = 10.;  // facteur pour deduire la somme des carres
                                             // des ecarts de la loi a priori
const int PRIOR_PRECISION = 2;         // precision hyperparametres

const double MAX_NB_WORD = 1.e7;       // nombre maximum de mots

const double STATIONARY_PROBABILITY_THRESHOLD = 1.e-8;  // seuil pour le calcul des probabilites
                                                        // stationnaires d'un modele en equilibre
const int STATIONARY_PROBABILITY_LENGTH = 10000;  // longueur maximum pour le calcul des probabilites
                                                  // stationnaires d'un modele en equilibre
const double LEAVE_INCREMENT = 1.e-6;  // seuil pour stopper le calcul de la probabilite
                                       // de quitter un etat/observation sans possibilite d'y revenir

const double CTM_BIC_THRESHOLD = 6.;   // seuil pour elaguer les memoires
const double CTM_KT_THRESHOLD = 12.;   // seuil pour elaguer les memoires
const double LOCAL_BIC_THRESHOLD = 10.;  // seuil pour elaguer les memoires
const double CONTEXT_THRESHOLD = 5.;  // seuil pour elaguer les memoires

const int MIN_NB_STATE_SEQUENCE = 1;   // nombre de sequences d'etats 1ere iteration de l'algorithme MCEM
const int MAX_NB_STATE_SEQUENCE = 10;  // nombre de sequences d'etats maximum pour l'algorithme MCEM
const double NB_STATE_SEQUENCE_PARAMETER = 1.;  // parametre nombre de sequences d'etats pour l'algorithme MCEM

const double OCCUPANCY_THRESHOLD = 0.99999;  // seuil sur la fonction de repartition
                                             // pour borner une loi d'occupation d'un etat
const double OCCUPANCY_MEAN = 10.;     // temps moyen d'occupation d'un etat

const int POSTERIOR_PROBABILITY_NB_SEQUENCE = 300;  // nombre maximum de sequences pour la sortie des probabilites
                                                    // a posteriori des sequences d'etats les plus probables

const int NB_STATE_SEQUENCE = 10;      // nombre de sequences d'etats calculees

const int COUNT_MAX_LENGTH = 10000;    // longueur maximum des sequences pour l'extraction des comptages
const int COUNTING_MAX_LENGTH = 500;   // longueur maximum des sequences pour le calcul des lois de comptage

const int NB_SEQUENCE = 100000;        // nombre maximum de sequences simulees
const int MAX_LENGTH = 1000000;        // longueur maximum des sequences simulees
const int CUMUL_LENGTH = 1000000;      // longueur maximum cumulee des sequences simulees

enum {
  SEQUENCE ,
  TREND ,
  SUBTRACTION_RESIDUAL ,
  DIVISION_RESIDUAL ,
  STANDARDIZED_RESIDUAL ,
  SEGMENTATION_ENTROPY ,
  SEGMENTATION_DIVERGENCE
};

enum {
  ADAPTATIVE ,
  FIXED
};

enum {
  DELETION ,
  BEGIN_END_DELETION ,
  INSERTION ,
  BEGIN_END_INSERTION ,
  MATCH ,
  SUBSTITUTION ,
  TRANSPOSITION
};

enum {
  DATA ,
  GAP ,
  BEGIN_END_GAP
};

const int NB_ALIGNMENT = 1000000;      // nombre maximum d'alignements
const int DISPLAY_NB_ALIGNMENT = 30;   // nombre maximum d'alignements
                                       // pour la sortie detaillee ecran
const int FILE_NB_ALIGNMENT = 300;     // nombre maximum d'alignements
                                       // pour la sortie detaillee fichier

const double INDEL_FACTOR_1 = 0.51;    // facteur pour deduire le cout d'elision/insertion -
                                       // alignement simple
const double INDEL_FACTOR_N = 0.51;    // facteur pour deduire le cout d'elision/insertion -
                                       // alignement multiple
const double TRANSPOSITION_FACTOR = 0.;  // facteur pour deduire le cout de transposition
const double INDEL_DISTANCE = 1.0;     // cout d'elision/insertion

enum {
  CHANGE_POINT ,
  SEGMENT
};

const double ROUNDOFF_ERROR = 1.e-10;  // erreur sur une somme de doubles
const int NB_SEGMENTATION = 10;        // nombre de segmentations calculees

enum {
  APPROXIMATED ,
  EXACT
};

const double FREQUENCY_RATIO = 0.1;    // rapport des frequences pour stopper le calcul
                                       // de la fonction de correlation
// const int CORRELATION_MIN_FREQUENCY = 20;  frequence minimum pour stopper le calcul
                                           // de la fonction de correlation

const int MAX_DIFFERENCING_ORDER = 3;  // ordre maximum de differenciation
const int POINTWISE_AVERAGE_NB_SEQUENCE = 250;  // nombre maximum de sequences ecrites
                                                // pour la sortie fichier
const int ABSORBING_RUN_LENGTH = 5;    // longueur par defaut de la serie finale absorbante
const int MAX_ABSORBING_RUN_LENGTH = 20;  // longueur maximum de la serie finale absorbante
const double ABSORBING_RUN_STANDARD_DEVIATION_FACTOR = 10.;  // facteur pour definir l'ecart-type
                                                             // des series finales absorbantes reelles



/****************************************************************
 *
 *  Definition des classes :
 */


class CategoricalSequenceProcess;
class VariableOrderMarkov;
class VariableOrderMarkovData;
class VariableOrderMarkovIterator;
class HiddenVariableOrderMarkov;
class SemiMarkov;
class SemiMarkovData;
class SemiMarkovIterator;
class HiddenSemiMarkov;
class NonhomogeneousMarkov;
class NonhomogeneousMarkovData;
class Sequences;
class MarkovianSequences;
class SequenceCharacteristics;
class Vectors;

class Switching_sequence;  // ajout par Florence Chaubert


class CategoricalSequenceProcess : public CategoricalProcess {  // processus d'observation categoriel
                                                                // pour des sequences
public :

    Distribution *length;   // loi des longueurs des sequences
    Curves *index_value;    // probabilites des valeurs en fonction de l'index
    double *no_occurrence;  // probabilite de ne pas observer une valeur
    Distribution **first_occurrence;  // lois du temps avant la premiere occurrence d'une valeur
    double *leave;          // probabilite de quitter une valeur
    Distribution **recurrence_time;  // lois du temps de retour dans une valeur
    double *absorption;     // probabilite d'etre absorbe par une valeur
    DiscreteParametric **sojourn_time;  // lois du temps de sejour dans un valeur
    Distribution **nb_run;  // lois du nombre de series par sequence des valeurs
    Distribution **nb_occurrence;  // lois du nombre d'occurrences par sequence des valeurs

    void create_characteristic(const Distribution &ilength , bool* homogeneity ,
                               bool counting_flag = true);
    void create_characteristic(const Distribution &ilength , bool sojourn_time_flag = true ,
                               bool counting_flag = true);
    void copy(const CategoricalSequenceProcess &process , bool characteristic_flag = true);
    void init_occupancy(const CategoricalSequenceProcess &process , int occupancy_nb_value);
    void remove();

    CategoricalSequenceProcess(int inb_state = 0 , int inb_value = 0 ,
                                 int observation_flag = false);
    CategoricalSequenceProcess(int inb_state , DiscreteParametric **occupancy);
    CategoricalSequenceProcess(const CategoricalProcess &process);
    CategoricalSequenceProcess(const CategoricalSequenceProcess &process ,
                                 char manip = 'c' , int param = true);
    ~CategoricalSequenceProcess();
    CategoricalSequenceProcess& operator=(const CategoricalSequenceProcess &process);

    Distribution* weight_computation() const;

    std::ostream& ascii_print(std::ostream &os , int process ,
                              FrequencyDistribution **empirical_observation ,
                              FrequencyDistribution *marginal_distribution ,
                              const SequenceCharacteristics *characteristics ,
                              bool exhaustive , bool file_flag , Forward **forward = NULL) const;
    std::ostream& spreadsheet_print(std::ostream &os , int process ,
                                    FrequencyDistribution **empirical_observation = NULL ,
                                    FrequencyDistribution *marginal_distribution = NULL ,
                                    const SequenceCharacteristics *characteristics = NULL ,
                                    Forward **forward = NULL) const;
    bool plot_print(const char *prefix , const char *title , int process ,
                    FrequencyDistribution **empirical_observation = NULL ,
                    FrequencyDistribution *marginal_distribution = NULL ,
                    const SequenceCharacteristics *characteristics = NULL ,
                    const FrequencyDistribution *hlength = NULL , Forward **forward = NULL) const;
    void plotable_write(MultiPlotSet &plot , int &index , int process ,
                        FrequencyDistribution **empirical_observation = NULL ,
                        FrequencyDistribution *marginal_distribution = NULL ,
                        const SequenceCharacteristics *characteristics = NULL ,
                        const FrequencyDistribution *hlength = NULL , Forward **forward = NULL) const;
};


CategoricalSequenceProcess* occupancy_parsing(StatError &error , ifstream &in_file ,
                                                int &line , const Chain &chain ,
                                                double cumul_threshold = CUMUL_THRESHOLD);
bool test_hidden(int nb_output_process , CategoricalSequenceProcess **process);



class Correlation : public StatInterface , public Curves {  // fonctions de correlation

    friend class VariableOrderMarkovChain;
    friend class VariableOrderMarkov;
    friend class Sequences;

    friend std::ostream& operator<<(std::ostream &os , const Correlation &correlation)
    { return correlation.ascii_write(os); }

private :

    int type;               // type de coefficient (PEARSON/SPEARMAN/KENDALL)
    int *variable_type;     // types de variables (OBSERVED/THEORETICAL STATE/OUTPUT)
    int *variable1;         // 1ere variables
    int *variable2;         // 2eme variables
    double *white_noise;    // fonction theorique pour un bruit blanc filtre

    void copy(const Correlation &correl);
    void remove();
    bool plot_print(const char *path , double *confidence_limit) const;

public :

    Correlation();
    Correlation(int itype , int max_lag , int ivariable1 , int ivariable2);
    Correlation(int inb_curve , int ilength , bool frequency_flag , int itype);
    Correlation(const Correlation &correl)
    :Curves(correl) { copy(correl); }
    ~Correlation();
    Correlation& operator=(const Correlation &correl);

    Correlation* merge(StatError &error , int nb_correl , const Correlation **icorrel) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    bool white_noise_correlation(StatError &error , int nb_point , double *filter ,
                                 int residual = true);
    bool white_noise_correlation(StatError &error , const Distribution &dist);
    bool white_noise_correlation(StatError &error , int order);

    // acces membres de la classe

    int get_type() const { return type; }
    int get_variable_type(int index) const { return variable_type[index]; }
    int get_variable1(int index) const { return variable1[index]; }
    int get_variable2(int index) const { return variable2[index]; }
    double get_white_noise(int lag) const { return white_noise[lag]; }
};



class DistanceMatrix;
class RenewalData;
class Tops;
class VectorDistance;
class Vectors;

class Sequences : public StatInterface {  // sequences

    friend class Vectors;

    friend Sequences* sequences_ascii_read(StatError &error , const char *path ,
                                           bool old_format);
    friend std::ostream& operator<<(std::ostream &os , const Sequences &seq)
    { return seq.ascii_write(os); }

protected :

    int nb_sequence;        // nombre de sequences
    int *identifier;        // identificateurs des sequences
    int max_length;         // longueur maximum des sequences
    int cumul_length;       // longueur cumulee des sequences
    int *length;            // longueurs des sequences
    FrequencyDistribution *length_distribution;  // loi empirique des longueurs des sequences
    int **vertex_identifier;  // identificateurs des vertex d'un MTG associe
    int index_parameter_type;  // type du parametre d'index (TIME/POSITION)
    FrequencyDistribution *index_parameter_distribution;  // loi empirique des parametres d'index explicites
    FrequencyDistribution *index_interval;  // intervalles entre parametres d'index explicites
    int **index_parameter;  // parametres d'index explicites
    int nb_variable;        // nombre de variables
    int *type;              // type de chaque variable (INT_VALUE/REAL_VALUE/STATE/NB_INTERNODE)
    double *min_value;      // valeurs minimums
    double *max_value;      // valeurs maximums
    FrequencyDistribution **marginal_distribution;  // lois marginales empiriques
    Histogram **marginal_histogram;  // histogrammes marginaux
    int ***int_sequence;    // sequences, variables entieres
    double ***real_sequence;  // sequences, variables reelles

    void init(int inb_sequence , int *iidentifier , int *ilength , int **ivertex_identifier ,
              int iindex_parameter_type , int inb_variable , int *itype ,
              bool vertex_identifier_copy , bool init_flag);
    void init(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,
              bool init_flag);
    void copy(const Sequences &seq);
    void reverse(const Sequences &seq);
    void add_state_variable(const Sequences &seq);
    void remove_index_parameter(const Sequences &seq);
    void explicit_index_parameter(const Sequences &seq);
    void remove();

    bool increasing_index_parameter_checking(StatError &error , bool strict ,
                                             const char *pattern_label) const;
    bool increasing_sequence_checking(StatError &error , int variable , bool strict ,
                                      const char *pattern_label , const char *variable_label) const;

    void cluster(const Sequences &seq , int variable , int step , int mode);
    void transcode(const Sequences &seq , int ivariable , int min_symbol , int max_symbol ,
                   int *symbol , bool add_flag = false);
    void cluster(const Sequences &seq , int variable , int nb_class , double *limit);
    void select_variable(const Sequences &seq , int *variable);

    bool pointwise_average_ascii_print(StatError &error , const char *path , int *size ,
                                       bool standard_deviation , int output) const;
    bool pointwise_average_spreadsheet_print(StatError &error , const char *path , int *size ,
                                             bool standard_deviation , int output) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;
    std::ostream& ascii_print(std::ostream &os , char format , bool comment_flag ,
                              double *posterior_probability = NULL , double *entropy = NULL ,
                              double *nb_state_sequence = NULL , int line_nb_character = LINE_NB_CHARACTER) const;
    bool plot_print(const char *path , int ilength) const;

    void max_length_computation();
    void cumul_length_computation();
    void build_length_frequency_distribution();

    void index_parameter_computation();
    void build_index_parameter_frequency_distribution();
    void index_interval_computation();

    void min_value_computation(int variable);
    void max_value_computation(int variable);
    void build_marginal_frequency_distribution(int variable);
    void build_marginal_histogram(int variable , double step = D_DEFAULT ,
                                  double imin_value = D_INF);

    std::ostream& alignment_ascii_print(std::ostream &os , int width , int ref_index , int test_index ,
                                        const Sequences &alignment , int alignment_index) const;
    std::ostream& alignment_spreadsheet_print(std::ostream &os , int ref_index , int test_index ,
                                              const Sequences &alignment , int alignment_index) const;

    double indel_distance_computation(const VectorDistance &vector_dist ,
                                      double **rank , double **max_symbol_distance) const;
    double indel_distance_computation(const VectorDistance &vector_dist ,
                                      int index , int position , double **rank ,
                                      double **max_symbol_distance) const;
    double substitution_distance_computation(const VectorDistance &vector_dist , int ref_index ,
                                             int test_index , int ref_position , int test_position ,
                                             double **rank , const Sequences *test_seq = NULL) const;
    double substitution_distance_computation(int ref_index , int test_index , int ref_position,
                                             int test_position , double substitution_distance) const;

    std::ostream& multiple_alignment_ascii_print(std::ostream &os) const;
    bool multiple_alignment_ascii_print(StatError &error , const char *path) const;

    Sequences* multiple_alignment(const Sequences &test_seq , const VectorDistance &vector_dist ,
                                  double **rank , double **max_symbol_distance , bool begin_free ,
                                  bool end_free , int indel_cost , double indel_factor) const;

    void correlation_computation(Correlation &correl , int variable1 , int variable2 ,
                                 int normalization , bool individual_mean = false) const;

    std::ostream& profile_ascii_print(std::ostream &os , int index , int nb_segment ,
                                      double **profiles , const char *label ,
                                      double **piecewise_function = NULL , long double **change_point = NULL ,
                                      long double **begin_conditonal_entropy = NULL ,
                                      long double **end_conditional_entropy = NULL ,
                                      long double **change_point_entropy = NULL) const;
    std::ostream& profile_spreadsheet_print(std::ostream &os , int index , int nb_segment ,
                                            double **profiles , const char *label ,
                                            double **piecewise_function = NULL , long double **change_point = NULL ,
                                            long double **begin_conditonal_entropy = NULL ,
                                            long double **end_conditional_entropy = NULL ,
                                            long double **change_point_entropy = NULL) const;
    std::ostream& profile_plot_print(std::ostream &os , int index , int nb_segment ,
                                     double **profiles , double **piecewise_function = NULL ,
                                     long double **change_point = NULL ,
                                     long double **begin_conditonal_entropy = NULL ,
                                     long double **end_conditional_entropy = NULL ,
                                     long double **change_point_entropy = NULL) const;
    void change_point_profile_plotable_write(MultiPlot &plot , int index , int nb_segment ,
                                             long double **change_point) const;
    void entropy_profile_plotable_write(MultiPlot &plot , int index ,
                                        long double *begin_conditional_entropy ,
                                        long double *end_conditional_entropy ,
                                        long double *change_point_entropy) const;

    void gamma_hyperparameter_computation(int index , int variable ,
                                          double *hyperparam) const;
    void gaussian_gamma_hyperparameter_computation(int index , int variable ,
                                                   double *hyperparam) const;
    int nb_parameter_computation(int index , int nb_segment , int *model_type) const;
    double one_segment_likelihood(int index , int *model_type , double **rank) const;
    Sequences* segmentation_output(int *nb_segment , int *model_type , std::ostream &os ,
                                   int output = SEQUENCE , int *ichange_point = NULL);
    double segmentation(int *nb_segment , int *model_type , double **rank ,
                        double *isegmentation_likelihood = NULL , int *nb_parameter = NULL ,
                        double *segment_penalty = NULL);
    double forward_backward(int index , int nb_segment , int *model_type , double **rank ,
                            double *likelihood , long double *segmentation_entropy ,
                            long double *first_order_entropy ,
                            long double *change_point_entropy , double *uniform_entropy ,
                            long double *marginal_entropy) const;
    double forward_backward(int index , int nb_segment , int *model_type , double **rank ,
                            std::ostream *os , MultiPlotSet *plot_set , int output ,
                            char format) const;
    double forward_backward_sampling(int index , int nb_segment , int *model_type ,
                                     double **rank , std::ostream &os , char format ,
                                     int nb_segmentation) const;
    double N_segmentation(int index , int nb_segment , int *model_type , double **irank ,
                          std::ostream &os , char format , int inb_segmentation ,
                          double likelihood) const;
    double forward_backward_dynamic_programming(int index , int nb_segment , int *model_type ,
                                                double **rank , std::ostream *os ,
                                                MultiPlotSet *plot_set , int output ,
                                                char format , double likelihood = D_INF) const;

    std::ostream& profile_ascii_print(std::ostream &os , int index , int nb_state ,
                                      double **profiles , double *begin_conditional_entropy ,
                                      double *marginal_entropy , double *begin_partial_entropy ,
                                      double *end_conditional_entropy = NULL , double *end_partial_entropy = NULL) const;
    std::ostream& profile_spreadsheet_print(std::ostream &os , int index , int nb_state ,
                                            double **profiles , double *begin_conditional_entropy ,
                                            double *marginal_entropy , double *begin_partial_entropy ,
                                            double *end_conditional_entropy = NULL , double *end_partial_entropy = NULL) const;
    std::ostream& profile_plot_print(std::ostream &os , int index , int nb_state ,
                                     double **profiles , double *begin_conditional_entropy ,
                                     double *marginal_entropy , double *begin_partial_entropy ,
                                     double *end_conditional_entropy = NULL , double *end_partial_entropy = NULL) const;
    void profile_plotable_write(MultiPlot &plot , int index , int nb_state ,
                                double **profiles) const;
    void entropy_profile_plotable_write(MultiPlot &plot , int index , double *begin_entropy ,
                                        double *end_entropy = NULL , double *marginal_entropy = NULL) const;

public :

    Sequences();
    Sequences(int inb_sequence , int inb_variable);
    Sequences(int inb_sequence , int *iidentifier , int *ilength ,
              int **ivertex_identifier , int iindex_parameter_type , int inb_variable ,
              int *itype , bool vertex_identifier_copy = true , bool init_flag = false)
    { init(inb_sequence , iidentifier , ilength , ivertex_identifier ,
           iindex_parameter_type , inb_variable , itype ,
           vertex_identifier_copy , init_flag); }
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int iindex_parameter_type ,  // interface AML
              int inb_variable , int itype , int ***iint_sequence);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,  // interface AML
              double ***ireal_sequence);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int **ivertex_identifier ,
              int iindex_parameter_type , int **iindex_parameter , int inb_variable ,
              int *itype , int ***iint_sequence , double ***ireal_sequence);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,
              bool init_flag = false)
    { init(inb_sequence , iidentifier , ilength , inb_variable , init_flag); }
    Sequences(const FrequencyDistribution &ilength_distribution , int inb_variable ,
              int *itype , bool init_flag = false);
    Sequences(const RenewalData &timev);
    Sequences(const Sequences &seq , int variable , int itype);
    Sequences(const Sequences &seq , int inb_sequence , int *index);
    Sequences(const Sequences &seq , bool *auxiliary);
    Sequences(const Sequences &seq , char transform = 'c' , int param = UNCHANGED);
    ~Sequences();
    Sequences& operator=(const Sequences &seq);

    DiscreteDistributionData* extract(StatError &error , int variable) const;

    Vectors* build_vectors(bool index_variable) const;
    Vectors* extract_vectors(StatError &error , int feature_type , int variable = I_DEFAULT ,
                             int value = I_DEFAULT) const;

    MarkovianSequences* markovian_sequences(StatError &error) const;
    Tops* tops(StatError &error) const;

    bool check(StatError &error , const char *pattern_label);

    TimeEvents* extract_time_events(StatError &error , int variable ,
                                    int begin_date , int end_date ,
                                    int previous_date = I_DEFAULT , int next_date = I_DEFAULT) const;
    RenewalData* extract_renewal_data(StatError &error , int variable ,
                                      int begin_index , int end_index) const;

    Sequences* merge(StatError &error , int nb_sample , const Sequences **iseq) const;

    Sequences* shift(StatError &error , int variable , int shift_param) const;
    Sequences* shift(StatError &error , int variable , double shift_param) const;
    Sequences* thresholding(StatError &error , int variable , int threshold , int mode) const;
    Sequences* thresholding(StatError &error , int variable , double threshold , int mode) const;
    Sequences* cluster(StatError &error , int variable , int step ,
                       int mode = FLOOR) const;
    Sequences* transcode(StatError &error , int variable , int *symbol) const;
    Sequences* cluster(StatError &error , int variable , int nb_class ,
                       int *ilimit) const;
    Sequences* cluster(StatError &error , int variable , int nb_class ,
                       double *ilimit) const;
    Sequences* scaling(StatError &error , int variable , int scaling_coeff) const;
    Sequences* scaling(StatError &error , int variable , double scaling_coeff) const;
    Sequences* round(StatError &error , int variable = I_DEFAULT ,
                     int mode = ROUND) const;

    Sequences* index_parameter_select(StatError &error , std::ostream &os ,
                                      int min_index_parameter ,
                                      int max_index_parameter , bool keep) const;
    Sequences* value_select(StatError &error , std::ostream &os , int variable ,
                            int imin_value , int imax_value , bool keep = true) const;
    Sequences* value_select(StatError &error , std::ostream &os , int variable ,
                            double imin_value , double imax_value , bool keep = true) const;
    Sequences* select_individual(StatError &error , int inb_sequence , int *iidentifier ,
                                 bool keep = true) const;

    Sequences* remove_index_parameter(StatError &error) const;
    Sequences* explicit_index_parameter(StatError &error) const;
    Sequences* select_variable(StatError &error , int inb_variable , int *ivariable ,
                               bool keep = true) const;
    Sequences* merge_variable(StatError &error , int nb_sample , const Sequences **iseq ,
                              int ref_sample = I_DEFAULT) const;
    Sequences* shift_variable(StatError &error , int variable , int lag) const;

    Sequences* reverse(StatError &error) const;
    Sequences* length_select(StatError &error , std::ostream &os , int min_length ,
                             int imax_length , bool keep = true) const;
    Sequences* remove_run(StatError &error , int variable , int ivalue ,
                          char position , int max_run_length = I_DEFAULT) const;
    Sequences* index_parameter_extract(StatError &error , int min_parameter_index ,
                                       int max_parameter_index = I_DEFAULT) const;
    Sequences* segmentation_extract(StatError &error , int variable , int nb_value ,
                                    int *ivalue , bool keep = true ,
                                    bool concatenation = false) const;

    Sequences* cumulate(StatError &error , int variable = I_DEFAULT) const;
    Sequences* difference(StatError &error , int variable = I_DEFAULT ,
                          bool first_element = false , bool normalization = false) const;
    Sequences* sequence_normalization(StatError &error , int variable = I_DEFAULT) const;
    Sequences* moving_average(StatError &error , int nb_point , double *filter ,
                              int variable = I_DEFAULT , bool begin_end = false ,
                              int output = TREND) const;
    Sequences* moving_average(StatError &error , const Distribution &dist ,
                              int variable = I_DEFAULT , bool begin_end = false ,
                              int output = TREND) const;

    Sequences* pointwise_average(StatError &error , bool circular = false ,
                                 bool standard_deviation = false , int output = SEQUENCE ,
                                 const char *path = NULL , char format = 'a') const;

    Sequences* recurrence_time_sequences(StatError &error , int variable , int value) const;
    Sequences* sojourn_time_sequences(StatError &error , int variable) const;

    Sequences* transform_position(StatError &error , int step) const;

    Sequences* cross(StatError &error) const;

    std::ostream& line_write(std::ostream &os) const;

    virtual std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                           bool exhaustive = false) const;
    virtual bool ascii_data_write(StatError &error , const char *path ,
                                  char format = 'c' , bool exhaustive = false) const;
    bool plot_data_write(StatError &error , const char *prefix ,
                         const char *title = NULL) const;
    MultiPlotSet* get_plotable_data(StatError &error) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    int min_index_parameter_computation() const;
    int max_index_parameter_computation(bool last_position = false) const;

    void marginal_frequency_distribution_computation(int variable);
    bool select_step(StatError &error , int variable , double step ,
                     double imin_value = D_INF);

    double mean_computation(int variable) const;
    double variance_computation(int variable , double mean) const;
    double mean_absolute_deviation_computation(int variable , double mean) const;
    double mean_absolute_difference_computation(int variable) const;
    double skewness_computation(int variable , double mean , double variance) const;
    double kurtosis_computation(int variable , double mean , double variance) const;
    double* mean_direction_computation(int variable , int unit) const;

    FrequencyDistribution* value_index_interval_computation(StatError &error ,
                                                            int variable , int value) const;

    Correlation* correlation_computation(StatError &error , int variable1 , int variable2 ,
                                         int itype = PEARSON , int max_lag = I_DEFAULT ,
                                         int normalization = EXACT , bool individual_mean = false) const;
    Correlation* partial_autocorrelation_computation(StatError &error , int variable ,
                                                     int itype = PEARSON , int max_lag = I_DEFAULT) const;

    DistanceMatrix* alignment(StatError &error , std::ostream *os , const VectorDistance &ivector_dist ,
                              int ref_identifier = I_DEFAULT , int test_identifier = I_DEFAULT ,
                              bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,
                              double indel_factor = INDEL_FACTOR_1 , bool transposition_flag = false ,
                              double transposition_factor = TRANSPOSITION_FACTOR ,
                              const char *result_path = NULL , char result_format = 'a' ,
                              const char *alignment_path = NULL , char alignment_format = 'a') const;
    DistanceMatrix* alignment(StatError &error , std::ostream *os , int ref_identifier = I_DEFAULT ,
                              int test_identifier = I_DEFAULT , bool begin_free = false , bool end_free = false ,
                              const char *result_path = NULL , char result_format = 'a' ,
                              const char *alignment_path = NULL , char alignment_format = 'a') const;

    Sequences* multiple_alignment(StatError &error , std::ostream &os , const VectorDistance &ivector_dist ,
                                  bool begin_free = false , bool end_free = false , int indel_cost = ADAPTATIVE ,
                                  double indel_factor = INDEL_FACTOR_N , int algorithm = AGGLOMERATIVE ,
                                  const char *path = NULL) const;

    Sequences* segmentation(StatError &error , std::ostream &os , int iidentifier ,
                            int nb_segment , int *ichange_point , int *model_type ,
                            int output = SEQUENCE) const;
    Sequences* segmentation(StatError &error , std::ostream &os , int *nb_segment ,
                            int *model_type , int iidentifier = I_DEFAULT ,
                            int output = SEQUENCE) const;
    Sequences* segmentation(StatError &error , std::ostream &os , int iidentifier ,
                            int max_nb_segment , int *model_type ,
                            int output = SEQUENCE) const;

    Sequences* hierarchical_segmentation(StatError &error , std::ostream &os , int iidentifier ,
                                         int max_nb_segment , int *model_type) const;

    bool segment_profile_write(StatError &error , std::ostream &os , int iidentifier ,
                               int nb_segment , int *model_type , int output = SEGMENT ,
                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;
    bool segment_profile_write(StatError &error , const char *path , int iidentifier ,
                               int nb_segment , int *model_type , int output = SEGMENT ,
                               char format = 'a' , int segmentation = FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;
    bool segment_profile_plot_write(StatError &error , const char *prefix ,
                                    int iidentifier , int nb_segment , int *model_type ,
                                    int output = SEGMENT , const char *title = NULL) const;
    MultiPlotSet* segment_profile_plotable_write(StatError &error , int iidentifier ,
                                                 int nb_segment , int *model_type ,
                                                 int output = SEGMENT) const;

    // acces membres de la classe

    int get_nb_sequence() const { return nb_sequence; }
    int get_identifier(int iseq) const { return identifier[iseq]; }
    int get_max_length() const { return max_length; }
    int get_cumul_length() const { return cumul_length; }
    int get_length(int index_seq) const { return length[index_seq]; }
    FrequencyDistribution* get_length_distribution() const { return length_distribution; }
    int get_vertex_identifier(int iseq , int index) const
    { return vertex_identifier[iseq][index]; }
    int get_index_parameter_type() const { return index_parameter_type; }
    FrequencyDistribution* get_index_parameter_distribution() const { return index_parameter_distribution; }
    FrequencyDistribution* get_index_interval() const { return index_interval; }
    int get_index_parameter(int iseq , int index) const
    { return index_parameter[iseq][index]; }
    int get_nb_variable() const { return nb_variable; }
    int get_type(int variable) const { return type[variable]; }
    double get_min_value(int variable) const { return min_value[variable]; }
    double get_max_value(int variable) const { return max_value[variable]; }
    FrequencyDistribution* get_marginal_distribution(int variable) const
    { return marginal_distribution[variable]; }
    Histogram* get_marginal_histogram(int variable) const
    { return marginal_histogram[variable]; }
    int get_int_sequence(int iseq , int variable , int index) const
    { return int_sequence[iseq][variable][index]; }
    double get_real_sequence(int iseq , int variable , int index) const
    { return real_sequence[iseq][variable][index]; }
};


Sequences* sequences_ascii_read(StatError &error , const char *path ,
                                bool old_format = false);



class SequenceCharacteristics {  // caracteristiques des sequences pour une variable

public :

    int nb_value;           // nombre de valeurs a partir de 0
    Curves *index_value;    // probabilites des valeurs en fonction de l'index
    Curves *explicit_index_value;    // probabilites des valeurs en fonction de l'index explicite
    FrequencyDistribution **first_occurrence;  // lois empiriques du temps avant la premiere
                                               // occurrence d'une valeur
    FrequencyDistribution **recurrence_time;  // lois empiriques du temps de retour dans une valeur
    FrequencyDistribution **sojourn_time;  // lois empiriques du temps de sejour dans une valeur
                                           // (series completes)
    FrequencyDistribution **initial_run;  // lois empiriques du temps de sejour dans une valeur
                                          // (series initiales censurees a gauche)
    FrequencyDistribution **final_run;  // lois empiriques du temps de sejour dans une valeur
                                        // (series finales censurees a droite)
    FrequencyDistribution **nb_run;  // lois empiriques du nombre de series 
                                     //par sequence des valeurs
    FrequencyDistribution **nb_occurrence;  // lois empiriques du nombre d'occurrences
                                            // par sequence des valeurs

    void copy(const SequenceCharacteristics &characteristics);
    void reverse(const SequenceCharacteristics &characteristics);
    void remove();

    void create_sojourn_time_frequency_distribution(int max_length , int initial_run_flag = false);

    std::ostream& ascii_print(std::ostream &os , int type ,
                              const FrequencyDistribution &length_distribution ,
                              bool exhaustive , bool comment_flag) const;
    std::ostream& spreadsheet_print(std::ostream &os , int type ,
                                    const FrequencyDistribution &length_distribution) const;
    bool plot_print(const char *prefix , const char *title , int variable ,
                    int nb_variable , int type , const FrequencyDistribution &length_distribution) const;
    void plotable_write(MultiPlotSet &plot , int &index , int variable ,
                        int type , const FrequencyDistribution &length_distribution) const;

    SequenceCharacteristics(int inb_value = I_DEFAULT);
    SequenceCharacteristics(const SequenceCharacteristics &characteristics ,
                            bool initial_run_flag);
    SequenceCharacteristics(const SequenceCharacteristics &characteristics ,
                            char transform = 'c');
    ~SequenceCharacteristics();
    SequenceCharacteristics& operator=(const SequenceCharacteristics &characteristics);
};



class Function;

class SelfTransition : public Curves {  // probabilites de rester dans un etat en fonction de l'index

public :

    SelfTransition(int ilength)
    :Curves(1 , ilength , true , false) {}

    Function* monomolecular_regression() const;
    Function* logistic_regression() const;
};



class VariableOrderMarkovChain;
class VariableOrderMarkovChainData;

class MarkovianSequences : public Sequences {  // trajectoires correspondant a
                                               // un processus markovien
    friend class VariableOrderMarkovChain;
    friend class VariableOrderMarkov;
    friend class HiddenVariableOrderMarkov;
    friend class SemiMarkov;
    friend class HiddenSemiMarkov;
    friend class NonhomogeneousMarkov;

    friend std::ostream& operator<<(std::ostream &os , const MarkovianSequences &seq)
    { return seq.ascii_write(os); }

protected :

    double *min_interval;   // intervalles minimums entre 2 valeurs
    SelfTransition **self_transition;  // probabilites de rester dans un etat
                                       // en fonction de l'index
    FrequencyDistribution ***observation_distribution;  // lois empiriques d'observation
    Histogram ***observation_histogram;  // histogrammes d'observation
    SequenceCharacteristics **characteristics;  // caracteristiques pour une variable donnee

    void init();
    void copy(const MarkovianSequences &seq , int param = UNCHANGED);
    void reverse(const MarkovianSequences &seq);
    void add_state_variable(const MarkovianSequences &seq , int param);
    void remove();

    MarkovianSequences* transcode(StatError &error ,
                                  const CategoricalSequenceProcess *process) const;
    MarkovianSequences* build_auxiliary_variable(DiscreteParametricProcess **discrete_process ,
                                                 ContinuousParametricProcess **continuous_process) const;

    MarkovianSequences* remove_variable_1() const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;
    bool plot_print(const char *prefix , const char *title , int variable ,
                    int nb_variable) const;
    void plotable_write(MultiPlotSet &plot , int &index , int variable) const;

    void state_variable_init(int itype = STATE);

    void min_interval_computation(int variable);

    void transition_count_computation(const VariableOrderMarkovChainData &chain_data ,
                                      const VariableOrderMarkovChain &markov ,
                                      bool begin = true , bool non_terminal = false) const;
    void transition_count_computation(const ChainData &chain_data ,
                                      const SemiMarkov *smarkov = NULL) const;

    void self_transition_computation(int state);
    Distribution* weight_computation() const;
    void observation_frequency_distribution_computation(int variable , int nb_state);
    bool test_hidden(int variable) const;

    void build_index_value(int variable);
    void build_explicit_index_value(int variable);
    void build_first_occurrence_frequency_distribution(int variable);
    void build_recurrence_time_frequency_distribution(int variable);
    void build_sojourn_time_frequency_distribution(int variable , int initial_run_flag = false);
    void build_nb_run_frequency_distribution(int variable);
    void build_nb_occurrence_frequency_distribution(int variable);

    void censored_sojourn_time_frequency_distribution_computation(FrequencyDistribution **initial_run ,
                                                                  FrequencyDistribution **final_run ,
                                                                  FrequencyDistribution **single_run) const;

    std::ostream& linear_model_spreadsheet_print(std::ostream &os , int variable ,
                                                 ContinuousParametricProcess *process) const;
    bool linear_model_plot_print(const char *prefix , const char *title , int variable ,
                                 ContinuousParametricProcess *process) const;
    void linear_model_plotable_write(MultiPlotSet &plot , int &index , int variable ,
                                     ContinuousParametricProcess *process) const;

    template <typename Type>
    void gamma_estimation(Type ***state_sequence_count , int variable ,
                          ContinuousParametricProcess *process , int iter) const;
    template <typename Type>
    void zero_inflated_gamma_estimation(Type ***state_sequence_count , int variable ,
                                        ContinuousParametricProcess *process , int iter) const;
    template <typename Type>
    void gaussian_estimation(Type ***state_sequence_count , int variable ,
                             ContinuousParametricProcess *process) const;
    template <typename Type>
    void von_mises_estimation(Type ***state_sequence_count , int variable ,
                              ContinuousParametricProcess *process) const;
    template <typename Type>
    void  linear_model_estimation(Type ***state_sequence_count , int variable ,
                                  ContinuousParametricProcess *process) const;

    std::ostream& likelihood_write(std::ostream &os , int nb_model , double **likelihood ,
                                   const char *label , bool exhaustive = false ,
                                   char algorithm = 'd') const;
    bool likelihood_write(StatError &error , const char *path , int nb_model ,
                          double **likelihood , const char *label , char algorithm = 'd') const;

public :

    MarkovianSequences();
    MarkovianSequences(int inb_sequence , int *iidentifier , int *ilength ,
                       int **ivertex_identifier , int iindex_parameter_type , int inb_variable ,
                       int *itype , bool vertex_identifier_copy = true , bool init_flag = false)
    :Sequences(inb_sequence , iidentifier , ilength , ivertex_identifier ,
               iindex_parameter_type , inb_variable , itype ,
               vertex_identifier_copy , init_flag) { init(); }
    MarkovianSequences(const FrequencyDistribution &ilength_distribution , int inb_variable ,
                       int *itype , bool init_flag = false)
    :Sequences(ilength_distribution , inb_variable , itype , init_flag) { init(); }
    MarkovianSequences(const MarkovianSequences &seq , int variable , int itype)
    :Sequences(seq , variable , itype) { init(); }
    MarkovianSequences(const Sequences &seq);
    MarkovianSequences(const MarkovianSequences &seq , bool *auxiliary);
    MarkovianSequences(const MarkovianSequences &seq , char transform = 'c' ,
                       int param = UNCHANGED);
    ~MarkovianSequences();
    MarkovianSequences& operator=(const MarkovianSequences &seq);

    DiscreteDistributionData* extract(StatError &error , int type ,
                                      int variable , int value) const;

    MarkovianSequences* merge(StatError &error , int nb_sample ,
                              const MarkovianSequences **iseq) const;

    MarkovianSequences* cluster(StatError &error , int variable , int step ,
                                int mode = FLOOR) const;
    MarkovianSequences* transcode(StatError &error , int ivariable , int *symbol ,
                                  bool add_flag = false) const;
    MarkovianSequences* consecutive_values(StatError &error , std::ostream &os ,
                                           int ivariable , bool add_flag = false) const;
    MarkovianSequences* cluster(StatError &error , int ivariable , int nb_class ,
                                int *ilimit , bool add_flag = false) const;
    MarkovianSequences* cluster(StatError &error , int variable , int nb_class ,
                                double *ilimit) const;

    MarkovianSequences* remove_index_parameter(StatError &error) const;
    MarkovianSequences* explicit_index_parameter(StatError &error) const;
    MarkovianSequences* select_variable(StatError &error , int inb_variable ,
                                        int *ivariable , bool keep = true) const;
    MarkovianSequences* merge_variable(StatError &error , int nb_sample ,
                                       const MarkovianSequences **iseq , int ref_sample = I_DEFAULT) const;

    MarkovianSequences* initial_run_computation(StatError &error) const;
    MarkovianSequences* add_absorbing_run(StatError &error ,
                                          int sequence_length = I_DEFAULT ,
                                          int run_length = I_DEFAULT) const;

    MarkovianSequences* split(StatError &error , int step) const;

    std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                   bool exhaustive = false) const;
    bool ascii_data_write(StatError &error , const char *path ,
                          char format = 'c' , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    bool transition_count(StatError &error , std::ostream &os , int max_order ,
                          bool begin = false , int estimator = MAXIMUM_LIKELIHOOD ,
                          const char *path = NULL) const;
    bool word_count(StatError &error , std::ostream &os , int variable , int word_length ,
                    int begin_state = I_DEFAULT , int end_state = I_DEFAULT ,
                    int min_frequency = 1) const;
    bool mtg_write(StatError &error , const char *path , int *itype) const;

    int cumulative_distribution_function_computation(int variable , double **cdf) const;
    int cumulative_distribution_function_computation(int variable , int state ,
                                                     double **cdf) const;

    double iid_information_computation() const;

    void self_transition_computation();
    void self_transition_computation(bool *homogeneity);
    void sojourn_time_frequency_distribution_computation(int variable);

    void build_observation_frequency_distribution(int nb_state);
    void build_observation_histogram(int variable , int nb_state , double step = D_DEFAULT);
    void build_observation_histogram(int nb_state);
    bool select_step(StatError &error , int variable , double step ,
                     double imin_value = D_INF);

    void build_characteristic(int variable = I_DEFAULT , bool sojourn_time_flag = true ,
                              bool initial_run_flag = false);

    NonhomogeneousMarkov* nonhomogeneous_markov_estimation(StatError &error , int *ident ,
                                                           bool counting_flag = true) const;

    VariableOrderMarkov* variable_order_markov_estimation(StatError &error , std::ostream &os ,
                                                          char model_type , int min_order = 0 ,
                                                          int max_order = ORDER ,
                                                          int algorithm = LOCAL_BIC ,
                                                          double threshold = LOCAL_BIC_THRESHOLD ,
                                                          int estimator = LAPLACE ,
                                                          bool global_initial_transition = true ,
                                                          bool global_sample = true ,
                                                          bool counting_flag = true) const;
    VariableOrderMarkov* variable_order_markov_estimation(StatError &error ,
                                                          const VariableOrderMarkov &imarkov ,
                                                          bool global_initial_transition = true ,
                                                          bool counting_flag = true) const;
    VariableOrderMarkov* variable_order_markov_estimation(StatError &error ,
                                                          char model_type , int order = 1 ,
                                                          bool global_initial_transition = true ,
                                                          bool counting_flag = true) const;

    VariableOrderMarkov* lumpability_estimation(StatError &error , std::ostream &os , int *symbol ,
                                                int penalty_type = BIC , int order = 1 ,
                                                bool counting_flag = true) const;

    SemiMarkov* semi_markov_estimation(StatError &error , std::ostream &os , char model_type ,
                                       int estimator = COMPLETE_LIKELIHOOD ,
                                       bool counting_flag = true , int nb_iter = I_DEFAULT ,
                                       int mean_computation_method = COMPUTED) const;

    HiddenVariableOrderMarkov* hidden_variable_order_markov_estimation(StatError &error , std::ostream &os ,
                                                                       const HiddenVariableOrderMarkov &ihmarkov ,
                                                                       bool global_initial_transition = true ,
                                                                       bool common_dispersion = false ,
                                                                       bool counting_flag = true ,
                                                                       bool state_sequence = true ,
                                                                       int nb_iter = I_DEFAULT) const;
    HiddenVariableOrderMarkov* hidden_variable_order_markov_stochastic_estimation(StatError &error , std::ostream &os ,
                                                                                  const HiddenVariableOrderMarkov &ihmarkov ,
                                                                                  bool global_initial_transition = true ,
                                                                                  bool common_dispersion = false ,
                                                                                  int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                                                  int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                                                  double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                                                  bool counting_flag = true ,
                                                                                  bool state_sequence = true ,
                                                                                  int nb_iter = I_DEFAULT) const;

    HiddenSemiMarkov* hidden_semi_markov_estimation(StatError &error , std::ostream &os ,
                                                    const HiddenSemiMarkov &ihsmarkov ,
                                                    bool common_dispersion = false ,
                                                    int estimator = COMPLETE_LIKELIHOOD ,
                                                    bool counting_flag = true ,
                                                    bool state_sequence = true ,
                                                    int nb_iter = I_DEFAULT ,
                                                    int mean_computation_method = COMPUTED) const;
    HiddenSemiMarkov* hidden_semi_markov_estimation(StatError &error , std::ostream &os ,
                                                    char model_type , int nb_state , bool left_right ,
                                                    double occupancy_mean = D_DEFAULT ,
                                                    bool common_dispersion = false ,
                                                    int estimator = COMPLETE_LIKELIHOOD ,
                                                    bool counting_flag = true ,
                                                    bool state_sequence = true ,
                                                    int nb_iter = I_DEFAULT ,
                                                    int mean_computation_method = COMPUTED) const;
    HiddenSemiMarkov* hidden_semi_markov_stochastic_estimation(StatError &error , std::ostream &os ,
                                                               const HiddenSemiMarkov &ihsmarkov ,
                                                               bool common_dispersion = false ,
                                                               int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                               int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                               double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                               int estimator = COMPLETE_LIKELIHOOD ,
                                                               bool counting_flag = true ,
                                                               bool state_sequence = true ,
                                                               int nb_iter = I_DEFAULT) const;
    HiddenSemiMarkov* hidden_semi_markov_stochastic_estimation(StatError &error , std::ostream &os ,
                                                               char model_type , int nb_state , bool left_right ,
                                                               double occupancy_mean = D_DEFAULT ,
                                                               bool common_dispersion = false ,
                                                               int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                               int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                               double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                               int estimator = COMPLETE_LIKELIHOOD ,
                                                               bool counting_flag = true ,
                                                               bool state_sequence = true ,
                                                               int nb_iter = I_DEFAULT) const;

    bool lumpability_test(StatError &error , std::ostream &os , int *symbol , int order = 1) const;

    bool comparison(StatError &error , std::ostream &os , int nb_model ,
                    const VariableOrderMarkov **imarkov , const char *path = NULL) const;

    bool comparison(StatError &error , std::ostream &os , int nb_model ,
                    const SemiMarkov **ismarkov , const char *path = NULL) const;

    bool comparison(StatError &error , std::ostream &os , int nb_model ,
                    const HiddenVariableOrderMarkov **ihmarkov ,
                    int algorithm = FORWARD , const char *path = NULL) const;

    bool comparison(StatError &error , std::ostream &os , int nb_model ,
                    const HiddenSemiMarkov **ihsmarkov ,
                    int algorithm = FORWARD , const char *path = NULL) const;

    // acces membres de la classe

    Curves* get_self_transition(int state) const { return self_transition[state]; }
    FrequencyDistribution*** get_observation_distribution() const
    { return observation_distribution; }
    FrequencyDistribution** get_observation_distribution(int variable) const
    { return observation_distribution[variable]; }
    FrequencyDistribution* get_observation_distribution(int variable , int state) const
    { return observation_distribution[variable][state]; }
    Histogram*** get_observation_histogram() const { return observation_histogram; }
    Histogram** get_observation_histogram(int variable) const
    { return observation_histogram[variable]; }
    Histogram* get_observation_histogram(int variable , int state) const
    { return observation_histogram[variable][state]; }
    SequenceCharacteristics* get_characteristics(int variable) const
    { return characteristics[variable]; }
};



#include "continuous_parametric_sequence_estimation.hpp"



#endif
