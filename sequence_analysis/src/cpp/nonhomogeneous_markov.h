/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: nonhomogeneous_markov.h 3257 2007-06-06 12:56:12Z dufourko $
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



#ifndef NONHOMOGENEOUS_MARKOV_H
#define NONHOMOGENEOUS_MARKOV_H



namespace sequence_analysis {



/****************************************************************
 *
 *  Constantes :
 */


  const double START_RATIO = 0.03;       // proportion de l'echantillon pour l'initialisation
                                         // des parametres (debut)
  const double END_RATIO = 0.1;          // proportion de l'echantillon pour l'initialisation
                                         // des parametres (fin)
  const int REGRESSION_NB_ELEMENT = 100;  // taille minimum de l'echantillon pour
                                          // la regression non-lineare
  const double GRADIENT_DESCENT_COEFF = 1.;  // coefficient algorithme de gradient
  const double RESIDUAL_SQUARE_SUM_DIFF = 1.e-6;  // seuil pour stopper les iterations
                                                  // de l'algorithme de gradient
  const int REGRESSION_NB_ITER = 1000;   // nombre d'iterations pour la regression non-lineaire



/****************************************************************
 *
 *  Definition des classes :
 */


  class Function : public stat_tool::RegressionKernel {  // fonction d'evolution des probabilites
                                                         // de rester dans un etat
  public :

    double *residual;       // residus
    int *frequency;         // effectifs correspondant a chaque index

    void copy(const Function&);
    void remove();

    Function();
    Function(int iident , int length , double *iparameter);
    Function(int iident , int length);
    Function(const Function &function);
    ~Function();
    Function& operator=(const Function &function);

    std::ostream& ascii_print(std::ostream &os , bool exhaustive , bool file_flag ,
                              const stat_tool::Curves *curves = NULL) const;
    std::ostream& spreadsheet_print(std::ostream &os , const stat_tool::Curves *curves = NULL) const;
    bool plot_print(const char *path , double residual_standard_deviation = D_DEFAULT) const;

    double regression_square_sum_computation(double self_transition_mean) const;
    void residual_computation(const SelfTransition &self_transition);
    double residual_mean_computation() const;
    double residual_variance_computation(double residual_mean) const;
    double residual_square_sum_computation() const;
  };


  Function* function_parsing(stat_tool::StatError &error , std::ifstream &in_file , int &line ,
                             int length , double min = 0. , double max = 1.);



  class NonhomogeneousMarkovData;


  class NonhomogeneousMarkov : public stat_tool::StatInterface , protected stat_tool::Chain {  // chaine de Markov non-homogene

    friend class MarkovianSequences;
    friend class NonhomogeneousMarkovData;

    friend NonhomogeneousMarkov* nonhomogeneous_markov_ascii_read(stat_tool::StatError &error ,
                                                                  const char *path , int length);
    friend std::ostream& operator<<(std::ostream &os , const NonhomogeneousMarkov &markov)
    { return markov.ascii_write(os , markov.markov_data); }

  protected :

    NonhomogeneousMarkovData *markov_data;  // pointeur sur un objet NonhomogeneousMarkovData
    bool *homogeneity;      // homogeneite des etats
    Function **self_transition;  // fonction d'evolution des probabilites
                                 // de rester dans un etat
    CategoricalSequenceProcess *process;

    void copy(const NonhomogeneousMarkov &markov , bool data_flag = true ,
              bool characteristic_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const NonhomogeneousMarkovData *seq ,
                              bool exhaustive = false , bool file_flag  = false) const;
    std::ostream& spreadsheet_write(std::ostream &os , const NonhomogeneousMarkovData *seq) const;
    bool plot_write(const char *prefix , const char *title ,
                    const NonhomogeneousMarkovData *seq) const;
    stat_tool::MultiPlotSet* get_plotable(const NonhomogeneousMarkovData *seq) const;

    int nb_parameter_computation() const;

    void transition_update(int state , int index , stat_tool::Chain &index_chain) const;
    void index_state_distribution();
    void state_no_occurrence_probability(int state , double increment = LEAVE_INCREMENT);
    void state_first_occurrence_distribution(int state , int min_nb_value = 1 ,
                                             double cumul_threshold = CUMUL_THRESHOLD);
    void state_nb_pattern_mixture(int state , char pattern);

  public :

    NonhomogeneousMarkov();
    NonhomogeneousMarkov(int inb_state , int *ident);
    NonhomogeneousMarkov(const stat_tool::Chain *pchain , const Function **pself_transition , int length);
    NonhomogeneousMarkov(const NonhomogeneousMarkov &markov , bool data_flag = true ,
                         bool characteristic_flag = true)
    :stat_tool::Chain(markov) { copy(markov , data_flag , characteristic_flag); }
    ~NonhomogeneousMarkov();
    NonhomogeneousMarkov& operator=(const NonhomogeneousMarkov &markov);

    DiscreteParametricModel* extract(stat_tool::StatError &error , int type , int state) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const char *path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    void characteristic_computation(int length , bool counting_flag);
    void characteristic_computation(const NonhomogeneousMarkovData &seq , bool counting_flag ,
                                    bool length_flag = true);

    double likelihood_computation(const MarkovianSequences &seq ,
                                  int index = I_DEFAULT) const;

    NonhomogeneousMarkovData* simulation(stat_tool::StatError &error ,
                                         const stat_tool::FrequencyDistribution &hlength ,
                                         bool counting_flag = true) const;
    NonhomogeneousMarkovData* simulation(stat_tool::StatError &error , int nb_sequence ,
                                         int length , bool counting_flag = true) const;
    NonhomogeneousMarkovData* simulation(stat_tool::StatError &error , int nb_sequence ,
                                         const MarkovianSequences &iseq ,
                                         bool counting_flag = true) const;

    // acces membres de la classe

    NonhomogeneousMarkovData* get_markov_data() const { return markov_data; }
    bool get_homogeneity(int state) const { return homogeneity[state]; }
    Function* get_self_transition(int state) const { return self_transition[state]; }
    CategoricalSequenceProcess* get_process() const { return process; }
  };


  NonhomogeneousMarkov* nonhomogeneous_markov_ascii_read(stat_tool::StatError &error , const char *path ,
                                                         int length = DEFAULT_LENGTH);



  class NonhomogeneousMarkovData : public MarkovianSequences {  // structure de donnees correspondant
                                                                // a une chaine de Markov non-homogene

    friend class MarkovianSequences;
    friend class NonhomogeneousMarkov;

    friend std::ostream& operator<<(std::ostream &os , const NonhomogeneousMarkovData &seq)
    { return seq.ascii_write(os , false); }

  private :

    NonhomogeneousMarkov *markov;  // pointeur sur un objet NonhomogeneousMarkov
    stat_tool::ChainData *chain_data;  // etats initaux et transitions
    double likelihood;      // vraisemblance des sequences

    void copy(const NonhomogeneousMarkovData &seq , bool model_flag = true);

  public :

    NonhomogeneousMarkovData();
    NonhomogeneousMarkovData(const stat_tool::FrequencyDistribution &ihlength);
    NonhomogeneousMarkovData(const MarkovianSequences &seq);
    NonhomogeneousMarkovData(const NonhomogeneousMarkovData &seq , bool model_flag = true ,
                             char transform = 'c')
    :MarkovianSequences(seq , transform) { copy(seq , model_flag); }
    ~NonhomogeneousMarkovData();
    NonhomogeneousMarkovData& operator=(const NonhomogeneousMarkovData &seq);

    stat_tool::DiscreteDistributionData* extract(stat_tool::StatError &error , int type , int state) const;
    NonhomogeneousMarkovData* remove_index_parameter(stat_tool::StatError &error) const;
    NonhomogeneousMarkovData* explicit_index_parameter(stat_tool::StatError &error) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const char *path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    void build_transition_count();

    // acces membres de la classe

    NonhomogeneousMarkov* get_markov() const { return markov; }
    stat_tool::ChainData* get_chain_data() const { return chain_data; }
    double get_likelihood() const { return likelihood; }
  };


};  // namespace sequence_analysis



#endif
