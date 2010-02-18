/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): J.-B. Durand and Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: mv_mixture.h 5336 2008-07-24 15:40:10Z guedon $
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

#ifndef MVMIXTURE_H
#define MVMIXTURE_H

#include "markovian.h"

/****************************************************************
 *
 *  Constantes :
 */

const double MVMIXTURE_LIKELIHOOD_DIFF = 1.e-6;

/****************************************************************
 *
 *  Definition des classes :
 */

#include "vectors.h"

class FrequencyDistribution;
class Mv_Mixture_data;
class DiscreteParametricProcess;
class NonparametricProcess;

// melange de lois multivariees
// a variables independantes

class Mv_Mixture : public StatInterface {

    friend class FrequencyDistribution;
    friend class Vectors;
    friend class Mv_Mixture_data;

    friend Mv_Mixture* mv_mixture_building(StatError &error , int nb_component , double *weight ,
                                           const DiscreteParametricProcess **component);
    friend Mv_Mixture* mv_mixture_ascii_read(StatError &error , const char *path ,
                                             double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Mv_Mixture &mixt)
    { return mixt.ascii_write(os , mixt.mixture_data , false , false); }

private :

    Mv_Mixture_data *mixture_data;  // pointeur sur un objet Mv_Mixture_data
    int nb_component;       // nombre de composantes
    int nb_var;       // dimension
    DiscreteParametric *weight;     // poids de chaque composante
    DiscreteParametricProcess **pcomponent; // composantes parametriques
    NonparametricProcess **npcomponent; // composantes non parametriques

    void copy(const Mv_Mixture &mixt , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const Mv_Mixture_data *mixt_data ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Mv_Mixture_data *mixt_data) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Mv_Mixture_data *mixt_data) const;
    plotable::MultiPlotSet* get_plotable(const Mv_Mixture_data *mixt_data) const;

    int nb_parameter_computation(double min_probability) const;
    double penalty_computation() const;

    /** Conditional density of observations */
    void get_output_conditional_distribution(const Vectors &mixt_data,
                         double** &output_cond,
                         bool log_computation=false) const;

    /** Conditional distribution of states */
    void get_posterior_distribution(const Vectors &mixt_data,
                    double** output_cond,
                    double** &posterior_dist) const;

    /** MAP algorithm */
    std::vector<int>* state_computation(StatError &error, const Vectors &vec,
                                        int algorithm=VITERBI, int index=I_DEFAULT,
                                        double** posterior_dist=NULL) const;

    /** Initialization of EM algorithm */
    void init();

public :

    Mv_Mixture();
    Mv_Mixture(int inb_component , double *pweight , int inb_variable,
           DiscreteParametricProcess **ppcomponent, NonparametricProcess **pnpcomponent);
    Mv_Mixture(int inb_component , int inb_variable, const DiscreteParametricProcess **ppcomponent,
           const NonparametricProcess **pnpcomponent);
    Mv_Mixture(const Mv_Mixture &mixt , bool *variable_flag , int inb_variable);
    Mv_Mixture(int inb_component, int inb_variable, int *nb_value, bool *force_param=NULL);
    Mv_Mixture(const Mv_Mixture &mixt , bool data_flag = true)
    { copy(mixt , data_flag); }
    ~Mv_Mixture();
    Mv_Mixture& operator=(const Mv_Mixture &mixt);

    /** extract parametric component */
    DiscreteParametricModel* extract_parametric_model(StatError &error , int ivariable,
                           int index) const;
    /** extract nonparametric component */
    Distribution* extract_nonparametric_model(StatError &error , int ivariable,
                          int index) const;
    /** extract marginal mixture distribution */
    Distribution* extract_distribution(StatError &error , int ivariable) const;
    Mv_Mixture_data* extract_data(StatError &error) const;

   /** Permutation of the states of \e self */
   void state_permutation(StatError& error, int* perm) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = 0) const;
    plotable::MultiPlotSet* get_plotable() const;

    double likelihood_computation(const Vectors &mixt_data,
                  bool log_computation=false) const;

    Mv_Mixture_data* simulation(StatError &error , int nb_element) const;

    /** add restored states to Vectors */
    Mv_Mixture_data* cluster(StatError &error,  const Vectors &vec,
                             int algorithm=VITERBI) const;

    /** return "true" if process ivariable is parametric */
    bool is_parametric(int ivariable) const;

    // acces membres de la classe

    Mv_Mixture_data* get_mixture_data() const { return mixture_data; }
    int get_nb_component() const { return nb_component; }
    int get_nb_variable() const { return nb_var; }
    DiscreteParametric* get_weight() const { return weight; }
    DiscreteParametricProcess* get_parametric_process(int variable) const;
    NonparametricProcess* get_nonparametric_process(int variable) const;
    DiscreteParametric* get_parametric_component(int variable, int index) const;
    Distribution* get_nonparametric_component(int variable, int index) const;
};


Mv_Mixture* mv_mixture_building(StatError &error , int nb_component ,
                int nb_variable, double *weight,
                DiscreteParametricProcess **ppcomponent,
                NonparametricProcess **pnpcomponent);
Mv_Mixture* mv_mixture_ascii_read(StatError &error , const char *path ,
                  double cumul_threshold = CUMUL_THRESHOLD);



// structure de donnees correspondant
 // a un melange
class Mv_Mixture_data : public Vectors {

    friend class FrequencyDistribution;
    friend class Mv_Mixture;
    friend class Vectors;

    friend std::ostream& operator<<(std::ostream &os , const Mv_Mixture_data &mixt_data)
    { return mixt_data.ascii_write(os , false); }

private :

    Mv_Mixture *mixture;       // pointeur sur un objet Mv_Mixture
    int nb_component;          // nombre de composantes
    FrequencyDistribution *weight;      // loi empirique des poids
    /// component[variable][state]
    FrequencyDistribution ***component;  // composantes empiriques pour chaque variable

    void copy(const Mv_Mixture_data &mixt_data , bool model_flag = true);
    void remove();

    /** Permutation of the states of \e self.*/
    void state_permutation(int *perm);

public :

    Mv_Mixture_data();
    Mv_Mixture_data(const Vectors &vec, int inb_component);
    Mv_Mixture_data(const Mv_Mixture &mixt);
    Mv_Mixture_data(const Mv_Mixture_data &mixt_data , bool model_flag = true)
      :Vectors(mixt_data) { copy(mixt_data , model_flag); }
    virtual ~Mv_Mixture_data();
    Mv_Mixture_data& operator=(const Mv_Mixture_data &mixt_data);

    DiscreteDistributionData* extract(StatError &error , int variable, int index) const;
    DiscreteDistributionData* extract_marginal(StatError &error , int variable) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = 0) const;
    plotable::MultiPlotSet* get_plotable() const;

    double information_computation() const;

    // acces membres de la classe

    Mv_Mixture* get_mixture() const { return mixture; }
    int get_nb_component() const { return nb_component; }
    FrequencyDistribution* get_weight() const { return weight; }
    FrequencyDistribution* get_component(int variable, int index) const { return component[variable][index]; }
};


#endif
