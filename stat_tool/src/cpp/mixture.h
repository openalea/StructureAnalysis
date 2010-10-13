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



#ifndef MIXTURE_H
#define MIXTURE_H



/****************************************************************
 *
 *  Constantes :
 */


const int MIXTURE_NB_COMPONENT = 100;   // nombre maximum de composantes

const double MIXTURE_PARAMETER = 20.;  // parametre initial pour une loi binomiale negative
const double MIN_WEIGHT_STEP = 0.1;    // pas minimum d'initialisation des poids
const double MAX_WEIGHT_STEP = 0.5;    // pas maximum d'initialisation des poids
const int MIXTURE_COEFF = 2;           // coefficient arrondi estimateur
const double MIXTURE_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
const int MIXTURE_NB_ITER = 500;        // nombre maximum d'iterations EM



/****************************************************************
 *
 *  Definition des classes :
 */


class FrequencyDistribution;
class MixtureData;


// class Mixture : public StatInterface , protected Distribution {
class Mixture : public StatInterface , public Distribution {  // melange de lois discretes

    friend class FrequencyDistribution;
    friend class MixtureData;

    friend Mixture* mixture_building(StatError &error , int nb_component , double *weight ,
                                     const DiscreteParametric **component);
    friend Mixture* mixture_ascii_read(StatError &error , const char *path ,
                                       double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Mixture &mixt)
    { return mixt.ascii_write(os , mixt.mixture_data , false , false); }

private :

    MixtureData *mixture_data;  // pointeur sur un objet MixtureData
    int nb_component;       // nombre de composantes
    DiscreteParametric *weight;  // poids de chaque composante
    DiscreteParametric **component; // composantes

    void copy(const Mixture &mixt , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const MixtureData *mixt_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const MixtureData *mixt_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const MixtureData *mixt_histo) const;
    MultiPlotSet* get_plotable(const MixtureData *mixt_histo) const;

    int nb_parameter_computation() const;
    double penalty_computation() const;

    void init(const FrequencyDistribution &histo , bool *estimate , int min_inf_bound ,
              bool component_flag);
    void expectation_step(MixtureData *mixt_histo , int nb_element) const;
    void variance_correction(MixtureData *mixt_histo , bool *estimate ,
                             int min_inf_bound) const;
    bool component_order_test() const;

public :

    Mixture();
    Mixture(int inb_component , double *pweight , const DiscreteParametric **pcomponent);
    Mixture(const Mixture &mixt , bool *component_flag , int inb_value);
    Mixture(int inb_component , const DiscreteParametric **pcomponent);
    Mixture(const Mixture &mixt , bool data_flag = true)
    :Distribution(mixt) { copy(mixt , data_flag); }
    ~Mixture();
    Mixture& operator=(const Mixture &mixt);

    DiscreteParametricModel* extract(StatError &error , int index) const;
    MixtureData* extract_data(StatError &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void computation(int min_nb_value = 1 , double cumul_threshold = CUMUL_THRESHOLD ,
                     bool component_flag = true);
    double likelihood_computation(const MixtureData &mixt_histo) const;
    MixtureData* simulation(StatError &error , int nb_element) const;

    // acces membres de la classe

    MixtureData* get_mixture_data() const { return mixture_data; }
    int get_nb_component() const { return nb_component; }
    DiscreteParametric* get_weight() const { return weight; }
    DiscreteParametric* get_component(int index) const { return component[index]; }
};


Mixture* mixture_building(StatError &error , int nb_component , double *weight ,
                          const DiscreteParametric **component);
Mixture* mixture_ascii_read(StatError &error , const char *path ,
                            double cumul_threshold = CUMUL_THRESHOLD);



// class MixtureData : public StatInterface , protected FrequencyDistribution {
class MixtureData : public StatInterface , public FrequencyDistribution {  // structure de donnees correspondant
                                                                           // a un melange de lois discretes
    friend class FrequencyDistribution;
    friend class Mixture;

    friend std::ostream& operator<<(std::ostream &os , const MixtureData &mixt_histo)
    { return mixt_histo.ascii_write(os , false); }

private :

    Mixture *mixture;       // pointeur sur un objet Mixture
    int nb_component;       // nombre de composantes
    FrequencyDistribution *weight;      // loi empirique des poids
    FrequencyDistribution **component;  // composantes empiriques

    void copy(const MixtureData &mixt_histo , bool model_flag = true);
    void remove();

public :

    MixtureData();
    MixtureData(const FrequencyDistribution &histo , int inb_component);
    MixtureData(const Mixture &mixt);
    MixtureData(const MixtureData &mixt_histo , bool model_flag = true)
    :FrequencyDistribution(mixt_histo) { copy(mixt_histo , model_flag); }
    ~MixtureData();
    MixtureData& operator=(const MixtureData &mixt_histo);

    DiscreteDistributionData* extract(StatError &error , int index) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    double information_computation() const;

    // acces membres de la classe

    Mixture* get_mixture() const { return mixture; }
    int get_nb_component() const { return nb_component; }
    FrequencyDistribution* get_weight() const { return weight; }
    FrequencyDistribution* get_component(int index) const { return component[index]; }
};



#endif
