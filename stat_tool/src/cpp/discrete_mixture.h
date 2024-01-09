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
 *       $Id: discrete_mixture.h 17991 2015-04-23 06:46:12Z guedon $
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



#ifndef DISCRETE_MIXTURE_H
#define DISCRETE_MIXTURE_H



namespace stat_tool {



/****************************************************************
 *
 *  Constantes :
 */


  const int DISCRETE_MIXTURE_NB_COMPONENT = 100;   // nombre maximum de composantes

  const double NEGATIVE_BINOMIAL_PARAMETER = 20.;  // parametre initial pour une loi binomiale negative
  const double MIN_WEIGHT_STEP = 0.1;    // pas minimum d'initialisation des poids
  const double MAX_WEIGHT_STEP = 0.5;    // pas maximum d'initialisation des poids
  const int MIXTURE_COEFF = 2;           // coefficient arrondi estimateur
  const double DISCRETE_MIXTURE_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
  const int DISCRETE_MIXTURE_NB_ITER = 500;        // nombre maximum d'iterations EM



/****************************************************************
 *
 *  Definition des classes :
 */


  class DiscreteMixtureData;


  class DiscreteMixture : public StatInterface , public Distribution {  // melange de lois discretes

    friend class FrequencyDistribution;
    friend class DiscreteMixtureData;

    friend DiscreteMixture* discrete_mixture_building(StatError &error , int nb_component , double *weight ,
                                                      const DiscreteParametric **component);
    friend DiscreteMixture* discrete_mixture_ascii_read(StatError &error , const char *path ,
                                                        double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const DiscreteMixture &mixt)
    { return mixt.ascii_write(os , mixt.mixture_data , false , false); }

  private :

    DiscreteMixtureData *mixture_data;  // pointeur sur un objet DiscreteMixtureData
    int nb_component;       // nombre de composantes
    DiscreteParametric *weight;  // poids de chaque composante
    DiscreteParametric **component; // composantes

    void copy(const DiscreteMixture &mixt , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const DiscreteMixtureData *mixt_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const DiscreteMixtureData *mixt_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const DiscreteMixtureData *mixt_histo) const;
    MultiPlotSet* get_plotable(const DiscreteMixtureData *mixt_histo) const;

    int nb_parameter_computation() const;
    double penalty_computation() const;

    void init(const FrequencyDistribution &histo , bool *estimate , int min_inf_bound ,
              bool component_flag);
    void expectation_step(DiscreteMixtureData *mixt_histo , int nb_element) const;
    void variance_correction(DiscreteMixtureData *mixt_histo , bool *estimate ,
                             int min_inf_bound) const;
    bool component_order_test() const;

  public :

    DiscreteMixture();
    DiscreteMixture(int inb_component , double *pweight , const DiscreteParametric **pcomponent);
    DiscreteMixture(const DiscreteMixture &mixt , bool *component_flag , int inb_value);
    DiscreteMixture(int inb_component , const DiscreteParametric **pcomponent);
    DiscreteMixture(const DiscreteMixture &mixt , bool data_flag = true)
    :Distribution(mixt) { copy(mixt , data_flag); }
    ~DiscreteMixture();
    DiscreteMixture& operator=(const DiscreteMixture &mixt);

    DiscreteParametricModel* extract(StatError &error , int index) const;
    DiscreteMixtureData* extract_data(StatError &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void computation(int min_nb_value = 1 , double cumul_threshold = CUMUL_THRESHOLD ,
                     bool component_flag = true);
    double likelihood_computation(const DiscreteMixtureData &mixt_histo) const;
    DiscreteMixtureData* simulation(StatError &error , int nb_element) const;

    // acces membres de la classe

    DiscreteMixtureData* get_mixture_data() const { return mixture_data; }
    int get_nb_component() const { return nb_component; }
    DiscreteParametric* get_weight() const { return weight; }
    DiscreteParametric* get_component(int index) const { return component[index]; }
  };


  DiscreteMixture* discrete_mixture_building(StatError &error , int nb_component , double *weight ,
                                             const DiscreteParametric **component);
  DiscreteMixture* discrete_mixture_ascii_read(StatError &error , const char *path ,
                                               double cumul_threshold = CUMUL_THRESHOLD);



  class DiscreteMixtureData : public StatInterface , public FrequencyDistribution {  // structure de donnees correspondant
                                                                                     // a un melange de lois discretes
    friend class FrequencyDistribution;
    friend class DiscreteMixture;

    friend std::ostream& operator<<(std::ostream &os , const DiscreteMixtureData &mixt_histo)
    { return mixt_histo.ascii_write(os , false); }

  private :

    DiscreteMixture *mixture;  // pointeur sur un objet DiscreteMixture
    int nb_component;       // nombre de composantes
    FrequencyDistribution *weight;  // loi empirique des poids
    FrequencyDistribution **component;  // composantes empiriques

    void copy(const DiscreteMixtureData &mixt_histo , bool model_flag = true);
    void remove();

  public :

    DiscreteMixtureData();
    DiscreteMixtureData(const FrequencyDistribution &histo , int inb_component);
    DiscreteMixtureData(const FrequencyDistribution &histo , const DiscreteMixture *pmixture);
    DiscreteMixtureData(const DiscreteMixture &mixt);
    DiscreteMixtureData(const DiscreteMixtureData &mixt_histo , bool model_flag = true)
    :FrequencyDistribution(mixt_histo) { copy(mixt_histo , model_flag); }
    ~DiscreteMixtureData();
    DiscreteMixtureData& operator=(const DiscreteMixtureData &mixt_histo);

    DiscreteDistributionData* extract(StatError &error , int index) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    double information_computation() const;

    // acces membres de la classe

    DiscreteMixture* get_mixture() const { return mixture; }
    int get_nb_component() const { return nb_component; }
    FrequencyDistribution* get_weight() const { return weight; }
    FrequencyDistribution* get_component(int index) const { return component[index]; }
  };


};  // namespace stat_tool



#endif
