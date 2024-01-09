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
 *       $Id: compound.h 17979 2015-04-23 06:37:30Z guedon $
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



#ifndef COMPOUND_H
#define COMPOUND_H



namespace stat_tool {



/****************************************************************
 *
 *  Constantes :
 */


  const double COMPOUND_THRESHOLD = 0.99999;  // seuil sur la fonction de repartition
                                              // pour borner une loi

  const double COMPOUND_INIT_PROBABILITY = 0.001;  // seuil pour l'initialisation de la probabilite
  const double COMPOUND_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
  const int COMPOUND_NB_ITER = 10000;    // nombre maximum d'iterations EM
  const double COMPOUND_DIFFERENCE_WEIGHT = 0.5;  // poids par defaut de la penalisation
                                                  // (cas des differences 1ere ou 2nde)
  const double COMPOUND_ENTROPY_WEIGHT = 0.1;  // poids par defaut de la penalisation (cas de l'entropie)
  const int COMPOUND_COEFF = 10;         // coefficient arrondi estimateur



/****************************************************************
 *
 *  Definition des classes :
 */


  class CompoundData;


  class Compound : public StatInterface , public Distribution {  // loi composee

    friend class FrequencyDistribution;
    friend class CompoundData;

    friend Compound* compound_ascii_read(StatError &error , const char *path ,
                                         double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Compound &compound)
    { return compound.ascii_write(os , compound.compound_data , false , false); }

  private :

    CompoundData *compound_data;  // pointeur sur un objet CompoundData
    DiscreteParametric *sum_distribution;  // loi de la somme
    DiscreteParametric *distribution;  // loi elementaire

    void copy(const Compound &compound , bool data_flag = true);

    std::ostream& ascii_write(std::ostream &os , const CompoundData *compound_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const CompoundData *compound_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const CompoundData *compound_histo) const;
    MultiPlotSet* get_plotable(const CompoundData *compound_histo) const;

    void computation(DiscreteParametric **power_dist , int min_nb_value ,
                     double cumul_threshold , bool sum_flag , bool dist_flag);
    void expectation_step(const FrequencyDistribution &histo , DiscreteParametric **power_dist ,
                          Reestimation<double> *sum_reestim , Reestimation<double> *reestim) const;

  public :

    Compound();
    Compound(const DiscreteParametric &sum_dist , const DiscreteParametric &dist ,
             double cumul_threshold = COMPOUND_THRESHOLD);
    Compound(const DiscreteParametric &sum_dist , const DiscreteParametric &dist , char type);
    Compound(const Compound &compound , bool data_flag = true)
    :Distribution(compound) { copy(compound , data_flag); }
    ~Compound();
    Compound& operator=(const Compound &compound);

    CompoundData* extract_data(StatError &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void computation(int min_nb_value = 1 ,
                     double cumul_threshold = COMPOUND_THRESHOLD ,
                     bool sum_flag = true , bool dist_flag = true);
    CompoundData* simulation(StatError &error , int nb_element) const;

    // acces membres de la classe

    CompoundData* get_compound_data() const { return compound_data; }
    DiscreteParametric* get_sum_distribution() const { return sum_distribution; }
    DiscreteParametric* get_distribution() const { return distribution; }
  };


  Compound* compound_ascii_read(StatError &error , const char *path ,
                                double cumul_threshold = COMPOUND_THRESHOLD);



  class CompoundData : public StatInterface , public FrequencyDistribution {  // structure de donnees correspondant
                                                                              // a une loi composee
    friend class FrequencyDistribution;
    friend class Compound;

    friend std::ostream& operator<<(std::ostream &os , const CompoundData &compound_histo)
    { return compound_histo.ascii_write(os , false); }

  private :

    Compound *compound;     // pointeur sur un objet Compound
    FrequencyDistribution *sum_frequency_distribution;   // loi empirique de la somme
    FrequencyDistribution *frequency_distribution;   // loi empirique elementaire

    void copy(const CompoundData &compound_histo , bool model_flag = true);

  public :

    CompoundData();
    CompoundData(const FrequencyDistribution &histo , const Compound &icompound);
    CompoundData(const Compound &icompound);
    CompoundData(const CompoundData &compound_histo , bool model_flag = true)
    :FrequencyDistribution(compound_histo) { copy(compound_histo , model_flag); }
    ~CompoundData();
    CompoundData& operator=(const CompoundData &compound_histo);

    DiscreteDistributionData* extract(StatError &error , char type) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    // acces membres de la classe

    Compound* get_compound() const { return compound; }
    FrequencyDistribution* get_sum_frequency_distribution() const { return sum_frequency_distribution; }
    FrequencyDistribution* get_frequency_distribution() const { return frequency_distribution; }
  };


};  // namespace stat_tool



#endif
