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
 *       $Id: convolution.h 17985 2015-04-23 06:42:55Z guedon $
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



#ifndef CONVOLUTION_H
#define CONVOLUTION_H



namespace stat_tool {



/****************************************************************
 *
 *  Constantes :
 */


  const int CONVOLUTION_NB_DISTRIBUTION = 10;  // nombre maximum de lois elementaires
  const double CONVOLUTION_THRESHOLD = 0.9999;  // seuil sur la fonction de repartition
                                              // pour borner une loi

  const double CONVOLUTION_INIT_PROBABILITY = 0.001;  // seuil pour l'initialisation de la probabilite
  const double CONVOLUTION_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
  const int CONVOLUTION_NB_ITER = 10000;  // nombre maximum d'iterations EM
  const double CONVOLUTION_DIFFERENCE_WEIGHT = 0.5;  // poids par defaut de la penalisation
                                                   // (cas des differences 1ere ou 2nde)
  const double CONVOLUTION_ENTROPY_WEIGHT = 0.1;  // poids par defaut de la penalisation (cas de l'entropie)
  const int CONVOLUTION_COEFF = 10;      // coefficient arrondi estimateur



/****************************************************************
 *
 *  Definition des classes :
 */


  class ConvolutionData;


  class Convolution : public StatInterface , public Distribution {  // produit de convolution
                                                                    // de lois discretes
    friend class FrequencyDistribution;
    friend class ConvolutionData;

    friend Convolution* convolution_building(StatError &error , int nb_dist ,
                                             const DiscreteParametric **dist);
    friend Convolution* convolution_ascii_read(StatError &error , const char *path ,
                                               double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Convolution &convol)
    { return convol.ascii_write(os , convol.convolution_data , false , false); }

  private :

    ConvolutionData *convolution_data;  // pointeur sur un objet ConvolutionData
    int nb_distribution;    // nombre de lois elementaires
    DiscreteParametric **distribution;  // lois elementaires

    void copy(const Convolution &convol , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const ConvolutionData *convol_histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const ConvolutionData *convol_histo) const;
    bool plot_write(const char *prefix , const char *title ,
                    const ConvolutionData *convol_histo) const;
    MultiPlotSet* get_plotable(const ConvolutionData *convol_histo) const;

    void expectation_step(const FrequencyDistribution &histo , const Distribution **partial_convol ,
                          Reestimation<double> **reestim) const;

  public :

    Convolution();
    Convolution(int nb_dist , const DiscreteParametric **pdist);
    Convolution(const DiscreteParametric &known_dist , const DiscreteParametric &unknown_dist);
    Convolution(const Convolution &convol , bool data_flag = true)
    :Distribution(convol) { copy(convol , data_flag); }
    ~Convolution();
    Convolution& operator=(const Convolution &convol);

    DiscreteParametricModel* extract(StatError &error , int index) const;
    ConvolutionData* extract_data(StatError &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void computation(int min_nb_value = 1 , double cumul_threshold = CONVOLUTION_THRESHOLD ,
                     bool *dist_flag = NULL);
    ConvolutionData* simulation(StatError &error , int nb_element) const;

    // acces membres de la classe

    ConvolutionData* get_convolution_data() const { return convolution_data; }
    int get_nb_distribution() const { return nb_distribution; }
    DiscreteParametric* get_distribution(int index) const { return distribution[index]; }
  };


  Convolution* convolution_building(StatError &error , int nb_dist ,
                                    const DiscreteParametric **dist);
  Convolution* convolution_ascii_read(StatError &error , const char *path ,
                                      double cumul_threshold = CONVOLUTION_THRESHOLD);



  class ConvolutionData : public StatInterface , public FrequencyDistribution {  // structure de donnees correspondant
                                                                                 // a un produit de convolution de lois discretes
    friend class FrequencyDistribution;
    friend class Convolution;

    friend std::ostream& operator<<(std::ostream &os , const ConvolutionData &convol_histo)
    { return convol_histo.ascii_write(os , false); }

  private :

    Convolution *convolution;  // pointeur sur un objet Convolution
    int nb_distribution;       // nombre de lois elementaires empiriques
    FrequencyDistribution **frequency_distribution;  // lois elementaires empiriques

    void copy(const ConvolutionData &convol_histo , bool model_flag = true);
    void remove();

  public :

    ConvolutionData();
    ConvolutionData(const FrequencyDistribution &histo , int nb_dist);
    ConvolutionData(const Convolution &convol);
    ConvolutionData(const ConvolutionData &convol_histo , bool model_flag = true)
    :FrequencyDistribution(convol_histo) { copy(convol_histo , model_flag); }
    ~ConvolutionData();
    ConvolutionData& operator=(const ConvolutionData &convol_histo);

    DiscreteDistributionData* extract(StatError &error , int index) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    // acces membres de la classe

    Convolution* get_convolution() const { return convolution; }
    int get_nb_distribution() const { return nb_distribution; }
    FrequencyDistribution* get_frequency_distribution(int index) const { return frequency_distribution[index]; }
  };


};  // namespace stat_tool



#endif
