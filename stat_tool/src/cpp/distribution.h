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
 *       $Id: distribution.h 17998 2015-04-23 06:56:11Z guedon $
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



#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H



namespace stat_tool {


  class DiscreteDistributionData;


  class DiscreteParametricModel : public StatInterface , public DiscreteParametric {  // loi discrete parametrique

    friend class Distribution;  // Hack pour Windows
    friend class FrequencyDistribution;
    friend class DiscreteDistributionData;

    friend DiscreteParametricModel* discrete_parametric_ascii_read(StatError &error ,
                                                                   const char *path ,
                                                                   double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const DiscreteParametricModel &dist)
    { return dist.ascii_write(os , dist.frequency_distribution , false , false); }

  private :

    DiscreteDistributionData *frequency_distribution;  // pointeur sur un objet DiscreteDistributionData

    std::ostream& ascii_write(std::ostream &os , const DiscreteDistributionData *histo ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const DiscreteDistributionData *histo) const;
    bool plot_write(StatError &error , const char *prefix , const char *title ,
                    const DiscreteDistributionData *histo) const;
    MultiPlotSet* get_plotable(const DiscreteDistributionData *histo) const;

  public :

    DiscreteParametricModel(int inb_value = 0 , int iident = CATEGORICAL ,
                     int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
                     double iparameter = D_DEFAULT, double iprobability = D_DEFAULT)
    :DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability)
    { frequency_distribution = NULL; }
    DiscreteParametricModel(int iident , int iinf_bound , int isup_bound , double iparameter ,
                            double iprobability , double cumul_threshold = CUMUL_THRESHOLD)
    :DiscreteParametric(iident , iinf_bound , isup_bound , iparameter , iprobability , cumul_threshold)
    { frequency_distribution = NULL; }
    DiscreteParametricModel(const FrequencyDistribution &histo);
    DiscreteParametricModel(const Distribution &dist)
    :DiscreteParametric(dist) { frequency_distribution = NULL; }
    DiscreteParametricModel(const DiscreteParametric &dist)
    :DiscreteParametric(dist) { frequency_distribution = NULL; }
    DiscreteParametricModel(const Distribution &dist , const FrequencyDistribution *histo);
    DiscreteParametricModel(const DiscreteParametric &dist , const FrequencyDistribution *histo);
    DiscreteParametricModel(const DiscreteParametricModel &dist , bool data_flag = true);
    ~DiscreteParametricModel();
    DiscreteParametricModel& operator=(const DiscreteParametricModel &dist);

    DiscreteDistributionData* extract_data(StatError &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    DiscreteDistributionData* simulation(StatError &error , int nb_element) const;

    DiscreteDistributionData* get_frequency_distribution() const { return frequency_distribution; }
  };


  DiscreteParametricModel* discrete_parametric_ascii_read(StatError &error , const char *path ,
                                                        double cumul_threshold = CUMUL_THRESHOLD);



  class DiscreteDistributionData : public StatInterface , public FrequencyDistribution {  // loi discrete empirique

    friend class DiscreteParametricModel;

    friend std::ostream& operator<<(std::ostream &os , const DiscreteDistributionData &histo)
    { return histo.ascii_write(os); }

  private :

    DiscreteParametricModel *distribution;  // pointeur sur un objet DiscreteParametricModel

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;

  public :

    DiscreteDistributionData(int inb_value = 0)
    :FrequencyDistribution(inb_value) { distribution = NULL; }
    DiscreteDistributionData(const Distribution &dist)
    :FrequencyDistribution(dist) { distribution = NULL; }
    DiscreteDistributionData(const FrequencyDistribution &histo)
    :FrequencyDistribution(histo) { distribution = NULL; }
    DiscreteDistributionData(int inb_element , int *pelement)
    :FrequencyDistribution(inb_element , pelement) { distribution = NULL; }
    DiscreteDistributionData(const FrequencyDistribution &histo , char transform ,
                             int param , int mode = FLOOR)
    :FrequencyDistribution(histo , transform , param , mode) { distribution = NULL; }
    DiscreteDistributionData(int nb_histo , const FrequencyDistribution **phisto)
    :FrequencyDistribution(nb_histo , phisto) { distribution = NULL; }
    DiscreteDistributionData(const FrequencyDistribution &histo , const Distribution *dist);
    DiscreteDistributionData(const FrequencyDistribution &histo , const DiscreteParametric *dist);
    DiscreteDistributionData(const DiscreteDistributionData &histo , bool model_flag = true);
    ~DiscreteDistributionData();
    DiscreteDistributionData& operator=(const DiscreteDistributionData &histo);

    DiscreteParametricModel* extract_model(StatError &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    DiscreteParametric* get_distribution() const { return distribution; }
  };


};  // namespace stat_tool



#endif
