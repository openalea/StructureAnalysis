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



#ifndef TOPS_H
#define TOPS_H



/****************************************************************
 *
 *  Constantes :
 */


const double TOP_MIN_PROBABILITY = 0.05;  // probabilite minimum
const double MIN_RHYTHM_RATIO = 0.33;  // rapport de rythme minimum
const int DEFAULT_MAX_POSITION = 20;   // position maximum par defaut
const int MAX_POSITION = 100;          // position maximum
const int PLOT_NB_AXILLARY = 10;       // nombre maximum de lois du nombre d'entrenoeuds
                                       // axes portes affichees (sortie Gnuplot)

const int MAX_NEIGHBORHOOD = 5;        // seuil sur le voisinage
const int NB_NEIGHBOR = 30;            // nombre minimum de voisins
const double PROBABILITY_DIFF = 0.05;  // ecart maximum pour considerer les probabilites
                                       // de croissance egales

const int NB_TOP = 10000;              // nombre maximum de cimes
const int NB_TRIAL = 1000;             // nombre maximum d'epreuves axe porteur
const int NB_AXILLARY = 4;             // nombre maximum d'axillaires par noeud
const int TOP_SIZE = 2000000;          // taille memoire maximum (en int) d'un ensemble de cimes



/****************************************************************
 *
 *  Definition des classes :
 */


class TopParameters : public StatInterface {  // parametres de croissance d'une cime

    friend class Tops;

    friend TopParameters* top_parameters_ascii_read(StatError &error , const char *path ,
                                                    int max_position);
    friend std::ostream& operator<<(std::ostream &os , const TopParameters &parameters)
    { return parameters.ascii_write(os , parameters.tops , false , false); }

private :

    Tops *tops;             // pointeur sur un objet Tops
    double probability;     // probabilite de croissance axe porteur
    double axillary_probability;  // probabilite de croissance axes portes
    double rhythm_ratio;    // rapport de rythme (axes portes / axe porteur)
    int max_position;       // position maximum
    Distribution **axillary_nb_internode;  // lois du nombre d'entrenoeuds axes portes

    void copy(const TopParameters &parameters , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const Tops *itops ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Tops *itops) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Tops *itops) const;
    MultiPlotSet* get_plotable(const Tops *itops) const;

    void position_nb_internode_computation(int position);

public :

    TopParameters(double iprobability = D_DEFAULT , double iaxillary_probability = D_DEFAULT ,
                  double irhythm_ratio = D_DEFAULT , int imax_position = 0);
    TopParameters(const TopParameters &parameters , bool data_flag = true)
    { copy(parameters , data_flag); }
    virtual ~TopParameters();
    TopParameters& operator=(const TopParameters&);

    DiscreteParametricModel* extract(StatError &error , int position) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void axillary_nb_internode_computation(int imax_position);

    Tops* simulation(StatError &error , int nb_top ,
                     const Distribution &nb_trial ,
                     const Distribution &nb_axillary) const;
    Tops* simulation(StatError &error , int nb_top ,
                     int nb_trial , int nb_axillary = 1) const;

    // acces membres de la classe

    Tops* get_tops() const { return tops; }
    double get_probability() const { return probability; }
    double get_axillary_probability() const { return axillary_probability; }
    double get_rhythm_ratio() const { return rhythm_ratio; }
    int get_max_position() const { return max_position; }
    Distribution* get_axillary_nb_internode(int position) const
    { return axillary_nb_internode[position]; }
};


TopParameters* top_parameters_ascii_read(StatError &error , const char *path ,
                                         int max_position = DEFAULT_MAX_POSITION);



class Tops : public Sequences {  // ensemble de cimes

    friend class TopParameters;

    friend Tops* tops_ascii_read(StatError &error , const char *path , bool old_format);
    friend std::ostream& operator<<(std::ostream &os , const Tops &tops)
    { return tops.ascii_write(os); }

private :

    TopParameters *top_parameters;  // pointeur sur un objet TopParameters
    FrequencyDistribution *nb_internode;  // loi empirique du nombre d'entrenoeuds axe porteur
    int max_position;        // position maximum
    FrequencyDistribution **axillary_nb_internode;  // lois empiriques du nombre d'entrenoeuds axes portes

    void copy(const Tops &tops , bool model_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;

    void max_position_computation();
    int min_position_computation() const;

public :

    Tops();
    Tops(int nb_top , int *iidentifier , int *nb_position , bool init_flag = true);
    Tops(const Sequences &seq);
    Tops(const Tops &tops , int inb_sequence , int *index);
    Tops(const Tops &tops , bool model_flag = true , bool reverse_flag = false);
    Tops(int nb_sample , const Tops **ptops);
    ~Tops();
    Tops& operator=(const Tops &tops);

    DiscreteDistributionData* extract(StatError &error , int position) const;

    Tops* shift(StatError &error , int inb_internode) const;
    Tops* select_individual(StatError &error , int inb_sequence ,
                            int *iidentifier , bool keep = true) const;
    Tops* reverse(StatError &error) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_data_write(std::ostream &os , char format = 'c' , bool exhaustive = false) const;
    bool ascii_data_write(StatError &error , const char *path ,
                          char format = 'c' , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void build_nb_internode_frequency_distribution();

    TopParameters* estimation(StatError &error , int imin_position , int imax_position ,
                              int neighborhood , bool equal_probability = false) const;
    TopParameters* estimation(StatError &error , int neighborhood , bool equal_probability = false) const
    { return estimation(error , 1 , max_position , neighborhood , equal_probability); }

    // acces membres de la classe

    TopParameters* get_top_parameters() const { return top_parameters; }
    FrequencyDistribution* get_nb_internode() const { return nb_internode; }
    int get_max_position() const { return max_position; }
    FrequencyDistribution* get_axillary_nb_internode(int position) const
    { return axillary_nb_internode[position]; }
};


Tops* tops_ascii_read(StatError &error , const char *path , bool old_format = false);



#endif
