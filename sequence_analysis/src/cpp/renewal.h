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



#ifndef RENEWAL_H
#define RENEWAL_H



/****************************************************************
 *
 *  Constantes :
 */


const int DEFAULT_TIME = 20;           // temps d'observation par defaut
const int MAX_TIME = 500;              // temps d'observation maximum
const int PLOT_NEVENT_TIME = 10;       // nombre maximum de lois du temps avant
                                       // le n-eme evenement affichees (sortie Gnuplot)
const int PLOT_NB_TIME = 5;            // nombre maximum de lois du nombre d'evenements
                                       // affichees dans le cas de melange (sortie Gnuplot)

const double RENEWAL_THRESHOLD = 0.99999;  // seuil sur la fonction de repartition
                                           // pour borner une loi
const double RB_THRESHOLD = 2000.;     // seuil pour utiliser le calcul rapide de la loi du
                                       // nombre d'ev correspondant a une loi binomiale
const double RNB_THRESHOLD = 2000.;    // seuil pour utiliser le calcul rapide de la loi du
                                       // nombre d'ev correspondant a une loi binomiale negative

enum {
  INTER_EVENT ,
  WITHIN_OBSERVATION_PERIOD ,
  LENGTH_BIAS ,
  BACKWARD_RECURRENCE_TIME ,
  FORWARD_RECURRENCE_TIME ,
  NB_EVENT ,
  MIXTURE
};

const double MIN_NB_EVENT = 0.4;       // nombre d'evenements moyen minimum
const double MIN_INTER_EVENT = 1.;     // temps moyen minimum entre 2 evenements
const double RENEWAL_INIT_PROBABILITY = 0.001;  // seuil pour l'initialisation de la probabilite
const double RENEWAL_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
const int RENEWAL_NB_ITER = 10000;     // nombre maximum d'iterations EM
const double RENEWAL_DIFFERENCE_WEIGHT = 0.5;  // poids par defaut de la penalisation
                                               // (cas des differences 1ere ou 2nde)
const double RENEWAL_ENTROPY_WEIGHT = 0.05;  // poids par defaut de la penalisation (cas de l'entropie)
const int RENEWAL_COEFF = 10;          // coefficient arrondi estimateur

const int NB_COMPLETE_INTERVAL = 3;    // nombre minimum d'intervalles de temps complets
const double MEAN_COEFF = 2.;          // coefficient sur la moyenne pour compenser le biais par la longueur
const double MAX_VALUE_COEFF = 10.;    // coefficient pour deduire la valeur maximum de la loi inter-evenement

const int RENEWAL_NB_ELEMENT = 1000000;  // taille maximum de l'echantillon pour la simulation



/****************************************************************
 *
 *  Definition des classes :
 */


class LengthBias : public DiscreteParametric {  // loi biaisee par la longueur

/*    friend class Renewal;
    friend class RenewalData; */

public :

    LengthBias(int inb_value = 0 , int iident = NONPARAMETRIC ,
               int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
               double iparameter = D_DEFAULT, double iprobability = D_DEFAULT)
    :DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    LengthBias(const DiscreteParametric &inter_event)
    :DiscreteParametric(inter_event) { computation(inter_event); }
    LengthBias(const LengthBias &length_bias)
    :DiscreteParametric((DiscreteParametric&)length_bias) {}

    void computation(const DiscreteParametric&);
};



class Backward : public DiscreteParametric {  // loi de l'intervalle de temps apres le dernier evenement

/*    friend class Renewal;
    friend class RenewalData; */

public :

    Backward(int inb_value = 0 , int iident = NONPARAMETRIC ,
             int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
             double iparameter = D_DEFAULT, double iprobability = D_DEFAULT)
    :DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Backward(const Backward &dist , int ialloc_nb_value = I_DEFAULT)
    :DiscreteParametric(dist , 'c' , ialloc_nb_value) {}

    void computation(const DiscreteParametric &inter_event , const Distribution &time);
};



class NbEvent : public DiscreteParametric {  // loi du nombre d'evenements

//    friend class Renewal;

// private :
public :

    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int time;               // temps d'observation

    void binomial_computation();
    void negative_binomial_computation();

// public :

    NbEvent(char itype = 'v' , int itime = 0 , int inb_value = 0 , int iident = NONPARAMETRIC ,
            int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
            double iparameter = D_DEFAULT , double iprobability = D_DEFAULT);
    NbEvent(char itype , int itime , DiscreteParametric &inter_event);
    NbEvent(const NbEvent &nb_event , int ialloc_nb_value = I_DEFAULT);

    void ordinary_computation(DiscreteParametric &inter_event);
    void computation(DiscreteParametric &inter_event);

    // acces membres de la classe

/*    char get_type() const { return type; }
    int get_time() const { return time; } */
};



class RenewalIterator;
class TimeEvents;
class RenewalData;

class Renewal : public StatInterface {  // processus de renouvellement

    friend class RenewalIterator;
    friend class TimeEvents;
    friend class RenewalData;

    friend Renewal* renewal_building(StatError &error , const DiscreteParametric &inter_event ,
                                     char type, int time);
    friend Renewal* renewal_ascii_read(StatError& error , const char *path ,
                                       char type, int time,
                                       double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Renewal &renew)
    { return renew.ascii_write(os , renew.renewal_data , false , false); }

private :

    int nb_iterator;        // nombre d'iterateurs pointant sur l'objet
    RenewalData *renewal_data;  // pointeur sur un objet RenewalData
    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int nb_event_max;       // borne max sur le nombre d'evenements
    Distribution *time;     // loi du temps d'observation
    DiscreteParametric *inter_event;  // loi inter-evenement
    LengthBias *length_bias;  // loi biaisee par la longueur
    Backward *backward;     // loi de l'intervalle de temps apres le dernier evenement
    Forward *forward;       // loi de l'intervalle de temps residuel
    DiscreteParametric **nevent_time;  // lois du temps avant le n-eme evenement
    NbEvent **nb_event;    // lois du nombre d'evenements
                           // pour un temps d'observation donne
    Distribution *mixture;  // melange de lois du nombre d'evenements
    Curves *index_event;    // probabilites de non-evenement/evenement fonction du temps

    void init(int inf_bound , int sup_bound , double parameter , double probability);
    void init(int ident , int inf_bound , int sup_bound , double parameter , double probability);
    void copy(const Renewal &renew , bool data_flag = true);
    void remove();
    void type_init(int itype);

    std::ostream& ascii_write(std::ostream &os , const RenewalData *timev ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os , const RenewalData *timev) const;
    bool plot_write(const char *prefix , const char *title ,
                    const RenewalData *timev) const;
    MultiPlotSet* get_plotable(const RenewalData *timev) const;

    void index_event_computation();

    void expectation_step(const TimeEvents &timev , Reestimation<double> *reestim) const;
    void expectation_step(const TimeEvents &timev , Reestimation<double> *inter_event_reestim ,
                          Reestimation<double> *length_bias_reestim , int estimator ,
                          bool combination = false , int mean_computation = COMPUTED) const;

public :

    Renewal();
    Renewal(char itype , const FrequencyDistribution &htime ,
            const DiscreteParametric &iinter_event);
    Renewal(char itype , const Distribution &itime ,
            const DiscreteParametric &iinter_event);
    Renewal(const RenewalData &irenewal_data ,
            const DiscreteParametric &iinter_event);
    Renewal(const Renewal &renew , bool data_flag = true)
    { copy(renew , data_flag); }
    virtual ~Renewal();
    void conditional_delete();
    Renewal& operator=(const Renewal &renew);

    DiscreteParametricModel* extract(StatError &error , int dist_type ,
                                     int itime = I_DEFAULT) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void computation(bool inter_event_flag = true , char itype = 'v' ,
                     const Distribution *dtime = NULL);

    double likelihood_computation(const TimeEvents &timev) const;

    RenewalData* simulation(StatError &error , char itype ,
                            const FrequencyDistribution &ihtime) const;
    RenewalData* simulation(StatError &error , char itype ,
                            int nb_element , int itime) const;
    RenewalData* simulation(StatError &error , char itype ,
                            int nb_element , const TimeEvents &itimev) const;

    // acces membres de la classe

    int get_nb_iterator() const { return nb_iterator; }
    RenewalData* get_renewal_data() const { return renewal_data; }
    char get_type() const { return type; }
    Distribution* get_time() const { return time; }
    DiscreteParametric* get_inter_event() const { return inter_event; }
    LengthBias* get_length_bias() const { return length_bias; }
    Backward* get_backward() const { return backward; }
    Forward* get_forward() const { return forward; }
    DiscreteParametric* get_nevent_time(int inb_event) const { return nevent_time[inb_event]; }
    NbEvent* get_nb_event(int itime) const { return nb_event[itime]; }
    Distribution* get_mixture() const { return mixture; }
    Curves* get_index_event() const { return index_event; }
};


Renewal* renewal_building(StatError &error , const DiscreteParametric &inter_event ,
                          char type = 'e' , int time = DEFAULT_TIME);
Renewal* renewal_ascii_read(StatError& error , const char *path ,
                            char type = 'e' , int time = DEFAULT_TIME ,
                            double cumul_threshold = RENEWAL_THRESHOLD);



class RenewalIterator {  // iterateur processus de renouvellement

private :

    Renewal *renewal;       // pointeur sur un objet Renewal
    int interval;           // intervalle de temps
    int counter;            // compteur
    int length;             // longueur de la sequence
    int *sequence;          // sequence

    void copy(const RenewalIterator &iterator);

public :

    RenewalIterator(Renewal *irenewal , int ilength = 1);
    RenewalIterator(const RenewalIterator &iterator)
    { copy(iterator); }
    ~RenewalIterator();
    RenewalIterator& operator=(const RenewalIterator &iterator);

    void simulation(int ilength = 1 , char type = 'v');

    // acces membres de la classe

    Renewal* get_renewal() const { return renewal; }
    int get_interval() const { return interval; }
    int get_counter() const { return counter; }
    int get_length() const { return length; }
    int get_sequence(int index) const { return sequence[index]; }
};



class TimeEvents : public StatInterface {  // echantillons {temps, nombre d'evenements, effectif}

    friend class FrequencyDistribution;
    friend class Renewal;

    friend TimeEvents* time_events_ascii_read(StatError &error , const char *path);
    friend TimeEvents* old_time_events_ascii_read(StatError &error , const char *path);
    friend std::ostream& operator<<(std::ostream &os , const TimeEvents &timev)
    { return timev.ascii_write(os , true); }

protected :

    int nb_element;         // effectif total
    int nb_class;           // nombre de classes
    int *time;              // temps d'observation
    int *nb_event;          // nombre d'evenements
    int *frequency;         // effectif de chacune des classes
                            // {temps, nombre d'evenements}
    FrequencyDistribution *htime;  // loi empirique du temps d'observation
    FrequencyDistribution **hnb_event;  // lois empiriques du nombre d'evenements
                                        // pour un temps d'observation donne
    FrequencyDistribution *mixture;  // loi empirique du nombre d'evenements

    void build_frequency_distribution();
    void build_sample();
    void build(int inb_element , int *itime , int *inb_event);
    void copy(const TimeEvents&);
    void merge(int nb_sample , const TimeEvents **ptimev);
    void remove();

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , char type) const;
    std::ostream& ascii_file_write(std::ostream &os , bool exhaustive , char type = 'v') const;
    std::ostream& spreadsheet_write(std::ostream &os , char type = 'v') const;

    void nb_element_computation();
    double min_inter_event_computation() const;

public :

    TimeEvents(int inb_class = 0);
    TimeEvents(int inb_element , int *itime , int *inb_event)
    { build(inb_element , itime , inb_event); }
    TimeEvents(int nb_sample , const TimeEvents **ptimev) { merge(nb_sample , ptimev); }
    TimeEvents(const TimeEvents &timev) { copy(timev); }
    virtual ~TimeEvents();
    TimeEvents& operator=(const TimeEvents &timev);

    DiscreteDistributionData* extract(StatError &error , int histo_type ,
                                      int itime = I_DEFAULT) const;

    TimeEvents* time_scaling(StatError &error , int scaling_coeff) const;
    TimeEvents* time_select(StatError &error , int min_time ,
                            int max_time) const;
    TimeEvents* nb_event_select(StatError &error , int min_nb_event ,
                                int max_nb_event) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = true) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = true) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    double information_computation() const;

    Renewal* estimation(StatError &error , std::ostream &os , char type ,
                        const DiscreteParametric &iinter_event , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int equilibrium_estimator = COMPLETE_LIKELIHOOD ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Renewal* estimation(StatError &error , std::ostream &os , char type , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int equilibrium_estimator = COMPLETE_LIKELIHOOD ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;

    // acces membres de la classe

    int get_nb_element() const { return nb_element; }
    int get_nb_class() const { return nb_class; }
    FrequencyDistribution* get_htime() const { return htime; }
    FrequencyDistribution* get_hnb_event(int itime) const { return hnb_event[itime]; }
    FrequencyDistribution* get_mixture() const { return mixture; }
};


TimeEvents* time_events_ascii_read(StatError &error , const char *path);
TimeEvents* old_time_events_ascii_read(StatError &error , const char *path);



class RenewalData : public TimeEvents {  // donnees correspondant a
                                           // un processus de renouvellement
    friend class Renewal;
    friend class Sequences;

    friend std::ostream& operator<<(std::ostream &os , RenewalData &timev)
    { return timev.ascii_write(os , false); }

private :

    Renewal *renewal;       // pointeur sur un objet Renewal
    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int *length;            // longueurs des sequences
    int **sequence;         // sequences
    FrequencyDistribution *inter_event; // intervalles de temps entre 2 evenements
    FrequencyDistribution *within;  // intervalles de temps entre 2 evenements a l'interieur
                                    // de la periode d'observation
    FrequencyDistribution *length_bias;  // intervalles de temps entre 2 evenements
                             // recouvrant une date d'observation
    FrequencyDistribution *backward;  // intervalles de temps apres le dernier evenement
    FrequencyDistribution *forward;  // intervalles de temps residuel
    Curves *index_event;    // probabilites de non-evenement/evenement fonction du temps

    void copy(const RenewalData &timev , bool model_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_write(std::ostream &os) const;

    void build_index_event(int offset = 1);

public :

    RenewalData();
    RenewalData(int nb_element , int itime);
    RenewalData(int itype , const Renewal &renew);
    RenewalData(const TimeEvents &timev , int itype);
    RenewalData(int nb_sample , const RenewalData **itimev);
    RenewalData(const RenewalData &timev , bool model_flag = true)
    :TimeEvents(timev) { copy(timev , model_flag); }
    ~RenewalData();
    RenewalData& operator=(const RenewalData&);

    RenewalData* merge(StatError &error , int nb_sample , const RenewalData **itimev) const;
    DiscreteDistributionData* extract(StatError &error , int histo_type ,
                                      int itime = I_DEFAULT) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    Renewal* estimation(StatError &error , std::ostream &os , const DiscreteParametric &iinter_event ,
                        int estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                        int mean_computation = COMPUTED , double weight = D_DEFAULT ,
                        int penalty_type = SECOND_DIFFERENCE , int outside = ZERO) const;
    Renewal* estimation(StatError &error , std::ostream &os , int estimator = LIKELIHOOD ,
                        int nb_iter = I_DEFAULT , int mean_computation = COMPUTED ,
                        double weight = D_DEFAULT , int penalty_type = SECOND_DIFFERENCE ,
                        int outside = ZERO) const;

    // acces membres de la classe

    Renewal* get_renewal() const { return renewal; }
    char get_type() const { return type; }
    int get_length(int index_seq) const { return length[index_seq]; }
    int get_sequence(int index_seq , int index) const
    { return sequence[index_seq][index]; }
    FrequencyDistribution* get_inter_event() const { return inter_event; }
    FrequencyDistribution* get_within() const { return within; }
    FrequencyDistribution* get_length_bias() const { return length_bias; }
    FrequencyDistribution* get_backward() const { return backward; }
    FrequencyDistribution* get_forward() const { return forward; }
    Curves* get_index_event() const { return index_event; }
};



#endif
