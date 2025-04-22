#ifndef CLASSE_HPP_INCLUDED
#define CLASSE_HPP_INCLUDED

//#include "sequence_analysis/sequences.h"


/*######################################
Objet contenant la séquence observée (Y_t)_t, et les covariables.
######################################*/

class MixedHSMMData //public //SemiMarkovData
{

    private :
    // Attributs





    public :
    // Méthodes   




};


/*######################################
Objet contenant les covariables (X^q)_q [non dynamique], ou bien ((X^q_t)_t)_q [dynamique]
######################################*/

class CovariateData
{

    private :
    // Attributs

    bool m_is_dynamic;      // (true : dynamique, false : non dynamique)
    int m_nb_individual;    // (nombre d'individus/séquences, c'est à dire Q)
    int m_nb_covariate;     // (nombre de covariables, c'est à dire dimension de X^q ou X^q_t)

    double m_data[5];






    public :
    // Méthodes

    CovariateData (bool is_dynamic, int nb_individual, int nb_covariate, double data[5]);  // Constructeur 



};


#endif