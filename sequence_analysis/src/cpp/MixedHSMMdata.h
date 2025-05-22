#ifndef MIXEDHSMMDATA_H
#define MIXEDHSMMDATA_H

#include "semi_markov.h"

namespace sequence_analysis {


/*
Attributs+Méthodes de [Sequences] dans "sequences.h" ligne 456-1027
Attributs+Méthodes de [MarkovianSequences] dans "sequences.h" ligne 1094-fin
Attributs+Méthodes de [SemiMarkovData] dans "semi_markov.h" ligne 303-fin


Méthodes détaillées de [Sequences] dans "sequences1.cpp", "sequences2.cpp", "sequences3.cpp"
Méthodes détaillées de [MarkovianSequences] dans "sequences1.cpp" ligne 2688-2760 et "markovian_sequences1.cpp", "markovian_sequences2.cpp"
Méthodes détaillées de [SemiMarkovData] dans "semi_markov.cpp" ligne 2797-fin
*/



/*##################################################################################################################
Objet contenant :
- les séquences observées ((Y^q_t)_t)_q, qui correspondent à un objet SemiMarkovData via l'héritage
- les covariables (X^q)_q [non dynamique], ou bien ((X^q_t)_t)_q [dynamique], qui sont contenues dans un objet Sequences
##################################################################################################################*/

class MixedHSMMdata : public SemiMarkovData
{

    private :
    /*************************************** Attributs *****************************************/

    bool covariate_are_dynamic;            // (true : dynamique, false : non dynamique)
    Sequences* covariate;                  // pointeur vers un objet Sequences contenant les covariables

    /*
    ----- Attributs "utiles" hérités de Sequences :

    int *identifier;        /// sequence identifiers  (????)
    int nb_sequence;        /// number of sequences
    int *length;            /// sequence lengths
    int max_length;         /// maximum sequence length
    int nb_variable;        /// number of variables
    stat_tool::variable_nature *type;  /// variable types (INT_VALUE/REAL_VALUE/STATE)

    int ***int_sequence;    /// sequences, integer-valued variables  (pointeur=NULL si REAL_VALUE)
    double ***real_sequence;  /// sequences, real-valued variables   (pointeur=NULL si INT_VALUE)

    ----- Attributs "utiles" hérités de MarkovianSequences :
    rien ?

    ----- Attributs "utiles" hérités de SemiMarkovData :
    double likelihood;        /// log-likelihood for the observed sequences
    SemiMarkov *semi_markov;  /// pointer on a SemiMarkov object  (pour le modèle)

    */
   /*******************************************************************************************/


    public :
    /*************************************** Méthodes ****************************************/

    // Constructeur nul par défaul (0 séquence produite)
    MixedHSMMdata ();

    // Constructeur par défaul
    MixedHSMMdata (bool _covariate_are_dynamic, const MarkovianSequences& observed_data, const Sequences& covariate_data);

    // Destructeur
    ~MixedHSMMdata ();



    // class member access
    bool get_covariate_are_dynamic() const { return covariate_are_dynamic; }
    Sequences* get_covariate() const { return new Sequences (*covariate); }

    /*****************************************************************************************/




};


}

#endif