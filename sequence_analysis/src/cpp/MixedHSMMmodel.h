#ifndef MIXEDHSMMMODEL_H
#define MIXEDHSMMMODEL_H

#include "hidden_semi_markov.h"

namespace sequence_analysis {

/*
Attributs+Méthodes de [SemiMarkovChain] dans "semi_markov.h" ligne 78-109
Attributs+Méthodes de [SemiMarkov] dans "semi_markov.h" ligne 127-264
Attributs+Méthodes de [HiddenSemiMarkov] dans "hidden_semi_markov.h" ligne 79-fin


Méthodes détaillées de [SemiMarkovChain] dans "semi_markov.cpp" ligne début-311
Méthodes détaillées de [SemiMarkov] dans "semi_markov.cpp" ligne 314-2792
Méthodes détaillées de [HiddenSemiMarkov] dans "semi_markov.cpp" ligne 2795-fin
*/



/*##################################################################################################################
Objet contenant :

##################################################################################################################*/

class MixedHSMMmodel: public HiddenSemiMarkov
{

    private :
    /*************************************** Attributs *****************************************/

    

   /*******************************************************************************************/


    public :
    /*************************************** Méthodes ****************************************/

    // Constructeur nul par défaul (0 séquence produite)
    MixedHSMMmodel ();

    // Constructeur par défaul
    //MixedHSMMmodel ();



    // class member access


    /*****************************************************************************************/




};

}


#endif