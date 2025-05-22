#include "MixedHSMMdata.h"

namespace sequence_analysis {


// Constructeur nul par défaul (0 séquence produite)
MixedHSMMdata::MixedHSMMdata () :
SemiMarkovData::SemiMarkovData()
{
    covariate_are_dynamic = false;
    covariate = NULL;
}


// Il faudra effectuer un test en python pour "longueur observed_data == longueur covariate_data"
// Ce constructeur copie observed_data et covariate_data, donc pas optimal

// Constructeur à partir d'un objet SemiMarkovData pour les données observées, et d'un objet Sequences pour les covariables
MixedHSMMdata::MixedHSMMdata (bool _covariate_are_dynamic, const MarkovianSequences& observed_data, const Sequences& covariate_data) :
SemiMarkovData::SemiMarkovData(observed_data)
{

    covariate_are_dynamic = _covariate_are_dynamic;

    covariate = new Sequences (covariate_data);  // ça pose problème si l'objet
}


// Destructeur
MixedHSMMdata::~MixedHSMMdata ()
{
    if (covariate != NULL) {
        delete covariate;
        covariate = NULL;
    }
}

}




