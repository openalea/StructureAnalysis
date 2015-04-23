/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *        File author(s): Yann Guedon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_base.cpp 6703 2009-08-04 15:58:19Z cokelaer $
 *
 *-----------------------------------------------------------------------------*/



#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequence_analysis/sequences.h"
// #include "sequence_analysis/tops.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/hidden_variable_order_markov.h"
#include "sequence_analysis/nonhomogeneous_markov.h"
#include "sequence_analysis/renewal.h"
#include "sequence_analysis/sequence_label.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"


using namespace boost::python;
using namespace boost;
using namespace stat_tool;
using namespace sequence_analysis;



void class_constant_sequence()
{
   // sequences
  scope().attr("DEFAULT_LENGTH") = DEFAULT_LENGTH;
  scope().attr("SEQUENCE_NB_VARIABLE") = SEQUENCE_NB_VARIABLE;

  //const int PLOT_NB_SEQUENCE = 200;
  //const int PLOT_TITLE_NB_SEQUENCE = 15;

  enum_<sequence_analysis::wrap_util::UniqueInt<5, 100> > ("IndexParameterType")
     .value("IMPLICIT_TYPE", IMPLICIT_TYPE)
     .value("TIME", TIME)
     .value("TIME_INTERVAL", TIME_INTERVAL)
     .value("POSITION", POSITION)
     .value("POSITION_INTERVAL", POSITION_INTERVAL)
     .export_values();

  // OBERSEVER VALUES enumerate todo:

  enum_<sequence_analysis::wrap_util::UniqueInt<11, 101> > ("MarkovianSequenceType")
      .value("OBSERVATION", OBSERVATION)
      .value("INTENSITY", INTENSITY)
      .value("SELF_TRANSITION", SELF_TRANSITION)
      .value("COUNTING", COUNTING)
      .value("FIRST_OCCURRENCE", FIRST_OCCURRENCE)
      .value("RECURRENCE_TIME",  RECURRENCE_TIME)
      .value("SOJOURN_TIME", SOJOURN_TIME)
      .value("INITIAL_RUN", INITIAL_RUN)
      .value("FINAL_RUN", FINAL_RUN)
      .value("NB_RUN", NB_RUN)
      .value("NB_OCCURRENCE", NB_OCCURRENCE)
      .value("LENGTH", LENGTH)
      .value("SEQUENCE_CUMUL", SEQUENCE_CUMUL)
      .value("SEQUENCE_MEAN", SEQUENCE_MEAN)
      .export_values();

  // deFAULT/REVERSE enumerate todo:

  enum_<sequence_analysis::wrap_util::UniqueInt<4, 102> > ("Algorithm")
      .value("CTM_BIC", CTM_BIC)
      .value("CTM_KT", CTM_KT)
      .value("LOCAL_BIC",LOCAL_BIC)
      .value("CONTEXT", CONTEXT)
      .export_values();


  enum_<sequence_analysis::wrap_util::UniqueInt<4, 103> > ("Estimator")
    .value("MAXIMUM_LIKELIHOOD", MAXIMUM_LIKELIHOOD)
    .value("LAPLACE", LAPLACE)
    .value("ADAPTATIVE_LAPLACE", ADAPTATIVE_LAPLACE)
    .value("UNIFORM_SUBSET",UNIFORM_SUBSET)
    .value("UNIFORM_CARDINALITY", UNIFORM_CARDINALITY)
    .export_values();


  enum_<sequence_analysis::wrap_util::UniqueInt<6, 104> > ("ChangeType")
    .value("CATEGORICAL_CHANGE", CATEGORICAL_CHANGE)
    .value("POISSON_CHANGE", POISSON_CHANGE)
    .value("ORDINAL_GAUSSIAN_CHANGE", ORDINAL_GAUSSIAN_CHANGE)
    .value("GAUSSIAN_CHANGE", GAUSSIAN_CHANGE)
    .value("MEAN_CHANGE", MEAN_CHANGE)
    .value("VARIANCE_CHANGE", VARIANCE_CHANGE)
    .value("MEAN_VARIANCE_CHANGE", MEAN_VARIANCE_CHANGE)
    .export_values();

/*


  const double MAX_NB_WORD = 1.e7;       // nombre maximum de mots

  const double STATIONARY_PROBABILITY_THRESHOLD = 1.e-8;  // seuil pour le calcul des probabilites
                                                          // stationnaires d'un modele en equilibre
  const int STATIONARY_PROBABILITY_LENGTH = 10000;  // longueur maximum pour le calcul des probabilites
                                                    // stationnaires d'un modele en equilibre
  const double LEAVE_INCREMENT = 1.e-6;  // seuil pour stopper le calcul de la probabilite
                                         // de quitter un etat/observation sans possibilite d'y revenir

*/
  scope().attr("CTM_BIC_THRESHOLD") = CTM_BIC_THRESHOLD;
  scope().attr("CTM_KT_THRESHOLD") = CTM_KT_THRESHOLD;
  scope().attr("LOCAL_BIC_THRESHOLD") = LOCAL_BIC_THRESHOLD;
  scope().attr("CONTEXT_THRESHOLD") = CONTEXT_THRESHOLD;

  scope().attr("LOCAL_BIC_THRESHOLD") = LOCAL_BIC_THRESHOLD;
  scope().attr("OCCUPANCY_THRESHOLD") = OCCUPANCY_THRESHOLD;

  /*                                             // pour borner une loi d'occupation d'un etat
  const double OCCUPANCY_MEAN = 10.;     // temps moyen d'occupation d'un etat
*/


  scope().attr("MIN_NB_STATE_SEQUENCE") = MIN_NB_STATE_SEQUENCE;
  scope().attr("MAX_NB_STATE_SEQUENCE") = MAX_NB_STATE_SEQUENCE;
  scope().attr("NB_STATE_SEQUENCE_PARAMETER") = NB_STATE_SEQUENCE_PARAMETER;
  scope().attr("NB_STATE_SEQUENCE") = NB_STATE_SEQUENCE;


/*
  const int POSTERIOR_PROBABILITY_NB_SEQUENCE = 300; // nombre maximum de sequences pour la sortie des probabilites
                                                     // a posteriori des sequences d'etats les plus probables


  const int COUNT_MAX_LENGTH = 10000;    // longueur maximum des sequences pour l'extraction des comptages
  const int COUNTING_MAX_LENGTH = 500;   // longueur maximum des sequences pour le calcul des lois de comptage

  const int NB_SEQUENCE = 100000;        // nombre maximum de sequences simulees
  const int MAX_LENGTH = 1000000;        // longueur maximum des sequences simulees
  const int CUMUL_LENGTH = 1000000;      // longueur maximum cumulee des sequences simulees

*/

  enum_<sequence_analysis::wrap_util::UniqueInt<5, 105> >("OutputType")
    .value("SEQUENCE", SEQUENCE)
    .value("TREND", TREND)
    .value("SUBTRACTION_RESIDUAL", SUBTRACTION_RESIDUAL)
    .value("DIVISION_RESIDUAL", DIVISION_RESIDUAL)
    .value("STANDARDIZED_RESIDUAL",STANDARDIZED_RESIDUAL  )
    .export_values();


  enum_<sequence_analysis::wrap_util::UniqueInt<2, 106> > ("IndelCost")
     .value("ADAPTATIVE", ADAPTATIVE)
     .value("FIXED", FIXED)
     .export_values();
  /*

  enum {
    DELETION ,
    BEGIN_END_DELETION ,
    INSERTION ,
    BEGIN_END_INSERTION ,
    MATCH ,
    SUBSTITUTION ,
    TRANSPOSITION
  };
*/
/*enum {
    DATA ,
    GAP ,
    BEGIN_END_GAP
  };

  const int NB_ALIGNMENT = 1000000;      // nombre maximum d'alignements
  const int DISPLAY_NB_ALIGNMENT = 30;   // nombre maximum d'alignements
                                         // pour la sortie detaillee ecran
  const int FILE_NB_ALIGNMENT = 300;     // nombre maximum d'alignements
                                         // pour la sortie detaillee fichier
*/
  scope().attr("INDEL_FACTOR_1") = INDEL_FACTOR_1;
  scope().attr("INDEL_FACTOR_N") = INDEL_FACTOR_N;
  scope().attr("TRANSPOSITION_FACTOR") = TRANSPOSITION_FACTOR;

/*
  const double INDEL_DISTANCE = 1.0;     // cout d'elision/insertion
*/

  enum_<sequence_analysis::wrap_util::UniqueInt<2, 107> > ("ChangePointType")
     .value("CHANGE_POINT", CHANGE_POINT)
     .value("SEGMENT", SEGMENT)
     .export_values();
  
/*
    const double ROUNDOFF_ERROR = 1.e-10;  // erreur sur une somme de doubles
*/  
    scope().attr("NB_SEGMENTATION") = NB_SEGMENTATION;


  enum_<sequence_analysis::wrap_util::UniqueInt<2, 108> > ("NormType")
    .value("APPROXIMATED", APPROXIMATED)
    .value("EXACT", EXACT)
    .export_values();

/*
const double FREQUENCY_RATIO = 0.1;    // rapport des frequences pour stopper le calcul
                                       // de la fonction de correlation
// const int CORRELATION_MIN_FREQUENCY = 20;  frequence minimum pour stopper le calcul
                                           // de la fonction de correlation

const int MAX_DIFFERENCING_ORDER = 3;  // ordre maximum de differenciation
const int POINTWISE_AVERAGE_NB_SEQUENCE = 250;  // nombre maximum de sequences ecrites
                                                // pour la sortie fichier
const int ABSORBING_RUN_LENGTH = 5;    // longueur par defaut de la serie finale absorbante
const int MAX_ABSORBING_RUN_LENGTH = 20;  // longueur maximum de la serie finale absorbante


*/




  // others is this the same as in stat_tool, is it used ?
   enum_<sequence_analysis::wrap_util::UniqueInt<6, 109> > ("Type")
     .value("INT_VALUE", INT_VALUE)
     .value("REAL_VALUE", REAL_VALUE)
     .value("STATE", STATE)
     .value("OLD_INT_VALUE", OLD_INT_VALUE)
//     .value("NB_INTERNODE", NB_INTERNODE)
     .value("AUXILIARY", AUXILIARY)
     .export_values();







   // hidden semi markov
   /*const double SEMI_MARKOV_LIKELIHOOD_DIFF = 1.e-6;  // seuil pour stopper les iterations EM
   const int EXPLORATION_NB_ITER = 10;    // nombre d'iterations de la phase d'exploration
   const int STOCHASTIC_EXPLORATION_NB_ITER = 5;  // nombre d'iterations de la phase d'exploration
   const int SEMI_MARKOV_NB_ITER = 500;   // nombre maximum d'iterations EM

   const double MIN_SMOOTHED_PROBABILITY = 1.e-3;  // seuil sur les probabilites lissees
    */

   enum_<sequence_analysis::wrap_util::UniqueInt<3, 110> >("SemiMarkovState")
       .value("SSTATE", SSTATE)
       .value("IN_STATE", IN_STATE)
       .value("OUT_STATE", OUT_STATE)
       .export_values();

   // hidden variable order
   //const double VARIABLE_ORDER_MARKOV_LIKELIHOOD_DIFF = 1.e-6;  // seuil pour stopper les iterations EM
   //const int VARIABLE_ORDER_MARKOV_NB_ITER = 100;  // nombre maximum d'iterations EM


   // non homogeneous_markov
   /*
    * const double START_RATIO = 0.03;       // proportion de l'echantillon pour l'initialisation
                                       // des parametres (debut)
  const double END_RATIO = 0.1;          // proportion de l'echantillon pour l'initialisation
                                       // des parametres (fin)
  const int REGRESSION_NB_ELEMENT = 100;  // taille minimum de l'echantillon pour
                                        // la regression non-lineare
  const double GRADIENT_DESCENT_COEFF = 1.;  // coefficient algorithme de gradient
  const double RESIDUAL_SQUARE_SUM_DIFF = 1.e-6;  // seuil pour stopper les iterations
                                                  // de l'algorithme de gradient
  const int REGRESSION_NB_ITER = 1
    */


   // vom
   /*

*/

   scope().attr("MAX_LAG") = MAX_LAG;
  

/*                                        // d'autocorrelation
   const int MEMORY_MIN_COUNT = 10;       // effectif minimum pour comparer
                                          // une memoire a ses fils
   const double LAPLACE_COEFF = 1.;       // coefficient pour l'estimateur de Laplace

   enum {
     NON_TERMINAL ,
     TERMINAL ,
     COMPLETION
   };
    */

   // tops

   /*
const double MIN_RHYTHM_RATIO = 0.33;  // rapport de rythme minimum
*/
//   scope().attr("MIN_RHYTHM_RATIO") = MIN_RHYTHM_RATIO;
//   scope().attr("TOP_MIN_PROBABILITY") = TOP_MIN_PROBABILITY;
//   scope().attr("DEFAULT_MAX_POSITION") = DEFAULT_MAX_POSITION;
//   scope().attr("MAX_POSITION") = MAX_POSITION;
/*
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


    * */


   // renewal
/*
*/
   scope().attr("DEFAULT_TIME") = DEFAULT_TIME;
/*
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

   */

   enum_<sequence_analysis::wrap_util::UniqueInt<7, 111> >("RenewalType")
     .value("INTER_EVENT", INTER_EVENT)
     .value("WITHIN_OBSERVATION_PERIOD", WITHIN_OBSERVATION_PERIOD)
     .value("LENGTH_BIAS", LENGTH_BIAS)
     .value("BACKWARD_RECURRENCE_TIME", BACKWARD_RECURRENCE_TIME)
     .value("FORWARD_RECURRENCE_TIME", FORWARD_RECURRENCE_TIME)
     .value("NB_EVENT", NB_EVENT)
     .value("MIXTURE", MIXTURE)
     .export_values();
   /*

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


*/
    // semi markov
/*
   const int LEAVE_LENGTH = 10000;         // longueur maximum pour le calcul de
                                           // la probabilite de quitter un etat

   const double OCCUPANCY_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
   const int OCCUPANCY_NB_ITER = 10000;   // nombre maximum d'iterations EM
   const int OCCUPANCY_COEFF = 10;        // coefficient arrondi estimateur pour les lois
                                          // d'occupation des etats

   enum {
     MARKOVIAN ,
     SEMI_MARKOVIAN
   };
*/

}
