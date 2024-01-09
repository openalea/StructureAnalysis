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
 *       $Id: renewal_algorithms.cpp 18065 2015-04-23 10:50:49Z guedon $
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



#include <math.h>
#include <iomanip>

#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "renewal.h"
#include "sequence_label.h"

using namespace std;



namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Calcul de la vraisemblance des donnees d'intervalles de temps.
 *
 *  arguments : loi de l'intervalle de temps residuel, lois empiriques des intervalles
 *              de temps complets, censures a gauche, a droite et
 *              de la longueur de la periode d'observation dans le cas 0 evenement.
 *
 *--------------------------------------------------------------*/

double DiscreteParametric::renewal_likelihood_computation(const Forward &forward_dist ,
                                                          const FrequencyDistribution &within ,
                                                          const FrequencyDistribution &backward ,
                                                          const FrequencyDistribution &forward ,
                                                          const FrequencyDistribution *no_event) const

{
  double likelihood , buff;
  FrequencyDistribution *histo;


  likelihood = likelihood_computation(within);

  if (likelihood != D_INF) {
    histo = new FrequencyDistribution(backward , 's' , 1);
    buff = survivor_likelihood_computation(*histo);
    delete histo;

    if (buff != D_INF) {
      likelihood += buff;
      buff = forward_dist.likelihood_computation(forward);

      if (buff != D_INF) {
        likelihood += buff;

        if (no_event) {
          histo = new FrequencyDistribution(*no_event , 's' , 1);
          buff = forward_dist.survivor_likelihood_computation(*histo);
          delete histo;

          if (buff != D_INF) {
            likelihood += buff;
          }
          else {
            likelihood = D_INF;
          }
        }
      }

      else {
        likelihood = D_INF;
      }
    }

    else {
      likelihood = D_INF;
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation correspondant a la loi inter-evenement
 *  (estimateur EM d'un processus de renouvellement en equilibre a partir
 *   de donnees d'intervalles de temps).
 *
 *  arguments : lois empiriques des intervalles de temps complets, censures a gauche,
 *              a droite et de la longueur de la periode d'observation dans le cas 0 evenement,
 *              pointeurs sur les quantites de reestimation de la loi inter-evenement et
 *              de la loi biaisee par la longueur.
 *
 *--------------------------------------------------------------*/

void DiscreteParametric::expectation_step(const FrequencyDistribution &within ,
                                          const FrequencyDistribution &backward ,
                                          const FrequencyDistribution &forward ,
                                          const FrequencyDistribution *no_event ,
                                          Reestimation<double> *inter_event_reestim ,
                                          Reestimation<double> *length_bias_reestim , int iter) const

{
  register int i , j;
  int reestim_offset , reestim_nb_value , *pfrequency;
  double sum , *ifrequency , *lfrequency , *pmass , *pcumul , *norm;


  // initialisations

  ifrequency = inter_event_reestim->frequency;
  lfrequency = length_bias_reestim->frequency;
  for (i = 0;i < inter_event_reestim->alloc_nb_value;i++) {
    *ifrequency++ = 0.;
    *lfrequency++ = 0.;
  }

  // calcul des quantites de reestimation de la loi inter-evenement

  ifrequency = inter_event_reestim->frequency + within.offset;
  pfrequency = within.frequency + within.offset;
  for (i = within.offset;i < within.nb_value;i++) {
    *ifrequency++ += *pfrequency++;
  }

  pfrequency = backward.frequency + backward.offset;
  pcumul = cumul + backward.offset;
  ifrequency = inter_event_reestim->frequency + backward.offset + 1;
  pmass = mass + backward.offset + 1;
  sum = 0.;

  for (i = backward.offset;i < backward.nb_value;i++) {
    sum += *pfrequency++ / (1. - *pcumul++);
    *ifrequency++ += *pmass++ * sum;
  }
  for (i = backward.nb_value;i < nb_value - 1;i++) {
    *ifrequency++ += *pmass++ * sum;
  }

  // calcul des quantites de reestimation de la loi biaisee par la longueur

  pfrequency = forward.frequency + forward.offset;
  pcumul = cumul + forward.offset - 1;
  lfrequency = length_bias_reestim->frequency + forward.offset;
  pmass = mass + forward.offset;
  sum = 0.;

  for (i = forward.offset;i < forward.nb_value;i++) {
    sum += *pfrequency++ / (1. - *pcumul++);
    *lfrequency++ += *pmass++ * sum;
  }
  for (i = forward.nb_value;i < nb_value;i++) {
    *lfrequency++ += *pmass++ * sum;
  }

  if (no_event) {
    norm = new double[no_event->nb_value];

    for (i = no_event->offset;i < no_event->nb_value;i++) {
      if (no_event->frequency[i] > 0) {
        pmass = mass + i + 1;
        norm[i] = 0.;
        for (j = i + 1;j < nb_value;j++) {
          norm[i] += (j - i) * *pmass++;
        }
      }
    }

    lfrequency = length_bias_reestim->frequency + no_event->offset + 1;
    pmass = mass + no_event->offset + 1;
    for (i = no_event->offset + 1;i < nb_value;i++) {
      pfrequency = no_event->frequency + no_event->offset;
      sum = 0.;
      for (j = no_event->offset;j < MIN(i , no_event->nb_value);j++) {
        if ((*pfrequency > 0) && (norm[j] > 0.)) {
          sum += *pfrequency * (i - j) / norm[j];
        }
        pfrequency++;
      }

      *lfrequency++ += *pmass++ * sum;
    }

    delete [] norm;
  }

  reestim_offset = 1;
  reestim_nb_value = inter_event_reestim->alloc_nb_value;

  ifrequency = inter_event_reestim->frequency + inter_event_reestim->alloc_nb_value;
  lfrequency = length_bias_reestim->frequency + inter_event_reestim->alloc_nb_value;
  while ((*--ifrequency == 0) && (*--lfrequency == 0) && (reestim_nb_value > 2)) {
    reestim_nb_value--;
  }
  inter_event_reestim->nb_value = reestim_nb_value;
  length_bias_reestim->nb_value = reestim_nb_value;

  ifrequency = inter_event_reestim->frequency + reestim_offset;
  lfrequency = length_bias_reestim->frequency + reestim_offset;
  while ((*ifrequency++ == 0) && (*lfrequency++ == 0) && (reestim_offset < reestim_nb_value - 1)) {
    reestim_offset++;
  }
  inter_event_reestim->offset = reestim_offset;
  length_bias_reestim->offset = reestim_offset;

  inter_event_reestim->nb_element_computation();
  length_bias_reestim->nb_element_computation();

# ifdef DEBUG
  if ((iter < 10) || ((iter < 100) && (iter % 10 == 0)) ||
      ((iter < 1000) && (iter % 100 == 0)) || (iter % 1000 == 0)) {
    inter_event_reestim->max_computation();
    inter_event_reestim->mean_computation();
    inter_event_reestim->variance_computation();

    length_bias_reestim->max_computation();
    length_bias_reestim->mean_computation();
    length_bias_reestim->variance_computation();

    cout << "\nquantites de reestimation loi inter_evenement :" << *inter_event_reestim << endl;
    cout << "\nquantites de reestimation loi biaisee par la longueur :" << *length_bias_reestim << endl;
  }
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul de la moyenne d'une loi par bissection d'intervalle.
 *
 *  arguments : pointeurs sur les quantites de reestimation de la loi et
 *              de la loi biaisee par la longueur.
 *
 *--------------------------------------------------------------*/

double interval_bisection(Reestimation<double> *distribution_reestim ,
                          Reestimation<double> *length_bias_reestim)

{
  register int i;
  double ratio , inf_ratio , sup_ratio , mean , inf_mean , sup_mean , *dfrequency , *lfrequency;

# ifdef DEBUG
  int iter = 0;
# endif


  // initialisations : calculs des 2 premieres valeurs

  dfrequency = distribution_reestim->frequency + distribution_reestim->offset;
  lfrequency = length_bias_reestim->frequency + distribution_reestim->offset;
  inf_ratio = 0.;
  sup_ratio = 0.;
  inf_mean = distribution_reestim->offset;
  sup_mean = distribution_reestim->nb_value - 1;

  for (i = distribution_reestim->offset;i < distribution_reestim->nb_value;i++) {
    inf_ratio += (*dfrequency + *lfrequency) * i /
                 (distribution_reestim->nb_element * sup_mean + length_bias_reestim->nb_element * i);
    sup_ratio += (*dfrequency++ + *lfrequency++) * i /
                 (distribution_reestim->nb_element * inf_mean + length_bias_reestim->nb_element * i);
  }

  do {
    dfrequency = distribution_reestim->frequency + distribution_reestim->offset;
    lfrequency = length_bias_reestim->frequency + distribution_reestim->offset;
    ratio = 0.;
    mean = (inf_mean + sup_mean) / 2.;

    for (i = distribution_reestim->offset;i < distribution_reestim->nb_value;i++) {
      ratio += (*dfrequency++ + *lfrequency++) * i /
               (distribution_reestim->nb_element * mean + length_bias_reestim->nb_element * i);
    }

#   ifdef DEBUG
    cout << STAT_label[STATL_ITERATION] << " " << iter++ << ": " << mean << " " << ratio << endl;
#   endif

    if (ratio < 1.) {
      inf_ratio = ratio;
      sup_mean = mean;
    }
    else {
      sup_ratio = ratio;
      inf_mean = mean;
    }
  }
  while (sup_ratio - inf_ratio > BISECTION_RATIO_THRESHOLD);

  mean = (inf_mean + sup_mean) / 2.;

# ifdef DEBUG
  cout << STAT_label[STATL_MEAN] << ": " << mean << " " << ratio << endl;
# endif

  return mean;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un processus de renouvellement en equilibre
 *  par l'algorithme EM a partir de donnees d'intervalles de temps.
 *
 *  arguments : reference sur un objet StatError, stream, lois empiriques des intervalles
 *              de temps censures a gauche, a droite et de la longueur
 *              de la periode d'observation dans le cas 0 evenement,
 *              reference sur la loi inter-evenement initiale, type d'estimateur
 *              (vraisemblance ou vraisemblance penalisee), nombre d'iterations,
 *              methode de calcul de la moyenne de la loi inter-evenement,
 *              poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi), moyenne de la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::estimation(StatError &error , ostream &os ,
                                                           const FrequencyDistribution &backward ,
                                                           const FrequencyDistribution &forward ,
                                                           const FrequencyDistribution *no_event ,
                                                           const DiscreteParametric &iinter_event ,
                                                           int estimator , int nb_iter ,
                                                           int mean_computation_method , double weight ,
                                                           int penalty_type , int outside ,
                                                           double iinter_event_mean) const

{
  using namespace sequence_analysis;

  bool status = true;
  register int i;
  int inb_value , max_nb_value;
  double likelihood , previous_likelihood , inter_event_mean , *penalty;
  DiscreteParametricModel *inter_event;
  Forward *forward_dist;
  Reestimation<double> *inter_event_reestim , *length_bias_reestim;
  FrequencyDistribution *backward_forward;
  const FrequencyDistribution *phisto[2];


  inter_event = NULL;
  error.init();

  if (nb_element < NB_COMPLETE_INTERVAL) {
    status = false;
    error.update(SEQ_error[SEQR_NB_COMPLETE_INTERVAL_TOO_SMALL]);
  }

  if (offset == 0) {
    status = false;
    error.update(SEQ_error[SEQR_COMPLETE_MIN_VALUE]);
  }
  if (forward.offset == 0) {
    status = false;
    error.update(SEQ_error[SEQR_FORWARD_MIN_VALUE]);
  }
  if ((no_event) && (no_event->offset == 0)) {
    status = false;
    error.update(SEQ_error[SEQR_NO_EVENT_MIN_VALUE]);
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if ((weight != D_DEFAULT) && (weight <= 0.)) {
    status = false;
    error.update(STAT_error[STATR_PENALTY_WEIGHT]);
  }

  if ((mean_computation_method == ESTIMATED) && (iinter_event_mean == D_DEFAULT)) {
    status = false;
    error.update(SEQ_error[SEQR_MEAN_COMPUTATION_METHOD]);
  }

  inb_value = nb_value;
  if (backward.nb_value + 1 > inb_value) {
    inb_value = backward.nb_value + 1;
  }
  if (forward.nb_value > inb_value) {
    inb_value = forward.nb_value;
  }

  if ((no_event) && (no_event->nb_value > inb_value)) {
    max_nb_value = no_event->nb_value;
  }
  else {
    max_nb_value = inb_value;
  }

  if ((iinter_event.offset > offset) || (iinter_event.nb_value < max_nb_value)) {
    status = false;
    error.update(SEQ_error[SEQR_INTER_EVENT_SUPPORT]);
  }

  if (status) {
    phisto[0] = new FrequencyDistribution(backward , 's' , 1);
    phisto[1] = &forward;
    backward_forward = new FrequencyDistribution(2 , phisto);
    delete phisto[0];

#   ifdef MESSAGE
    {
      int max_nb_element , width[2];
      long old_adjust;


      old_adjust = os.setf(ios::right , ios::adjustfield);

      width[0] = column_width(max_nb_value - 1);

      max_nb_element = nb_element;
      if (backward_forward->nb_element > max_nb_element) {
        max_nb_element = backward_forward->nb_element;
      }
      if ((no_event) && (no_event->nb_element > max_nb_element)) {
        max_nb_element = no_event->nb_element;
      }
      width[1] = column_width(max_nb_element) + ASCII_SPACE;

      os << "\n   | " << SEQ_label[SEQL_OBSERVATION_INTER_EVENT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
         << " | " << SEQ_label[SEQL_BACKWARD] << "/" << SEQ_label[SEQL_FORWARD]
         << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      if (no_event) {
        os << " | no-event " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      os << endl;

      for (i = 0;i < max_nb_value;i++) {
        os << setw(width[0]) << i;

        if (i < nb_value) {
          os << setw(width[1]) << frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }

        if (i < backward_forward->nb_value) {
          os << setw(width[1]) << backward_forward->frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }

        if (no_event) {
          if (i < no_event->nb_value) {
            os << setw(width[1]) << no_event->frequency[i];
          }
          else {
            os << setw(width[1]) << " ";
          }
        }

        os << "    |  ";
        if (i < backward.nb_value) {
          os << setw(width[1]) << backward.frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }

        if (i < forward.nb_value) {
          os << setw(width[1]) << forward.frequency[i];
        }
        else {
          os << setw(width[1]) << " ";
        }

        os << endl;
      }
      os << endl;

      os << setw(width[0]) << " "
         << setw(width[1]) << nb_element
         << setw(width[1]) << backward_forward->nb_element;
      if (no_event) {
        os << setw(width[1]) << no_event->nb_element;
      }
      os << "    |  " << setw(width[1]) << backward.nb_element
         << setw(width[1]) << forward.nb_element << "\n" << endl;

      os.setf((FMTFLAGS)old_adjust , ios::adjustfield);
    }
#   endif

    // creation de la loi inter-evenement

    inter_event = new DiscreteParametricModel(iinter_event , this);
    inter_event->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
    forward_dist = new Forward(*inter_event);

    if (estimator == PENALIZED_LIKELIHOOD) {
      penalty = new double[inter_event->nb_value];

      if (weight == D_DEFAULT) {
        if (penalty_type != ENTROPY) {
          weight = RENEWAL_DIFFERENCE_WEIGHT;
        }
        else {
          weight = RENEWAL_ENTROPY_WEIGHT;
        }
      }

      if (no_event) {
        weight *= (nb_element + backward.nb_element + forward.nb_element + no_event->nb_element);
      }
      else {
        weight *= (nb_element + backward.nb_element + forward.nb_element);
      }
    }

    inter_event_reestim = new Reestimation<double>(inter_event->nb_value);
    length_bias_reestim = new Reestimation<double>(inter_event->nb_value);

    likelihood = D_INF;
    i = 0;

    do {
      i++;

      inter_event->expectation_step(*this , backward , forward , no_event ,
                                    inter_event_reestim , length_bias_reestim , i);

      switch (estimator) {

      case LIKELIHOOD : {
        switch (mean_computation_method) {
        case ESTIMATED :
          inter_event_mean = iinter_event_mean;
          break;
        case COMPUTED :
          inter_event_mean = interval_bisection(inter_event_reestim , length_bias_reestim);
          break;
        case ONE_STEP_LATE :
          inter_event_mean = inter_event->mean;
          break;
        }

        inter_event_reestim->equilibrium_process_estimation(length_bias_reestim , inter_event ,
                                                            inter_event_mean);
        break;
      }

      case PENALIZED_LIKELIHOOD : {
        switch (mean_computation_method) {
        case ESTIMATED :
          inter_event_mean = iinter_event_mean;
          break;
        case ONE_STEP_LATE :
          inter_event_mean = inter_event->mean;
          break;
        }

        inter_event_reestim->penalized_likelihood_equilibrium_process_estimation(length_bias_reestim ,
                                                                                 inter_event , inter_event_mean ,
                                                                                 weight , penalty_type , penalty ,
                                                                                 outside);
        break;
      }
      }

      forward_dist->computation(*inter_event);
      previous_likelihood = likelihood;
      likelihood = inter_event->renewal_likelihood_computation(*forward_dist , *this , backward ,
                                                               forward , no_event);

#     ifdef MESSAGE
      if ((i < 10) || ((i < 100) && (i % 10 == 0)) || ((i < 1000) && (i % 100 == 0)) || (i % 1000 == 0)) {
        os << STAT_label[STATL_ITERATION] << " " << i << "   "
           << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
           << STAT_label[STATL_SMOOTHNESS] << ": " << inter_event->second_difference_norm_computation();
        if (estimator == PENALIZED_LIKELIHOOD) {
          os << "   cumul: " << inter_event->cumul[inter_event->nb_value - 1];
        }

        if ((no_event) && (no_event->offset + 1 == no_event->nb_value) && (backward_forward->nb_value > nb_value) &&
            ((forward.nb_element + no_event->nb_element) * (1. - inter_event->cumul[inb_value - 2]) > 0.)) {
          if (mean_computation_method == ESTIMATED) {
            inter_event_mean = iinter_event_mean;
          }
          else {
            inter_event_mean = inter_event->mean;
          }

          os << "   smaller upper bound: "
             << inb_value - 1 + (no_event->nb_element * inter_event_mean) /
                                ((forward.nb_element + no_event->nb_element) * (1. - inter_event->cumul[inb_value - 2]));
        }

/*        if (backward_forward->nb_value > nb_value) {
          register int j;
          double term = forward.nb_element * (backward_forward->nb_value - 1) /
                        inter_event->mean + nb_element + backward.nb_element;
          if (no_event) {
            term += no_event->nb_element * (backward_forward->nb_value - 1) / inter_event->mean;
          }
          for (j = backward_forward->offset;j < backward_forward->nb_value;j++) {
            term -= backward_forward->frequency[j] / (1. - inter_event->cumul[j - 1]);
          }

          os << " |   " << term;
        } */

        os << endl;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (i < RENEWAL_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > RENEWAL_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (i < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
         << STAT_label[STATL_SMOOTHNESS] << ": " << inter_event->second_difference_norm_computation();
      if (estimator == PENALIZED_LIKELIHOOD) {
        os << "   cumul: " << inter_event->cumul[inter_event->nb_value - 1];
      }

      if ((no_event) && (no_event->offset + 1 == no_event->nb_value) && (backward_forward->nb_value > nb_value) &&
          ((forward.nb_element + no_event->nb_element) * (1. - inter_event->cumul[inb_value - 2]) > 0.)) {
        if (mean_computation_method == ESTIMATED) {
          inter_event_mean = iinter_event_mean;
        }
        else {
          inter_event_mean = inter_event->mean;
        }

        os << "   smaller upper bound: "
           << inb_value - 1 + (no_event->nb_element * inter_event_mean) /
                              ((forward.nb_element + no_event->nb_element) * (1. - inter_event->cumul[inb_value - 2]));
      }
      os << endl;
#     endif

    }

    else {
      delete inter_event;
      inter_event = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    delete backward_forward;

    delete forward_dist;
    if (estimator == PENALIZED_LIKELIHOOD) {
      delete [] penalty;
    }

    delete inter_event_reestim;
    delete length_bias_reestim;
  }

  return inter_event;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un processus de renouvellement en equilibre
 *  par l'algorithme EM a partir de donnees d'intervalles de temps.
 *
 *  arguments : reference sur un objet StatError, stream, lois empiriques des intervalles
 *              de temps censures a gauche, a droite et de la longueur de la periode
 *              d'observation dans le cas 0 evenement, type d'estimateur
 *              (vraisemblance ou vraisemblance penalisee), nombre d'iterations,
 *              methode de calcul de la moyenne de la loi inter-evenement,
 *              poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* FrequencyDistribution::estimation(StatError &error , ostream &os ,
                                                           const FrequencyDistribution &backward ,
                                                           const FrequencyDistribution &forward ,
                                                           const FrequencyDistribution *no_event ,
                                                           int estimator , int nb_iter ,
                                                           int mean_computation_method , double weight ,
                                                           int penalty_type , int outside) const

{
  using namespace sequence_analysis;

  register int i;
  int nb_histo , *pfrequency;
  double *pmass;
  DiscreteParametric *iinter_event;
  DiscreteParametricModel *inter_event;
  FrequencyDistribution *interval;
  const FrequencyDistribution *phisto[4];


  nb_histo = 3;
  phisto[0] = this;
  phisto[1] = new FrequencyDistribution(backward , 's' , 1);
  phisto[2] = &forward;
  if (no_event) {
    nb_histo++;
    phisto[3] = new FrequencyDistribution(*no_event , 's' , 1);
  }

  interval = new FrequencyDistribution(nb_histo , phisto);
  delete phisto[1];
  if (no_event) {
    delete phisto[3];
  }

  iinter_event = new DiscreteParametric((int)(interval->nb_value * MAX_VALUE_COEFF));

  iinter_event->offset = interval->offset;

  pmass = iinter_event->mass;
  for (i = 0;i < interval->offset;i++) {
    *pmass++ = 0.;
  }

  pfrequency = interval->frequency + interval->offset;
  for (i = interval->offset;i < interval->nb_value;i++) {
    *pmass++ = (double)*pfrequency++ / (double)(interval->nb_element + 1);
  }

  for (i = interval->nb_value;i < iinter_event->nb_value - 1;i++) {
    *pmass++ = 0.;
  }
  *pmass = 1. / (double)(interval->nb_element + 1);

  iinter_event->cumul_computation();

  iinter_event->max = (double)max / (double)(interval->nb_element + 1);
  iinter_event->mean_computation();
  iinter_event->variance_computation();

  delete interval;

# ifdef DEBUG
  iinter_event->ascii_print(cout);
# endif

  inter_event = estimation(error , os , backward , forward , no_event ,
                           *iinter_event , estimator , nb_iter , mean_computation_method ,
                           weight , penalty_type , outside);
  delete iinter_event;

  return inter_event;
}


};  // namespace stat_tool



using namespace stat_tool;


namespace sequence_analysis {



/*--------------------------------------------------------------*
 *
 *  Calcul de la quantite d'information d'un objet TimeEvents.
 *
 *--------------------------------------------------------------*/

double TimeEvents::information_computation() const

{
  register int i;
  double information , buff;


  information = htime->information_computation();

  if (information != D_INF) {
    for (i = htime->offset;i < htime->nb_value;i++) {
      if (htime->frequency[i] > 0) {
        buff = hnb_event[i]->information_computation();

        if (buff != D_INF) {
          information += buff;
        }
        else {
          information = D_INF;
          break;
        }
      }
    }
  }

  return information;
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la  vraisemblance d'echantillons {temps, nombre d'evenements}
 *  pour une melange de lois du nombre d'evenements donne.
 *
 *  argument : reference sur un objet TimeEvents.
 *
 *--------------------------------------------------------------*/

double Renewal::likelihood_computation(const TimeEvents &timev) const

{
  register int i;
  double likelihood , buff;


  likelihood = time->likelihood_computation(*(timev.htime));

  if (likelihood != D_INF) {
    for (i = timev.htime->offset;i < timev.htime->nb_value;i++) {
      if (timev.htime->frequency[i] > 0) {
        buff = nb_event[i]->likelihood_computation(*(timev.hnb_event[i]));

        if (buff != D_INF) {
          likelihood += buff;
        }
        else {
          likelihood = D_INF;
          break;
        }
      }
    }
  }

  return likelihood;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation correspondant a la loi inter-evenement
 *  (estimateur EM d'un processus de renouvellement ordinaire a partir de donnees de comptage).
 *
 *  arguments : reference sur les echantillons {temps, nombre d'evenements},
 *              pointeur sur les quantites de reestimation.
 *
 *--------------------------------------------------------------*/

void Renewal::expectation_step(const TimeEvents &timev ,
                               Reestimation<double> *inter_event_reestim) const

{
  register int i , j;
  int min_time , max_time , *ptime , *pnb_event , *pfrequency;
  double num , denom , *ifrequency , *pmass;


  // initialisation

  ifrequency = inter_event_reestim->frequency;
  for (i = 0;i < inter_event_reestim->alloc_nb_value;i++) {
    *ifrequency++ = 0.;
  }

  ptime = timev.time;
  pnb_event = timev.nb_event;
  pfrequency = timev.frequency;

  for (i = 0;i < timev.nb_class;i++) {

    // cas pas d'evenement

    if (*pnb_event == 0) {
      min_time = MAX(inter_event->offset , *ptime + 1);

      if (min_time < inter_event->nb_value) {
        denom = nb_event[*ptime]->mass[*pnb_event];

        if (denom > 0.) {
          ifrequency = inter_event_reestim->frequency + min_time;
          pmass = inter_event->mass + min_time;
          for (j = min_time;j < inter_event->nb_value;j++) {
            *ifrequency++ += *pfrequency * *pmass++ / denom;
          }
        }
      }
    }

    // cas au moins 1 evenement

    else {
      denom = 0.;
      if (*pnb_event < nb_event[*ptime]->nb_value) {
        denom = nb_event[*ptime]->mass[*pnb_event];
      }

      if (denom > 0.) {
        max_time = MIN(nevent_time[*pnb_event]->nb_value - 1 , *ptime);

        ifrequency = inter_event_reestim->frequency + inter_event->offset;
        pmass = inter_event->mass + inter_event->offset;

        for (j = inter_event->offset;j < inter_event->nb_value;j++) {
          num = 0.;

          // cas 1 evenement : intervalle complet

          if (*pnb_event == 1) {
            if ((*ptime - j >= 0) && (*ptime - j < nevent_time[*pnb_event]->nb_value)) {
              num = 1. - nevent_time[*pnb_event]->cumul[*ptime - j];
            }
          }

          // cas plus de 1 evenement : intervalles complets

          else {
            if (*ptime - j >= nevent_time[*pnb_event - 1]->offset) {
              if (*ptime - j < nevent_time[*pnb_event - 1]->nb_value) {
                num = *pnb_event * (nevent_time[*pnb_event - 1]->cumul[*ptime - j] -
                                    nevent_time[*pnb_event]->cumul[*ptime - j]);
              }
              else if (*ptime - j < nevent_time[*pnb_event]->nb_value) {
                num = *pnb_event * (1. - nevent_time[*pnb_event]->cumul[*ptime - j]);
              }
            }
          }

          // intervalle tronque

          if ((max_time > *ptime - j) && (max_time >= nevent_time[*pnb_event]->offset)) {
            num += nevent_time[*pnb_event]->cumul[max_time];
            if (*ptime - j >= nevent_time[*pnb_event]->offset) {
              num -= nevent_time[*pnb_event]->cumul[*ptime - j];
            }
          }

          *ifrequency++ += *pfrequency * *pmass++ * num / denom;
        }
      }
    }

    ptime++;
    pnb_event++;
    pfrequency++;
  }

  inter_event_reestim->nb_value_computation();
  inter_event_reestim->offset_computation();
  inter_event_reestim->nb_element_computation();
  inter_event_reestim->max_computation();
  inter_event_reestim->mean_computation();
  inter_event_reestim->variance_computation();

# ifdef DEBUG
  cout << "\n" << (timev.mixture->mean + 1) * timev.nb_element << " | "
       << " " << inter_event_reestim->nb_element << endl;

  cout << "\nquantites de reestimation loi inter_evenement :" << *inter_event_reestim << endl;
# endif

}


/*--------------------------------------------------------------*
 *
 *  Calcul des quantites de reestimation correspondant a la loi inter-evenement
 *  (estimateur EM d'un processus de renouvellement en equilibre a partir de donnees de comptage).
 *
 *  arguments : reference sur les echantillons {temps, nombre d'evenements},
 *              pointeurs sur les quantites de reestimation de la loi inter-evenement et
 *              de la loi biaisee par la longueur, type d'estimateur,
 *              combinaison ou non des quantites de reestimation,
 *              methode de calcul de la moyenne de la loi inter-evenement.
 *
 *--------------------------------------------------------------*/

void Renewal::expectation_step(const TimeEvents &timev ,
                               Reestimation<double> *inter_event_reestim ,
                               Reestimation<double> *length_bias_reestim , int estimator ,
                               bool combination , int mean_computation_method) const

{
  register int i , j;
  int max_time , offset , nb_value , *ptime , *pnb_event , *pfrequency;
  double complete_num , censored_num , denom , inter_event_mean , *ifrequency ,
         *lfrequency , *pmass;


  // initialisations

  ifrequency = inter_event_reestim->frequency;
  for (i = 0;i < inter_event_reestim->alloc_nb_value;i++) {
    *ifrequency++ = 0.;
  }

  if (estimator == COMPLETE_LIKELIHOOD) {
    lfrequency = length_bias_reestim->frequency;
    for (i = 0;i < length_bias_reestim->alloc_nb_value;i++) {
      *lfrequency++ = 0.;
    }
  }

  ptime = timev.time;
  pnb_event = timev.nb_event;
  pfrequency = timev.frequency;

  for (i = 0;i < timev.nb_class;i++) {

    // cas pas d'evenement

    if (*pnb_event == 0) {
      if ((estimator == COMPLETE_LIKELIHOOD) && (*ptime + 1 < inter_event->nb_value)) {
        denom = nb_event[*ptime]->mass[*pnb_event] * inter_event->mean;

        if (denom > 0.) {
          lfrequency = length_bias_reestim->frequency + *ptime + 1;
          pmass = inter_event->mass + *ptime + 1;
          for (j = *ptime + 1;j < inter_event->nb_value;j++) {
            *lfrequency++ += *pfrequency * (j - *ptime) * *pmass++ / denom;
          }
        }
      }
    }

    // cas au moins 1 evenement

    else {
      denom = 0.;
      if (*pnb_event < nb_event[*ptime]->nb_value) {
        denom = nb_event[*ptime]->mass[*pnb_event];
      }

      if (denom > 0.) {

        // intervalle tronque initial

/*        if (estimator == COMPLETE_LIKELIHOOD) {
          lfrequency = length_bias_reestim->frequency + inter_event->offset;
          pmass = inter_event->mass + inter_event->offset;
          num = 0.;

          for (j = 1;j < inter_event->nb_value;j++) {
            if (j <= *ptime) {

              // cas 1 evenement

              if (*pnb_event == 1) {
                if (*ptime - j < aux_nevent_time[*pnb_event]->nb_value) {
                  num += 1. - aux_nevent_time[*pnb_event]->cumul[*ptime - j];
                }
              }

              // cas plus de 1 evenement

              else {
                if (*ptime - j >= aux_nevent_time[*pnb_event - 1]->offset) {
                  if (*ptime - j < aux_nevent_time[*pnb_event - 1]->nb_value) {
                    num += aux_nevent_time[*pnb_event - 1]->cumul[*ptime - j] -
                           aux_nevent_time[*pnb_event]->cumul[*ptime - j];
                  }
                  else if (*ptime - j < aux_nevent_time[*pnb_event]->nb_value) {
                    num += 1. - aux_nevent_time[*pnb_event]->cumul[*ptime - j];
                  }
                }
              }
            }

            if (j >= inter_event->offset) {
              *lfrequency++ += *pfrequency * *pmass++ * num / (denom * inter_event->mean);
            }
          }
        } */

        max_time = MIN(nevent_time[*pnb_event]->nb_value - 1 , *ptime);

        ifrequency = inter_event_reestim->frequency + inter_event->offset;
        if (estimator == COMPLETE_LIKELIHOOD) {
          lfrequency = length_bias_reestim->frequency + inter_event->offset;
        }
        pmass = inter_event->mass + inter_event->offset;

        for (j = inter_event->offset;j < inter_event->nb_value;j++) {
          complete_num = 0.;

          // cas plus de 1 evenement : intervalles complets

          if (*pnb_event > 1) {
            if (*ptime - j >= nevent_time[*pnb_event - 1]->offset) {
              if (*ptime - j < nevent_time[*pnb_event - 1]->nb_value) {
                complete_num = (*pnb_event - 1) * (nevent_time[*pnb_event - 1]->cumul[*ptime - j] -
                                                   nevent_time[*pnb_event]->cumul[*ptime - j]);
              }
              else if (*ptime - j < nevent_time[*pnb_event]->nb_value) {
                complete_num = (*pnb_event - 1) * (1. - nevent_time[*pnb_event]->cumul[*ptime - j]);
              }
            }
          }

          // intervalles tronques initial et final

          censored_num = 0.;
          if ((max_time > *ptime - j) && (max_time >= nevent_time[*pnb_event]->offset)) {
            censored_num += nevent_time[*pnb_event]->cumul[max_time];
            if (*ptime - j >= nevent_time[*pnb_event]->offset) {
              censored_num -= nevent_time[*pnb_event]->cumul[*ptime - j];
            }
          }

          *ifrequency++ += *pfrequency * *pmass * (complete_num + censored_num) / denom;
          if (estimator == COMPLETE_LIKELIHOOD) {
            *lfrequency++ += *pfrequency * *pmass * censored_num / denom;
          }
          pmass++;
        }
      }
    }

    ptime++;
    pnb_event++;
    pfrequency++;
  }

  switch (estimator) {

  case PARTIAL_LIKELIHOOD : {
    inter_event_reestim->nb_value_computation();
    inter_event_reestim->offset_computation();
    inter_event_reestim->nb_element_computation();
    break;
  }

  case COMPLETE_LIKELIHOOD : {
    offset = 1;
    nb_value = inter_event_reestim->alloc_nb_value;

    ifrequency = inter_event_reestim->frequency + inter_event_reestim->alloc_nb_value;
    lfrequency = length_bias_reestim->frequency + inter_event_reestim->alloc_nb_value;
    while ((*--ifrequency == 0) && (*--lfrequency == 0) && (nb_value > 2)) {
      nb_value--;
    }
    inter_event_reestim->nb_value = nb_value;
    length_bias_reestim->nb_value = nb_value;

    ifrequency = inter_event_reestim->frequency + offset;
    lfrequency = length_bias_reestim->frequency + offset;
    while ((*ifrequency++ == 0) && (*lfrequency++ == 0) && (offset < nb_value - 1)) {
      offset++;
    }
    inter_event_reestim->offset = offset;
    length_bias_reestim->offset = offset;

    inter_event_reestim->nb_element_computation();
    length_bias_reestim->nb_element_computation();
    break;

#   ifdef DEBUG
    inter_event_reestim->max_computation();
    inter_event_reestim->mean_computation();
    inter_event_reestim->variance_computation();

    length_bias_reestim->max_computation();
    length_bias_reestim->mean_computation();
    length_bias_reestim->variance_computation();

    cout << "\n" << timev.nb_element << " | "  << length_bias_reestim->nb_element << " || "
         << timev.mixture->mean * timev.nb_element << " | " << " " << inter_event_reestim->nb_element << endl;

    cout << "\nquantites de reestimation loi inter_evenement :" << *inter_event_reestim << endl;
    cout << "\nquantites de reestimation loi biaisee par la longueur :" << *length_bias_reestim << endl;
#   endif

  }
  }

  if ((estimator == COMPLETE_LIKELIHOOD) && (combination)) {
    switch (mean_computation_method) {
    case ESTIMATED :
      inter_event_mean = timev.htime->mean / timev.mixture->mean;
      break;
    case COMPUTED :
      inter_event_mean = interval_bisection(inter_event_reestim , length_bias_reestim);
      break;
    case ONE_STEP_LATE :
      inter_event_mean = inter_event->mean;
      break;
    }

#   ifdef DEBUG
    if (mean_computation_method != ESTIMATED) {
      cout << SEQ_label[SEQL_INTER_EVENT] << " " << STAT_label[STATL_MEAN] << ": "
           << inter_event_mean << " (" << timev.htime->mean / timev.mixture->mean << ") | ";
    }
#   endif

    inter_event_reestim->equilibrium_process_combination(length_bias_reestim , inter_event_mean);
  }

  else {
    inter_event_reestim->max_computation();
    inter_event_reestim->mean_computation();
    inter_event_reestim->variance_computation();
  }

# ifdef DEBUG
  cout << "\nquantites de reestimation loi inter_evenement :" << *inter_event_reestim << endl;
# endif

}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un processus de renouvellement par l'algorithme EM
 *  a partir d'echantillons {temps d'observation, nombre d'evenements}.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              type de processus ('o' : ordinaire, 'e' : en equilibre),
 *              reference sur la loi inter-evenement initiale,
 *              type d'estimateur (vraisemblance, vraisemblance penalisee ou
 *              estimation d'une loi parametrique), nombre d'iterations,
 *              type d'estimateur dans le cas d'un processus de renouvellement en equilibre,
 *              methode de calcul de la moyenne de la loi inter-evenement,
 *              poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

Renewal* TimeEvents::estimation(StatError &error , ostream &os , char type ,
                                const DiscreteParametric &iinter_event , int estimator ,
                                int nb_iter , int equilibrium_estimator , int mean_computation_method ,
                                double weight , int penalty_type , int outside) const

{
  bool status = true;
  register int i;
  int nb_likelihood_decrease;
  double likelihood , previous_likelihood , information , hlikelihood , inter_event_mean , *penalty;
  DiscreteParametric *pinter_ev;
  Reestimation<double> *inter_event_reestim , *length_bias_reestim;
  FrequencyDistribution *hreestim;
  Renewal *renew;


  renew = NULL;
  error.init();

  if (mixture->nb_value <= 2) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_NB_EVENT_TOO_SMALL]);
  }
  if (mixture->mean < MIN_NB_EVENT) {
    status = false;
    error.update(SEQ_error[SEQR_NB_EVENT_TOO_SMALL]);
  }

  if (min_inter_event_computation() < MIN_INTER_EVENT) {
    status = false;
    error.update(SEQ_error[SEQR_TIME_UNIT]);
  }

  if ((nb_iter != I_DEFAULT) && (nb_iter < 1)) {
    status = false;
    error.update(STAT_error[STATR_NB_ITERATION]);
  }

  if ((weight != D_DEFAULT) && (weight <= 0.)) {
    status = false;
    error.update(STAT_error[STATR_PENALTY_WEIGHT]);
  }

  if (status) {
    information = information_computation();

    // creation du processus de renouvellement

    renew = new Renewal(type , *htime , iinter_event);
    renew->renewal_data = new RenewalData(*this , type);

    pinter_ev = renew->inter_event;

    if (estimator == PENALIZED_LIKELIHOOD) {
      penalty = new double[pinter_ev->alloc_nb_value];

      if (weight == D_DEFAULT) {
        if (penalty_type != ENTROPY) {
          weight = RENEWAL_DIFFERENCE_WEIGHT;
        }
        else {
          weight = RENEWAL_ENTROPY_WEIGHT;
        }
      }

      if (equilibrium_estimator == PARTIAL_LIKELIHOOD) {
        weight *= mixture->mean * nb_element;
      }
      else {
        weight *= (mixture->mean + 1) * nb_element;
      }
    }

    inter_event_reestim = new Reestimation<double>(pinter_ev->alloc_nb_value);

    if (type == 'e') {
      if (equilibrium_estimator == COMPLETE_LIKELIHOOD) {
        length_bias_reestim = new Reestimation<double>(pinter_ev->alloc_nb_value);
      }
    }

    renew->computation();

    renew->init(CATEGORICAL , I_DEFAULT , I_DEFAULT , D_DEFAULT , D_DEFAULT);
    likelihood = D_INF;
    i = 0;

    do {
      i++;

      // calcul des quantites de reestimation

      switch (type) {
      case 'o' :
        renew->expectation_step(*this , inter_event_reestim);
        break;
      case 'e' :
        renew->expectation_step(*this , inter_event_reestim ,
                                length_bias_reestim , equilibrium_estimator);
        break;
      }

      if (estimator != PENALIZED_LIKELIHOOD) {
        if ((type == 'o') || (equilibrium_estimator == PARTIAL_LIKELIHOOD)) {
          inter_event_reestim->distribution_estimation(pinter_ev);
        }

        else {
          switch (mean_computation_method) {
          case ESTIMATED :
            inter_event_mean = htime->mean / mixture->mean;
            break;
          case COMPUTED :
            inter_event_mean = interval_bisection(inter_event_reestim , length_bias_reestim);
            break;
          case ONE_STEP_LATE :
            inter_event_mean = pinter_ev->mean;
            break;
          }

          inter_event_reestim->equilibrium_process_estimation(length_bias_reestim , pinter_ev ,
                                                              inter_event_mean);
        }
      }

      else {
        if ((type == 'o') || (equilibrium_estimator == PARTIAL_LIKELIHOOD)) {
          inter_event_reestim->penalized_likelihood_estimation(pinter_ev , weight ,
                                                               penalty_type , penalty ,
                                                               outside);
        }

        else {
          switch (mean_computation_method) {
          case ESTIMATED :
            inter_event_mean = htime->mean / mixture->mean;
            break;
          case ONE_STEP_LATE :
            inter_event_mean = pinter_ev->mean;
            break;
          }

          inter_event_reestim->penalized_likelihood_equilibrium_process_estimation(length_bias_reestim ,
                                                                                   pinter_ev , inter_event_mean ,
                                                                                   weight , penalty_type ,
                                                                                   penalty , outside);
        }
      }

      // calcul du melange de lois du nombre d'evenements et
      // de la log-vraisemblance correspondante

      renew->computation();

      previous_likelihood = likelihood;
      likelihood = renew->likelihood_computation(*this);

#     ifdef MESSAGE
      if ((i < 10) || ((i < 100) && (i % 10 == 0)) || ((i < 1000) && (i % 100 == 0)) || (i % 1000 == 0)) {
        os << STAT_label[STATL_ITERATION] << " " << i << "   "
           << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
           << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << "   "
           << STAT_label[STATL_SMOOTHNESS] << ": "  << pinter_ev->second_difference_norm_computation();
        if (estimator == PENALIZED_LIKELIHOOD) {
          os << "   cumul: " << pinter_ev->cumul[pinter_ev->nb_value - 1];
        }
        os << endl;
      }
#     endif

    }
    while ((likelihood != D_INF) && (((nb_iter == I_DEFAULT) && (i < RENEWAL_NB_ITER) &&
             ((likelihood - previous_likelihood) / -likelihood > RENEWAL_LIKELIHOOD_DIFF)) ||
            ((nb_iter != I_DEFAULT) && (i < nb_iter))));

    if (likelihood != D_INF) {

#     ifdef MESSAGE
      os << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
         << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
         << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << "   "
         << STAT_label[STATL_SMOOTHNESS] << ": " << pinter_ev->second_difference_norm_computation();
      if (estimator == PENALIZED_LIKELIHOOD) {
        os << "   cumul: " << pinter_ev->cumul[pinter_ev->nb_value - 1];
      }
      os << endl;
#     endif

      if (estimator == PARAMETRIC_REGULARIZATION) {
        hreestim = new FrequencyDistribution(pinter_ev->alloc_nb_value);

        likelihood = D_INF;
        nb_likelihood_decrease = 0;

        i = 0;
        do {
          i++;

          // calcul des quantites de reestimation

          switch (type) {
          case 'o' :
            renew->expectation_step(*this , inter_event_reestim);
            break;
          case 'e' :
            renew->expectation_step(*this , inter_event_reestim ,
                                    length_bias_reestim , equilibrium_estimator ,
                                    true , mean_computation_method);
            break;
          }

          hreestim->update(inter_event_reestim , (int)(inter_event_reestim->nb_element *
                           MAX(sqrt(inter_event_reestim->variance) , 1.) * RENEWAL_COEFF));
          hlikelihood = hreestim->Reestimation<int>::type_parametric_estimation(pinter_ev , 1 , true ,
                                                                                RENEWAL_THRESHOLD);

          if (hlikelihood == D_INF) {
            likelihood = D_INF;
          }

          // calcul du melange de lois du nombre d'evenements et
          // de la log-vraisemblance correspondante

          else {
            renew->init(pinter_ev->ident , pinter_ev->inf_bound , pinter_ev->sup_bound ,
                        pinter_ev->parameter , pinter_ev->probability);
            renew->computation();

            previous_likelihood = likelihood;
            likelihood = renew->likelihood_computation(*this);

            if (likelihood < previous_likelihood) {
              nb_likelihood_decrease++;
            }
            else {
              nb_likelihood_decrease = 0;
            }

#           ifdef DEBUG
            if ((i < 10) || (i % 10 == 0)) {
              os << STAT_label[STATL_ITERATION] << " " << i << "   "
                 << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
                 << STAT_label[STATL_SMOOTHNESS] << ": " << pinter_ev->second_difference_norm_computation() << endl;
            }
#           endif

          }
        }
        while ((likelihood != D_INF) && (i < RENEWAL_NB_ITER) &&
               (((likelihood - previous_likelihood) / -likelihood > RENEWAL_LIKELIHOOD_DIFF) ||
                (hlikelihood == D_INF) || (nb_likelihood_decrease == 1)));

        delete hreestim;

#       ifdef MESSAGE
        if (likelihood != D_INF) {
          os << "\n" << i << " " << STAT_label[STATL_ITERATIONS] << "   "
             << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   "
             << STAT_label[STATL_SMOOTHNESS] << ": " << pinter_ev->second_difference_norm_computation() << endl;
        }
#       endif

      }
    }

    if (estimator == PENALIZED_LIKELIHOOD) {
      delete [] penalty;
    }

    delete inter_event_reestim;

    if (type == 'e') {
      if (equilibrium_estimator == COMPLETE_LIKELIHOOD) {
        delete length_bias_reestim;
      }
    }

    if (likelihood != D_INF) {

      // mise a jour du nombre de parametres inconnus

      pinter_ev->nb_parameter_update();
      for (i = renew->time->offset;i < renew->time->nb_value;i++) {
        if (renew->time->mass[i] > 0.) {
          renew->nb_event[i]->nb_parameter = pinter_ev->nb_parameter;
        }
      }
      renew->mixture->nb_parameter = pinter_ev->nb_parameter;
    }

    else {
      delete renew;
      renew = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }
  }

  return renew;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un processus de renouvellement par l'algorithme EM
 *  a partir d'echantillons {temps d'observation, nombre d'evenements}.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              type de processus ('o' : ordinaire, 'e' : en equilibre),
 *              type d'estimateur (vraisemblance, vraisemblance penalisee ou
 *              estimation d'une loi parametrique), nombre d'iterations,
 *              type d'estimateur dans le cas d'un processus de renouvellement en equilibre,
 *              methode de calcul de la moyenne de la loi inter-evenement,
 *              poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

Renewal* TimeEvents::estimation(StatError &error , ostream &os , char type ,
                                int estimator , int nb_iter , int equilibrium_estimator ,
                                int mean_computation_method , double weight ,
                                int penalty_type , int outside) const

{
  double proba;
  DiscreteParametric *iinter_event;
  Renewal *renew;


  proba = mixture->mean / htime->mean;
  if (proba > 1. - RENEWAL_INIT_PROBABILITY) {
    proba = 1. - RENEWAL_INIT_PROBABILITY;
  }
  else if (proba < RENEWAL_INIT_PROBABILITY) {
    proba = RENEWAL_INIT_PROBABILITY;
  }

  iinter_event = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                        proba , RENEWAL_THRESHOLD);

# ifdef DEBUG
  iinter_event->ascii_print(cout);
# endif

  renew = estimation(error , os , type , *iinter_event , estimator , nb_iter ,
                     equilibrium_estimator , mean_computation_method , weight ,
                     penalty_type , outside);
  delete iinter_event;

  return renew;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un processus de renouvellement en equilibre
 *  par l'algorithme EM a partir de donnees d'intervalles de temps.
 *
 *  arguments : reference sur un objet StatError, stream,
 *              reference sur la loi inter-evenement initiale, type d'estimateur
 *              (vraisemblance ou vraisemblance penalisee), nombre d'iterations,
 *              methode de calcul de la moyenne de la loi inter-evenement,
 *              poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

Renewal* RenewalData::estimation(StatError &error , ostream &os ,
                                 const DiscreteParametric &iinter_event ,
                                 int estimator , int nb_iter ,
                                 int mean_computation_method , double weight ,
                                 int penalty_type , int outside) const

{
  register int i , j;
  int *psequence;
  DiscreteParametricModel *inter_event;
  Renewal *renew;
  FrequencyDistribution *within_backward , *within_forward , *no_event;


  within_backward = new FrequencyDistribution(MIN(backward->nb_value , htime->nb_value - 1));
  within_forward = new FrequencyDistribution(MIN(forward->nb_value , htime->nb_value));
  no_event = new FrequencyDistribution(htime->nb_value);

  for (i = 0;i < nb_element;i++) {
    psequence = sequence[i] + length[i];
    for (j = 0;j < length[i];j++) {
      if (*--psequence == 1) {
        (within_backward->frequency[j])++;
        break;
      }
    }

    psequence = sequence[i];
    for (j = 0;j < length[i];j++) {
      if (*psequence++ == 1) {
        (within_forward->frequency[j + 1])++;
        break;
      }
    }

    if (j == length[i]) {
      (no_event->frequency[j])++;
    }
  }

  within_backward->nb_value_computation();
  within_backward->offset_computation();
  within_backward->nb_element_computation();

# ifdef DEBUG
  within_backward->max_computation();
  within_backward->mean_computation();
  within_backward->variance_computation();

  cout << *within_backward;
  cout << *backward;
# endif

  within_forward->nb_value_computation();
  within_forward->offset_computation();
  within_forward->nb_element_computation();

# ifdef DEBUG
  within_forward->max_computation();
  within_forward->mean_computation();
  within_forward->variance_computation();

  cout << *within_forward;
  cout << *forward;
# endif

  no_event->nb_value_computation();
  no_event->offset_computation();
  no_event->nb_element_computation();

# ifdef DEBUG
  no_event->max_computation();
  no_event->mean_computation();
  no_event->variance_computation();

  cout << *no_event;
# endif

  if (no_event->nb_element == 0) {
    delete no_event;
    no_event = NULL;
  }

  inter_event = within->estimation(error , os , *within_backward , *within_forward , no_event ,
                                   iinter_event , estimator , nb_iter , mean_computation_method ,
                                   weight , penalty_type , outside , htime->mean / mixture->mean);

  delete within_backward;
  delete within_forward;
  delete no_event;

  if (inter_event) {
    renew = new Renewal(*this , *((DiscreteParametric*)inter_event));
  }
  else {
    renew = NULL;
  }

  return renew;
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'un processus de renouvellement en equilibre
 *  par l'algorithme EM a partir de donnees d'intervalles de temps.
 *
 *  arguments : reference sur un objet StatError, stream, type d'estimateur
 *              (vraisemblance ou vraisemblance penalisee), nombre d'iterations,
 *              methode de calcul de la moyenne de la loi inter-evenement,
 *              poids de la penalisation, type de penalisation,
 *              type de gestion des effets de bord (zero a l'exterieur du support ou
 *              prolongation de la loi).
 *
 *--------------------------------------------------------------*/

Renewal* RenewalData::estimation(StatError &error , ostream &os , int estimator ,
                                 int nb_iter , int mean_computation_method , double weight ,
                                 int penalty_type , int outside) const

{
  double proba;
  DiscreteParametric *iinter_event;
  Renewal *renew;


  proba = mixture->mean / htime->mean;
  if (proba > 1. - RENEWAL_INIT_PROBABILITY) {
    proba = 1. - RENEWAL_INIT_PROBABILITY;
  }
  else if (proba < RENEWAL_INIT_PROBABILITY) {
    proba = RENEWAL_INIT_PROBABILITY;
  }

  iinter_event = new DiscreteParametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                        proba , RENEWAL_THRESHOLD);

# ifdef DEBUG
  iinter_event->ascii_print(cout);
# endif

  renew = estimation(error , os , *iinter_event , estimator , nb_iter ,
                     mean_computation_method , weight , penalty_type , outside);
  delete iinter_event;

  return renew;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un processus de renouvellement.
 *
 *  arguments : reference sur un objet StatError,
 *              type de simulation ('o' : ordinaire, 'e' : en equilibre),
 *              loi empirique du temps d'observation.
 *
 *--------------------------------------------------------------*/

RenewalData* Renewal::simulation(StatError &error , char itype ,
                                 const FrequencyDistribution &ihtime) const

{
  bool status = true , compute;
  register int i , j , k , m;
  int offset , time_interval , cumul_time , *ptime , *pnb_event , *psequence;
  Distribution *dtime;
  Renewal *renew;
  RenewalData *timev;


  timev = NULL;
  error.init();

  if ((ihtime.nb_element < 1) || (ihtime.nb_element > RENEWAL_NB_ELEMENT)) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }
  if (ihtime.offset < MAX(inter_event->offset , 2)) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_OBSERVATION_TIME]);
  }
  if (ihtime.nb_value - 1 > MAX_TIME) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
  }

  if (status) {
    dtime = new Distribution(ihtime);

    if ((itype != type) || (*dtime != *time)) {
      compute = true;
    }
    else {
      compute = false;
    }

    // creation d'un objet RenewalData

    timev = new RenewalData(itype , *this);

    timev->renewal = new Renewal(*this , false);
    renew = timev->renewal;

    timev->length = new int[ihtime.nb_element];
    timev->sequence = new int*[ihtime.nb_element];

    switch (itype) {
    case 'o' :
      offset = 0;
      break;
    case 'e' :
      offset = 1;
      break;
    }

    // 1er au n-eme evenement

    ptime = new int[ihtime.nb_element];
    pnb_event = new int[ihtime.nb_element];

    i = 0;
    for (j = ihtime.offset;j < ihtime.nb_value;j++) {
      for (k = 0;k < ihtime.frequency[j];k++) {

        // temps avant le 1er evenement (processus de renouvellement en equilibre)

        if (itype == 'e') {
          if (i == 0) {
            cumul_time = renew->forward->simulation();
          }
          else {
            cumul_time -= *(ptime - 1);
          }
          (timev->forward->frequency[cumul_time])++;

          time_interval = cumul_time;
        }

        // temps d'observation

        *ptime = j;

        timev->length[i] = *ptime + 1 - offset;
        timev->sequence[i] = new int[timev->length[i]];

        // temps avant le 1er evenement (processus de renouvellement ordinaire)

        if (itype == 'o') {
          time_interval = renew->inter_event->simulation();
          cumul_time = time_interval;
          (timev->inter_event->frequency[time_interval])++;
          if (time_interval <= *ptime) {
            (timev->within->frequency[time_interval])++;
          }
          else {
            (timev->length_bias->frequency[time_interval])++;
          }
        }

        psequence = timev->sequence[i] - 1;
        if (itype == 'o') {
          *++psequence = 1;
        }
        for (m = 1;m < MIN(time_interval , *ptime + 1);m++) {
          *++psequence = 0;
        }
        if (time_interval <= *ptime) {
          *++psequence = 1;
        }

        *pnb_event = 0;
        while (cumul_time <= *ptime) {
          (*pnb_event)++;
          time_interval = renew->inter_event->simulation();
          cumul_time += time_interval;

          for (m = cumul_time - time_interval + 1;m < MIN(cumul_time , *ptime + 1);m++) {
            *++psequence = 0;
          }

          (timev->inter_event->frequency[time_interval])++;
          if (cumul_time <= *ptime) {
            (timev->within->frequency[time_interval])++;
            *++psequence = 1;
          }
          else {
            (timev->length_bias->frequency[time_interval])++;
          }
        }

        (timev->backward->frequency[*ptime - (cumul_time - time_interval)])++;
        if (itype == 'o') {
          (timev->forward->frequency[cumul_time - *ptime])++;
        }

        ptime++;
        pnb_event++;
        i++;
      }
    }
    ptime -= ihtime.nb_element;
    pnb_event -= ihtime.nb_element;

    // construction des echantillons {temps, nombre d'evenements, frequence} et
    // des lois empiriques du temps d'observation et du nombre d'evenements

    timev->build(ihtime.nb_element , ptime , pnb_event);
    delete [] ptime;
    delete [] pnb_event;

    // extraction des caracteristiques des lois empiriques des intervalles de temps entre 2 evenements,
    // des intervalles de temps entre 2 evenements a l'interieur de la periode d'observation,
    // des intervalles de temps recouvrant une date d'observation,
    // des intervalles de temps apres le dernier evenement, des intervalles de temps residuel

    timev->inter_event->nb_value_computation();
    timev->inter_event->offset_computation();
    timev->inter_event->nb_element_computation();
    timev->inter_event->max_computation();
    timev->inter_event->mean_computation();
    timev->inter_event->variance_computation();

    timev->within->nb_value_computation();
    timev->within->offset_computation();
    timev->within->nb_element_computation();
    timev->within->max_computation();
    timev->within->mean_computation();
    timev->within->variance_computation();

    timev->length_bias->nb_value_computation();
    timev->length_bias->offset_computation();
    timev->length_bias->nb_element_computation();
    timev->length_bias->max_computation();
    timev->length_bias->mean_computation();
    timev->length_bias->variance_computation();

    timev->backward->nb_value_computation();
    timev->backward->offset_computation();
    timev->backward->nb_element_computation();
    timev->backward->max_computation();
    timev->backward->mean_computation();
    timev->backward->variance_computation();

    timev->forward->nb_value_computation();
    timev->forward->offset_computation();
    timev->forward->nb_element_computation();
    timev->forward->max_computation();
    timev->forward->mean_computation();
    timev->forward->variance_computation();

    // extraction des probabilites de non-evenement/evenement fonction du temps

    timev->build_index_event(offset);

    if (compute) {
      renew->computation(false , itype , dtime);
    }
    delete dtime;
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un processus de renouvellement.
 *
 *  arguments : reference sur un objet StatError,
 *              type de simulation ('o' : ordinaire, 'e' : en equilibre),
 *              nombre d'echantillons, temps d'observation.
 *
 *--------------------------------------------------------------*/

RenewalData* Renewal::simulation(StatError &error , char itype ,
                                 int nb_element , int itime) const

{
  bool status = true;
  RenewalData *timev;


  timev = NULL;
  error.init();

  if ((nb_element < 1) || (nb_element > RENEWAL_NB_ELEMENT)) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }
  if (itime < MAX(inter_event->offset , 2)) {
    status = false;
    error.update(SEQ_error[SEQR_SHORT_OBSERVATION_TIME]);
  }
  if (itime > MAX_TIME) {
    status = false;
    error.update(SEQ_error[SEQR_LONG_OBSERVATION_TIME]);
  }

  if (status) {
    FrequencyDistribution htime(itime + 1);

    htime.nb_element = nb_element;
    htime.offset = itime;
    htime.max = nb_element;
    htime.mean = itime;
    htime.variance = 0.;
    htime.frequency[itime] = nb_element;

    timev = simulation(error , itype , htime);
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un processus de renouvellement.
 *
 *  arguments : reference sur un objet StatError,
 *              type de simulation ('o' : ordinaire, 'e' : en equilibre),
 *              nombre d'echantillons, reference sur un objet TimeEvents.
 *
 *--------------------------------------------------------------*/

RenewalData* Renewal::simulation(StatError &error , char itype ,
                                 int nb_element , const TimeEvents &itimev) const

{
  FrequencyDistribution *htime;
  RenewalData *timev;


  error.init();

  if ((nb_element < 1) || (nb_element > RENEWAL_NB_ELEMENT)) {
    timev = NULL;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  else {
    htime = itimev.htime->frequency_scale(nb_element);

    timev = simulation(error , itype , *htime);
    delete htime;
  }

  return timev;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe RenewalIterator.
 *
 *  arguments : pointeur sur un objet Renewal, longueur de la sequence.
 *
 *--------------------------------------------------------------*/

RenewalIterator::RenewalIterator(Renewal *irenewal , int ilength)

{
  renewal = irenewal;
  (renewal->nb_iterator)++;

  interval = 0;
  counter = 0;

  length = ilength;
  sequence = new int[length];
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet RenewalIterator.
 *
 *  argument : reference sur un objet RenewalIterator.
 *
 *--------------------------------------------------------------*/

void RenewalIterator::copy(const RenewalIterator &iter)

{
  register int i;
  int *psequence , *isequence;


  renewal = iter.renewal;
  (renewal->nb_iterator)++;

  interval = iter.interval;
  counter = iter.counter;
  length = iter.length;

  sequence = new int[length];

  psequence = sequence;
  isequence = iter.sequence;
  for (i = 0;i < length;i++) {
    *psequence++ = *isequence++;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe RenewalIterator.
 *
 *--------------------------------------------------------------*/

RenewalIterator::~RenewalIterator()

{
  (renewal->nb_iterator)--;
  delete [] sequence;
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe RenewalIterator.
 *
 *  argument : reference sur un objet RenewalIterator.
 *
 *--------------------------------------------------------------*/

RenewalIterator& RenewalIterator::operator=(const RenewalIterator &iter)

{
  if (&iter != this) {
    (renewal->nb_iterator)--;
    delete [] sequence;

    copy(iter);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Simulation par un processus de renouvellement.
 *
 *  arguments : longueur de la sequence,
 *              type d'initialisation ('o' : ordinaire, 'e' : en equilibre).
 *
 *--------------------------------------------------------------*/

void RenewalIterator::simulation(int ilength , char type)

{
  register int i;
  int offset , *psequence;


  switch (type) {
  case 'o' :
    offset = 1;
    break;
  default :
    offset = 0;
    break;
  }

  if (ilength + offset != length) {
    length = ilength + offset;
    delete [] sequence;
    sequence = new int[length];
  }

  psequence = sequence;

  switch (type) {
  case 'o' :
    interval = renewal->inter_event->simulation();
    *psequence++ = 1;
    counter = 1;
    break;
  case 'e' :
    interval = renewal->forward->simulation();
    counter = 1;
    break;
  }

  for (i = offset;i < length;i++) {
    if (counter < interval) {
      *psequence++ = 0;
      counter++;
    }

    else {
      interval = renewal->inter_event->simulation();
      *psequence++ = 1;
      counter = 1;
    }
  }
}


};  // namespace sequence_analysis
