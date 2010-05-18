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



#include <math.h>
#include "stat_tool/stat_tools.h"
#include "stat_tool/compound.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequences.h"
#include "tops.h"
#include "sequence_label.h"

using namespace std;



/*--------------------------------------------------------------*
 *
 *  Calcul de la position maximum.
 *
 *--------------------------------------------------------------*/

void Tops::max_position_computation()

{
  if (max_position == 0) {
    register int i;


    for (i = 0;i < nb_sequence;i++) {
      if (index_parameter[i][length[i] - 1] > max_position) {
        max_position = index_parameter[i][length[i] - 1];
      }
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la position minimum.
 *
 *--------------------------------------------------------------*/

int Tops::min_position_computation() const

{
  register int i;
  int min_position = 0;


  if (max_position > 0) {
    min_position = max_position;

    for (i = 0;i < nb_sequence;i++) {
      if (index_parameter[i][0] < min_position) {
        min_position = index_parameter[i][0];
      }
    }
  }

  return min_position;
}


/*--------------------------------------------------------------*
 *
 *  Construction des lois empiriques du nombre d'entrenoeuds
 *  axe porteur/axes portes.
 *
 *--------------------------------------------------------------*/

void Tops::build_nb_internode_frequency_distribution()

{
  register int i , j;
  int max_nb_internode , *bmax_length , *pposition , *plength;


  // creation de la loi empirique du nombre d'entrenoeuds de l'axe porteur

  max_nb_internode = 0;
  for (i = 0;i < nb_sequence;i++) {
    if (index_parameter[i][length[i]] > max_nb_internode) {
      max_nb_internode = index_parameter[i][length[i]];
    }
  }

  nb_internode = new FrequencyDistribution(max_nb_internode + 1);

  // constitution de la loi empirique du nombre d'entrenoeuds de l'axe porteur et
  // creation des lois empiriques du nombre d'entrenoeuds des axes portes

  max_position_computation();
  axillary_nb_internode = new FrequencyDistribution*[max_position + 1];

  bmax_length = new int[max_position + 1];
  for (i = 0;i <= max_position;i++) {
    bmax_length[i] = I_DEFAULT;
  }

  for (i = 0;i < nb_sequence;i++) {
    (nb_internode->frequency[index_parameter[i][length[i]]])++;

    pposition = index_parameter[i];
    plength = int_sequence[i][0];

    for (j = 0;j < length[i];j++) {
      if (*plength > bmax_length[*pposition]) {
        bmax_length[*pposition] = *plength;
      }
      pposition++;
      plength++;
    }
  }

  axillary_nb_internode[0] = NULL;
  for (i = 1;i <= max_position;i++) {
    if (bmax_length[i] != I_DEFAULT) {
      axillary_nb_internode[i] = new FrequencyDistribution(bmax_length[i] + 1);
    }
    else {
      axillary_nb_internode[i] = NULL;
    }
  }

  // constitution des lois empiriques du nombre d'entrenoeuds des axes portes

  for (i = 0;i < nb_sequence;i++) {
    pposition = index_parameter[i];
    plength = int_sequence[i][0];
    for (j = 0;j < length[i];j++) {
      (axillary_nb_internode[*pposition++]->frequency[*plength++])++;
    }
  }

  delete [] bmax_length;

  // calcul des caracteristiques des lois empiriques

  nb_internode->offset_computation();
  nb_internode->nb_element = nb_sequence;
  nb_internode->max_computation();
  nb_internode->mean_computation();
  nb_internode->variance_computation();

  for (i = 1;i <= max_position;i++) {
    if (axillary_nb_internode[i]) {
      axillary_nb_internode[i]->offset_computation();
      axillary_nb_internode[i]->nb_element_computation();
      axillary_nb_internode[i]->max_computation();
      axillary_nb_internode[i]->mean_computation();
      axillary_nb_internode[i]->variance_computation();
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la loi de croissance d'un axe porte a une position
 *  donne en fonction des parametres de la cime.
 *
 *  argument : position.
 *
 *--------------------------------------------------------------*/

void TopParameters::position_nb_internode_computation(int position)

{
  DiscreteParametric *nb_trial , *scale_nb_trial , *elementary;
  Compound *cdist;


  nb_trial = new DiscreteParametric(NEGATIVE_BINOMIAL , position , I_DEFAULT ,
                                    (double)position , probability);

  scale_nb_trial = new DiscreteParametric(*((Distribution*)nb_trial) , rhythm_ratio);

#  ifdef DEBUG
//   cout << "\ntest : " << nb_trial->variance / scale_nb_trial->variance << " | "
//        << (nb_trial->mean * nb_trial->mean) / (scale_nb_trial->mean * scale_nb_trial->mean) << endl;
#  endif

  delete nb_trial;

  elementary = new DiscreteParametric(BINOMIAL , 0 , 1 , D_DEFAULT , axillary_probability);

  cdist = new Compound(*scale_nb_trial , *elementary);
  delete scale_nb_trial;
  delete elementary;

  axillary_nb_internode[position] = new Distribution(*((Distribution*)cdist));
  delete cdist;
}


/*--------------------------------------------------------------*
 *
 *  Calcul des lois de croissance des axes portes.
 *
 *  argument : position maximum.
 *
 *--------------------------------------------------------------*/

void TopParameters::axillary_nb_internode_computation(int imax_position)

{
  if (imax_position != max_position) {
    register int i;


    if (max_position > 0) {
      for (i = 1;i <= max_position;i++) {
        delete axillary_nb_internode[i];
      }
      delete [] axillary_nb_internode;
    }

    max_position = imax_position;
    axillary_nb_internode = new Distribution*[max_position + 1];

    axillary_nb_internode[0] = NULL;
    for (i = 1;i <= max_position;i++) {
      position_nb_internode_computation(i);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Estimation des parametres d'une cime.
 *
 *  arguments : reference sur un objet StatError, positions minimum et
 *              maximum, voisinage, flag probabilites egales.
 *
 *--------------------------------------------------------------*/

TopParameters* Tops::estimation(StatError &error , int imin_position , int imax_position ,
                                int neighborhood , bool equal_probability) const

{
  bool status = true;
  register int i , j , k;
  int nb_neighbor , diff , total , variance_total , *pposition , *plength;
  double mixt_position , mean , mixt_mean2 , moment2 , variance;
  TopParameters *parameters;


  parameters = NULL;
  error.init();

  if ((imin_position < 1) || (imin_position > imax_position)) {
    status = false;
    error.update(SEQ_error[SEQR_MIN_POSITION]);
  }
  if ((imax_position < min_position_computation()) || (imax_position < imin_position)) {
    status = false;
    error.update(SEQ_error[SEQR_MAX_POSITION]);
  }

  if ((neighborhood < 0) || (neighborhood > MAX_NEIGHBORHOOD)) {
    status = false;
    error.update(SEQ_error[SEQR_NEIGHBORHOOD]);
  }

  if (status) {

    // verification du voisinage

    nb_neighbor = 0;

    for (i = 0;i < nb_sequence;i++) {
      pposition = index_parameter[i];

      for (j = 0;j < length[i] - 1;j++) {
        if (*pposition >= imin_position) {
          if (*pposition > imax_position) {
            break;
          }

          k = 1;
          while ((j + k < length[i]) && (*(pposition + k) - *pposition < neighborhood)) {
            k++;
          }

          while ((j + k < length[i]) && (*(pposition + k) - *pposition == neighborhood)) {
            nb_neighbor++;
            k++;
          }
        }

        pposition++;
      }
    }

    if (nb_neighbor < NB_NEIGHBOR) {
      status = false;
      error.update(SEQ_error[SEQR_NEIGHBORS]);
    }
  }

  if (status) {

    // creation de l'objet TopParameters

    parameters = new TopParameters();
    parameters->tops = new Tops(*this , false);

    // accumulation position moyenne, moyenne, moment d'ordre 2 et
    // variations de la longueur des axes portes

    total = 0;
    variance_total = 0;
    mixt_position = 0.;
    mean = 0.;
    moment2 = 0.;
    variance = 0.;

    for (i = 0;i < nb_sequence;i++) {
      pposition = index_parameter[i];
      plength = int_sequence[i][0];

      for (j = 0;j < length[i];j++) {
        if (*pposition >= imin_position) {
          if (*pposition > imax_position) {
            break;
          }

          total++;
          mean += *plength;
          moment2 += *plength * *plength;
          mixt_position += *pposition;

          k = 1;
          while ((j + k < length[i]) && (*(pposition + k) - *pposition < neighborhood)) {
            k++;
          }

          while ((j + k < length[i]) && (*(pposition + k) - *pposition == neighborhood)) {
            variance_total++;
            diff = *(plength + k) - *plength;
            variance += diff * diff;
            k++;
          }
        }

        pposition++;
        plength++;
      }
    }

    mixt_position /= total;
    mean /= total;
    moment2 /= total;
    variance /= (2 * variance_total);

    // accumulation melange des moyennes aux carres des longueurs
    // des axes portes pour chaque position

    mixt_mean2 = 0.;
    for (i = imin_position;i <= MIN(imax_position , max_position);i++) {
      if (axillary_nb_internode[i]) {
        mixt_mean2 += axillary_nb_internode[i]->nb_element *
                      axillary_nb_internode[i]->mean * axillary_nb_internode[i]->mean;
      }
    }
    mixt_mean2 /= total;

#   ifdef DEBUG
    cout << "\n" << nb_neighbor << " " << neighborhood << " | " << max_position << "\n\n"
         << mixt_position << " " << mean << " " << moment2 << " (" << total << ")\n"
         << variance << " (" << variance_total << ")\n"
         << mixt_mean2 << endl;
#   endif

    // extraction des parametres de croissance

    parameters->probability = 1. - ((double)neighborhood / 2. + mixt_position) * (moment2 - mixt_mean2) /
                              (mean * mean) - (neighborhood * neighborhood) / (2. * mixt_position) +
                              mixt_position * variance / (mean * mean);
    parameters->axillary_probability = 1. + (neighborhood / (2. * mixt_position)) *
                                       (moment2 - mixt_mean2) / mean + (neighborhood * neighborhood * mean) /
                                       (2. * mixt_position * mixt_position) - variance / mean;

#   ifdef DEBUG
    cout << *parameters;
#   endif

    if ((parameters->probability < TOP_MIN_PROBABILITY) || (parameters->probability > 1.) ||
        (parameters->axillary_probability < TOP_MIN_PROBABILITY) || (parameters->axillary_probability > 1.)) {
      delete parameters;
      parameters = NULL;
      error.update(STAT_error[STATR_ESTIMATION_FAILURE]);
    }

    else {
      if (equal_probability) {
        if (fabs(parameters->probability - parameters->axillary_probability) > PROBABILITY_DIFF) {
          delete parameters;
          parameters = NULL;
          error.update(SEQ_error[SEQR_EQUAL_PROBABILITY]);
        }

        else {
          parameters->probability = (parameters->probability + parameters->axillary_probability) / 2.;
          parameters->axillary_probability = parameters->probability;
        }
      }
    }

    if (parameters) {
      parameters->rhythm_ratio = (mean * parameters->probability) /
                                 (mixt_position * parameters->axillary_probability);
      parameters->axillary_nb_internode_computation(max_position);
    }
  }

  return parameters;
}


/*--------------------------------------------------------------*
 *
 *  Simulation d'un ensemble de cimes.
 *
 *  arguments : reference sur un objet StatError, nombre de cimes,
 *              loi du nombre d'epreuves axes porteurs,
 *              loi du nombre d'axes portes par noeud.
 *
 *--------------------------------------------------------------*/

Tops* TopParameters::simulation(StatError &error , int nb_top ,
                                const Distribution &nb_trial ,
                                const Distribution &nb_axillary) const

{
  bool status = true;
  register int i , j , k , m;
  int nb_trial_value , value , nb_axillary_value , axillary_nb_trial , *nb_position ,
      *pposition , *plength;
  DiscreteParametric *growth , *axillary_growth;
  TopParameters *parameters;
  Tops *tops;


  tops = NULL;
  error.init();

  if ((nb_top < 1) || (nb_top > NB_TOP)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_TOP]);
  }
  if ((nb_trial.offset < 2) || (nb_trial.nb_value > NB_TRIAL)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_TRIAL]);
  }
  if ((nb_axillary.nb_value < 2) || (nb_axillary.nb_value > NB_AXILLARY)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_AXILLARY_SHOOT]);
  }
  if ((int)(nb_top * nb_trial.mean * probability * nb_axillary.mean * 2) > TOP_SIZE) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  if (status) {

    // creation de l'objet Tops

    nb_position = new int[nb_top];
    for (i = 0;i < nb_top;i++) {
      nb_position[i] = (nb_trial.nb_value - 1) * (nb_axillary.nb_value - 1);
    }

    tops = new Tops(nb_top , NULL , nb_position);
    tops->top_parameters = new TopParameters(*this , false);
    parameters = tops->top_parameters;

    // construction de la loi de croissance de l'axe porteur et
    // creation de la loi de croissance des axes portes

    growth = new DiscreteParametric(BINOMIAL , 0 , 1 , D_DEFAULT , parameters->probability);
    axillary_growth = new DiscreteParametric((int)round(nb_trial.nb_value * parameters->rhythm_ratio) , BINOMIAL);

    for (i = 0;i < nb_top;i++) {
      nb_trial_value = nb_trial.simulation();

      pposition = tops->index_parameter[i];
      plength = tops->int_sequence[i][0];
      j = 0;
      k = 0;
      nb_position[i] = 0;

      while (j < nb_trial_value) {
        value = growth->simulation();
        j++;

        if (value > 0) {
          k++;
          nb_axillary_value = nb_axillary.simulation();

          if (nb_axillary_value > 0) {
            nb_position[i] += nb_axillary_value;
            for (m = 0;m < nb_axillary_value;m++) {
              *pposition++ = k;
            }

            axillary_nb_trial = (int)round(j * parameters->rhythm_ratio);

            if (axillary_nb_trial > 0) {
              axillary_growth->init(0 , axillary_nb_trial , D_DEFAULT , parameters->axillary_probability);
              axillary_growth->computation();

              for (m = 0;m < nb_axillary_value;m++) {
                *plength++ = axillary_growth->simulation();
              }
            }

            else {
              for (m = 0;m < nb_axillary_value;m++) {
                *plength++ = 0;
              }
            }
          }
        }
      }

      *pposition = k;
      tops->length[i] = nb_position[i];
    }

    tops->max_length_computation();
    tops->cumul_length_computation();
    delete tops->hlength;
    tops->build_length_frequency_distribution();

    tops->build_index_parameter_frequency_distribution();

    tops->min_value_computation(0);
    tops->max_value_computation(0);
    tops->build_marginal_frequency_distribution(0);

    tops->build_nb_internode_frequency_distribution();
    parameters->axillary_nb_internode_computation(tops->max_position);

    delete [] nb_position;
    delete growth;
    delete axillary_growth;
  }

  return tops;
}


/*--------------------------------------------------------------*
 *
 *  Simulation d'un ensemble de cimes.
 *
 *  arguments : reference sur un objet StatError, nombre de cimes,
 *              nombre d'epreuves axes porteurs,
 *              nombre d'axes portes par noeud.
 *
 *--------------------------------------------------------------*/

Tops* TopParameters::simulation(StatError &error , int nb_top ,
                                int nb_trial , int nb_axillary) const

{
  bool status = true;
  DiscreteParametric *nb_trial_dist , *nb_axillary_dist;
  Tops *tops;


  tops = NULL;
  error.init();

  if ((nb_top < 1) || (nb_top > NB_TOP)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_TOP]);
  }
  if ((nb_trial < 2) || (nb_trial > NB_TRIAL)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_TRIAL]);
  }
  if ((nb_axillary < 1) || (nb_axillary > NB_AXILLARY)) {
    status = false;
    error.update(SEQ_error[SEQR_NB_AXILLARY_SHOOT]);
  }
  if ((int)(nb_top * nb_trial * probability * nb_axillary * 2) > TOP_SIZE) {
    status = false;
    error.update(STAT_error[STATR_SAMPLE_SIZE]);
  }

  if (status) {
    nb_trial_dist = new DiscreteParametric(UNIFORM , nb_trial , nb_trial ,
                                           D_DEFAULT , D_DEFAULT);
    nb_axillary_dist = new DiscreteParametric(UNIFORM , nb_axillary , nb_axillary ,
                                              D_DEFAULT , D_DEFAULT);

    tops = simulation(error , nb_top , *nb_trial_dist , *nb_axillary_dist);

    delete nb_trial_dist;
    delete nb_axillary_dist;
  }

  return tops;
}
