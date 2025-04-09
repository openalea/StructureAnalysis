/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *
 *       Copyright 2005-2009 UMR DAP
 *
 *       File author(s): F. Chaubert-Pereira (chaubert@cirad.fr)
 *
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

#include<math.h>
#include<iostream>
#include<iomanip>
#include<limits.h>
#include "stat_tool/stat_label.h"
#include "tool/util_math.h"
#include "stat_tool/stat_tools.h"
#include "tool/config.h"
#include "stat_tool/curves.h"
#include "sequence_analysis/renewal.h"
#include "stat_tool/markovian.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/sequence_label.h"
#include "stat_tool/vectors.h"
#include "sequence_analysis/tops.h"
#include "continuous_histo.h"
#include "switching_sequence.h"
#include "switching_process.h"
#include "markov_switching.h"

using namespace std;


/*--------------------------------------------------------------
 *
 * Constructeur par defaut de la classe Markov_switching_data.
 *
 *--------------------------------------------------------------*/

Markov_switching_data::Markov_switching_data()
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}


/*--------------------------------------------------------------
 *
 * Constructeur de la classe Markov_switching_data. 
 *
 * arguments : nombre de variables, histogramme des longueurs des séquences, 
 *             nombre de covariables, nombre d'effets aléatoires, 
 *             nombre d'états, présence/absence de constante,
 *             flag initialisation
 *
 *--------------------------------------------------------------*/

Markov_switching_data::Markov_switching_data(int inb_variable , const Histogram &ihlength ,
					     int inb_covariable, int inb_random, int inb_val, int iconstant, bool init_flag)
  :Switching_sequence(inb_variable , ihlength , inb_covariable, inb_random, inb_val, iconstant, init_flag)
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}

/*--------------------------------------------------------------
 *
 * Constructeur de la classe Markov_switching_data. 
 *
 * arguments : nombre de variables, histogramme des longueurs des séquences, 
 *             nombre de covariables, nombre d'effets aléatoires, 
 *             nombre d'états,nombre de temps, présence/absence de constante,
 *             flag initialisation
 *
 *--------------------------------------------------------------*/

Markov_switching_data::Markov_switching_data(int inb_variable , const Histogram &ihlength ,
					     int inb_covariable, int inb_random, int inb_val, int T, int iconstant, bool init_flag)
  :Switching_sequence(inb_variable , ihlength , inb_covariable, inb_random, inb_val, T, iconstant, init_flag)
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}


/*--------------------------------------------------------------
 *
 * Constructeur de la classe Markov_switching_data. 
 *
 * arguments : nombre de séquences, nombre de covariables,  
 *             nombre d'effets aléatoires, nombre d'états, longueur
 *             des séquences, index, séquences, covariables, effets 
 *             aléatoires individuels, effets aléatoires temporels, 
 *             présence/absence de constante, identificateur des séquences,
 *             référence sur un objet Markov_switching
 *
 *--------------------------------------------------------------*/

Markov_switching_data::Markov_switching_data(int inb_sequence , int inb_covariable, int inb_random, int inb_val,
					     int *ilength , int **iindex, double ***ireal_sequence , 
					     double ***icovar, double **ieffect,double *iyear_effect,
					     int iconstant, int *iidentifier, const Markov_switching &isw_markov)
  :Switching_sequence(2 , inb_sequence , inb_covariable, inb_random, inb_val,
		      ilength , iindex, ireal_sequence, icovar, ieffect, iyear_effect, iidentifier, iconstant)
{
  register int i;

  type[0] = STATE;
  sw_markov = new Markov_switching(isw_markov, false);

  // extraction des caracteristiques des sequences simulees

  build_transition_count(*sw_markov);

  sw_markov->characteristic_computation(*this , true);

  // calcul de la vraisemblance

  likelihood = sw_markov->likelihood_computation(*this);

  posterior_probability = new double[inb_sequence];
  for (i = 0; i < inb_sequence; i++){
    posterior_probability[i] = 0.;
  }

}


/*--------------------------------------------------------------
 *
 * Construction d'un objet Markov_switching_data a partir
 * d'un objet Switching_sequence.
 *
 * arguments : référence sur un objet Switching_sequence, flag
 *             d'initialisation
 *
 *--------------------------------------------------------------*/

Markov_switching_data::Markov_switching_data(const Switching_sequence &isw_seq ,
                                                       bool initial_run_flag)
:Switching_sequence(isw_seq , 'c' , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}


/*--------------------------------------------------------------
 *
 * Construction d'un objet Markov_switching_data a partir
 * d'un objet Switching_sequence avec ajout d'une variable.
 *
 * arguments : référence sur un objet Switching_sequence
 *             numéro de variable
 *
 *--------------------------------------------------------------*/

Markov_switching_data::Markov_switching_data(const Switching_sequence &isw_seq ,
                                                       int variable)
:Switching_sequence(isw_seq , 'a' , variable , DEFAULT)
{
  sw_markov = 0;
  chain_data = 0;
  likelihood = D_INF;
  hidden_likelihood = D_INF;
  posterior_probability = 0;
}


/*--------------------------------------------------------------
 *
 * Copie d'un objet Markov_switching_data.
 *
 * arguments : référence sur un objet Markov_switching_data, flag modèle
 *
 *---------------------------------------------------------------*/

void Markov_switching_data::copy(const Markov_switching_data &isw_markov_data ,
                                      bool model_flag)
{
  register int i;

  if ((model_flag) && (isw_markov_data.sw_markov)) {
    sw_markov = new Markov_switching(*(isw_markov_data.sw_markov) , false);
  }
  else {
    sw_markov = 0;
  }

  if (isw_markov_data.chain_data) {
    chain_data = new Chain_data(*(isw_markov_data.chain_data));
  }
  else {
    chain_data = 0;
  }

  likelihood = isw_markov_data.likelihood;
  hidden_likelihood = isw_markov_data.hidden_likelihood;

  if (isw_markov_data.posterior_probability){
    posterior_probability = new double[nb_sequence];
    for (i = 0; i < nb_sequence; i++){
      posterior_probability[i] = isw_markov_data.posterior_probability[i];
    }
  }
  else {
    posterior_probability = 0;
  }

}


/*--------------------------------------------------------------
 *
 * Destructeur de la classe Markov_switching_data.
 *
 *--------------------------------------------------------------*/

Markov_switching_data::~Markov_switching_data()
{
  if (sw_markov){
    delete sw_markov;
  }

  if(chain_data){
    delete chain_data;
  }

  delete [] posterior_probability;

}


/*--------------------------------------------------------------
 *
 * Operateur d'assignement de la classe Markov_switching_data.
 *
 * arguments : référence sur un objet Markov_switching_data
 *
 *--------------------------------------------------------------*/

Markov_switching_data& Markov_switching_data::operator=(const Markov_switching_data &isw_markov_data)
{
  if (&isw_markov_data != this) {
    delete sw_markov;
    delete chain_data;

    remove();
    Sequences::remove();

    Sequences::copy(isw_markov_data);
    Switching_sequence::copy(isw_markov_data);
    copy(isw_markov_data);
  }

  return *this;
}


/*--------------------------------------------------------------
 *
 * Ecriture d'un objet Markov_switching_data
 *
 * arguments : stream, flag niveau de détail
 *
 *--------------------------------------------------------------*/

ostream& Markov_switching_data::ascii_write(ostream &os , bool exhaustive) const
{
  
  Switching_sequence::ascii_write(os, exhaustive);

  if (sw_markov) {
    sw_markov->ascii_write(os , this , exhaustive , 
			   false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov_switching_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Markov_switching_data::ascii_write(Format_error &error , const char *path ,
                                             bool exhaustive) const
{
  bool status = false;
  
  if (sw_markov) {
    ofstream out_file(path);
    error.init();
    
    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }
    
    else {
      status = true;
      sw_markov->ascii_write(out_file , this , exhaustive , true );
    }
  }
  
  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Markov_switching_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path
 *
 *--------------------------------------------------------------*/

bool Markov_switching_data::spreadsheet_write(Format_error &error , const char *path) const
{
  bool status = false;

  if (sw_markov) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      sw_markov->spreadsheet_write(out_file , this ,false); 
    }
  }
  
  return status;
}


/*--------------------------------------------------------------
 *
 * Comptage des etats initiaux et des transitions.
 *
 * arguments : référence sur un objet Chain_data, ordre, flag début
 *
 *--------------------------------------------------------------*/

void Switching_sequence::transition_count_computation(const Chain_data &chain_data ,
                                                       int order , bool begin) const
{
  register int i , j , k;
  int **ttransition , *pstate , *mstate , power[ORDER];


  for (i = 0;i < chain_data.nb_state;i++) {
    chain_data.initial[i] = 0;
  }

  for (i = 0;i < chain_data.nb_row;i++) {
    for (j = 0;j < chain_data.nb_state;j++) {
      chain_data.transition[i][j] = 0;
    }
  }

  // extraction des etats initiaux et des transitions

  i = 1;
  for (j = 0;j < order;j++) {
    power[j] = i;
    i *= chain_data.nb_state;
  }

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    (chain_data.initial[*pstate])++;

    if (!begin) {
       pstate += order - 1;
    }
    for (j = (begin ? 1 : order);j < length[i];j++) {
      ttransition = chain_data.transition;
      mstate = pstate;
      for (k = 0;k < MIN(j - 1 , order);k++) {
        ttransition += *mstate-- * power[order - 1 - k];
      }
      for (k = j - 1;k < order;k++) {
        ttransition += *mstate * power[order - 1 - k];
      }

      (*(*ttransition + *++pstate))++;
    }
  }

  ttransition = NULL;
  pstate = NULL;
  mstate = NULL;

  delete ttransition;
  delete pstate;
  delete mstate;


}


/*--------------------------------------------------------------
 *
 *Construction des comptages des etats initiaux et des transitions. 
 *
 * arguments : référence sur un objet Markov_switching, flag début
 *
 *--------------------------------------------------------------*/

void Markov_switching_data::build_transition_count(const Markov_switching &isw_markov, bool begin)
{
  chain_data = new Chain_data('o' , marginal[0]->nb_value ,
                              (int)pow((double)marginal[0]->nb_value , 1));
  transition_count_computation(*chain_data , 1 , begin);
}


/*--------------------------------------------------------------
 *
 * Calcul des variances de l'effet aléatoire - heterogénéité
 *
 * arguments : index, variances aléatoires, effets aléatoires individuels, 
 *             état
 *--------------------------------------------------------------*/

double Markov_switching_data::random_hetero_effect_variance(int index, double *random_variance, double **ieffect, int random)
{
  register int i;
  double sum, sum1;
  double variance;

  sum = 0.;
  sum1 = 0.;

  for (i = 0; i < nb_sequence; i++){
    sum = sum + ieffect[i][index] * ieffect[i][index];
    sum1 = sum1 + ieffect[i][index];
  } 

  sum = sum / nb_sequence;
  sum1 = sum1 / nb_sequence;

  variance = random_variance[random] * (sum - sum1*sum1);

  return variance;

}


/*--------------------------------------------------------------
 *
 * Calcul des variances de l'effet aléatoire - effet année
 *
 * arguments : variances aléatoires, effets aléatoires année, état
 *
 *--------------------------------------------------------------*/
double Markov_switching_data::random_year_effect_variance(double *random_variance, double *iyear_effect, int random)
{
  register int i;
  int nb_annee;
  double sum, sum1;
  double variance;

  sum = 0.;
  sum1 = 0.;

  nb_annee = nb_year(); //longueur de year_effect

  for (i = 0; i < nb_annee; i++){
    sum = sum + iyear_effect[i]* iyear_effect[i];
    sum1 = sum1 + iyear_effect[i];
  } 

  sum = sum / nb_annee;
  sum1 = sum1 / nb_annee;

  variance = random_variance[random] * (sum - sum1*sum1);

  return variance;

}


/*--------------------------------------------------------------
 *
 * Calcul des moyennes de l'effet aléatoire - hétérogénéité
 *
 * arguments : index, variances aléatoires, effets aléatoires
 *             individuels, état
 *
 *--------------------------------------------------------------*/

double Markov_switching_data::random_hetero_effect_mean(int index, double *random_variance, double **ieffect, int random)
{
  register int i;
  double sum1;
  double mean;

  sum1 = 0.;

  for (i = 0; i < nb_sequence; i++){
    sum1 = sum1 + ieffect[i][index];
  } 

  sum1 = sum1 / nb_sequence;

  mean = sqrt(random_variance[random]) * sum1;

  return mean;

}

/*--------------------------------------------------------------
 *
 * Calcul des moyennes de l'effet aléatoire - effet année
 *
 * arguments : variances aléatoires, effets aléatoires années, état
 *
 *--------------------------------------------------------------*/

double Markov_switching_data::random_year_effect_mean(double *random_variance, double *iyear_effect, int random)
{
  register int i;
  int nb_annee;
  double sum1 = 0.;
  double mean;

  nb_annee = nb_year();
  for (i = 0; i < nb_annee; i++){
    sum1 = sum1 + iyear_effect[i]*iyear_effect[i];
  } 

  sum1 = sum1 / nb_annee;

  mean = sqrt(random_variance[random]) * sum1;

  return mean;

}
