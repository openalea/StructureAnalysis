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

#include <math.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
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
#include "continuous_distribution.h"
#include "switching_sequence.h"

using namespace std;


/*--------------------------------------------------------------
 *
 * Constructeur par défaut de la classe Switching_sequence
 *
 *--------------------------------------------------------------*/

Switching_sequence::Switching_sequence()
{
  constant = 0;
  nb_covariable = 0;
  nb_random = 0;
  nb_val = 0;
  observation = 0;
  covar = 0;
  index = 0;
  effect = 0;
  year_effect = 0;
  self_transition = 0;
  characteristics = 0;
}


/*-------------------------------------------------------------
 * 
 * Initialisation par défaut des champs de la classe Switching_sequence
 *
 *-------------------------------------------------------------*/

void Switching_sequence::init()
{  
  register int i;

  self_transition = 0;
  observation = 0;

  characteristics = new Sequence_characteristics*[nb_variable];
  for (i = 0;i < nb_variable;i++) {
    characteristics[i] = 0;
  }
}


/*-------------------------------------------------------------
 *
 * Calcul du nombre d'années distinctes et donc de la longueur 
 * de year_effect 
 *
 *-------------------------------------------------------------*/

int Switching_sequence::nb_year()
{  
  register int i, j;
  int nb_annee = 1;

  for (i = 0; i < nb_sequence; i++){
    for ( j = 0; j < length[i]; j++){
      if ( index[i][j] > nb_annee){
	nb_annee = index[i][j];
      }
    }
  }

  return nb_annee;
  
}


/*-------------------------------------------------------------
 *
 * Constructeur de la classe Switching_sequence
 *
 * arguments : nombre de variables, type de chaque variable,
 *             nombre de séquences, identificateur des séquences,
 *             longueurs des séquences, nombre de covariables,
 *             nombre d'effets aléatoires, nombre de valeurs,
 *             présence/absence d'une constante, flag initialisation
 *
 *-------------------------------------------------------------*/

Switching_sequence::Switching_sequence(int inb_variable, int *itype, int inb_sequence, int *iidentifier, int *ilength, 
				       int inb_covariable, int inb_random, int inb_val, int iconstant, bool init_flag)
  :Sequences(inb_sequence, iidentifier, ilength, inb_variable, itype, init_flag)
// Sequences( inb_variable, itype, inb_sequence, iidentifier, ilength, init_flag)
{
  register int i, j;

  init();

  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }

  nb_random = inb_random;
  nb_val = inb_val;

  covar = new double**[nb_sequence];
  for (i = 0; i < nb_sequence; i++) {
    covar[i] = new double*[length[i]];
    for (j = 0; j < length[i]; j++) {
      covar[i][j] = 0;
    }
  }

  if ((nb_val > 0) && (nb_random !=0)){
    effect = new double*[nb_sequence];
    for (i = 0; i < nb_sequence; i++) {
      effect[i] = 0;
    }
  }
  else {
    effect = 0;
  }

  if ((nb_val == 0) && (nb_random != 0)){
    index = new int*[nb_sequence];
    for ( i = 0; i < nb_sequence; i++){
      index[i] = 0;
    }
  }
  else {
    index = 0;
  }

  year_effect = 0;
  
}


/*-------------------------------------------------------------
 *
 * Constructeur de la classe Switching_sequence
 *
 * arguments : nombre de variables, nombre de séquences, 
 *             identificateur des séquences,
 *             longueurs des séquences, nombre de covariables,
 *             nombre d'effets aléatoires, nombre de valeurs,
 *             présence/absence d'une constante, flag initialisation
 *
 *-------------------------------------------------------------*/

Switching_sequence:: Switching_sequence(int inb_variable, int inb_sequence, int *iidentifier, int *ilength, 
					int inb_covariable, int inb_random, 
					int inb_val, int iconstant, bool init_flag)
{
  register int i, j, k;

  double *pisequence;

  nb_sequence = inb_sequence;

  identifier = new int[nb_sequence];
  if (iidentifier) {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = iidentifier[i];
    }
  }
  else {
    for (i = 0;i < nb_sequence;i++) {
      identifier[i] = i + 1;
    }
  }

  length = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    length[i] = ilength[i];
  }

  max_length_computation();
  cumul_length_computation();
  build_length_histogram();

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = REAL_VALUE;
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal[i] = 0;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_sequence[i][j] = 0; 
      real_sequence[i][j] = new double[length[i]];

      if (init_flag) {
        pisequence = real_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = 0;
        }
      }
    }
  }


  init();
  
  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }
  
  nb_random = inb_random;
  nb_val = inb_val;

  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = 0;
    }
  }

  if ((nb_val > 0) && (nb_random !=0)){
    effect = new double*[nb_sequence];
    for (i = 0; i < nb_sequence; i++) {
      effect[i] = 0;
    }
  }
  else {
    effect = 0;
  }

  if ((nb_val == 0) && (nb_random != 0)){
    index = new int*[nb_sequence];
    for ( i = 0; i < nb_sequence; i++){
      index[i] = 0;
    }
  }
  else {
    index = 0;
  }

  year_effect = 0;

}


/*-------------------------------------------------------------
 *
 * Constructeur de la classe Switching_sequence
 *
 * arguments : nombre de variables, histogramme des longueurs des séquences
 *             nombre de covariables, nombre d'effets aléatoires, nombre de valeurs,
 *             présence/absence d'une constante, flag initialisation
 *
 *-------------------------------------------------------------*/

Switching_sequence::Switching_sequence(int inb_variable, const Histogram &ihlength, int inb_covariable, 
				       int inb_random, int inb_val, int iconstant, bool init_flag)
{
  register int i, j, k;
  int nb_annee = 1;

  double *pcovar, *peffect, *pyear_effect;
  int *pindex;

  int *plength;
  double *pisequence;

  nb_sequence = ihlength.nb_element;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  plength = length;
  for (i = ihlength.offset;i < ihlength.nb_value;i++) {
    for (j = 0;j < ihlength.frequency[i];j++) {
      *plength++ = i;
    }
  }

  max_length = ihlength.nb_value - 1;
  cumul_length_computation();
  hlength = new Histogram(ihlength);

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = REAL_VALUE;
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal[i] = 0;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_sequence[i][j] = new int[length[i]];
      real_sequence[i][j] = new double[length[i]];

      if (init_flag) {
        pisequence = real_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = 0;
        }
      }
    }
  }


  init();

  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }
  
  nb_random = inb_random;
  nb_val = inb_val;

  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = new double[nb_covariable];

      if (init_flag) {
        pcovar = covar[i][j];
        for (k = 0;k < inb_covariable;k++) {
          *pcovar++ = 0;
        }
      }
    }
  }

  if ((nb_val > 0) && (nb_random != 0)) {
    effect = new double*[nb_sequence];
    for (i = 0; i < nb_sequence; i++) {
      effect[i] = new double[nb_val];
      if (init_flag) {
	peffect = effect[i];
	for (j = 0; j < nb_val; j++){
	  *peffect++ = 0;
	}
      }
    }
  }
  else {
    effect = 0;
  }

  if ((nb_val == 0) && (nb_random != 0)){
    index = new int*[nb_sequence];
    for ( i = 0; i < nb_sequence; i++){
      index[i] = new int[length[i]];
      if (init_flag) {
	pindex = index[i];
	for (j = 0; j < length[i]; j++){
	  *pindex++ = 0;
	}
      }
    }
  }
  else {
    index = 0;
  }

  if (index) {
    nb_annee = nb_year();
    year_effect = new double[nb_annee];
    if (init_flag){
      pyear_effect = year_effect;
      for (i = 0; i < nb_annee; i++) {
	*pyear_effect++ = 0. ;
      }
    }
  }
  else {
    year_effect = 0;
  }

}


/*-------------------------------------------------------------
 *
 * Constructeur de la classe Switching_sequence
 *
 * arguments : nombre de variables, histogramme des longueurs des séquences
 *             nombre de covariables, nombre d'effets aléatoires, nombre de valeurs,
 *             nombre d'années différentes,
 *             présence/absence d'une constante, flag initialisation
 *
 *-------------------------------------------------------------*/

Switching_sequence::Switching_sequence(int inb_variable, const Histogram &ihlength, int inb_covariable, 
				       int inb_random, int inb_val, int T, int iconstant, bool init_flag)
//:Sequences(ihlength, inb_variable, init_flag)
//Sequences( inb_variable, ihlength, init_flag)
{
  register int i, j, k;

  double *pcovar, *peffect, *pyear_effect;
  int *pindex;

  int *plength;
  double *pisequence;

  nb_sequence = ihlength.nb_element;

  identifier = new int[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    identifier[i] = i + 1;
  }

  length = new int[nb_sequence];
  plength = length;
  for (i = ihlength.offset;i < ihlength.nb_value;i++) {
    for (j = 0;j < ihlength.frequency[i];j++) {
      *plength++ = i;
    }
  }

  max_length = ihlength.nb_value - 1;
  cumul_length_computation();
  hlength = new Histogram(ihlength);

  index_parameter_type = IMPLICIT_TYPE;
  hindex_parameter = 0;
  index_interval = 0;
  index_parameter = 0;

  nb_variable = inb_variable;

  type = new int[nb_variable];
  min_value = new double[nb_variable];
  max_value = new double[nb_variable];
  marginal = new Histogram*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    type[i] = REAL_VALUE;
    min_value[i] = 0.;
    max_value[i] = 0.;
    marginal[i] = 0;
  }

  int_sequence = new int**[nb_sequence];
  real_sequence = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    int_sequence[i] = new int*[nb_variable];
    real_sequence[i] = new double*[nb_variable];
    for (j = 0;j < nb_variable;j++) {
      int_sequence[i][j] = new int[length[i]];
      real_sequence[i][j] = new double[length[i]];

      if (init_flag) {
        pisequence = real_sequence[i][j];
        for (k = 0;k < length[i];k++) {
          *pisequence++ = 0;
        }
      }
    }
  }

  init();


  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }
  
  nb_random = inb_random;
  nb_val = inb_val;

  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = new double[nb_covariable];

      if (init_flag) {
        pcovar = covar[i][j];
        for (k = 0;k < inb_covariable;k++) {
          *pcovar++ = 0;
        }
      }
    }
  }

  if ((nb_val > 0) && (nb_random != 0)) {
    effect = new double*[nb_sequence];
    for (i = 0; i < nb_sequence; i++) {
      effect[i] = new double[nb_val];
      if (init_flag) {
	peffect = effect[i];
	for (j = 0; j < nb_val; j++){
	  *peffect++ = 0;
	}
      }
    }
  }
  else {
    effect = 0;
  }

  if ((nb_val == 0) && (nb_random != 0)) {
    index = new int*[nb_sequence];
    for ( i = 0; i < nb_sequence; i++){
      index[i] = new int[length[i]];
      if (init_flag) {
	pindex = index[i];
	for (j = 0; j < length[i]; j++){
	  *pindex++ = 0;
	}
      }
    }
  }
  else {
    index = 0;
  }
  
  if (index) {
    year_effect = new double[T];
    if (init_flag){
      pyear_effect = year_effect;
      for (i = 0; i < T; i++) {
	*pyear_effect++ = 0. ;
      }
    }
  }
  else {
    year_effect = 0;
  }

}


/*-------------------------------------------------------------
 *
 * Constructeur de la classe Switching_sequence
 *
 * arguments : nombre de variables, nombre de séquences,
 *             nombre de covariables, nombre d'effets aléatoires, 
 *             nombre de valeurs, longueurs des séquences, paramètre d'index,
 *             séquences, covariables, effet aléatoire individuel, 
 *             effet aléatoire temporel,identificateurs des séquences,
 *             présence/absence d'une constante
 *
 *-------------------------------------------------------------*/

Switching_sequence:: Switching_sequence(int inb_variable, int inb_sequence, int inb_covariable, int inb_random, int inb_val,
					int *ilength, int **iindex, double ***ireal_sequence, double ***icovar, double **ieffect, 
					double *iyear_effect, int *iidentifier, int iconstant)
  :Sequences(inb_sequence, iidentifier, ilength, inb_variable, ireal_sequence)
//Sequences(inb_variable , inb_sequence , ilength , isequence)
{
  register int i, j, k;
  int nb_annee;

  double *pcovar , *ccovar, *peffect, *ceffect, *pyear_effect, *cyear_effect;
  int *pindex, *cindex;

  init();

  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }

  nb_random = inb_random;
  nb_val = inb_val;

  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = new double[nb_covariable];
      
      pcovar = covar[i][j];
      ccovar = icovar[i][j];
      

      if (!constant) {
	*pcovar++ = 1;
      }
      if((nb_covariable > 1) || (constant)){
	for (k = 0;k < inb_covariable; k++) {
	  *pcovar++ = *ccovar++;
	}
      }
    }
  }

  if ((nb_val > 0) && (nb_random != 0)){
    effect = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      effect[i] = new double[nb_val];
      
      peffect = effect[i];
      ceffect = ieffect[i];
      for (j = 0;j < nb_val;j++) {

	if (nb_random > 0){
	  *peffect++ = *ceffect++;
	}
      }
    }
  }
  else {
    effect = 0;
  }

  if ((nb_val == 0) && (nb_random != 0)){
    index = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      index[i] = new int[length[i]];
      
      pindex = index[i];
      cindex = iindex[i];
      for (j = 0;j < length[i];j++) {
	*pindex++ = *cindex++;
      }
    }
  }
  else {
    index = 0 ;
  }

  
  if ((nb_val == 0) && (nb_random != 0)){
    nb_annee = nb_year();
    year_effect = new double[nb_annee];
    pyear_effect = year_effect;
    cyear_effect = iyear_effect;
    for (i = 0;i < nb_annee;i++) {
      if (nb_random > 0){
	*pyear_effect++ = *cyear_effect++;
      }
    }
  }
  else {
    year_effect = 0;
  }
  
  //   build_characteristic(); // en commentaire suite aux changements faits par Yann dans la classe sequences
}


/*-------------------------------------------------------------
 *
 * Constructeur de la classe Switching_sequence
 *
 * arguments : objet de type Sequences, nombre de covariables, 
 *             nombre d'effets aléatoires, nombre de valeurs, 
 *             paramètre d'index, covariables, effet aléatoire individuel, 
 *             effet aléatoire temporel, présence/absence d'une constante
 *
 *-------------------------------------------------------------*/

Switching_sequence::Switching_sequence(const Sequences &seq, int inb_covariable, int inb_random, int inb_val, int **iindex, 
				       double ***icovar, double **ieffect, double *iyear_effect, int iconstant)
  :Sequences(seq)
{
  register int i, j, k;
  int nb_annee;

  double *pcovar , *ccovar, *peffect, *ceffect, *pyear_effect, *cyear_effect;
  int *pindex, *cindex;
  
  init();
  
  constant = iconstant;
  if (!constant){
    nb_covariable = inb_covariable + 1;
  }
  else {
    nb_covariable = inb_covariable;
  }

  nb_random = inb_random;
  nb_val = inb_val;


  covar = new double**[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    covar[i] = new double*[length[i]];
    for (j = 0;j < length[i];j++) {
      covar[i][j] = new double[nb_covariable];
      
      pcovar = covar[i][j];
      ccovar = icovar[i][j];
      if (!constant) {
	*pcovar++ = 1;
      }

      if((nb_covariable > 1) || (constant)) {
	for (k = 0;k < inb_covariable;k++) {
	  *pcovar++ = *ccovar++;
	}
      }
    }
  }
 

  if ((nb_val > 0) && (nb_random != 0)){
    effect = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      effect[i] = new double[nb_val];
      
      peffect = effect[i];
      ceffect = ieffect[i];
      for (j = 0;j < nb_val;j++) {

	if (nb_random > 0){
	  *peffect++ = *ceffect++;
	}
      }
    }
  }
  else {
    effect = 0;
  }


  if ((nb_val == 0) && (nb_random != 0)){
    index = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      index[i] = new int[length[i]];
      
      pindex = index[i];
      cindex = iindex[i];
      for (j = 0;j < length[i];j++) {
	*pindex++ = *cindex++;
      }
    }
  }
  else {
    index = 0;
  }

  if ((nb_val == 0) && (nb_random != 0)){
    nb_annee = nb_year();
    year_effect = new double[nb_annee];
    pyear_effect = year_effect;
    cyear_effect = iyear_effect;
    for (i = 0;i < nb_annee; i++) {
      if (nb_random > 0){
	*pyear_effect++ = *cyear_effect++;
      }
    }
  }
  else {
    year_effect = 0;
  }

  build_characteristic();
}


/*--------------------------------------------------------------
 *
 * Destructeur des champs de la classe Switching_sequence.
 *
 *--------------------------------------------------------------*/

void Switching_sequence::remove()
{
  register int i , j;

  if (self_transition) {
    for (i = 0;i < marginal[0]->nb_value;i++) {
      delete self_transition[i];
    }
    delete [] self_transition;
  }

  if (observation) {
    for (i = 1;i < nb_variable;i++) {
      for (j = 0;j < marginal[0]->nb_value;j++) {
        delete observation[i][j];
      }
      delete [] observation [i];
    }
    delete [] observation;
  }

  if (characteristics) {
    for (i = 0;i < nb_variable;i++) {
      delete characteristics[i];
    }
    delete [] characteristics;
  }

  if (covar) {
    for (i = 0;i < nb_sequence;i++) {
      for (j = 0;j < length[i];j++) {
        delete [] covar[i][j];
      }
      delete [] covar[i];
    }
    delete [] covar;
  }

  if (effect) {
    for (i = 0;i < nb_sequence;i++) {
      delete [] effect[i];
    }
    delete [] effect;
  }

  if (index) {
    for (i = 0;i < nb_sequence;i++) {
      delete [] index[i];
    }
    delete [] index;
  }

  if(year_effect){
    delete [] year_effect;
  }

}


/*--------------------------------------------------------------
 *
 * Destructeur de la classe Switching_sequence.
 *
 *--------------------------------------------------------------*/

Switching_sequence::~Switching_sequence()
{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Switching_sequence.
 *
 *  arguments : reference sur un objet Switching_sequence, flag inversion.
 *
 *--------------------------------------------------------------*/

void Switching_sequence::copy(const Switching_sequence &sw_seq , int param)
{
  bool initial_run_flag;

  register int i , j, k;
  int nb_annee;

  double *pcovar, *ccovar, *peffect, *ceffect, *pyear_effect, *cyear_effect ;
  int *pindex, *cindex;  

  constant = sw_seq.constant;
  nb_covariable = sw_seq.nb_covariable;
  nb_random = sw_seq.nb_random;
  nb_val = sw_seq.nb_val;

  if ((sw_seq.self_transition) && (param != REVERSE)) {
    self_transition = new Self_transition*[marginal[0]->nb_value];
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (sw_seq.self_transition[i]) {
        self_transition[i] = new Self_transition(*(sw_seq.self_transition[i]));
      }
      else {
        self_transition[i] = 0;
      }
    }
  }
  
  else {
    self_transition = 0;
  }

  if (sw_seq.observation) {
    observation = new Continuous_histo**[nb_variable];
    observation[0] = 0;

    for (i = 1;i < nb_variable;i++) {
      observation[i] = new Continuous_histo*[marginal[0]->nb_value];
      for (j = 0;j < marginal[0]->nb_value;j++) {
        observation[i][j] = new Continuous_histo(*(sw_seq.observation[i][j]));
      }
    }
  }

  else {
    observation = 0;
  }

  if (sw_seq.covar) {
    covar = new double**[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      covar[i] = new double*[sw_seq.length[i]];
      for (j = 0;j < sw_seq.length[i];j++) {
	covar[i][j] = new double[nb_covariable];
	pcovar = covar[i][j];
	
	ccovar = sw_seq.covar[i][j];
	if (sw_seq.covar[i][j]){
	  for (k = 0;k < nb_covariable;k++) {
	    *pcovar++ = *ccovar++;
	  }
	}
	else {
	  pcovar = ccovar;
	}
      }
    }
  }
  else {
    covar = 0;
  }

  if ((sw_seq.nb_val > 0) && (sw_seq.nb_random != 0)){
    effect = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      effect[i] = new double[sw_seq.nb_val];
      
      peffect = effect[i];      
      ceffect = sw_seq.effect[i];
      if (sw_seq.effect[i]){
	for (j = 0;j < sw_seq.nb_val;j++) {
	  *peffect++ = *ceffect++;
	}
      }
      else {
	peffect = ceffect;
      }
    }
  }
  else {
    effect = 0;
  }
   
  if ((!effect) && (sw_seq.nb_random != 0)){
    if (sw_seq.index) {
      index = new int*[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
	index[i] = new int[sw_seq.length[i]];
	
	pindex = index[i];
	cindex = sw_seq.index[i];
	if (sw_seq.index[i]){
	  for (j = 0;j < sw_seq.length[i];j++) {
	    *pindex++ = *cindex++;
	  }
	}
	else {
	  pindex = cindex;
	}
      }
    }
  }
  else {
    index = 0;
  }
  
  if ((!effect) && (sw_seq.nb_random != 0)){
    if (sw_seq.year_effect){

      nb_annee = nb_year();
      year_effect = new double[nb_annee];
	
      pyear_effect = year_effect;
      cyear_effect = sw_seq.year_effect;
      for (i = 0;i < nb_annee; i++) {
	*pyear_effect++ = *cyear_effect++;
      }
    }
    else {
      pyear_effect = cyear_effect;
    }
  }
  else {
    year_effect = 0;
  }

  characteristics = new Sequence_characteristics*[nb_variable];

  for (i = 0;i < nb_variable;i++) {
    if (sw_seq.characteristics[i]) {
      if ((param == ADD_INITIAL_RUN) || (param == REMOVE_INITIAL_RUN)) {
        switch (param) {
        case ADD_INITIAL_RUN :
          initial_run_flag = true;
          break;
        case REMOVE_INITIAL_RUN :
          initial_run_flag = false;
          break;
        }

        characteristics[i] = new Sequence_characteristics(*(sw_seq.characteristics[i]) , initial_run_flag);

        if (((sw_seq.characteristics[i]->initial_run) && (!initial_run_flag)) ||
            ((!(sw_seq.characteristics[i]->initial_run)) && (initial_run_flag))) {
	  build_sojourn_time_histogram(i , initial_run_flag);
        }
      }

      else if (param == REVERSE) {
        characteristics[i] = new Sequence_characteristics(*(sw_seq.characteristics[i]), 'r');

        build_index_value(i);
	build_first_occurrence_histogram(i);
	
	if (!(sw_seq.characteristics[i]->initial_run)) {
	  build_sojourn_time_histogram(i);
	}
      }

      else {
        characteristics[i] = new Sequence_characteristics(*(sw_seq.characteristics[i]));
      }
    }

    else {
      characteristics[i] = 0;
    }
  }

}


/*--------------------------------------------------------------*
 *
 *  Constructeur avec ajout par copie de la classe Switching_sequence.
 *
 *  arguments : reference sur un objet Switching_sequence, variable à ajouter,
 *              inversion.
 *
 *--------------------------------------------------------------*/

void Switching_sequence::add_variable(const Switching_sequence &sw_seq ,
				      int variable , int param)
{
  bool initial_run_flag;

  register int i , j, k;
  int nb_annee;

  double *pcovar, *ccovar, *peffect, *ceffect, *pyear_effect, *cyear_effect;
  int *pindex, *cindex;

  self_transition = 0;
  observation = 0;

  characteristics = new Sequence_characteristics*[nb_variable];

  constant = sw_seq.constant;
  nb_covariable = sw_seq.nb_covariable;
  nb_random = sw_seq.nb_random;
  nb_val = sw_seq.nb_val;

  if (sw_seq.covar) {
    covar = new double**[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      covar[i] = new double*[sw_seq.length[i]];
      for (j = 0;j < sw_seq.length[i];j++) {
	covar[i][j] = new double[nb_covariable];
	pcovar = covar[i][j];
	
	ccovar = sw_seq.covar[i][j];
	if (sw_seq.covar[i][j]){
	  for (k = 0;k < nb_covariable;k++) {
	    *pcovar++ = *ccovar++;
	  }
	}
	else {
	  pcovar = ccovar;
	}
      }
    }
  }
  else {
    covar = 0;
  }
 
  if ((sw_seq.nb_val > 0) && (sw_seq.nb_random != 0)) {
    effect = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      effect[i] = new double[sw_seq.nb_val];
      
      peffect = effect[i];      
      ceffect = sw_seq.effect[i];
      if (sw_seq.effect[i]){
	for (j = 0;j < sw_seq.nb_val;j++) {
	  *peffect++ = *ceffect++;
	}
      }
      else {
	peffect = ceffect;
      }
    }
  }
  else {
    effect = 0;
  }
   
  if ((!effect) && (sw_seq.nb_random != 0)){
    if (sw_seq.index) {
      index = new int*[nb_sequence];
      for (i = 0;i < nb_sequence;i++) {
	index[i] = new int[sw_seq.length[i]];
	pindex = index[i];
	cindex = sw_seq.index[i];
	
	if (sw_seq.index[i]){
	  for (j = 0;j < sw_seq.length[i];j++) {
	    *pindex++ = *cindex++;
	  }
	}
	else {
	  pindex = cindex;
	}
      }
    }
  }
  else {
    index = 0;
  }
  
  if ((!effect) && (sw_seq.nb_random != 0)){
    if (sw_seq.year_effect){
      nb_annee = nb_year();
      year_effect = new double[nb_annee];
	
      pyear_effect = year_effect;
      cyear_effect = sw_seq.year_effect;
      for (i = 0;i < nb_annee; i++) {
	*pyear_effect++ = *cyear_effect++;
      }
    }
    else {
      pyear_effect = cyear_effect;
    }
  }  
  else {
    year_effect = 0;
  }
  
  i = 0;
  for (j = 0;j < nb_variable;j++) {
    if (j != variable) {
      if (sw_seq.characteristics[i]) {
        if ((param == ADD_INITIAL_RUN) || (param == REMOVE_INITIAL_RUN)) {
          switch (param) {
          case ADD_INITIAL_RUN :
            initial_run_flag = true;
            break;
          case REMOVE_INITIAL_RUN :
            initial_run_flag = false;
            break;
          }

          characteristics[j] = new Sequence_characteristics(*(sw_seq.characteristics[i]) , initial_run_flag);

          if (((sw_seq.characteristics[i]->initial_run) && (!initial_run_flag)) ||
              ((!(sw_seq.characteristics[i]->initial_run)) && (initial_run_flag))) {
	    build_sojourn_time_histogram(j , initial_run_flag);
          }
        }

        else {
          characteristics[j] = new Sequence_characteristics(*(sw_seq.characteristics[i]));
        }
      }
      
      else {
        characteristics[j] = 0;
      }
      
      i++;
    }
    
    else {
      characteristics[j] = 0;
    }
  }

}


/*--------------------------------------------------------------*
 *
 *  Constructeur par copie de la classe Switching_sequence.
 *
 *  arguments : reference sur un objet Switching_sequence, 
 *              type de transformation ('c' : copie, 'a' : addition 
 *              d'une variable, 'r' : suppression du parametre d'index),
 *              inversion, ajout.
 *
 *--------------------------------------------------------------*/

Switching_sequence::Switching_sequence(const Switching_sequence &sw_seq , char transform ,
				       int param1 , int param2)
{

  switch (transform) {
  case 'c' :
    Sequences::copy(sw_seq , (param1 == REVERSE ? true : false)); 
    copy(sw_seq , param1);
    break;
  case 'a' :
    Sequences::add_state_variable(sw_seq);
    add_variable(sw_seq , param1 , param2);
    break;
  default :
    Sequences::copy(sw_seq);
    copy(sw_seq);
    break;
  }
}

/*--------------------------------------------------------------
 *
 * Operateur d'assignement de la classe Switching_sequence.
 *
 * arguments : reference sur un objet Switching_sequence
 *
 *--------------------------------------------------------------*/

Switching_sequence& Switching_sequence::operator=(const Switching_sequence &sw_seq)
{

  if (&sw_seq != this) {
    remove();  
    Sequences::remove();
    Sequences::copy(sw_seq);
    copy(sw_seq);
  }
  
  return *this;
}

/*--------------------------------------------------------------
 *
 * Initialisation de la 1ere variable.
 *
 * arguments : type de la variable
 *
 *--------------------------------------------------------------*/

void Switching_sequence::state_variable_init(int itype)
{
  register int i , j;

  if (itype != type[0]) {
    if (type[0] == STATE) {
      if (self_transition) {
        for (i = 0;i < marginal[0]->nb_value;i++) {
          delete self_transition[i];
        }
        delete [] self_transition;

        self_transition = 0;
      }

      if (observation) {
        for (i = 1;i < nb_variable;i++) {
          for (j = 0;j < marginal[0]->nb_value;j++) {
            delete observation[i][j];
          }
          delete [] observation[i];
        }
        delete [] observation;

        observation = 0;
      }
    }

    type[0] = itype;
  }

  for (i = 1;i < nb_variable;i++) {
    type[i] = REAL_VALUE;
  }
}


/*--------------------------------------------------------------
 *
 * Suppression de la 1ere variable.
 *
 *--------------------------------------------------------------*/

Switching_sequence* Switching_sequence::remove_variable_1() const
{
  register int i, j, k;
  int nb_annee;

  int *variable , *itype, *pindex, *cindex;
  double *pcovar, *ccovar, *peffect, *ceffect, *pyear_effect, *cyear_effect;

  Switching_sequence *isw_seq;

  variable = new int[nb_variable - 1];
  itype = new int[nb_variable - 1];
  for (i = 0;i < nb_variable - 1;i++) {
    variable[i] = i + 1;
    itype[i] = type[i + 1];
  }

  if(!constant){
    isw_seq = new Switching_sequence(nb_variable - 1 , itype , nb_sequence ,
				     identifier , length , nb_covariable-1, nb_random, nb_val, constant,  false);
  }
  else {
    isw_seq = new Switching_sequence(nb_variable - 1 , itype , nb_sequence ,
				     identifier , length , nb_covariable, nb_random, nb_val, constant,  false);
  }

  isw_seq->Sequences::select_variable(*this , variable);

  if( covar ){
    isw_seq->covar = new double**[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      isw_seq->covar[i] = new double*[length[i]];
      for (j = 0;j < length[i];j++) {
	isw_seq->covar[i][j] = new double[nb_covariable];
	pcovar = isw_seq->covar[i][j];
	
	ccovar = covar[i][j];
	for (k = 0;k < nb_covariable;k++) {
	  *pcovar++ = *ccovar++;
	}
      }
    }
  }
  else {
    isw_seq->covar = 0;
  }

  if (effect) {
    isw_seq->effect = new double*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      isw_seq->effect[i] = new double[nb_val];
      
      peffect = isw_seq->effect[i];
      ceffect = effect[i];
      
      for (j = 0;j < nb_val;j++) {
	*peffect++ = *ceffect++;
      }
    }
  }
  else {
    isw_seq->effect = 0;
  }
  
  if (index){
    isw_seq->index = new int*[nb_sequence];
    for (i = 0;i < nb_sequence;i++) {
      isw_seq->index[i] = new int[length[i]];
      
      pindex = isw_seq->index[i];
      cindex = index[i];
      for (j = 0;j < length[i];j++) {
	*pindex++ = *cindex++;
      }
    }
  }
  else {
    isw_seq->index = 0;
  }
  
  if (year_effect){
    nb_annee = isw_seq->nb_year();
    isw_seq->year_effect = new double[nb_annee];
    
    pyear_effect = isw_seq->year_effect;
    cyear_effect = year_effect;
    
    for (i = 0;i < nb_annee; i++) {
      *pyear_effect++ = *cyear_effect++;
    }
  }
  else {
    isw_seq->year_effect = 0;
  }


  for (i = 0;i < isw_seq->nb_variable;i++) {
    if (characteristics[i + 1]) {
      isw_seq->characteristics[i] = new Sequence_characteristics(*(characteristics[i + 1]));
    }
  }

  return isw_seq;
}


/*--------------------------------------------------------------
 *
 * Calcul de la quantite d'information en faisant l'hypothese
 * de variables aleatoires independantes et equidistribuees.
 *
 *--------------------------------------------------------------*/

double Switching_sequence::iid_information_computation() const
{
  register int i;
  double information = 0.;

  for (i = (((type[0] != STATE) || (nb_variable == 1)) ? 0 : 1);i < nb_variable;i++) {
    information += marginal[i]->information_computation();
  }

  return information;
}


/*-------------------------------------------------------------- 
 *
 * Extraction des probabilites de chaque valeur en fonction de l'index
 * (pour une variable donnee).
 *
 * arguments : numéro de variable
 *
 *--------------------------------------------------------------*/

void Switching_sequence::build_index_value(int variable)
{
  register int i , j;
  int total , *frequency;

  // creation d'un objet Curves

  characteristics[variable]->index_value = new Curves(marginal[variable]->nb_value ,
                                                      max_length , true , false);
  frequency = new int[marginal[variable]->nb_value];

  // calcul des probabilites de chaque valeur en fonction de l'index

  for (i = 0;i < max_length;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      frequency[j] = 0;
    }

    for (j = 0;j < nb_sequence;j++) {
      if (i < length[j]) {
        frequency[int_sequence[j][variable][i]]++;
      }
    }

    total = 0;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      total += frequency[j];
    }
    characteristics[variable]->index_value->frequency[i] = total;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      characteristics[variable]->index_value->point[j][i] = (double)frequency[j] / (double)total;
    }
  }

  delete [] frequency;
}


/*--------------------------------------------------------------
 *
 * Construction des histogrammes du temps avant la premiere occurrence
 * d'une valeur (pour une variable donnee).
 *
 * arguments : numéro de variable
 *
 *--------------------------------------------------------------*/

void Switching_sequence::build_first_occurrence_histogram(int variable)
{
  bool *occurrence;
  register int i , j;
  int nb_value , *psequence;
  Histogram **first_occurrence;

  // creation des histogrammes

  first_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    first_occurrence[i] = new Histogram(max_length);
  }

  // mise a jour des histogrammes

  occurrence = new bool[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    nb_value = 0;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      occurrence[j] = false;
    }

    psequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      if (!occurrence[*psequence]) {
        occurrence[*psequence] = true;
        (first_occurrence[*psequence]->frequency[j])++;

        nb_value++;
        if (nb_value == marginal[variable]->nb_value) {
          break;
        }
      }

      psequence++;
    }
  }

  delete [] occurrence;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    first_occurrence[i]->nb_value_computation();
    first_occurrence[i]->offset_computation();
    first_occurrence[i]->nb_element_computation();
    first_occurrence[i]->max_computation();
    first_occurrence[i]->mean_computation();
    first_occurrence[i]->variance_computation();
  }

  characteristics[variable]->first_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->first_occurrence[i] = new Histogram(*(first_occurrence[i]));
    delete first_occurrence[i];
  }
  delete [] first_occurrence;

  psequence = NULL;
  delete psequence;

}


/*--------------------------------------------------------------
 *
 * Construction des histogrammes du temps de retour dans une valeur
 * (pour une variable donnee).
 *
 * arguments : numéro de variable
 *
 *--------------------------------------------------------------*/

void Switching_sequence::build_recurrence_time_histogram(int variable)
{
  register int i , j;
  int *index ;
  int *psequence;
  Histogram **recurrence_time;

  // creation des histogrammes

  recurrence_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    recurrence_time[i] = new Histogram(max_length);
  }

  // mise a jour des histogrammes

  index = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      index[j] = I_DEFAULT;
    }
    psequence = int_sequence[i][variable];

    for (j = 0;j < length[i];j++) {
      if (index[*psequence] != I_DEFAULT) {
        (recurrence_time[*psequence]->frequency[j - index[*psequence]])++;
      }
      index[*psequence++] = j;
    }
  }

  delete [] index;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    recurrence_time[i]->nb_value_computation();
    recurrence_time[i]->offset_computation();
    recurrence_time[i]->nb_element_computation();
    recurrence_time[i]->max_computation();
    recurrence_time[i]->mean_computation();
    recurrence_time[i]->variance_computation();
  }

  characteristics[variable]->recurrence_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->recurrence_time[i] = new Histogram(*(recurrence_time[i]));
    delete recurrence_time[i];
  }
  delete [] recurrence_time;

  psequence = NULL;
  delete psequence;

}


/*--------------------------------------------------------------
 *
 *  Accumulation du temps de sejour dans une valeur
 * (pour une variable donnee).
 *
 * arguments : numéro de variable
 *
 *--------------------------------------------------------------*/

void Switching_sequence::sojourn_time_histogram_computation(int variable)
{
  register int i , j;
  int run_length , *pfrequency , *psequence;

  // initialisation des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i]->offset = 1;
    characteristics[variable]->sojourn_time[i]->nb_value = characteristics[variable]->sojourn_time[i]->alloc_nb_value;

    pfrequency = characteristics[variable]->sojourn_time[i]->frequency;
    for (j = 0;j < characteristics[variable]->sojourn_time[i]->nb_value;j++) {
      *pfrequency++ = 0;
    }

    if (characteristics[variable]->initial_run) {
      characteristics[variable]->initial_run[i]->offset = 1;
      characteristics[variable]->initial_run[i]->nb_value = characteristics[variable]->initial_run[i]->alloc_nb_value;

      pfrequency = characteristics[variable]->initial_run[i]->frequency;
      for (j = 0;j < characteristics[variable]->initial_run[i]->nb_value;j++) {
        *pfrequency++ = 0;
      }
    }

    characteristics[variable]->final_run[i]->offset = 1;
    characteristics[variable]->final_run[i]->nb_value = characteristics[variable]->final_run[i]->alloc_nb_value;

    pfrequency = characteristics[variable]->final_run[i]->frequency;
    for (j = 0;j < characteristics[variable]->final_run[i]->nb_value;j++) {
      *pfrequency++ = 0;
    }
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    psequence = int_sequence[i][variable];
    run_length = 1;
    for (j = 1;j < length[i];j++) {
      if (*(psequence + 1) != *psequence) {
        if ((characteristics[variable]->initial_run) && (run_length == j)) {
          (characteristics[variable]->initial_run[*psequence]->frequency[run_length])++;
        }
        else {
          (characteristics[variable]->sojourn_time[*psequence]->frequency[run_length])++;
        }
        run_length = 0;
      }

      run_length++;
      psequence++;
    }

    if ((characteristics[variable]->initial_run) && (run_length == length[i])) {
      (characteristics[variable]->initial_run[*psequence]->frequency[run_length])++;
    }
    (characteristics[variable]->final_run[*psequence]->frequency[run_length])++;
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i]->nb_value_computation();
    characteristics[variable]->sojourn_time[i]->offset_computation();
    characteristics[variable]->sojourn_time[i]->nb_element_computation();
    characteristics[variable]->sojourn_time[i]->max_computation();
    characteristics[variable]->sojourn_time[i]->mean_computation();
    characteristics[variable]->sojourn_time[i]->variance_computation();
  }

  if (characteristics[variable]->initial_run) {
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i]->nb_value_computation();
      characteristics[variable]->initial_run[i]->offset_computation();
      characteristics[variable]->initial_run[i]->nb_element_computation();
      characteristics[variable]->initial_run[i]->max_computation();
      characteristics[variable]->initial_run[i]->mean_computation();
      characteristics[variable]->initial_run[i]->variance_computation();
    }
  }

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->final_run[i]->nb_value_computation();
    characteristics[variable]->final_run[i]->offset_computation();
    characteristics[variable]->final_run[i]->nb_element_computation();
    characteristics[variable]->final_run[i]->max_computation();
    characteristics[variable]->final_run[i]->mean_computation();
    characteristics[variable]->final_run[i]->variance_computation();
  }

  pfrequency = NULL;
  psequence = NULL;

  delete pfrequency;
  delete psequence;


}


/*--------------------------------------------------------------
 *
 * Construction des histogrammes du temps de sejour dans une valeur
 * (pour une variable donnee).
 *
 * arguments : numéro de variable, flag d'initialisation
 *
 *--------------------------------------------------------------*/

void Switching_sequence::build_sojourn_time_histogram(int variable , int initial_run_flag)
{
  register int i , j;
  int run_length , *psequence;
  Histogram **sojourn_time , **initial_run , **final_run;

  // creation des histogrammes

  sojourn_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    sojourn_time[i] = new Histogram(max_length + 1);
  }

  if (initial_run_flag) {
    initial_run = new Histogram*[marginal[variable]->nb_value];
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      initial_run[i] = new Histogram(max_length + 1);
    }
  }

  final_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    final_run[i] = new Histogram(max_length + 1);
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    psequence = int_sequence[i][variable];
    run_length = 1;
    for (j = 1;j < length[i];j++) {
      if (*(psequence + 1) != *psequence) {
        if ((initial_run_flag) && (run_length == j)) {
          (initial_run[*psequence]->frequency[run_length])++;
        }
        else {
          (sojourn_time[*psequence]->frequency[run_length])++;
        }
        run_length = 0;
      }

      run_length++;
      psequence++;
    }

    if ((initial_run_flag) && (run_length == length[i])) {
      (initial_run[*psequence]->frequency[run_length])++;
    }
    (final_run[*psequence]->frequency[run_length])++;
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    sojourn_time[i]->nb_value_computation();
    sojourn_time[i]->offset_computation();
    sojourn_time[i]->nb_element_computation();
    sojourn_time[i]->max_computation();
    sojourn_time[i]->mean_computation();
    sojourn_time[i]->variance_computation();
  }

  if (initial_run_flag) {
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      initial_run[i]->nb_value_computation();
      initial_run[i]->offset_computation();
      initial_run[i]->nb_element_computation();
      initial_run[i]->max_computation();
      initial_run[i]->mean_computation();
      initial_run[i]->variance_computation();
    }
  }

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    final_run[i]->nb_value_computation();
    final_run[i]->offset_computation();
    final_run[i]->nb_element_computation();
    final_run[i]->max_computation();
    final_run[i]->mean_computation();
    final_run[i]->variance_computation();
  }

  characteristics[variable]->sojourn_time = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->sojourn_time[i] = new Histogram(*(sojourn_time[i]));
    delete sojourn_time[i];
  }
  delete [] sojourn_time;

  if (initial_run_flag) {
    characteristics[variable]->initial_run = new Histogram*[marginal[variable]->nb_value];
    for (i = 0;i < marginal[variable]->nb_value;i++) {
      characteristics[variable]->initial_run[i] = new Histogram(*(initial_run[i]));
      delete initial_run[i];
    }
    delete [] initial_run;
  }

  characteristics[variable]->final_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->final_run[i] = new Histogram(*(final_run[i]));
    delete final_run[i];
  }
  delete [] final_run;

  psequence = NULL;
  delete psequence;

}


/*--------------------------------------------------------------
 *
 *  Construction des histogrammes du nombre de series d'une valeur par sequence
 * (pour une variable donnee).
 *
 * arguments : numéro de variable
 *
 *--------------------------------------------------------------*/

void Switching_sequence::build_nb_run_histogram(int variable)
{
  register int i , j;
  int *psequence , *count;
  Histogram **nb_run;

  // creation des histogrammes

  nb_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_run[i] = new Histogram((max_length % 2 == 0 ?
                               max_length / 2 : max_length / 2 + 1) + 1);
  }

  // mise a jour des histogrammes

  count = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      count[j] = 0;
    }

    psequence = int_sequence[i][variable];
    count[*psequence++]++;
    for (j = 1;j < length[i];j++) {
      if (*psequence != *(psequence - 1)) {
        count[*psequence]++;
      }
      psequence++;
    }

    for (j = 0;j < marginal[variable]->nb_value;j++) {
      (nb_run[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_run[i]->nb_value_computation();
    nb_run[i]->offset_computation();
    nb_run[i]->nb_element = nb_sequence;
    nb_run[i]->max_computation();
    nb_run[i]->mean_computation();
    nb_run[i]->variance_computation();
  }

  characteristics[variable]->nb_run = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_run[i] = new Histogram(*(nb_run[i]));
    delete nb_run[i];
  }
  delete [] nb_run;

  psequence = NULL;
  delete psequence;


}


/*--------------------------------------------------------------
 *
 * Construction des histogrammes du nombre d'occurrences
 * d'une valeur par sequence (pour une variable donnee).
 *
 * arguments : numéro de variable
 *
 *--------------------------------------------------------------*/

void Switching_sequence::build_nb_occurrence_histogram(int variable)
{
  register int i , j;
  int *psequence , *count;
  Histogram **nb_occurrence;

  // creation des histogrammes

  nb_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_occurrence[i] = new Histogram(max_length + 1);
  }

  // mise a jour des histogrammes

  count = new int[marginal[variable]->nb_value];

  for (i = 0;i < nb_sequence;i++) {
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      count[j] = 0;
    }

    psequence = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      count[*psequence++]++;
    }

    for (j = 0;j < marginal[variable]->nb_value;j++) {
      (nb_occurrence[j]->frequency[count[j]])++;
    }
  }

  delete [] count;

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[variable]->nb_value;i++) {
    nb_occurrence[i]->nb_value_computation();
    nb_occurrence[i]->offset_computation();
    nb_occurrence[i]->nb_element = nb_sequence;
    nb_occurrence[i]->max_computation();
    nb_occurrence[i]->mean_computation();
    nb_occurrence[i]->variance_computation();
  }

  characteristics[variable]->nb_occurrence = new Histogram*[marginal[variable]->nb_value];
  for (i = 0;i < marginal[variable]->nb_value;i++) {
    characteristics[variable]->nb_occurrence[i] = new Histogram(*(nb_occurrence[i]));
    delete nb_occurrence[i];
  }
  delete [] nb_occurrence;

  psequence = NULL;
  delete psequence;

}


/*--------------------------------------------------------------
 *
 * Extraction des caracteristiques d'un echantillon de sequences.
 * Calcul du nombre de valeurs et de l'histogramme des valeurs,
 * extraction des probabilites des valeurs en fonction de l'index,
 * construction des histogrammes du temps avant la premiere occurrence d'une valeur,
 * des histogrammes du temps de retour dans une valeur,
 * des histogrammes du temps de sejour dans une valeur,
 * des histogrammes du nombre de series d'une valeur par sequence et
 * des histogrammes du nombre d'occurrences d'une valeur par sequence.
 *
 * arguments : numéro de variable, flag temps de séjour, flag
 *             initialisation
 *--------------------------------------------------------------*/

void Switching_sequence::build_characteristic(int variable , bool sojourn_time_flag ,
					      bool initial_run_flag)
{
  register int i , j;
  bool build;

  for (i = 0;i < nb_variable;i++) {
    if ((variable == I_DEFAULT) || (i == variable)) {
      build = true;

      if (marginal[i]->nb_value > NB_OUTPUT) {
        build = false;
      }

      else {
        for (j = 0;j < marginal[i]->nb_value;j++) {
          if (marginal[i]->frequency[j] == 0) {
            build = false;
            break;
          }
        }
      }

      if (build) {
        if (sojourn_time_flag) {
          characteristics[i] = new Sequence_characteristics(marginal[i]->nb_value);
        }

        build_index_value(i);

	build_first_occurrence_histogram(i);
	build_recurrence_time_histogram(i);


        switch (sojourn_time_flag) {
        case false :
	  sojourn_time_histogram_computation(i);
          break;
        case true :
	  build_sojourn_time_histogram(i , initial_run_flag);
          break;
        }

        if (max_length <= COUNT_MAX_LENGTH) {
	  build_nb_run_histogram(i);
	  build_nb_occurrence_histogram(i);
        }
      }
    }
  }
}


/*--------------------------------------------------------------
 *
 * Ecriture d'un objet Switching_sequence.
 *
 * arguments : stream, flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Switching_sequence::ascii_write(ostream &os , bool exhaustive , bool comment_flag) const
{
  register int i;

  os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;
  os << nb_covariable << " " << "COVARIABLE(S)" << endl;
  if (nb_random > 0){
    os << nb_random << " " << "RANDOM(S)" << endl;
  }
  if(effect){
    os << "Inter individual heterogeneity effect" <<endl;
  }
  if(year_effect){
    os << "Common environment effect" << endl;
  }

  if (!constant) {
    os << "WITH_CONSTANT" << " " << endl;
  }
  else {
    os << "WITHOUT_CONSTANT" << " " << endl;
  }

  if ((self_transition) && (exhaustive)) {
    for (i = 0;i < marginal[0]->nb_value;i++) {
      if (self_transition[i]) {
        os << "\n";
        if (comment_flag) {
          os << "# ";
        }
        os << "   | " << STAT_label[STATL_STATE] << " " << i << " - "
           << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;

        self_transition[i]->ascii_print(os , comment_flag);
      }
    }
  }

  for (i = 0;i < nb_variable;i++) {
    os << "\n" << STAT_word[STATW_VARIABLE] << " " << i + 1 << " : "
       << STAT_variable_word[type[i]] << "   ";
    if (comment_flag) {
      os << "# ";
    }
    if (marginal[i])  {
      os << marginal[i]->nb_value << " ";
      
      switch (type[i]) {
      case STATE :
	os << STAT_label[marginal[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
	break;
      case INT_VALUE :
	os << STAT_label[marginal[i]->nb_value == 1 ? STATL_VALUE : STATL_VALUES] << endl;
	break;
      }

      os << "\n";
      if (comment_flag) {
	os << "# ";
      }
      os << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
	 << STAT_label[STATL_HISTOGRAM] << " - ";
      marginal[i]->ascii_characteristic_print(os , false , comment_flag);
      
      if ((marginal[i]->nb_value <= ASCII_NB_VALUE) || (exhaustive)) {
	os << "\n";
	if (comment_flag) {
	  os << "# ";
	}
	os << "   | " << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
	   << STAT_label[STATL_HISTOGRAM] << endl;
	marginal[i]->ascii_print(os , comment_flag);
      }
    }
 
    if (characteristics[i]) {
      characteristics[i]->ascii_print(os , type[i] , *hlength , exhaustive , comment_flag);
    }
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
  hlength->ascii_characteristic_print(os , false , comment_flag);

  if (exhaustive) {
    os << "\n";
    if (comment_flag) {
      os << "# ";
    }
    os << "   | " << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hlength->ascii_print(os , comment_flag);
  }

  os << "\n";
  if (comment_flag) {
    os << "# ";
  }
  os << SEQ_label[SEQL_CUMUL_LENGTH] << ": " << cumul_length << endl;

  return os;
}


/*--------------------------------------------------------------
 *
 * Ecriture d'un objet Switching_sequence.
 *
 * arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Switching_sequence::ascii_write(ostream &os , bool exhaustive) const
{
  return ascii_write(os , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Switching_sequence dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Switching_sequence::ascii_write(Format_error &error , const char *path ,
				     bool exhaustive) const

{
  bool status;
  ofstream out_file(path);


  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    ascii_write(out_file , exhaustive , false);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Switching_sequence.
 *
 *  arguments : stream, format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Switching_sequence::ascii_data_write(ostream &os , char format , bool exhaustive) const
{
  ascii_write(os , exhaustive , false);
  ascii_print(os , format , false);

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Switching_sequence dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              format (ligne/colonne), flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Switching_sequence::ascii_data_write(Format_error &error , const char *path ,
					  char format , bool exhaustive) const
{
  bool status = false;
  ofstream out_file(path);

  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    if (format != 'a') {
      ascii_write(out_file , exhaustive , true);
    }
    ascii_print(out_file , format , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 * Ecriture d'un objet Switching_sequence dans un fichier au 
 * format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Switching_sequence::spreadsheet_write(Format_error &error , const char *path) const
{
  bool status;
  register int i;
  Curves *smoothed_curves;
  ofstream out_file(path);

  error.init();

  if (!out_file) {
    status = false;
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;

    out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;
    out_file << nb_covariable << "\t" << "COVARIABLE(S)" << endl;
    if (!constant) {
      out_file << "WITH_CONSTANT" << "\t" << endl;
    }
    else {
      out_file << "WITHOUT_CONSTANT" << "\t" << endl;
    }
    if (self_transition) {
      for (i = 0;i < marginal[0]->nb_value;i++) {
        if (self_transition[i]) {
          out_file << "\n\t" << STAT_label[STATL_STATE] << " " << i << " - "
                   << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;
          self_transition[i]->spreadsheet_print(out_file);
	  
          smoothed_curves = new Curves(*(self_transition[i]) , 's');
	  
          out_file << "\n\t" << STAT_label[STATL_STATE] << " " << i << " - "
                   << SEQ_label[SEQL_SMOOTHED] << " " << SEQ_label[SEQL_OBSERVED] << " " << SEQ_label[SEQL_SELF_TRANSITION] << endl;
          smoothed_curves->spreadsheet_print(out_file);
	  
          delete smoothed_curves;
        }
      }
    }
    
    for (i = 0;i < nb_variable;i++) {
      if(marginal[i]) {
	out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << i + 1 << "\t\t"
		 << marginal[i]->nb_value << " ";
	
	switch (type[i]) {
	case STATE :
	  out_file << STAT_label[marginal[i]->nb_value == 1 ? STATL_STATE : STATL_STATES] << endl;
	  break;
	case INT_VALUE :
	  out_file << STAT_label[marginal[i]->nb_value == 1 ? STATL_VALUE : STATL_VALUES] << endl;
	  break;
	}
	
	out_file << "\n" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
		 << STAT_label[STATL_HISTOGRAM] << "\t";
	marginal[i]->spreadsheet_characteristic_print(out_file);
	
	out_file << "\n\t" << STAT_label[type[i] == STATE ? STATL_STATE : STATL_VALUE] << " "
		 << STAT_label[STATL_HISTOGRAM] << endl;
	marginal[i]->spreadsheet_print(out_file);
      }
      if (characteristics[i]) {
	characteristics[i]->spreadsheet_print(out_file , type[i] , *hlength);
      }
    }
    
    out_file << "\n" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    hlength->spreadsheet_characteristic_print(out_file);
    
    out_file << "\n\t" << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_HISTOGRAM] << endl;
    hlength->spreadsheet_print(out_file);
    
    out_file << "\n" << SEQ_label[SEQL_CUMUL_LENGTH] << "\t" << cumul_length << endl;
  }

  return status;
}


/*--------------------------------------------------------------
 *
 * Accumulation des observations (pour une variable donnee).
 *
 * arguments : numéro de variable
 *
 *--------------------------------------------------------------*/

void Switching_sequence::observation_histogram_computation(int variable)
{
  register int i , j;
  double *pfrequency;
  int *pstate ;
  int *poutput;// seulement pour les entiers 

  // initialisation des histogrammes

  for (i = 0;i < marginal[0]->nb_value;i++) {
    pfrequency = observation[variable][i]->frequency;
    for (j = 0;j < marginal[variable]->nb_value;j++) {
      *pfrequency++ = 0;
    }
  }

  // mise a jour des histogrammes

  for (i = 0;i < nb_sequence;i++) {
    pstate = int_sequence[i][0];
    poutput = int_sequence[i][variable];
    for (j = 0;j < length[i];j++) {
      (observation[variable][*pstate++]->frequency[*poutput++])++;
    }
  }

  // extraction des caracteristiques des histogrammes

  for (i = 0;i < marginal[0]->nb_value;i++) {
    if (!characteristics[variable]) {
      observation[variable][i]->nb_class_computation();
    }
    observation[variable][i]->nb_element_computation();

    if (!characteristics[variable]) {
      observation[variable][i]->mean_computation();
      observation[variable][i]->variance_computation();
    }
  }

  pfrequency = NULL;
  pstate = NULL;
  poutput = NULL;
  delete pfrequency;
  delete pstate;
  delete poutput;

}


/*--------------------------------------------------------------
 *
 * Accumulation des observations.
 *
 *--------------------------------------------------------------*/

void Switching_sequence::observation_histogram_computation()
{
  register int i;

  for (i = 1;i < nb_variable;i++) {
    observation_histogram_computation(i);
  }
}


/*--------------------------------------------------------------
 *
 * Construction des histogrammes correspondant aux probabilites 
 * d'observation.
 *
 * arguments : nombre d'états
 *
 *--------------------------------------------------------------*/

void Switching_sequence::create_observation_histogram(int nb_state)
{
  if ((nb_variable > 1) && (!observation)) {
    register int i , j;

    observation = new Continuous_histo**[nb_variable];
    observation[0] = 0;

    for (i = 1;i < nb_variable;i++) {
      observation[i] = new Continuous_histo*[nb_state];
      for (j = 0;j < nb_state;j++) {
        observation[i][j] = new Continuous_histo(marginal[i]->nb_value);
      }
    }
  }
}


/*--------------------------------------------------------------
 *
 * Construction des histogrammes correspondant aux probabilites
 * d'observation.
 *
 *--------------------------------------------------------------*/

void Switching_sequence::build_observation_histogram()
{
  create_observation_histogram(marginal[0]->nb_value); 
  observation_histogram_computation();
}


/*--------------------------------------------------------------
 *
 * Modification des effets aléatoires - heterogeneité
 *
 * arguments : effets aléatoires individuels
 *
 *--------------------------------------------------------------*/

void Switching_sequence::Set_effect(double **ieffect)
{
  register int i, j, k;
  double *peffect, *ceffect;

  effect = new double*[nb_sequence];
  for (i = 0;i < nb_sequence;i++) {
    effect[i] = new double[nb_val];
    
    peffect = effect[i];
    ceffect = ieffect[i];
    for (j = 0;j < nb_random;j++) {
      *peffect++ = *ceffect++;
    }
  }
  
}


/*--------------------------------------------------------------
 *
 * Modification des effets aléatoires - effet année
 *
 * arguments : effets aléatoires temporels
 *
 *--------------------------------------------------------------*/

void Switching_sequence::Set_year_effect(double *iyear_effect)
{
  register int i;
  int nb_annee;

  double *pyear_effect, *cyear_effect;

  nb_annee = nb_year();
  year_effect = new double[nb_annee];
  pyear_effect = year_effect;
  cyear_effect = iyear_effect;

  for (i = 0;i < nb_annee; i++) {
    *pyear_effect++ = *cyear_effect++;
  }
  
}




