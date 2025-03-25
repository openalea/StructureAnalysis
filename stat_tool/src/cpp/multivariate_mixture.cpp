/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr) and
 *                       Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: multivariate_mixture.cpp 5353 2008-07-29 12:33:46Z guedon $
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



#include <assert.h>

#include <sstream>
#include <iomanip>
#include <fstream>
//#include <string>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "stat_tools.h"
#include "distribution.h"
#include "discrete_mixture.h"
#include "markovian.h"
#include "vectors.h"
#include "stat_label.h"
// #include "mixture.h"
#include "multivariate_mixture.h"

using namespace std;
using namespace boost;


namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe MultivariateMixture.
 *
 *--------------------------------------------------------------*/

MultivariateMixture::MultivariateMixture()

{
  mixture_data = NULL;
  nb_component = 0;
  nb_var = 0;
  weight = NULL;
  pcomponent = NULL;
  npcomponent = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe MultivariateMixture.
 *
 *  arguments : nombre de composantes, poids, pointeurs sur les composantes parametriques
 *  et categorielles
 *
 *--------------------------------------------------------------*/

MultivariateMixture::MultivariateMixture(int inb_component , double *pweight , int inb_variable,
                                         DiscreteParametricProcess **ppcomponent,
                                         CategoricalProcess **pnpcomponent)

{
  register int var, i;

  mixture_data = NULL;
  nb_component = inb_component;
  nb_var = inb_variable;

  weight = new DiscreteParametric(nb_component);
  for (i = 0;i < nb_component;i++) {
    weight->mass[i] = *pweight++;
  }
  weight->cumul_computation();
  weight->max_computation();

  if (weight->ident != CATEGORICAL) {
    weight->computation(1 , CUMUL_THRESHOLD);
  }

  pcomponent = new DiscreteParametricProcess*[nb_var];
  npcomponent = new CategoricalProcess*[nb_var];

  for (var = 0;var < nb_var;var++) {
    if ((pnpcomponent != NULL) && (pnpcomponent[var] != NULL))
      {
         npcomponent[var]= new CategoricalProcess(*(pnpcomponent[var]));
         pcomponent[var]= NULL;
      }
      else
      {
         npcomponent[var]= NULL;
         pcomponent[var]= new DiscreteParametricProcess(*(ppcomponent[var]));
      }
   }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe MultivariateMixture.
 *
 *  arguments : reference sur un objet MultivariateMixture, flag sur les variables a copier,
 *              nombre de variables a copier.
 *
 *--------------------------------------------------------------*/

MultivariateMixture::MultivariateMixture(const MultivariateMixture &mixt ,
                                         bool *variable_flag, int inb_variable)

{
  register int var;


  mixture_data = NULL;
  nb_component = mixt.nb_component;

  weight = new DiscreteParametric(nb_component);

  nb_var = inb_variable;
  pcomponent = new DiscreteParametricProcess*[nb_var];
  npcomponent = new CategoricalProcess*[nb_var];

  for (var = 0;var < nb_var;var++) {
    if (variable_flag[var]) {
      if (mixt.pcomponent[var] != NULL) {
    pcomponent[var] = new DiscreteParametricProcess(*mixt.pcomponent[var], true);
    npcomponent[var] = NULL;
      }
      else {
    pcomponent[var] = NULL;
    npcomponent[var] = new CategoricalProcess(*mixt.npcomponent[var]);
      }
    }
  }
}

/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe MultivariateMixture.
 *
 *  arguments : nombre de composantes, nombre de variables,
 *  pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

MultivariateMixture::MultivariateMixture(int inb_component , int inb_variable,
                                         const DiscreteParametricProcess **ppcomponent,
                                         const CategoricalProcess **pnpcomponent)

{
  register int var;


  mixture_data = NULL;
  nb_component = inb_component;
  nb_var = inb_variable;

  weight = NULL;

  pcomponent = new DiscreteParametricProcess*[nb_var];
  npcomponent = new CategoricalProcess*[nb_var];

  for (var = 0;var < nb_var;var++) {
    if (ppcomponent[var] != NULL) {
    pcomponent[var] = new DiscreteParametricProcess(*ppcomponent[var], true);
    npcomponent[var] = NULL;
    }
    else {
       pcomponent[var] = NULL;
       npcomponent[var] = new CategoricalProcess(*pnpcomponent[var]);
    }
  }
}

/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe MultivariateMixture.
 *
 *  arguments : nombre de composantes, nombre de variables,
 *  nombre de valeurs pour chaque variable,
 *  choix de lois parametriques ou non
 *
 *--------------------------------------------------------------*/

MultivariateMixture::MultivariateMixture(int inb_component, int inb_variable,
                                         int *nb_value, bool *force_param) {


  register int var, i;
  double param;
  bool *fparam = NULL;

  DiscreteParametric *rand = new DiscreteParametric(UNIFORM, 0, 10 , D_DEFAULT , D_DEFAULT);

  mixture_data = NULL;
  nb_component = inb_component;
  nb_var = inb_variable;

  weight = NULL;

  pcomponent = new DiscreteParametricProcess*[nb_var];
  npcomponent = new CategoricalProcess*[nb_var];

  fparam= new bool[nb_var];
  if (force_param == NULL) {
    for (var = 0; var < nb_var; var++)
      fparam[var]= false;
  }
  else {
    for (var = 0; var < nb_var; var++)
      fparam[var]= force_param[var];
  }

  npcomponent= new CategoricalProcess*[nb_var];
  pcomponent= new DiscreteParametricProcess*[nb_var];

  for(var = 0; var < nb_var; var++) {
    if ((*nb_value <= NB_OUTPUT) && !(fparam[var])) {
      npcomponent[var] = new CategoricalProcess(nb_component, *nb_value++, true);
      pcomponent[var] = NULL;
    }
    else {
      npcomponent[var] = NULL;
      pcomponent[var] = new DiscreteParametricProcess(nb_component, (int)(*nb_value * SAMPLE_NB_VALUE_COEFF));
      for(i = 0; i < nb_component; i++) {
    delete pcomponent[var]->observation[i];
    param = (cumul_method(10, rand->cumul, 1.) + 1);
    // TODO: change default parameter values
    pcomponent[var]->observation[i] =
      new DiscreteParametric(NEGATIVE_BINOMIAL, 0, I_DEFAULT , 1., 1. / (double)((param * *nb_value)+1.));
      }
      nb_value++;
    }
  }
  delete [] fparam;
  fparam= NULL;
}

/*--------------------------------------------------------------*
 *
 *  Copie d'un objet MultivariateMixture.
 *
 *  arguments : reference sur un objet MultivariateMixture,
 *              flag copie de l'objet MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

void MultivariateMixture::copy(const MultivariateMixture &mixt , bool data_flag)

{
  register int var;


  if ((data_flag) && (mixt.mixture_data != NULL)) {
    mixture_data = new MultivariateMixtureData(*(mixt.mixture_data) , false);
  }
  else {
    mixture_data = NULL;
  }

  nb_component = mixt.nb_component;
  nb_var = mixt.nb_var;

  pcomponent = new DiscreteParametricProcess*[nb_var];
  npcomponent = new CategoricalProcess*[nb_var];

  weight = new DiscreteParametric(*(mixt.weight));

  for (var = 0;var < nb_var;var++) {
    if (mixt.pcomponent[var] != NULL) {
    pcomponent[var] = new DiscreteParametricProcess(*mixt.pcomponent[var], true);
    npcomponent[var] = NULL;
    }
    else {
       pcomponent[var] = NULL;
       npcomponent[var] = new CategoricalProcess(*mixt.npcomponent[var]);
    }
  }
  // init();
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet MultivariateMixture.
 *
 *--------------------------------------------------------------*/

void MultivariateMixture::remove()

{
  if (mixture_data != NULL) {
    delete mixture_data;
    mixture_data= NULL;
  }

  if (weight != NULL) {
    delete weight;
    weight = NULL;
  }

  if ((pcomponent != NULL) || (npcomponent != NULL)) {
    register int var;

    for (var = 0;var < nb_var;var++) {
      if (pcomponent[var] != NULL) {
    delete pcomponent[var];
    pcomponent[var] = NULL;
      }
      if (npcomponent[var] != NULL) {
    delete npcomponent[var];
    npcomponent[var] = NULL;
      }
    }
    delete [] pcomponent;
    pcomponent = NULL;
    delete [] npcomponent;
    npcomponent = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe MultivariateMixture.
 *
 *--------------------------------------------------------------*/

MultivariateMixture::~MultivariateMixture()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe MultivariateMixture.
 *
 *  argument : reference sur un objet MultivariateMixture.
 *
 *--------------------------------------------------------------*/

MultivariateMixture& MultivariateMixture::operator=(const MultivariateMixture &mixt)

{
  if (&mixt != this) {
    remove();
    copy(mixt);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une composante parametrique.
 *
 *  arguments : reference sur un objet StatError, variable,
 *  indice de la composante.
 *
 *--------------------------------------------------------------*/

DiscreteParametricModel* MultivariateMixture::extract_parametric_model(StatError &error ,
                                                                       int ivariable,
                                                                       int index) const

{
  bool status = true;
  DiscreteParametricModel *ppcomponent = NULL;
  FrequencyDistribution *hcomponent = NULL;

  if ((index < 1) || (index > nb_component)) {
    status = false;
    error.update(STAT_error[STATR_DISTRIBUTION_INDEX]);
  }

  if ((ivariable < 1) || (ivariable > nb_var)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  if (status) {
    if (pcomponent[ivariable] != NULL){
      index--;
      if (mixture_data != NULL) {
    if (mixture_data->type[0] == STATE) {
      if (mixture_data->component[ivariable+1] != NULL)
        hcomponent = mixture_data->component[ivariable+1][index];
    }
    else {
      if (mixture_data->component[ivariable] != NULL)
        hcomponent = mixture_data->component[ivariable][index];
    }
      }
      ppcomponent = new DiscreteParametricModel(*pcomponent[ivariable]->observation[index] ,
                     hcomponent);
    }
  }

  return ppcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Extraction d'une composante categorielle.
 *
 *  arguments : reference sur un objet StatError, variable,
 *  indice de la composante.
 *
 *--------------------------------------------------------------*/

Distribution* MultivariateMixture::extract_categorical_model(StatError &error ,
                                                               int ivariable, int index) const

{
  bool status = true;
  Distribution *pnpcomponent = NULL;


  if ((index < 1) || (index > nb_component)) {
    status = false;
    error.update(STAT_error[STATR_DISTRIBUTION_INDEX]);
  }

  if ((ivariable < 1) || (ivariable > nb_var)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  if (status) {
    if (npcomponent[ivariable] != NULL){
      index--;
      pnpcomponent = new Distribution(*npcomponent[ivariable]->observation[index]);
    }
  }

  return pnpcomponent;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi marginale du melange
 *
 *  argument : reference sur un objet StatError, index de la variable
 *
 *--------------------------------------------------------------*/

DiscreteMixture* MultivariateMixture::extract_distribution(StatError &error , int ivariable) const
{
  bool status = true, discrete_mixture = false;
  register int i , j, variable = ivariable - 1,  nb_value;
  double *pweight , *pmass;
  const DiscreteParametric **marginal_pcomponent = NULL;
  DiscreteMixture *pDistribution = NULL;
  DiscreteMixtureData *data = NULL;
  FrequencyDistribution *obs_marginal = NULL;

  if ((ivariable < 1) || (ivariable > nb_var)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  if (status) {
    marginal_pcomponent = new const DiscreteParametric*[nb_component];
    pweight = weight->mass;

    if (pcomponent[variable] != NULL) { // parametric distribution

      for (i = 0;i < nb_component;i++)
    	  marginal_pcomponent[i] = pcomponent[variable]->observation[i];
    }
    else { // npcomponent[variable] != NULL)
    	for (i = 0;i < nb_component;i++) {
    		nb_value = npcomponent[variable]->observation[i]->nb_value;
    		marginal_pcomponent[i] = new DiscreteParametric(nb_value, CATEGORICAL);
    		for (j = 0; j < nb_value; j++)
    			marginal_pcomponent[i]->mass[j] = npcomponent[variable]->observation[i]->mass[j];
    	}
    }
    pDistribution = new DiscreteMixture(nb_component, pweight, marginal_pcomponent);
	if (mixture_data != NULL) {
		obs_marginal = mixture_data->get_marginal_distribution(variable);
		data = new DiscreteMixtureData(*obs_marginal, pDistribution);

		pDistribution->mixture_data = data;
    }
    if (npcomponent[variable] != NULL)
    	for (i = 0;i < nb_component;i++)
    	    delete marginal_pcomponent[i];

    delete [] marginal_pcomponent;

    pDistribution->computation();
    pDistribution->cumul_computation();

    pDistribution->max_computation();
    pDistribution->mean_computation();
    pDistribution->variance_computation();
  }
  return pDistribution;
}

/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet MultivariateMixture.
 *
 *  argument : reference sur un objet StatError.
 *
 *--------------------------------------------------------------*/

MultivariateMixtureData* MultivariateMixture::extract_data(StatError &error) const

{
  MultivariateMixtureData *mixt_data;


  error.init();

  if (!mixture_data) {
    mixt_data = NULL;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  else {
    mixt_data = new MultivariateMixtureData(*mixture_data);
    mixt_data->mixture = new MultivariateMixture(*this , false);
  }

  return mixt_data;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet MultivariateMixture a partir de poids et de composantes.
 *
 *  arguments : reference sur un objet StatError, nombre de composantes,
 *              poids, pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

MultivariateMixture* multivariate_mixture_building(StatError &error , int nb_component ,
                                                   int nb_variable , double *weight,
                                                   DiscreteParametricProcess **ppcomponent,
                                                   CategoricalProcess **pnpcomponent)

{
  bool status;
  register int i;
  double cumul;
  MultivariateMixture *mixt;


  mixt = NULL;
  error.init();

  if ((nb_component < 2) || (nb_component > DISCRETE_MIXTURE_NB_COMPONENT)) {
    error.update(STAT_parsing[STATP_NB_DISTRIBUTION]);
  }

  else {
    status = true;
    cumul = 0.;

    for (i = 0;i < nb_component;i++) {
      if ((weight[i] <= 0.) || (weight[i] > 1. - cumul + DOUBLE_ERROR)) {
        status = false;
        cerr << "i="<<i << " and weight=" << weight[i]<< " cumul=" << cumul<<endl;
        std::flush(cerr);
        error.update(STAT_parsing[STATP_WEIGHT_VALUE]);
      }
      else {
        cumul += weight[i];
      }
    }

    if (status) {
      mixt = new MultivariateMixture(nb_component , weight , nb_variable ,
                                     ppcomponent , pnpcomponent);
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet MultivariateMixture a partir d'un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

MultivariateMixture* multivariate_mixture_ascii_read(StatError &error , const std::string & path ,
                                                     double cumul_threshold)

{
  string buffer;
  size_t position;
  typedef tokenizer<char_separator<char>> tokenizer;
  char_separator<char> separator(" \t");
  bool status , lstatus, categorical = false;
  int i , j;
  int line, nb_variable;
  long index , nb_component, value;
  double cumul , *weight = NULL;
  CategoricalProcess **np_observation;
  DiscreteParametricProcess **p_observation;
  MultivariateMixture *mixt;
  ifstream in_file(path.c_str());


  mixt = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    nb_component = 0;

    // read up to next line jump
    while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);

      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        switch (i) {

        // search for keyword MIXTURE

          case 0 : {
            if (*token != STAT_word[STATW_MIXTURE]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_MIXTURE] , line , i + 1);
            }
            break;
          }
        }

        i++;
      }

      if (i > 0) {
        if (i != 1) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      }
    } // while (getline(in_file , buffer))

    // read up to next line jump
  	while (getline(in_file , buffer)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.find('#');
      if (position != string::npos) {
        buffer.erase(position);
      }
      i = 0;

      tokenizer tok_buffer(buffer , separator);
      // separate line into words
      for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
        switch (i) {
        // search for the number of components
          case 0 : {
            nb_component = atoi(token->c_str());
            if ((lstatus) && ((nb_component < 2) || (nb_component > DISCRETE_MIXTURE_NB_COMPONENT))) {
              lstatus = false;
            }

            if (!lstatus) {
              status = false;
              error.update(STAT_parsing[STATP_NB_DISTRIBUTION] , line , i + 1);
            }
            break;
          }

          case 1 : {
            if (*token != STAT_word[STATW_DISTRIBUTIONS]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_KEYWORD] , STAT_word[STATW_DISTRIBUTIONS] , line , i + 1);
            }
            break;
          }
        }

        i++;
      } // for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++)


      if (i > 0) {
        if (i != 2) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      }
    } // while (getline(in_file , buffer))

    if (nb_component == 0) {
      status = false;
      error.update(STAT_parsing[STATP_FORMAT] , line);
    }
    else
      weight = new double[nb_component];

    if (status) {
      nb_variable = I_DEFAULT;

      np_observation= NULL;
      p_observation= NULL;

      // read weights

      index= 0;
      i = 0;

      while (getline(in_file , buffer)) {
        line++;

  #         ifdef DEBUG
        cout << line << "  " << buffer << endl;
  #         endif

        position = buffer.find('#');
        if (position != string::npos)
          buffer.erase(position);

        // i= 0;

        tokenizer tok_buffer(buffer , separator);
        // separate line into words
        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          if (i == 0) {
            if (*token != STAT_word[STATW_WEIGHTS]) {
              status = false;
              error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] ,
                          STAT_word[STATW_WEIGHTS] , line , i + 1);
            }
          }
          else {
            if (index < nb_component) {
                weight[index] = atoi(token->c_str());
                if (lstatus) {
                  if ((weight[index] <= 0.) || (weight[index] > 1. - cumul + DOUBLE_ERROR)) {
                    lstatus = false;
                  }
                  else {
                    cumul += weight[index];
                  }
                }

                if (!lstatus) {
                  status = false;
                  error.update(STAT_parsing[STATP_WEIGHT_VALUE] , line , i + 1);
                }
              index++;
            }
          }
        i++;
        } // for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++)
        if (i >= nb_component) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
      } // end : while (getline(in_file , buffer))

      if (index != nb_component) {
        status= false;
        error.update(STAT_parsing[STATP_FORMAT] , line);
      }

      if (cumul < CUMUL_THRESHOLD) {
        status = false;
        error.update(STAT_parsing[STATP_PROBABILITY_SUM]);
      }

      // read emission distributions
      // read up to next line jump
      while (getline(in_file , buffer)) {
        line++;

  #         ifdef DEBUG
        cout << line << "  " << buffer << endl;
  #         endif

        position = buffer.find('#');
        if (position != string::npos)
          buffer.erase(position);

        i= 0;
        index = 0;

        tokenizer tok_buffer(buffer , separator);
        // separate line into words
        for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
          switch (i) {
            // test for number of variables
              case 0 : {
                  value = atoi(token->c_str());
                  if (lstatus) {
                    if ((value < 1) || (value > NB_OUTPUT_PROCESS))
                        lstatus = false;
                    else
                        nb_variable = value;
                  }

                  if (!lstatus) {
                    status= false;
                    error.update(STAT_parsing[STATP_NB_OUTPUT_PROCESS],
                                  line, i+1);
                  }
                  break;
              }
              // test for keyword VARIABLE(S)
              case 1 : {
                  if (*token != STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES]) {
                    status = false;
                    error.correction_update(STAT_parsing[STATP_KEYWORD] ,
                                            STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES],
                                            line, i+1);
                  }
                  break;
              }
            } // end switch
            i++;
          } // end : for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++)

          if (i > 0) {
              if (i != 2)
              {
                status = false;
                error.update(STAT_parsing[STATP_FORMAT] , line);
              }
              break;
          }
        } // end : while (getline(in_file , buffer))

        if (nb_variable == I_DEFAULT) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }
        else {
          np_observation= new CategoricalProcess*[nb_variable];
          p_observation = new DiscreteParametricProcess*[nb_variable];
          for(i= 0; i < nb_variable; i++)
          {
              np_observation[i]= NULL;
              p_observation[i]= NULL;
          }

          // read up to next line jump
          while (getline(in_file , buffer)) {
            line++;

  #                  ifdef DEBUG
            cout << line << "  " << buffer << endl;
  #                  endif

            position = buffer.find('#');
            if (position != string::npos)
                      buffer.erase(position);

            i= 0;

            tokenizer tok_buffer(buffer , separator);
            // separate line into words
            for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++) {
                switch (i) {

                  // test for keyword VARIABLE
                case 0 : {
                    categorical= true;

                    if (*token == STAT_word[STATW_VARIABLE])
                      index++;
                    else {
                      status= false;
                      error.correction_update(STAT_parsing[STATP_KEYWORD], STAT_word[STATW_VARIABLE],
                                  line, i+1);
                    }
                    break;
                }
                // test for variable index
                case 1 : {
                    value = atoi(token->c_str());
                    if ((lstatus) && ((value != index) || (value > nb_variable)))
                      lstatus= false;

                    if (!lstatus) {
                      status= false;
                      error.update(STAT_parsing[STATP_VARIABLE_INDEX],
                                line, i+1);
                    }
                    break;
                }
                // test on separator
                case 2 : {
                    if (*token != ":") {
                      status= false;
                      error.update(STAT_parsing[STATP_SEPARATOR], line, i+1);
                    }
                    break;
                }
                // test on keywords CATEGORICAL / DISCRETE_PARAMETRIC
                case 3 : {
                    if ((*token == STAT_word[STATW_CATEGORICAL]) ||
                        (*token == STAT_word[STATW_NONPARAMETRIC]))
                      categorical = true;
                    else {
                    if ((*token == STAT_word[STATW_DISCRETE_PARAMETRIC]) ||
                        (*token == STAT_word[STATW_PARAMETRIC]))
                      categorical = false;
                    else {
                        status = false;
                        ostringstream correction_message;
                        correction_message << STAT_word[STATW_CATEGORICAL] << " or "
                                << STAT_word[STATW_DISCRETE_PARAMETRIC] << "(instead of "
                                << *token << ")";
                        error.correction_update(STAT_parsing[STATP_KEYWORD],
                                    (correction_message.str()).c_str(),
                                    line, i+1);
                      }
                      break;
                    }
                }
              } // end switch
              i++;
            } // for (tokenizer::iterator token = tok_buffer.begin();token != tok_buffer.end();token++)

            if (i > 0) {
              if (i != 4) {
                status= false;
                error.update(STAT_parsing[STATP_FORMAT], line);
              }

              switch (categorical) {

                  case true : {
                      np_observation[index-1]= CategoricalProcess::parsing(error, in_file, line,
                                                                                nb_component, MIXTURE, true);
                      if (np_observation[index-1] == NULL)
                        status= false;
                      break;
                  }

                  case false : {
                      p_observation[index-1]= DiscreteParametricProcess::parsing(error, in_file, line,
                                                                                  nb_component, MIXTURE,
                                                                                  cumul_threshold);
                      if (p_observation[index-1] == NULL)
                        status = false;
                      break;
                  }
                }
            }
          } // end : while (getline(in_file , buffer))

          if (index != nb_variable)  {
           status = false;
           error.update(STAT_parsing[STATP_FORMAT] , line);
          }


          if (cumul < CUMUL_THRESHOLD) {
            status = false;
            error.update(STAT_parsing[STATP_PROBABILITY_SUM]);
          }

          // read up to next line jump
          while (getline(in_file , buffer))  {
            line++;

#              ifdef DEBUG
         cout << line << "  " << buffer << endl;
#              endif

            position = buffer.find('#');
            if (position != string::npos) {
              buffer.erase(position);
            }

            if (!(trim_right_copy_if(buffer , is_any_of(" \t")).empty())) {
              status = false;
              error.update(STAT_parsing[STATP_FORMAT] , line);
            }
          } // while (getline(in_file , buffer))

          if (status) {
            mixt = new MultivariateMixture(nb_component , weight , nb_variable ,
                                          p_observation, np_observation);
          }

         if (nb_variable > 0) {
           for(i= 0; i < nb_variable; i++)
             {
               if (np_observation[i] != NULL)
             delete np_observation[i];
               if (p_observation[i] != NULL)
             delete p_observation[i];
             }

           delete [] np_observation;
           delete [] p_observation;
         }
      } // end if (nb_variable == I_DEFAULT)
    } // (status)
  } // (!in_file)
  return mixt;
}

/*****************************************************************
 *
 *  Permutation of the states of a mixture model based on
 *  a given permutation perm (which must be the same length as
 *  the number of states)
 *
 **/

void MultivariateMixture::state_permutation(StatError& error,
                                            int* perm) const
{
   register int i, j;
   bool status = true;
   // each element of check_perm must be used exactly once
   bool * const check_perm = new bool[nb_component];
   double* pmass = NULL;

   // check permutation
   error.init();

   for (i = 0; i < nb_component; i++)
      check_perm[i] = false; // indicates if ith element was used in perm

   for (i = 0; i < nb_component; i++)
   {
      if (check_perm[perm[i]])
         status = false;
      else
         check_perm[perm[i]] = true;
   }

   if (!status)
      error.update(STAT_error[STATR_NO_PERMUTATION]);
   else
   {
      // permutation of weights
      pmass = new double[nb_component];
      for (i = 0; i < nb_component; i++)
         pmass[perm[i]] = weight->mass[i];
      for (i= 0; i < nb_component; i++)
         weight->mass[i] = pmass[i];
      delete [] pmass;
      pmass = NULL;
      weight->cumul_computation();
      weight->max_computation();

      if (weight->ident != CATEGORICAL) {
         weight->computation(1 , CUMUL_THRESHOLD);
      }

      // permutation of observation distributions
      for (i = 0; i < nb_var; i++)
      {
         if (pcomponent[i] != NULL)
            pcomponent[i]->state_permutation(perm);
         if (npcomponent[i] != NULL)
            npcomponent[i]->state_permutation(perm);
      }

      if (mixture_data != NULL)
         mixture_data->state_permutation(perm);
   }
}

/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet MultivariateMixture.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& MultivariateMixture::line_write(ostream &os) const

{
  os << nb_component << " " << STAT_word[STATW_DISTRIBUTIONS] ;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un melange et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet MultivariateMixtureData,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& MultivariateMixture::ascii_write(ostream &os , const MultivariateMixtureData *mixt_data ,
                                          bool exhaustive , bool file_flag) const

{
  register int i, var,
    data_var; // index that corresponds to var in mixt_data
  int bnb_parameter;
  int *var_array = NULL;
  // double *scale;
// const Distribution *ptComponent[DISCRETE_MIXTURE_NB_COMPONENT];
  FrequencyDistribution **observation= NULL;
  StatError error;
  Vectors *vect_data = NULL;

  os << STAT_word[STATW_MIXTURE] << endl << endl;
  os << nb_component << " " << STAT_word[STATW_DISTRIBUTIONS] << endl << endl;
  // ascii_characteristic_print(os , false , file_flag);

  os << STAT_word[STATW_WEIGHTS] << endl;
  for (i = 0;i < nb_component;i++)
    os << weight->mass[i] << "  ";
  os << endl;


  if (nb_var > 0) {
    os << "\n" << nb_var << " "
       << STAT_word[nb_var == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    for(var= 1; var <= nb_var; var++)
      {
    os << "\n" << STAT_word[STATW_VARIABLE];
    os << " " << var;

    if (npcomponent[var-1] != NULL)
      os << " : " << STAT_word[STATW_CATEGORICAL];
    else
      os << " : " << STAT_word[STATW_DISCRETE_PARAMETRIC];

    os << endl;
    if (mixt_data != NULL) {
      if (mixt_data->type[0] == STATE)
        data_var = var;
      else
        data_var = var-1;
      os << "\n";
      if (file_flag) {
        os << "# ";
      }

      if (mixt_data->component != NULL)
        observation= mixt_data->component[data_var];
      else
        observation= NULL;
    }
    if (npcomponent[var-1] == NULL)
      pcomponent[var-1]->ascii_print(os, observation, NULL, exhaustive, file_flag);
    else
      npcomponent[var-1]->ascii_print(os, observation, NULL, exhaustive, file_flag);

    if (exhaustive) {
      /* for (i = 0;i < nb_component;i++) {
        if (pcomponent[var-1] != NULL)
          ptComponent[var-1] = pcomponent[var-1]->observation[i];
        else
          ptComponent[var-1] = npcomponent[var-1]->observation[i];

        if (mixt_data) {
          scale[i] = mixt_data->nb_vector * weight->mass[i];
        }
        else {
          scale[i] = weight->mass[i];
        }
        }*/

      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << "  ";
      if (mixt_data != NULL) {
        os << " | " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
      }
      os << " | " << STAT_label[STATL_MIXTURE];
      for (i = 0;i < nb_component;i++) {
        os << " | " << STAT_label[STATL_DISTRIBUTION] << " " << var;
      }
      if (mixt_data != NULL) {
        os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
           << STAT_label[STATL_FUNCTION];
      }
      os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
         << STAT_label[STATL_FUNCTION] << endl;

      // ascii_print(os , nb_component , ptComponent, scale , file_flag , true , mixt_data);
    }
      } // end for (var)

    if (mixt_data != NULL) {
      double likelihood , information;
      const register unsigned int dif_var = mixt_data->nb_variable - nb_var;
      MultivariateMixtureData * mixt_vec_data = NULL;
      
      bnb_parameter = nb_parameter_computation(MIN_PROBABILITY);
      if ((mixt_data->type[0] == STATE) 
          && (dif_var > 0)) {
	// 1ere variable = etat restaure
	// eventuellement, derniere variable = entropie
	var_array = new int[dif_var];
	var_array[0] = 1;
	if (dif_var > 1)
	  var_array[1] = mixt_data->nb_variable;

	vect_data = mixt_data->select_variable(error, dif_var, var_array, false);
	if (vect_data == NULL) {
	  cerr << error;
	  likelihood = D_INF;
	}
	else {
	  likelihood = likelihood_computation(*vect_data, true);
	  mixt_vec_data = new MultivariateMixtureData(*vect_data, this->nb_component);
	  mixt_vec_data->mixture = new MultivariateMixture(*this, false);
	  information = mixt_vec_data->information_computation(); 
	  delete mixt_vec_data->mixture;
	  mixt_vec_data->mixture = NULL;	  
	  delete mixt_vec_data;
	  mixt_vec_data = NULL;
	}
      }
      else {
	likelihood = likelihood_computation(*mixt_data, true);
	information = mixt_data->information_computation();
      }
      os << "\n";
      if (file_flag)
	os << "# ";

      os << STAT_label[STATL_INFORMATION] << ": " << information << " ("
         << information / mixt_data->nb_vector << ")" << endl;

      // print the likelihood

      if (likelihood != D_INF) {
	os << "\n";
	if (file_flag)
	  os << "# ";
	os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
	   << STAT_label[STATL_NORMALIZED] << ": " << likelihood / mixt_data->nb_vector 
	   << ")" << endl;
      }

      // print AIC, AICc and BIC
      if (likelihood != D_INF) {
	if (file_flag)
	  os << "# ";

	os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;
	
	bnb_parameter = nb_parameter_computation(MIN_PROBABILITY);
	
	os << "\n";
	if (file_flag)
	  os << "# ";

	os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
	   << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AIC] << "): "
	   << 2 * (likelihood - bnb_parameter) << endl;
	
	if (0 < bnb_parameter < mixt_data->nb_vector - 1) {
	  if (file_flag) 
	    os << "# ";

	  os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
	     << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AICc] << "): "
	     << 2 * (likelihood - (double)(bnb_parameter * mixt_data->nb_vector) /
		     (double)(mixt_data->nb_vector - bnb_parameter - 1)) << endl;
	}

	if (file_flag)
	  os << "# ";

	os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
	   << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BIC] << "): "
	   << 2 * likelihood - bnb_parameter * log((double)mixt_data->nb_vector) << endl;

      if (file_flag)
        os << "# ";

      os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
         << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BICc] << "): "
         << 2 * likelihood - penalty_computation() << endl;

      // calculer la vraisemblance classifiante ?
      // likelihood = likelihood_computation(*mixt_data);
      // information = mixt_data->information_computation();
      /*
      os << "\n";
      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << ": " << likelihood << "   ("
         << STAT_label[STATL_NORMALIZED] << ": " << likelihood / mixt_data->nb_vector << ")" << endl;

      if (file_flag) {
        os << "# ";
      }
      os << STAT_label[STATL_MAX_CLASSIFICATION_LIKELIHOOD] << ": " << information << "   ("
         << STAT_label[STATL_INFORMATION] << ": " << information / mixt_data->nb_vector << ")" << endl;
      */

      } // end if (likelihood > 0)
    } // end if (mixt_data != NULL)

    if (mixt_data != NULL) {
      if (exhaustive) {
	os << "\n";
	if (file_flag) 	  
	  os << "# ";

	os << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
	mixt_data->weight->ascii_characteristic_print(os , false , file_flag);
	
	os << "\n";
	if (file_flag)
	  os << "# ";
    
	os << "   | " << STAT_label[STATL_WEIGHT] 
	   << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
	   << " | " << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

	weight->Distribution::ascii_print(os , file_flag , false , false , mixt_data->weight);
      }
    }
  }

  if (var_array != NULL) {
    delete [] var_array;
    if (vect_data != NULL)
      delete vect_data;
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MultivariateMixture.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& MultivariateMixture::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , mixture_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MultivariateMixture dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool MultivariateMixture::ascii_write(StatError &error , std::string path ,
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
    ascii_write(out_file , mixture_data , exhaustive , true);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un melange et de la structure de donnees associee
 *  dans un fichier au format tableur. Les probabilites sont
 *  remises a l'echelle des effectifs dans les lois parametriques,
 *  et aussi pour les lois categorielles mais au niveau
 *  des histogrammes
 *
 *  arguments : stream, pointeur sur un objet MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

ostream& MultivariateMixture::spreadsheet_write(ostream &os ,
                                                const MultivariateMixtureData *mixt_data) const
{
  register int i, var,
    data_var; // index that corresponds to var in mixt_data
  int bnb_parameter;
  int *var_array = NULL;
  FrequencyDistribution **observation= NULL;
  StatError error;
  Vectors *vect_data = NULL;

  os << STAT_word[STATW_MIXTURE] << endl << endl;
  os << nb_component << "\t" << STAT_word[STATW_DISTRIBUTIONS] << endl << endl;

  os << STAT_word[STATW_WEIGHTS] << endl;
  for (i = 0;i < nb_component;i++)
    os << weight->mass[i] << "\t";
  os << endl;


  if (nb_var > 0) {
    os << "\n" << nb_var << "\t"
       << STAT_word[nb_var == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

    for(var= 1; var <= nb_var; var++)
      {
    os << "\n" << STAT_word[STATW_VARIABLE];
    os << "\t" << var;

    if (npcomponent[var-1] != NULL)
      os << "\t" << STAT_word[STATW_CATEGORICAL];
    else
      os << "\t" << STAT_word[STATW_DISCRETE_PARAMETRIC];

    os << endl;
    if (mixt_data != NULL) {
      if (mixt_data->type[0] == STATE)
        data_var = var;
      else
        data_var = var-1;
      os << "\n";

      if (mixt_data->component != NULL)
        observation= mixt_data->component[data_var];
      else
        observation= NULL;
    }
    if (npcomponent[var-1] == NULL)
      pcomponent[var-1]->spreadsheet_print(os, observation);
    else
      npcomponent[var-1]->spreadsheet_print(os, observation);

      } // end for (var)

    if (mixt_data != NULL) {
      double likelihood , information;

      bnb_parameter = nb_parameter_computation(MIN_PROBABILITY);
      if (mixt_data->type[0] == STATE) {
    var_array = new int[1];
    var_array[0] = 1;
    vect_data = mixt_data->select_variable(error, 1, var_array, false);
    if (vect_data == NULL) {
      cerr << error;
      likelihood = D_INF;
    }
    else
      likelihood = likelihood_computation(*vect_data, true);
      }
      else
    likelihood = likelihood_computation(*mixt_data, true);

      information = mixt_data->information_computation();
      os << "\n";

      os << STAT_label[STATL_INFORMATION] << "\t" << information << " ("
         << information / mixt_data->nb_vector << ")" << endl;

      // print the likelihood

      if (likelihood != D_INF)
    {
      os << "\n";
      os << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
         << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / mixt_data->nb_vector << "\t" << endl;
    }

      // print AIC, AICc and BIC
      if (likelihood != D_INF)
    {
      os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

      bnb_parameter = nb_parameter_computation(MIN_PROBABILITY);

      os << "\n";
      }
      os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t 2 * "
         << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AIC] << "): "
         << 2 * (likelihood - bnb_parameter) << endl;

      if (0 < bnb_parameter < mixt_data->nb_vector - 1) {
        os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t 2 * "
           << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AICc] << "): "
           << 2 * (likelihood - (double)(bnb_parameter * mixt_data->nb_vector) /
               (double)(mixt_data->nb_vector - bnb_parameter - 1)) << endl;
      }

      os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t 2 * "
         << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BIC] << "): "
         << 2 * likelihood - bnb_parameter * log((double)mixt_data->nb_vector) << endl;

      os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t 2 * "
         << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BICc] << "): "
         << 2 * likelihood - penalty_computation() << endl;


    } // end if (likelihood > 0)
  } // end if (mixt_data != NULL)

  if (mixt_data != NULL) {
    os << "\n";
    os << STAT_label[STATL_WEIGHT] << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
    mixt_data->weight->spreadsheet_print(os);

    os << "\n";
    os << "   | " << STAT_label[STATL_WEIGHT] << "\t" << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t"
       << " | " << STAT_label[STATL_WEIGHT] << "\t" << STAT_label[STATL_DISTRIBUTION] << endl;

    weight->Distribution::spreadsheet_print(os, false , false , false , mixt_data->weight);
  }

  if (var_array != NULL) {
    delete [] var_array;
    if (vect_data != NULL)
      delete vect_data;
  }

  return os;

}



/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MultivariateMixture dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool MultivariateMixture::spreadsheet_write(StatError &error , std::string path) const

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
    spreadsheet_write(out_file , mixture_data);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Gnuplot output of MultivariateMixture and associated data, if any.
 *
 *  Parameters: prefix for files, figure titles,
 *              pointer on MultivariateMixtureData object
 *
 *--------------------------------------------------------------*/

bool MultivariateMixture::plot_write(const char *prefix , const char *title ,
                                     const MultivariateMixtureData *mixt_data) const

{
  bool status = true;
  int var = 0, cvariable = 0; // variable index that corresponds to var in mixture_data
  DiscreteParametricModel *pparam = NULL;
  FrequencyDistribution **observation = NULL;
  StatError error;

  // affiche la loi des poids
  if (mixt_data != NULL) {
    pparam = new DiscreteParametricModel(*weight, mixt_data->weight);
    status= pparam->plot_write(error, prefix, title);
    delete pparam;
    pparam = NULL;
  }

  if (status) {
    for (var = 1; var <= nb_var; var++) {
      if (mixt_data != NULL) {
    switch (mixt_data->type[0])
      {
      case INT_VALUE :
        cvariable = var - 1;
        break;
      case STATE :
        cvariable = var;
        break;
      }

    if (mixt_data->component != NULL) {
      if (mixt_data->component[cvariable] != NULL) {
        observation = mixt_data->component[cvariable];
      }
      else
        observation = NULL;
    }
      }
      if (npcomponent[var-1] != NULL)
    npcomponent[var-1]->plot_print(prefix, title, var, observation);
      else
    pcomponent[var-1]->plot_print(prefix, title, var, observation);
    }
  }
  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet MultivariateMixture.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool MultivariateMixture::plot_write(StatError &error , const char *prefix ,
                                     const char *title) const

{
  bool status = plot_write(prefix , title , mixture_data);

  error.init();

  if (!status) {
    error.update(STAT_error[STATR_FILE_PREFIX]);
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Matplotlib output of MultivariateMixture and associated data, if any.
 *  Plot marginal mixture distributions.
 *
 *  Parameter: pointer on MultivariateMixtureData object
 *
 *--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/**
 *  \brief Plot of a Mixture object and the associated data structure.
 *
 *  \param[in] vec pointer on a MixtureData object.
 *
 *  \return        MultiPlotSet object.
 */
/*--------------------------------------------------------------*/
//TODO: split MultiPlotSet into an array of of MultiPlotSets?
MultiPlotSet* MultivariateMixture::get_plotable(const MultivariateMixtureData *vec) const

{
  int i , j;
  int nb_plot_set , index , variable;
  FrequencyDistribution *marginal_dist = NULL , **observation_dist = NULL;
  Histogram *marginal_histo = NULL , **observation_histo = NULL;
  FrequencyDistribution **component = NULL;
  MultiPlotSet *plot_set, *marginal_plot_set;
  DiscreteMixture *marginal_mixture = NULL;
  ostringstream title , legend;
  StatError error;


  // computation of the number of plots

  nb_plot_set = 1;

  for (i = 0;i < nb_var;i++) {
    if (vec) {
      switch (vec->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }
    } // end if (vec)

    if ((vec) && (vec->component)) {
      nb_plot_set += nb_component;
    }
    else {
      nb_plot_set++;
    }

    if ((npcomponent[i]) && (vec->marginal_distribution[variable])) {
      if ((npcomponent[i]->weight) &&
          (npcomponent[i]->mixture)) {
        nb_plot_set++;
      }
    }

    if ((pcomponent[i]) && (vec) && (vec->marginal_distribution[variable])) {
      if ((pcomponent[i]->weight) &&
          (pcomponent[i]->mixture)) {
        nb_plot_set += 2;
      }
    }
  }
  nb_plot_set += nb_component+2;

  plot_set = new MultiPlotSet(nb_plot_set , nb_var + 1);
  MultiPlotSet &plot = *plot_set;

  plot_set->border = "15 lw 0";

  index = 0;
  plot_set->variable_nb_viewpoint[0] = 0;

  plot[index].xrange = Range(0 , weight->nb_value - 1);

  title.str("");
  title << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION];
  if (vec) {

	 title << " and " << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    // weights

    plot[index].yrange = Range(0 , ceil(MAX(vec->marginal_distribution[0]->max ,
                                    weight->max * vec->marginal_distribution[0]->nb_element)
                                    * YSCALE));

    if (weight->nb_value - 1 < TIC_THRESHOLD) {
      plot[index].xtics = 1;
    }

    plot[index].resize(2);

    legend.str("");
    legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "impulses";

    vec->marginal_distribution[0]->plotable_frequency_write(plot[0][0]);

    legend.str("");
    legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION];
    plot[index][1].legend = legend.str();

    plot[index][1].style = "linespoints";

    weight->plotable_mass_write(plot[index][1] , vec->marginal_distribution[0]->nb_element);
  }

  else {
    plot[index].yrange = Range(0 , weight->max * YSCALE);
    legend.str("");
    legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION];
    plot[index][0].legend = legend.str();

    plot[index][0].style = "linespoints";

    weight->plotable_mass_write(plot[index][0] , 1);
    // weight->plotable_mass_write(plot[index][0] , vec->marginal_distribution[0]->nb_element);

  } // end if (vec)
  plot[index].title = title.str();
  index++;

  for (i = 0;i < nb_var;i++) {
    if (vec) {
      switch (vec->type[0]) {
      case STATE :
        variable = i + 1;
        break;
      default :
        variable = i;
        break;
      }


      marginal_dist = vec->marginal_distribution[variable];

      if (vec->component) {
    	  component = vec->component[variable];
      }
      marginal_histo = vec->marginal_histogram[variable];
      if ((vec != NULL) && (vec->marginal_distribution != NULL))
    	  observation_dist = vec->marginal_distribution;
    }

    if (npcomponent[i]) {
      plot_set->variable_nb_viewpoint[i] = 0;
      npcomponent[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                             marginal_dist , MIXTURE);
    }
    else
      pcomponent[i]->plotable_write(*plot_set , index , i + 1 , observation_dist ,
                                                     marginal_dist , MIXTURE);
    marginal_mixture = extract_distribution(error, i+1);
    marginal_plot_set = marginal_mixture->get_plotable();
    // TODO: investigate seg fault, may be create copy operator or plotable_write method.
	plot[index].resize(marginal_plot_set[0].size());
    for (j=0; j < marginal_plot_set[0].size(); j++)   {
    	plot[index][j] = (*marginal_plot_set)[0][j];
    }
    delete marginal_mixture;
    delete marginal_plot_set;
    index++;
  }


  return plot_set;
}

//
//MultiPlotSet* MultivariateMixture::get_plotable(const MultivariateMixtureData *mixt_data) const
//
//{
//  bool status = true;
//  int var = 0, cvariable = 0; // variable index that corresponds to var in mixture_data
//  DiscreteParametricModel *pparam = NULL;
//  FrequencyDistribution **observation = NULL;
//  StatError error;
//  MultiPlotSet *plotset=NULL;
///*
//  // affiche la loi des poids
//  if (mixt_data != NULL) {
//	pparam = new DiscreteParametricModel(*weight, mixt_data->weight);
//	status= pparam->plot_write(error, prefix, title);
//	delete pparam;
//	pparam = NULL;
//  }
//
//  if (status) {
//	for (var = 1; var <= nb_var; var++) {
//	  if (mixt_data != NULL) {
//	switch (mixt_data->type[0])
//	  {
//	  case INT_VALUE :
//		cvariable = var - 1;
//		break;
//	  case STATE :
//		cvariable = var;
//		break;
//	  }
//
//	if (mixt_data->component != NULL) {
//	  if (mixt_data->component[cvariable] != NULL) {
//		observation = mixt_data->component[cvariable];
//	  }
//	  else
//		observation = NULL;
//	}
//	  }
//	  if (npcomponent[var-1] != NULL)
//	npcomponent[var-1]->plot_print(prefix, title, var, observation);
//	  else
//	pcomponent[var-1]->plot_print(prefix, title, var, observation);
//	}
//  }*/
//   return plotset;
//}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet MultivariateMixture.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* MultivariateMixture::get_plotable() const

{
  return get_plotable(mixture_data);
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants.
 *
 *--------------------------------------------------------------*/

int MultivariateMixture::nb_parameter_computation(double min_probability) const

{
  register int i, var;
  int bnb_parameter = nb_component - 1;


  for (var = 0; var < nb_var; var ++) {
    if (pcomponent[var] != NULL) {
      for (i = 0;i < nb_component;i++)
    bnb_parameter += pcomponent[var]->observation[i]->nb_parameter_computation();
    }
    else
      bnb_parameter += npcomponent[var]->nb_parameter_computation(min_probability);
    }

  return bnb_parameter;
}


/*--------------------------------------------------------------*
 *
 *  Calcul d'une penalite adaptative.
 *
 *--------------------------------------------------------------*/

double MultivariateMixture::penalty_computation() const

{
  register int i, var;
  double penalty = 0.;


  if (mixture_data != NULL) {
    penalty += (nb_component - 1) * log((double)mixture_data->nb_vector);
    for (var = 0; var < nb_var; var ++) {
      if (pcomponent[var] != NULL) {
    for (i = 0;i < nb_component;i++)
      penalty += pcomponent[var]->observation[i]->nb_parameter_computation() *
        log(weight->mass[i] * mixture_data->nb_vector);
      }
      else
    for (i = 0;i < nb_component;i++)
      penalty += npcomponent[var]->nb_parameter_computation(MIN_PROBABILITY) *
        log(weight->mass[i] * mixture_data->nb_vector);
    }
  }

  return penalty;
}

/*--------------------------------------------------------------*
 *
 *  Lois d'observations parametriques pour une variable donnee
 *
 *--------------------------------------------------------------*/

DiscreteParametricProcess* MultivariateMixture::get_parametric_process(int variable) const
{
  DiscreteParametricProcess *ppcomponent = NULL;

  assert((variable >= 0) && (variable < nb_var));
  if (pcomponent[variable] != NULL)
    ppcomponent = new DiscreteParametricProcess(*pcomponent[variable]);

  return ppcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Lois d'observations parametriques pour une variable donnee
 *
 *--------------------------------------------------------------*/

CategoricalProcess* MultivariateMixture::get_categorical_process(int variable) const
{
  CategoricalProcess *pnpcomponent = NULL;

  assert((variable >= 0) && (variable < nb_var));
  if (npcomponent[variable] != NULL)
    pnpcomponent = new CategoricalProcess(*npcomponent[variable]);

  return pnpcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Loi d'observation parametrique pour une variable et un etat donnes
 *
 *--------------------------------------------------------------*/

DiscreteParametric* MultivariateMixture::get_parametric_component(int variable, int index) const
{
  DiscreteParametric *ppcomponent = NULL;

  assert((variable >= 0) && (variable < nb_var));
  if (pcomponent[variable] != NULL)
    ppcomponent = new DiscreteParametric(*pcomponent[variable]->observation[index]);

  return ppcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Loi d'observation categorielle pour une variable et un etat donnes
 *
 *--------------------------------------------------------------*/

Distribution* MultivariateMixture::get_categorical_component(int variable, int index) const
{
  Distribution *pnpcomponent = NULL;

  assert((variable >= 0) && (variable < nb_var));
  if (npcomponent[variable] != NULL)
    pnpcomponent = new Distribution(*npcomponent[variable]->observation[index]);

  return pnpcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Retourne vrai si la ieme variable est parametrique
 *
 *--------------------------------------------------------------*/

bool MultivariateMixture::is_parametric(int ivariable) const {

  assert((ivariable >= 0) && (ivariable < nb_var));
  return (pcomponent[ivariable] != NULL);
}

// ###################### MultivariateMixtureData #########################################



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

MultivariateMixtureData::MultivariateMixtureData()
{
  mixture = NULL;
  nb_component = 0;
  weight = NULL;
  component = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe MultivariateMixtureData.
 *
 *  arguments : reference sur un objet FrequencyDistribution,
 *              nombre de composantes.
 *
 *--------------------------------------------------------------*/

MultivariateMixtureData::MultivariateMixtureData(const Vectors &vec , int inb_component)
:Vectors(vec)

{
  register int i, var;
  int nb_val, nb_int_variable = 0;

  mixture = NULL;
  nb_component = inb_component;

  for (i = 0; i < nb_variable; i++)
    if ((vec.get_type(i) == STATE) || (vec.get_type(i) == INT_VALUE))
      nb_int_variable++;

  weight = new FrequencyDistribution(nb_component);

  if (nb_variable > 0) {
    component = new FrequencyDistribution**[nb_variable];
    for (var = 0; var < nb_int_variable; var++) {
      component[var] = new FrequencyDistribution*[nb_component];
      nb_val = (int)ceil(get_max_value(var))+1;
      for (i = 0; i < nb_component; i++) {
	component[var][i] = new FrequencyDistribution(nb_val);
      }
    }
    for (var = nb_int_variable; var < nb_variable; var++)
      component[var] = NULL;
  }
  else
    component = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe MultivariateMixtureData.
 *
 *  argument : reference sur un objet MultivariateMixture.
 *
 *--------------------------------------------------------------*/

MultivariateMixtureData::MultivariateMixtureData(const MultivariateMixture &mixt)
{
  register int i, var, nb_variables = mixt.nb_var;

  mixture = NULL;
  nb_component = mixt.nb_component;

  weight = new FrequencyDistribution(*(mixt.weight));

  component = new FrequencyDistribution**[nb_variables];
  for (var = 0; var < nb_variables; var++) {
    component[var] = new FrequencyDistribution*[nb_component];
    for (i = 0;i < nb_component;i++) {
      if (mixt.pcomponent[var] != NULL)
    component[var][i] = new FrequencyDistribution(*mixt.pcomponent[var]->observation[i]);
      else
    component[var][i] = new FrequencyDistribution(*mixt.npcomponent[var]->observation[i]);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet MultivariateMixtureData.
 *
 *  arguments : reference sur un objet MultivariateMixtureData,
 *              flag copie de l'objet MultivariateMixture.
 *
 *--------------------------------------------------------------*/

void MultivariateMixtureData::copy(const MultivariateMixtureData &mixt_data , bool model_flag)

{
  register int i, var;


  if ((model_flag) && (mixt_data.mixture != NULL)) {
    mixture = new MultivariateMixture(*(mixt_data.mixture) , false);
  }
  else {
    mixture = NULL;
  }

  nb_component = mixt_data.nb_component;

  weight = new FrequencyDistribution(*(mixt_data.weight));

  component = new FrequencyDistribution**[nb_variable];
  for (var = 0; var < nb_variable; var++) {
    if (mixt_data.component[var] != NULL) {
      component[var] = new FrequencyDistribution*[nb_component];
      for (i = 0;i < nb_component;i++)
    component[var][i] = new FrequencyDistribution(*mixt_data.component[var][i]);
    }
    else
      component[var] = NULL;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

void MultivariateMixtureData::remove()

{
  delete mixture;
  mixture = NULL;

  delete weight;
  weight = NULL;

  if (component != NULL) {
    register int i, var;

    for (var = 0; var < nb_variable; var++) {
      for (i = 0;i < nb_component;i++)
	if (component[var] != NULL) {
	  delete component[var][i];
	  component[var][i] = NULL;
	}
      delete [] component[var];
      component[var] = NULL;
    }
    delete [] component;
    component = NULL;
  }
}

/*****************************************************************
*
*  Permutation des etats d'un MultivariateMixtureData
*  a partir d'une permutation donnee.
*  La validite de la permutation doit etre verifiee avant l'appel
*
**/

void MultivariateMixtureData::state_permutation(int* perm)
{
   register int i, var, n;
   if ((component != NULL) && (weight != NULL))
   {
      FrequencyDistribution*** pobservation = new FrequencyDistribution**[nb_variable];
      double* pmass = new double[nb_component];

      for(var = 0; var < nb_variable; var++)
      {
         if (component[var] != NULL)
            pobservation[var] = new FrequencyDistribution*[nb_component];
         for (i = 0; i < nb_component; i++)
         {
            if (component[var] != NULL)
               pobservation[var][perm[i]] = component[var][i];
            pmass[perm[i]] = weight->frequency[i];
         }
         for (i = 0; i < nb_component; i++)
         {
            if (component[var] != NULL)
               component[var][i] = pobservation[var][i];
            weight->frequency[i] = pmass[i];
         }
         if (component[var] != NULL)
         {
            delete [] pobservation[var];
            pobservation[var] = NULL;
         }
      }
      delete [] pmass;
      pmass = NULL;

      weight->mean_computation();
      weight->variance_computation();

      if (type[0] == STATE)
      {
         for(n = 0; n < nb_vector; n++)
            int_vector[n][0] = perm[int_vector[n][0]];

      }
   }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

MultivariateMixtureData::~MultivariateMixtureData()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe MultivariateMixtureData.
 *
 *  argument : reference sur un objet MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

MultivariateMixtureData& MultivariateMixtureData::operator=(const MultivariateMixtureData &mixt_data)

{
  if (&mixt_data != this) {
    remove();

    Vectors::copy(mixt_data);
    copy(mixt_data);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une composante empirique.
 *
 *  arguments : reference sur un objet StatError, variable,
 *  indice de la composante.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* MultivariateMixtureData::extract(StatError &error , int ivariable,
                                                           int index) const

{
  bool status = true;
  DiscreteDistributionData *ppcomponent = NULL;

  error.init();

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }
  if ((index < 1) || (index > nb_component)) {
    status = false;
    error.update(STAT_error[STATR_FREQUENCY_DISTRIBUTION_INDEX]);
  }
  if (type[ivariable] == STATE) {
    status = false;
    error.update(STAT_error[STATR_FREQUENCY_DISTRIBUTION_INDEX]);
  }
  if (status) {
    index--;
    if (component[ivariable][index]->nb_element == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_SAMPLE]);
    }
  }
  if (status) {
    if (mixture != NULL) {
      if (mixture->pcomponent[ivariable] != NULL)
    ppcomponent = new DiscreteDistributionData(*component[ivariable][index] ,
                        mixture->pcomponent[ivariable]->observation[index]);
      else
    ppcomponent = new DiscreteDistributionData(*component[ivariable][index] ,
                        mixture->npcomponent[ivariable]->observation[index]);
    }
    else
      ppcomponent = new DiscreteDistributionData(*component[ivariable][index] , (Distribution*)NULL);
  }

  return ppcomponent;
}



/*--------------------------------------------------------------*
 *
 *  Extraction d'une composante empirique.
 *
 *  arguments : reference sur un objet StatError, variable,
 *  indice de la composante.
 *
 *--------------------------------------------------------------*/

DiscreteDistributionData* MultivariateMixtureData::extract_marginal(StatError &error ,
                                                                    int ivariable) const

{
  bool status = true;
  register int var;
  Distribution *pDistribution = NULL;
  FrequencyDistribution *pFrequencyDistribution = NULL;
  DiscreteDistributionData *marginal = NULL;

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }
  if (type[ivariable] == STATE) {
    status = false;
    error.update(STAT_error[STATR_FREQUENCY_DISTRIBUTION_INDEX]);
  }
  if (status) {

    var = ivariable-1;
    pFrequencyDistribution = get_marginal_distribution(var);
    if (mixture != NULL)
      pDistribution = mixture->extract_distribution(error, ivariable);

    marginal = new DiscreteDistributionData(*pFrequencyDistribution, pDistribution);
  }

  return marginal;
}

/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet MultivariateMixtureData.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& MultivariateMixtureData::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_vector << "   "
     << STAT_label[STATL_VARIABLE] << ": " << nb_variable << "   "
     << STAT_label[STATL_STATES] << ": " << nb_component;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MultivariateMixtureData.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& MultivariateMixtureData::ascii_write(ostream &os , bool exhaustive) const

{
  if (mixture != NULL) {
    mixture->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MultivariateMixtureData dans un fichier.
 *
 *  arguments : reference sur un objet StatError, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool MultivariateMixtureData::ascii_write(StatError &error , std::string path ,
                                          bool exhaustive) const

{
  bool status = false;


  if (mixture != NULL) {
    ofstream out_file(path.c_str());

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      mixture->ascii_write(out_file , this , exhaustive , true);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet MultivariateMixtureData dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet StatError, path.
 *
 *--------------------------------------------------------------*/

bool MultivariateMixtureData::spreadsheet_write(StatError &error , const char *path) const

{
  bool status = false;


  if (mixture) {
    ofstream out_file(path);

    error.init();

    if (!out_file) {
      status = false;
      error.update(STAT_error[STATR_FILE_NAME]);
    }

    else {
      status = true;
      mixture->spreadsheet_write(out_file , this);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie Gnuplot d'un objet MultivariateMixtureData.
 *
 *  arguments : reference sur un objet StatError, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool MultivariateMixtureData::plot_write(StatError &error , const char *prefix ,
                                         const char *title) const

{
  bool status = false;


  if (mixture) {
    status = mixture->plot_write(prefix , title , this);

    error.init();

    if (!status) {
      error.update(STAT_error[STATR_FILE_PREFIX]);
    }
  }

  return status;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet MultivariateMixtureData.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* MultivariateMixtureData::get_plotable() const

{
  MultiPlotSet *set;


  if (mixture) {
    set = mixture->get_plotable(this);
  }
  else {
    set = NULL;
  }

  return set;
}


};  // namespace stat_tool
