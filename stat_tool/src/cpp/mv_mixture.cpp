/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): J.-B. Durand and Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: mixture.cpp 5353 2008-07-29 12:33:46Z guedon $
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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



#include <sstream>
#include <assert.h>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
// #include <rw/vstream.h>
// #include <rw/rwfile.h>
#include "stat_tools.h"
#include "distribution.h"
#include "mixture.h"
#include "mv_mixture.h"
#include "markovian.h"
#include "stat_label.h"

using namespace std;


extern char* label(const char *file_name);



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Mv_Mixture.
 *
 *--------------------------------------------------------------*/

Mv_Mixture::Mv_Mixture()

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
 *  Constructeur de la classe Mv_Mixture.
 *
 *  arguments : nombre de composantes, poids, pointeurs sur les composantes parametriques
 *  et non parametriques
 *
 *--------------------------------------------------------------*/

Mv_Mixture::Mv_Mixture(int inb_component , double *pweight , int inb_variable,
		       Parametric_process **ppcomponent, 
		       Nonparametric_process **pnpcomponent)

{
  register int var, i;


  mixture_data = NULL;
  nb_component = inb_component;
  nb_var = inb_variable;

  weight = new Parametric(nb_component);
  for (i = 0;i < nb_component;i++) {
    weight->mass[i] = *pweight++;
  }
  weight->cumul_computation();
  weight->max_computation();
  
  if (weight->ident != NONPARAMETRIC) {
    weight->computation(1 , CUMUL_THRESHOLD);
  }

  pcomponent = new Parametric_process*[nb_var];
  npcomponent = new Nonparametric_process*[nb_var];

  for (var = 0;var < nb_var;var++) {
    if ((pnpcomponent != NULL) && (pnpcomponent[var] != NULL))
      {
         npcomponent[var]= new Nonparametric_process(*(pnpcomponent[var]));
         pcomponent[var]= NULL;
      }
      else
      {
         npcomponent[var]= NULL;
         pcomponent[var]= new Parametric_process(*(ppcomponent[var]));
      }
   }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mv_Mixture.
 *
 *  arguments : reference sur un objet Mv_Mixture, flag sur les variables a copier,
 *              nombre de variables a copier.
 *
 *--------------------------------------------------------------*/

Mv_Mixture::Mv_Mixture(const Mv_Mixture &mixt , bool *variable_flag, int inb_variable)

{
  register int var;


  mixture_data = NULL;
  nb_component = mixt.nb_component;

  weight = new Parametric(nb_component);

  nb_var = inb_variable;
  pcomponent = new Parametric_process*[nb_var];
  npcomponent = new Nonparametric_process*[nb_var];

  for (var = 0;var < nb_var;var++) {
    if (variable_flag[var]) {
      if (mixt.pcomponent[var] != NULL) {
	pcomponent[var] = new Parametric_process(*mixt.pcomponent[var]);
	npcomponent[var] = NULL;
      }
      else {
	pcomponent[var] = NULL;
	npcomponent[var] = new Nonparametric_process(*mixt.npcomponent[var]);
      }
    }
  }
}

/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mv_Mixture.
 *
 *  arguments : nombre de composantes, nombre de variables,
 *  pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

Mv_Mixture::Mv_Mixture(int inb_component , int inb_variable, 
		       const Parametric_process **ppcomponent, 
		       const Nonparametric_process **pnpcomponent)

{
  register int var;


  mixture_data = NULL;
  nb_component = inb_component;
  nb_var = inb_variable;

  weight = NULL;

  pcomponent = new Parametric_process*[nb_var];
  npcomponent = new Nonparametric_process*[nb_var];

  for (var = 0;var < nb_var;var++) {
    if (ppcomponent[var] != NULL) {
	pcomponent[var] = new Parametric_process(*ppcomponent[var]);
	npcomponent[var] = NULL;
    }
    else {
       pcomponent[var] = NULL;
       npcomponent[var] = new Nonparametric_process(*pnpcomponent[var]);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Mv_Mixture.
 *
 *  arguments : reference sur un objet Mv_Mixture,
 *              flag copie de l'objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

void Mv_Mixture::copy(const Mv_Mixture &mixt , bool data_flag)

{
  register int var;


  if ((data_flag) && (mixt.mixture_data != NULL)) {
    mixture_data = new Mv_Mixture_data(*(mixt.mixture_data) , false);
  }
  else {
    mixture_data = NULL;
  }

  nb_component = mixt.nb_component;
  nb_var = mixt.nb_var;

  pcomponent = new Parametric_process*[nb_var];
  npcomponent = new Nonparametric_process*[nb_var];

  weight = new Parametric(*(mixt.weight));

  for (var = 0;var < nb_var;var++) {
    if (mixt.pcomponent[var] != NULL) {
	pcomponent[var] = new Parametric_process(*mixt.pcomponent[var]);
	npcomponent[var] = NULL;
    }
    else {
       pcomponent[var] = NULL;
       npcomponent[var] = new Nonparametric_process(*mixt.npcomponent[var]);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Mv_Mixture.
 *
 *--------------------------------------------------------------*/

void Mv_Mixture::remove()

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
 *  Destructeur de la classe Mv_Mixture.
 *
 *--------------------------------------------------------------*/

Mv_Mixture::~Mv_Mixture()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Mv_Mixture.
 *
 *  argument : reference sur un objet Mv_Mixture.
 *
 *--------------------------------------------------------------*/

Mv_Mixture& Mv_Mixture::operator=(const Mv_Mixture &mixt)

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
 *  arguments : reference sur un objet Format_error, variable, 
 *  indice de la composante.
 *
 *--------------------------------------------------------------*/

Parametric_model* Mv_Mixture::extract_parametric_model(Format_error &error , int ivariable, 
						       int index) const

{
  bool status = true;
  Parametric_model *ppcomponent = NULL;


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
      ppcomponent = new Parametric_model(*pcomponent[ivariable]->observation[index] ,
	 				(mixture_data ? mixture_data->component[ivariable][index] : NULL));
    }
  }    

  return ppcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Extraction d'une composante non parametrique.
 *
 *  arguments : reference sur un objet Format_error, variable, 
 *  indice de la composante.
 *
 *--------------------------------------------------------------*/

Distribution* Mv_Mixture::extract_nonparametric_model(Format_error &error , 
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
      pnpcomponent = new Distribution(*npcomponent[ivariable]->get_observation(index));
    }
  }

  return pnpcomponent;
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi marginale du melange
 *
 *  argument : reference sur un objet Format_error, index de la variable
 *
 *--------------------------------------------------------------*/

Distribution* Mv_Mixture::extract_distribution(Format_error &error , int ivariable) const 
{
  bool status = true;
  register int i , j, variable = ivariable - 1;  
  double *pweight , *pmass;
  Distribution *pDistribution = NULL; 

  if ((ivariable < 1) || (ivariable > nb_var)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  if (status) {

    pDistribution = new Distribution();
    pDistribution->nb_value = 0;

    if (pcomponent[variable] != NULL) {
      for (i = 0;i < nb_component;i++) {
	if (pcomponent[variable]->observation[i]->nb_value > pDistribution->nb_value) {
	  pDistribution->nb_value = pcomponent[variable]->observation[i]->nb_value;
	}
      }

      pDistribution->offset = pDistribution->nb_value; // majorant
      for (i = 0;i < nb_component;i++) {
	if (pcomponent[variable]->observation[i]->offset < pDistribution->offset) {
	  pDistribution->offset = pcomponent[variable]->observation[i]->offset;
	}
      }

      pDistribution->mass = new double[pDistribution->nb_value];
      pDistribution->cumul = new double[pDistribution->nb_value];
      pmass = pDistribution->mass - 1;
      for (i = 0;i < pDistribution->nb_value;i++) {
	pweight = weight->mass;
	*++pmass = 0.;
	for (j = 0;j < nb_component;j++) {
	  if (i < pcomponent[variable]->observation[j]->nb_value) {
	    *pmass += *pweight * pcomponent[variable]->observation[j]->mass[i];
	  }
	  pweight++;
	}
      }
    }
    else { // npcomponent[variable] != NULL)
      for (i = 0;i < nb_component;i++) {
	if (npcomponent[variable]->get_observation(i)->nb_value > pDistribution->nb_value) {
	  pDistribution->nb_value = npcomponent[variable]->get_observation(i)->nb_value;
	}
      }

      pDistribution->offset = pDistribution->nb_value; // majorant
      for (i = 0;i < nb_component;i++) {
	if (npcomponent[variable]->get_observation(i)->offset < pDistribution->offset) {
	  pDistribution->offset = npcomponent[variable]->get_observation(i)->offset;
	}
      }

      pmass = pDistribution->mass - 1;
      for (i = 0;i < pDistribution->nb_value;i++) {
	pweight = weight->mass;
	*++pmass = 0.;
	for (j = 0;j < nb_component;j++) {
	  if (i < npcomponent[variable]->get_observation(j)->nb_value) {
	    *pmass += *pweight * npcomponent[variable]->get_observation(j)->mass[i];
	  }
	  pweight++;
	}
      }
    }

    pDistribution->cumul_computation();

    pDistribution->max_computation();
    pDistribution->mean_computation();
    pDistribution->variance_computation();
  }
  return pDistribution;
}

/*--------------------------------------------------------------*
 *
 *  Extraction de la partie "donnees" d'un objet Mv_Mixture.
 *
 *  argument : reference sur un objet Format_error.
 *
 *--------------------------------------------------------------*/

Mv_Mixture_data* Mv_Mixture::extract_data(Format_error &error) const

{
  Mv_Mixture_data *mixt_data;


  error.init();

  if (!mixture_data) {
    mixt_data = NULL;
    error.update(STAT_error[STATR_NO_DATA]);
  }

  else {
    mixt_data = new Mv_Mixture_data(*mixture_data);
    mixt_data->mixture = new Mv_Mixture(*this , false);
  }

  return mixt_data;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Mv_Mixture a partir de poids et de composantes.
 *
 *  arguments : reference sur un objet Format_error, nombre de composantes,
 *              poids, pointeurs sur les composantes.
 *
 *--------------------------------------------------------------*/

Mv_Mixture* mv_mixture_building(Format_error &error , int nb_component , int nb_variable , 
				double *weight, Parametric_process **ppcomponent,
				Nonparametric_process **pnpcomponent)

{
  bool status;
  register int i;
  double cumul;
  Mv_Mixture *mixt;


  mixt = NULL;
  error.init();

  if ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT)) {
    error.update(STAT_parsing[STATP_NB_DISTRIBUTION]);
  }

  else {
    status = true;
    cumul = 0.;

    for (i = 0;i < nb_component;i++) {
      if ((weight[i] <= 0.) || (weight[i] > 1. - cumul + DOUBLE_ERROR)) {
        status = false;
        error.update(STAT_parsing[STATP_WEIGHT_VALUE]);
      }
      else {
        cumul += weight[i];
      }
    }

    if (status) {
      mixt = new Mv_Mixture(nb_component , weight , nb_variable , 
			    ppcomponent , pnpcomponent);
    }
  }

  return mixt;
}


/*--------------------------------------------------------------*
 *
 *  Construction d'un objet Mv_Mixture a partir d'un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              seuil sur la fonction de repartition.
 *
 *--------------------------------------------------------------*/

Mv_Mixture* mv_mixture_ascii_read(Format_error &error , const char *path ,
				  double cumul_threshold)

{
  RWLocaleSnapshot locale("en");
  RWCString buffer , token;
  size_t position;
  bool status , lstatus, nonparametric = false;
  register int i , j;
  int line, nb_variable;
  long index , nb_component, value;
  double cumul , *weight = NULL;   
  Nonparametric_process **np_observation;
  Parametric_process **p_observation;
  Mv_Mixture *mixt;
  ifstream in_file(path);


  mixt = 0;
  error.init();

  if (!in_file) {
    error.update(STAT_error[STATR_FILE_NAME]);
  }

  else {
    status = true;
    line = 0;
    nb_component = 0;

    // lit jusqu'au prochain saut de ligne
    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
        switch (i) {

        // test mot cle MIXTURE

        case 0 : {
          if (token != STAT_word[STATW_MIXTURE]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_MIXTURE] , line , i + 1);
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

        break;
      }

    }

    // lit jusqu'au prochain saut de ligne
    while (buffer.readLine(in_file , false)) {
      line++;

#     ifdef DEBUG
      cout << line << "  " << buffer << endl;
#     endif

      position = buffer.first('#');
      if (position != RW_NPOS) {
        buffer.remove(position);
      }
      i = 0;

      RWCTokenizer next(buffer);

      while (!((token = next()).isNull())) {
        switch (i) {

        // test nombre de composantes

        case 0 : {
          lstatus = locale.stringToNum(token , &nb_component);
          if ((lstatus) && ((nb_component < 2) || (nb_component > MIXTURE_NB_COMPONENT))) {
            lstatus = false;
          }

          if (!lstatus) {
            status = false;
            error.update(STAT_parsing[STATP_NB_DISTRIBUTION] , line , i + 1);
          }
          break;
        }

        case 1 : {
          if (token != STAT_word[STATW_DISTRIBUTIONS]) {
            status = false;
            error.correction_update(STAT_parsing[STATP_KEY_WORD] , STAT_word[STATW_DISTRIBUTIONS] , line , i + 1);
          }
          break;
        }

        }

        i++;
      }


      if (i > 0) {
        if (i != 2) {
          status = false;
          error.update(STAT_parsing[STATP_FORMAT] , line);
        }

        break;
      }
    }

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

      // lecture des poids

      index= 0;
      i = 0;

      while (buffer.readLine(in_file , false))
	{
	  line++;
	  
#         ifdef DEBUG
	  cout << line << "  " << buffer << endl;
#         endif

	  position = buffer.first('#');
	  if (position != RW_NPOS)
	    buffer.remove(position);
	  
	  // i= 0;

	  RWCTokenizer next(buffer);

	  while (!((token = next()).isNull())) {
	    if (i == 0) {
	      if (token != STAT_word[STATW_WEIGHTS]) {
		status = false;
		error.correction_update(STAT_parsing[STATP_PARAMETER_NAME] , 
					STAT_word[STATW_WEIGHTS] , line , i + 1);
	      }
	    }
	    else {
	      if (index < nb_component) {
		lstatus = locale.stringToNum(token , weight + index);
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
	  }
	  if (i >= nb_component)
	    break;	   
	} // end : while (buffer.readLine(in_file , false))
	       
      if (index != nb_component) {
	status= false;
	error.update(STAT_parsing[STATP_FORMAT] , line);
      }
	       

      if (cumul < CUMUL_THRESHOLD) {
	status = false;
	error.update(STAT_parsing[STATP_PROBABILITY_SUM]);
      }    

      // lecture des lois d'observation
      while (buffer.readLine(in_file , false))
	{
	  line++;
	  
#         ifdef DEBUG
	  cout << line << "  " << buffer << endl;
#         endif

	  position = buffer.first('#');
	  if (position != RW_NPOS)
	    buffer.remove(position);
	  
	  i= 0;
	  index = 0;

	  RWCTokenizer next(buffer);

	  while (!((token = next()).isNull()))
	    {
	      switch (i)
		{
                     // test du nombre de variables

                     case 0 :
                     {
                        lstatus= locale.stringToNum(token, &value);
                        if (lstatus)
                        {
                           if ((value < 1) || (value > NB_OUTPUT_PROCESS))
                              lstatus = false;
                           else
                              nb_variable = value;
                        }

                        if (!lstatus)
                        {
                           status= false;
                           error.update(STAT_parsing[STATP_NB_OUTPUT_PROCESS],
                                        line, i+1);
                        }
                        break;
                     }

                     // test du mot-cle VARIABLE(S)

                     case 1 :
                     {
                        if (token != STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES])
                        {
                           status = false;
                           error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                                   STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES],
                                                   line, i+1);
                        }
                        break;
                     }
                  }
                  i++;
               } // end : while (!((token = next()).isNull()))

               if (i > 0)
               {
                  if (i != 2)
                  {
                     status = false;
                     error.update(STAT_parsing[STATP_FORMAT] , line);
                  }
                  break;
               }
            } // end : while (buffer.readLine(in_file , false))

            if (nb_variable == I_DEFAULT)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT] , line);
            }
            else
            {
               np_observation= new Nonparametric_process*[nb_variable];
               p_observation = new Parametric_process*[nb_variable];
               for(i= 0; i < nb_variable; i++)
               {
                  np_observation[i]= NULL;
                  p_observation[i]= NULL;
               }
	       
               while (buffer.readLine(in_file , false))
		 {
		   line++;
		   
#                  ifdef DEBUG
		   cout << line << "  " << buffer << endl;
#                  endif
		   
		   position = buffer.first('#');
		   if (position != RW_NPOS)
                     buffer.remove(position);
		   
		   i= 0;
		   
		   RWCTokenizer next(buffer);
		   
		   while (!((token = next()).isNull()))
		     {
		       switch (i)
			 {

			   // test du mot-cle VARIABLE

			 case 0 :
			   {
			     nonparametric= true;
			     
			     if (token == STAT_word[STATW_VARIABLE])
			       index++;
			     else
			       {
				 status= false;
				 error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_word[STATW_VARIABLE],
							 line, i+1);
			       }
			     break;
			   }
			   
			   // test sur l'index de la variable
			   
			 case 1 :
			   {
			     lstatus= locale.stringToNum(token, &value);
			     if ((lstatus) && ((value != index) || (value > nb_variable)))
			       lstatus= false;
			     
			     if (!lstatus)
			       {
				 status= false;
				 error.update(STAT_parsing[STATP_VARIABLE_INDEX],
					      line, i+1);
			       }
			     break;
			   }            

			   // test on separator

			 case 2 :
			   {
			     if (token != ":")
			       {
				 status= false;
				 error.update(STAT_parsing[STATP_SEPARATOR], line, i+1);
			       }
			     break;
			   }

			   // test sur les mots-cles NONPARAMETRIC / PARAMETRIC
			   
			 case 3 :
			   {
			     if (token == STAT_word[STATW_NONPARAMETRIC])
			       nonparametric = true;
			     else
			       {
				 if (token == STAT_word[STATW_PARAMETRIC])
				   nonparametric = false;
				 else
				   {
				     status = false;
				     ostringstream correction_message;
				     correction_message << STAT_word[STATW_NONPARAMETRIC] << " or "
							<< STAT_word[STATW_PARAMETRIC] << "(instead of "
							<< token << ")";
				     error.correction_update(STAT_parsing[STATP_KEY_WORD],
							     (correction_message.str()).c_str(),
							     line, i+1);
				   }
				 break;
			       }
			   }
			 }
		       
		       i++;
		     }
		   
		   if (i > 0)
		     {
		       if (i != 4)
			 {
			   status= false;
			   error.update(STAT_parsing[STATP_FORMAT], line);
			 }
		       
		       switch (nonparametric)
			 {
			   
			 case true :
			   {
			     np_observation[index-1]= observation_parsing(error, in_file, line,
									  nb_component, true);
			     if (np_observation[index-1] == NULL)
			       status= false;
			     break;
			   }

			 case false :
			   {
			     p_observation[index-1]= observation_parsing(error, in_file, line,
									 nb_component,
									 cumul_threshold);
			     if (p_observation[index-1] == NULL)
			       status = false;
			     break;
			   }
			 }
		     }
		 } // end : while (buffer.readLine(in_file , false))
	       
               if (index != nb_variable)
		 {
		   status= false;
		   error.update(STAT_parsing[STATP_FORMAT] , line);
		 }
	       

	       if (cumul < CUMUL_THRESHOLD) {
		 status = false;
		 error.update(STAT_parsing[STATP_PROBABILITY_SUM]);
	       }

	       while (buffer.readLine(in_file , false)) {
		 line++;
		 
#              ifdef DEBUG
		 cout << line << "  " << buffer << endl;
#              endif

		 position = buffer.first('#');
		 if (position != RW_NPOS) {
		   buffer.remove(position);
		 }
		 if (!(buffer.isNull())) {
		   status = false;
		   error.update(STAT_parsing[STATP_FORMAT] , line);
		 }
	       }
	       
	       if (status) {
		 mixt = new Mv_Mixture(nb_component , weight , nb_variable , 
				       p_observation, np_observation);
	       }

               if (nb_variable > 0)
		 {
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
    }
  }
  return mixt;
}

/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Mv_Mixture.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Mv_Mixture::line_write(ostream &os) const

{
  os << nb_component << " " << STAT_word[STATW_DISTRIBUTIONS] ;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un melange et de la structure de donnees associee.
 *
 *  arguments : stream, pointeur sur un objet Mv_Mixture_data,
 *              flag niveau de detail, flag fichier.
 *
 *--------------------------------------------------------------*/

ostream& Mv_Mixture::ascii_write(ostream &os , const Mv_Mixture_data *mixt_data ,
				 bool exhaustive , bool file_flag) const

{
  register int i, var;
  int bnb_parameter;
  // double *scale;
// const Distribution *ptComponent[MIXTURE_NB_COMPONENT];
  Histogram **observation= NULL;

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
	  os << " : " << STAT_word[STATW_NONPARAMETRIC];
	else
	  os << " : " << STAT_word[STATW_PARAMETRIC];

	os << endl;
	if (mixt_data != NULL) {
	  os << "\n";
	  if (file_flag) {
	    os << "# ";
	  }

	  if (mixt_data->component != NULL)
	    observation= mixt_data->component[var-1];
	  else
	    observation= NULL;
	}
	if (npcomponent[var-1] == NULL)
	  pcomponent[var-1]->ascii_print(os, observation, exhaustive, file_flag);
	else
	  npcomponent[var-1]->ascii_print(os, observation, exhaustive, file_flag);

	if (exhaustive) {
	  /* for (i = 0;i < nb_component;i++) {
	    if (pcomponent[var-1] != NULL)
	      ptComponent[var-1] = pcomponent[var-1]->observation[i];
	    else
	      ptComponent[var-1] = npcomponent[var-1]->get_observation(i);

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
	    os << " | " << STAT_label[STATL_HISTOGRAM];
	  }
	  os << " | " << STAT_label[STATL_MIXTURE];
	  for (i = 0;i < nb_component;i++) {
	    os << " | " << STAT_label[STATL_DISTRIBUTION] << " " << var;
	  }
	  if (mixt_data != NULL) {
	    os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
	       << STAT_label[STATL_FUNCTION];
	  }
	  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
	     << STAT_label[STATL_FUNCTION] << endl;

	  // ascii_print(os , nb_component , ptComponent, scale , file_flag , true , mixt_data);
	}
      } // end for (var)
    
    if (mixt_data != NULL) {
      double likelihood , information;

      bnb_parameter = nb_parameter_computation(MIN_PROBABILITY);
      likelihood = likelihood_computation(*mixt_data, true);
      information = mixt_data->information_computation();
      os << "\n";
      if (file_flag)
	os << "# ";

      os << STAT_label[STATL_INFORMATION] << ": " << information << " ("
         << information / mixt_data->nb_vector << ")" << endl;

      // print the likelihood

      if (likelihood != D_INF)
	{
	  os << "\n";
	  if (file_flag)
	    os << "# ";
	  os << STAT_label[STATL_LIKELIHOOD] << ": " << likelihood << "   ("
	     << STAT_label[STATL_NORMALIZED] << ": " << likelihood / mixt_data->nb_vector << ")" << endl;
	}

      // print AIC, AICc and BIC
      if (likelihood != D_INF)
	{
	  if (file_flag) {
	    os << "# ";
	  }
	  os << STAT_label[STATL_DEVIANCE] << ": " << 2 * (information - likelihood) << endl;

	  bnb_parameter = nb_parameter_computation(MIN_PROBABILITY);

	  os << "\n";
	  if (file_flag) {
	    os << "# ";
	  }
	  os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
	     << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AIC] << "): "
	     << 2 * (likelihood - bnb_parameter) << endl;

	  if (0 < bnb_parameter < mixt_data->nb_vector - 1) {
	    if (file_flag) {
	      os << "# ";
	    }
	    os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
	       << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AICc] << "): "
	       << 2 * (likelihood - (double)(bnb_parameter * mixt_data->nb_vector) /
		       (double)(mixt_data->nb_vector - bnb_parameter - 1)) << endl;
	  }

	  if (file_flag) {
	    os << "# ";
	  }
	  os << bnb_parameter << " " << STAT_label[STATL_FREE_PARAMETERS] << "   2 * "
	     << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BIC] << "): "
	     << 2 * likelihood - bnb_parameter * log((double)mixt_data->nb_vector) << endl;

	  if (file_flag) {
	    os << "# ";
	  }
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
	if (file_flag) {
	  os << "# ";
	}
	os << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
	mixt_data->weight->ascii_characteristic_print(os , false , file_flag);

	os << "\n";
	if (file_flag) {
	  os << "# ";
	}
	os << "   | " << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM] << " "
	   << " | " << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

	weight->Distribution::ascii_print(os , file_flag , false , false , mixt_data->weight);
      }
    }
  }
  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mv_Mixture.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Mv_Mixture::ascii_write(ostream &os , bool exhaustive) const

{
  return ascii_write(os , mixture_data , exhaustive , false);
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mv_Mixture dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Mv_Mixture::ascii_write(Format_error &error , const char *path ,
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
 *  dans un fichier au format tableur.
 *
 *  arguments : stream, pointeur sur un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

ostream& Mv_Mixture::spreadsheet_write(ostream &os , const Mv_Mixture_data *mixt_data) const

{
  /*  register int i;
  int bnb_parameter;
  double scale[MIXTURE_NB_COMPONENT];
  const Distribution *pcomponent[MIXTURE_NB_COMPONENT];


  os << STAT_word[STATW_MIXTURE] << "\t" << nb_component << "\t" << STAT_word[STATW_DISTRIBUTIONS] << endl;
  spreadsheet_characteristic_print(os);

  if (mixt_data) {
    double likelihood , information;
    Test test(CHI2);


    os << "\n" << STAT_label[STATL_HISTOGRAM] << "\t";
    mixt_data->spreadsheet_characteristic_print(os);

    likelihood = Distribution::likelihood_computation(*mixt_data);
    information = mixt_data->Histogram::information_computation();

    os << "\n" << STAT_label[STATL_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / mixt_data->nb_vector << endl;
    os << STAT_label[STATL_MAX_LIKELIHOOD] << "\t" << information << "\t"
       << STAT_label[STATL_INFORMATION] << "\t" << information / mixt_data->nb_vector << endl;
    os << STAT_label[STATL_DEVIANCE] << "\t" << 2 * (information - likelihood) << endl;

    bnb_parameter = nb_parameter_computation();

    os << "\n" << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
       << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AIC] << ")\t"
       << 2 * (likelihood - bnb_parameter) << endl;

    if (0 < bnb_parameter < mixt_data->nb_vector - 1) {
      os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
         << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[AICc] << ")\t"
         << 2 * (likelihood - (double)(bnb_parameter * mixt_data->nb_vector) /
            (double)(mixt_data->nb_vector - bnb_parameter - 1)) << endl;
    }

    os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
       << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BIC] << ")\t"
       << 2 * likelihood - bnb_parameter * log((double)mixt_data->nb_vector) << endl;

    os << bnb_parameter << "\t" << STAT_label[STATL_FREE_PARAMETERS] << "\t"
       << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " (" << STAT_criterion_word[BICc] << ")\t"
       << 2 * likelihood - penalty_computation() << endl;

    likelihood = likelihood_computation(*mixt_data);
    information = mixt_data->information_computation();

    os << "\n" << STAT_label[STATL_CLASSIFICATION_LIKELIHOOD] << "\t" << likelihood << "\t"
       << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / mixt_data->nb_vector << endl;
    os << STAT_label[STATL_MAX_CLASSIFICATION_LIKELIHOOD] << "\t" << information << "\t"
       << STAT_label[STATL_INFORMATION] << "\t" << information / mixt_data->nb_vector << endl;

    chi2_fit(*mixt_data , test);
    os << "\n";
    test.spreadsheet_print(os);
  }

  else {
    for (i = 0;i < nb_component;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << "\t"
         << STAT_word[STATW_WEIGHT] << "\t"  << weight->mass[i] << endl;
      component[i]->spreadsheet_print(os);
      component[i]->spreadsheet_characteristic_print(os , true);
    }
  }

  for (i = 0;i < nb_component;i++) {
    pcomponent[i] = component[i];

    if (mixt_data) {
      scale[i] = mixt_data->nb_vector * weight->mass[i];
    }
    else {
      scale[i] = weight->mass[i];
    }
  }

  os << "\n";
  if (mixt_data) {
    os << "\t" << STAT_label[STATL_HISTOGRAM];
  }
  os << "\t" << STAT_label[STATL_MIXTURE];
  for (i = 0;i < nb_component;i++) {
    os << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1;
  }
  if (mixt_data) {
    os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM] << " "
       << STAT_label[STATL_FUNCTION];
  }
  os << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
     << STAT_label[STATL_FUNCTION] << endl;

  spreadsheet_print(os , nb_component , pcomponent , scale , true , mixt_data);

  if (mixt_data) {
    os << "\n" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
    mixt_data->weight->spreadsheet_characteristic_print(os);

    os << "\n\t" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM] << " "
       << "\t" << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

    weight->Distribution::spreadsheet_print(os , false , false , false , mixt_data->weight);

    for (i = 0;i < nb_component;i++) {
      os << "\n" << STAT_word[STATW_DISTRIBUTION] << "\t" << i + 1 << "\t"
         << STAT_word[STATW_WEIGHT] << "\t"  << weight->mass[i] << endl;
      component[i]->spreadsheet_print(os);
      component[i]->spreadsheet_characteristic_print(os , true);

      os << "\n" << STAT_label[STATL_HISTOGRAM] << " " << i + 1 << "\t";
      mixt_data->component[i]->spreadsheet_characteristic_print(os , true);

      if (mixt_data->component[i]->nb_element > 0) {
        os << "\n\t" << STAT_label[STATL_HISTOGRAM] << " " << i + 1
           << "\t" << STAT_label[STATL_DISTRIBUTION] << " " << i + 1
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_HISTOGRAM]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION]
           << "\t" << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION]
           << " " << i + 1 << " " << STAT_label[STATL_FUNCTION] << endl;

        component[i]->Distribution::spreadsheet_print(os , true , false , false , mixt_data->component[i]);
      }
    }
  }
  */
  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mv_Mixture dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Mv_Mixture::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un melange et de la structure de donnees associee.
 *
 *  arguments : prefixe des fichiers, titre des figures,
 *              pointeur sur un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

bool Mv_Mixture::plot_write(const char *prefix , const char *title ,
			    const Mv_Mixture_data *mixt_data) const

{
  bool status = true;
  int var = 0, cvariable = 0; 
  Parametric_model *pparam = NULL;
  Histogram **observation = NULL;
  Format_error error;

  // affiche la loi des poids
  if (mixt_data != NULL) {
    pparam = new Parametric_model(*weight, mixt_data->weight);
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
	
	if (mixt_data->component != NULL)
	  observation = mixt_data->component[cvariable];
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
 *  Sortie Gnuplot d'un objet Mv_Mixture.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Mv_Mixture::plot_write(Format_error &error , const char *prefix ,
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
 *  Sortie graphique d'un melange et de la structure de donnees associee.
 *
 *  argument : pointeur sur un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Mv_Mixture::get_plotable(const Mv_Mixture_data *mixt_data) const

{
//   register int i , j;
//   int xmax;
//   double scale;
//   std::ostringstream legend;


//   // nombre de fenetres: nb_component + 2 si ajustement

   MultiPlotSet *plotset = new MultiPlotSet(mixt_data ? nb_component + 2 : 1);
//   MultiPlotSet &set = *plotset;

//   set.title = "Mv_Mixture fit";
//   set.border = "15 lw 0";

//   // 1ere vue : melange ajuste

//   if (nb_value - 1 < TIC_THRESHOLD) {
//     set[0].xtics = 1;
//   }

//   xmax = nb_value - 1;
//   if ((cumul[xmax] > 1. - DOUBLE_ERROR) &&
//       (mass[xmax] > PLOT_MASS_THRESHOLD)) {
//     xmax++;
//   }
//   set[0].xrange = Range(0 , xmax);

//   // definition du nombre de SinglePlot 

//   i = 0;

//   if (mixt_data) {
//     set[0].yrange = Range(0 , ceil(MAX(mixt_data->max , max * mixt_data->nb_element)
//                                    * YSCALE));

//     set[0].resize(nb_component + 2);
//     set[0][i].legend = STAT_label[STATL_HISTOGRAM];
//     set[0][i].style = "impulses";

//     mixt_data->plotable_frequency_write(set[0][i++]);

//     scale = mixt_data->nb_element;
//   }

//   else {
//     set[0].yrange = Range(0. , MIN(max * YSCALE , 1.));

//     set[0].resize(nb_component + 1);

//     scale = 1;
//   }

//   set[0][i].legend = STAT_label[STATL_MIXTURE];
//   set[0][i].style = "linespoints";

//   plotable_mass_write(set[0][i++] , scale);
  
//   for (j = 0;j < nb_component;j++) {
//     legend.str("");
//     legend << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
//     set[0][i + j].legend = legend.str();

//     set[0][i + j].style = "linespoints";

//     if (mixt_data) {
//       component[j]->plotable_mass_write(set[0][i + j] , weight->mass[j] * mixt_data->nb_element);
//     }
//     else {
//       component[j]->plotable_mass_write(set[0][i + j] , weight->mass[j]);
//     }
//   }

//   if (mixt_data) {

//     // 2eme vue : poids

//     if (weight->nb_value - 1 < TIC_THRESHOLD) {
//       set[1].xtics = 1;
//     }

//     set[1].xrange = Range(0 , weight->nb_value - 1);
//     set[1].yrange = Range(0 , ceil(MAX(mixt_data->weight->max ,
//                                    weight->max * mixt_data->weight->nb_element) 
//                                    * YSCALE));

//     set[1].resize(2);

//     legend.str("");
//     legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_HISTOGRAM];
//     set[1][0].legend = legend.str();
//     set[1][0].style = "impulses";

//     legend.str("");
//     legend << STAT_label[STATL_WEIGHT] << " " << STAT_label[STATL_DISTRIBUTION];
//     set[1][1].legend = legend.str();
//     set[1][1].style = "linespoints";

//     mixt_data->weight->plotable_frequency_write(set[1][0]);
//     weight->plotable_mass_write(set[1][1] , mixt_data->weight->nb_element);

//     i = 2;
//     for (j = 0;j < nb_component;j++) {
//       if (mixt_data->component[j]->nb_element > 0) {

//         // vues suivantes : composantes ajustees

//         if (component[j]->nb_value - 1 < TIC_THRESHOLD) {
//           set[i].xtics = 1;
//         }

//         xmax = component[j]->nb_value - 1;
//         if ((component[j]->cumul[xmax] > 1. - DOUBLE_ERROR) &&
//             (component[j]->mass[xmax] > PLOT_MASS_THRESHOLD)) {
//           xmax++;
//         }
//         set[0].xrange = Range(0 , xmax);

//         set[i].xrange = Range(0 , component[j]->nb_value - 1);
//         set[i].yrange = Range(0 , ceil(MAX(mixt_data->component[j]->max ,
//                                        component[j]->max * mixt_data->component[j]->nb_element)
//                                        * YSCALE));

//         set[i].resize(2);

//         legend.str("");
//         legend << STAT_label[STATL_HISTOGRAM] << " " << j + 1;
//         set[i][0].legend = legend.str();
//         set[i][0].style = "impulses";

//         legend.str("");
//         legend << STAT_label[STATL_DISTRIBUTION] << " " << j + 1;
//         set[i][1].legend = legend.str();
//         set[i][1].style = "linespoints";

//         mixt_data->component[j]->plotable_frequency_write(set[i][0]);
//         component[j]->plotable_mass_write(set[i][1] , mixt_data->component[j]->nb_element);
  
//         i++;
//       }
//     }
//   }

   return plotset;
}


/*--------------------------------------------------------------*
 *
 *  Sortie graphique d'un objet Mv_Mixture.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Mv_Mixture::get_plotable() const

{
  return get_plotable(mixture_data);
}


/*--------------------------------------------------------------*
 *
 *  Calcul du nombre de parametres independants.
 *
 *--------------------------------------------------------------*/

int Mv_Mixture::nb_parameter_computation(double min_probability) const

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

double Mv_Mixture::penalty_computation() const

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

Parametric_process* Mv_Mixture::get_parametric_process(int variable) const 
{
  Parametric_process *ppcomponent = NULL;

  assert((variable >= 0) && (variable < nb_var));
  if (pcomponent[variable] != NULL)
    ppcomponent = new Parametric_process(*pcomponent[variable]);

  return ppcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Lois d'observations parametriques pour une variable donnee
 *
 *--------------------------------------------------------------*/

Nonparametric_process* Mv_Mixture::get_nonparametric_process(int variable) const
{
  Nonparametric_process *pnpcomponent = NULL;

  assert((variable >= 0) && (variable < nb_var));
  if (npcomponent[variable] != NULL)
    pnpcomponent = new Nonparametric_process(*npcomponent[variable]);

  return pnpcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Loi d'observation parametrique pour une variable et un etat donnes
 *
 *--------------------------------------------------------------*/

Parametric* Mv_Mixture::get_parametric_component(int variable, int index) const
{
  Parametric *ppcomponent = NULL;

  assert((variable >= 0) && (variable < nb_var));
  if (pcomponent[variable] != NULL)
    ppcomponent = new Parametric(*pcomponent[variable]->observation[index]);

  return ppcomponent;
}

/*--------------------------------------------------------------*
 *
 *  Loi d'observation non parametrique pour une variable et un etat donnes
 *
 *--------------------------------------------------------------*/

Distribution* Mv_Mixture::get_nonparametric_component(int variable, int index) const
{
  Distribution *pnpcomponent = NULL;

  assert((variable >= 0) && (variable < nb_var));
  if (npcomponent[variable] != NULL)
    pnpcomponent = new Distribution(*npcomponent[variable]->get_observation(index));

  return pnpcomponent;
}

// ###################### Mv_Mixture_data #########################################



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut de la classe Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

Mv_Mixture_data::Mv_Mixture_data()
{
  mixture = NULL;
  nb_component = 0;
  weight = NULL;
  component = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mv_Mixture_data.
 *
 *  arguments : reference sur un objet Histogram,
 *              nombre d'histogrammes.
 *
 *--------------------------------------------------------------*/

Mv_Mixture_data::Mv_Mixture_data(const Vectors &vec , int inb_component)
:Vectors(vec)

{
  register int i, var;


  mixture = NULL;
  nb_component = inb_component;

  weight = new Histogram(nb_component);
  
  if (nb_variable > 0) {
    component = new Histogram**[nb_variable];
    for (var = 0;var < nb_variable;var++) {
      component[var] = new Histogram*[nb_component];
      for (i = 0;i < nb_component;i++) {
	component[var][i] = new Histogram((int)ceil(get_max_value(var))+1);
      }
    }
  }
  else
    component = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Mv_Mixture_data.
 *
 *  argument : reference sur un objet Mv_Mixture.
 *
 *--------------------------------------------------------------*/

Mv_Mixture_data::Mv_Mixture_data(const Mv_Mixture &mixt)
{
  register int i, var, nb_variables = mixt.nb_var;

  mixture = NULL;
  nb_component = mixt.nb_component;

  weight = new Histogram(*(mixt.weight));

  component = new Histogram**[nb_variables];
  for (var = 0; var < nb_variables; var++) {
    component[var] = new Histogram*[nb_component];
    for (i = 0;i < nb_component;i++) {
      if (mixt.pcomponent[var] != NULL)
	component[var][i] = new Histogram(*mixt.pcomponent[var]->observation[i]);
      else
	component[var][i] = new Histogram(*mixt.npcomponent[var]->get_observation(i));
    }}
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Mv_Mixture_data.
 *
 *  arguments : reference sur un objet Mv_Mixture_data,
 *              flag copie de l'objet Mv_Mixture.
 *
 *--------------------------------------------------------------*/

void Mv_Mixture_data::copy(const Mv_Mixture_data &mixt_data , bool model_flag)

{
  register int i, var;


  if ((model_flag) && (mixt_data.mixture != NULL)) {
    mixture = new Mv_Mixture(*(mixt_data.mixture) , false);
  }
  else {
    mixture = NULL;
  }

  nb_component = mixt_data.nb_component;

  weight = new Histogram(*(mixt_data.weight));

  component = new Histogram**[nb_variable];
  for (var = 0; var < nb_variable; var++) {
    component[var] = new Histogram*[nb_component];
    for (i = 0;i < nb_component;i++)
      component[var][i] = new Histogram(*mixt_data.component[var][i]);
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

void Mv_Mixture_data::remove()

{
  delete mixture;
  mixture = NULL;

  delete weight;
  weight = NULL;

  if (component != NULL) {
    register int i, var;

    for (var = 0; var < nb_variable; var++) {
      for (i = 0;i < nb_component;i++) {
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


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

Mv_Mixture_data::~Mv_Mixture_data()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Mv_Mixture_data.
 *
 *  argument : reference sur un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

Mv_Mixture_data& Mv_Mixture_data::operator=(const Mv_Mixture_data &mixt_data)

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
 *  Extraction d'un histogramme correspondant a une composante
 *
 *  arguments : reference sur un objet Format_error, variable,
 *  indice de l'histogramme.
 *
 *--------------------------------------------------------------*/

Distribution_data* Mv_Mixture_data::extract(Format_error &error , int ivariable, 
					    int index) const

{
  bool status = true;
  Distribution_data *ppcomponent = NULL;

  error.init();

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }
  if ((index < 1) || (index > nb_component)) {
    status = false;
    error.update(STAT_error[STATR_HISTOGRAM_INDEX]);
  }
  else {
    index--;
    if (component[ivariable][index]->nb_element == 0) {
      status = false;
      error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
    }
  }
  if (status) {
    if (mixture != NULL) {
      if (mixture->pcomponent[ivariable] != NULL)
	ppcomponent = new Distribution_data(*component[ivariable][index] ,
					    mixture->pcomponent[ivariable]->observation[index]);
      else
	ppcomponent = new Distribution_data(*component[ivariable][index] ,
					    mixture->npcomponent[ivariable]->get_observation(index));
    }
    else
      ppcomponent = new Distribution_data(*component[ivariable][index] , (Distribution*)NULL);
  }

  return ppcomponent;
}



/*--------------------------------------------------------------*
 *
 *  Extraction d'un sous-histogramme.
 *
 *  arguments : reference sur un objet Format_error, variable,
 *  indice de l'histogramme.
 *
 *--------------------------------------------------------------*/

Distribution_data* Mv_Mixture_data::extract_marginal(Format_error &error , int ivariable) const 

{
  bool status = true;
  register int var;  
  Distribution *pDistribution = NULL; 
  Histogram *pHistogram = NULL;
  Distribution_data *marginal = NULL;

  if ((ivariable < 1) || (ivariable > nb_variable)) {
    status = false;
    error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }

  if (status) {

    var = ivariable-1;
    pHistogram = get_marginal(var);
    if (mixture != NULL)
      pDistribution = mixture->extract_distribution(error, ivariable);
      
    marginal = new Distribution_data(*pHistogram, pDistribution);
  }

  return marginal;
}

/*--------------------------------------------------------------*
 *
 *  Ecriture sur une ligne d'un objet Mv_Mixture_data.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

ostream& Mv_Mixture_data::line_write(ostream &os) const

{
  os << STAT_label[STATL_SAMPLE_SIZE] << ": " << nb_vector << "   "
     << STAT_label[STATL_VARIABLE] << ": " << nb_variable << "   "
     << STAT_label[STATL_STATES] << ": " << nb_component;

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mv_Mixture_data.
 *
 *  arguments : stream, flag niveau de detail.
 *
 *--------------------------------------------------------------*/

ostream& Mv_Mixture_data::ascii_write(ostream &os , bool exhaustive) const

{
  if (mixture != NULL) {
    mixture->ascii_write(os , this , exhaustive , false);
  }

  return os;
}


/*--------------------------------------------------------------*
 *
 *  Ecriture d'un objet Mv_Mixture_data dans un fichier.
 *
 *  arguments : reference sur un objet Format_error, path,
 *              flag niveau de detail.
 *
 *--------------------------------------------------------------*/

bool Mv_Mixture_data::ascii_write(Format_error &error , const char *path ,
				  bool exhaustive) const

{
  bool status = false;


  if (mixture != NULL) {
    ofstream out_file(path);

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
 *  Ecriture d'un objet Mv_Mixture_data dans un fichier au format tableur.
 *
 *  arguments : reference sur un objet Format_error, path.
 *
 *--------------------------------------------------------------*/

bool Mv_Mixture_data::spreadsheet_write(Format_error &error , const char *path) const

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
 *  Sortie Gnuplot d'un objet Mv_Mixture_data.
 *
 *  arguments : reference sur un objet Format_error, prefixe des fichiers,
 *              titre des figures.
 *
 *--------------------------------------------------------------*/

bool Mv_Mixture_data::plot_write(Format_error &error , const char *prefix ,
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
 *  Sortie graphique d'un objet Mv_Mixture_data.
 *
 *--------------------------------------------------------------*/

MultiPlotSet* Mv_Mixture_data::get_plotable() const

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

