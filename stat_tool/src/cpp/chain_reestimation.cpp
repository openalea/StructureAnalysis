/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
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



#ifndef CHAIN_REESTIMATION_C
#define CHAIN_REESTIMATION_C



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut d'un objet Chain_reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
Chain_reestimation<Type>::Chain_reestimation()

{
  type = 'v';
  nb_state = 0;
  nb_row = 0;

  initial = 0;
  transition = 0;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet Chain_reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Chain_reestimation<Type>::init()

{
  register int i , j;
  Type *pinitial , *ptransition;


  if (initial) {
    pinitial = initial;
    for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
      *pinitial++ = 0;
    }
  }

  for (i = 0;i < nb_row;i++) {
    ptransition = transition[i];
    for (j = 0;j < nb_state;j++) {
      *ptransition++ = 0;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe Chain_reestimation.
 *
 *  arguments : type, nombre d'etats, nombre de lignes de la matrice
 *              des probabilites de transition, flag initialisation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
Chain_reestimation<Type>::Chain_reestimation(char itype , int inb_state ,
                                             int inb_row , bool init_flag)

{
  register int i;


  type = itype;
  nb_state = inb_state;
  nb_row = inb_row;

  if (type == 'o') {
    initial = new Type[nb_state];
  }
  else if (type == 'e') {
    initial = new Type[nb_row];
  }
  else {
    initial = 0;
  }

  transition = new Type*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new Type[nb_state];
  }

  if (init_flag) {
    init();
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Chain_reestimation.
 *
 *  argument : reference sur un objet Chain_reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Chain_reestimation<Type>::copy(const Chain_reestimation<Type> &chain_data)

{
  register int i , j;
  Type *pinitial , *cinitial , *ptransition , *ctransition;


  nb_state = chain_data.nb_state;
  nb_row = chain_data.nb_row;

  if (chain_data.initial) {
    initial = new Type[type == 'o' ? nb_state : nb_row];

    pinitial = initial;
    cinitial = chain_data.initial;
    for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
      *pinitial++ = *cinitial++;
    }
  }

  else {
    initial = 0;
  }

  transition = new Type*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new Type[nb_state];

    ptransition = transition[i];
    ctransition = chain_data.transition[i];
    for (j = 0;j < nb_state;j++) {
      *ptransition++ = *ctransition++;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet Chain_reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void Chain_reestimation<Type>::remove()

{
  delete [] initial;

  if (transition) {
    register int i;

    for (i = 0;i < nb_row;i++) {
      delete [] transition[i];
    }
    delete [] transition;
  }
}


/*--------------------------------------------------------------*
 *
 *  Destructeur de la classe Chain_reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
Chain_reestimation<Type>::~Chain_reestimation()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe Chain_reestimation.
 *
 *  argument : reference sur un objet Chain_reestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
Chain_reestimation<Type>& Chain_reestimation<Type>::operator=(const Chain_reestimation<Type> &chain_data)

{
  if (&chain_data != this) {
    remove();
    copy(chain_data);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Visualisation d'un objet Chain_reestimation.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ostream& Chain_reestimation<Type>::print(ostream &os) const

{
  register int i , j;


  os << nb_state << " " << STAT_label[STATL_STATES] << endl;

  if (initial) {
    if (type == 'o') {
      os << "\ninitial states" << endl;
    }
    else {
      os << "\nmemories" << endl;
    }

    for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
      os << initial[i] << " ";
    }
  }

  os << "\n\ntransition" << endl;
  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      os << transition[i][j] << " ";
    }
    os << endl;
  }
  os << endl;

  return os;
}



#endif
