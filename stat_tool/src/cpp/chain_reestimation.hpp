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
 *       $Id: chain_reestimation.hpp 17976 2015-04-23 06:34:55Z guedon $
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



#ifndef CHAIN_REESTIMATION_C
#define CHAIN_REESTIMATION_C



namespace stat_tool {



/*--------------------------------------------------------------*
 *
 *  Constructeur par defaut d'un objet ChainReestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ChainReestimation<Type>::ChainReestimation()

{
  type = 'v';
  nb_state = 0;
  nb_row = 0;

  initial = NULL;
  transition = NULL;
}


/*--------------------------------------------------------------*
 *
 *  Initialisation d'un objet ChainReestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void ChainReestimation<Type>::init()

{
  register int i , j;


  if (initial) {
    for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
      initial[i] = 0;
    }
  }

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      transition[i][j] = 0;
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Constructeur de la classe ChainReestimation.
 *
 *  arguments : type, nombre d'etats, nombre de lignes de la matrice
 *              des probabilites de transition, flag initialisation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ChainReestimation<Type>::ChainReestimation(char itype , int inb_state ,
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
    initial = NULL;
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
 *  Copie d'un objet ChainReestimation.
 *
 *  argument : reference sur un objet ChainReestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void ChainReestimation<Type>::copy(const ChainReestimation<Type> &chain_data)

{
  register int i , j;


  nb_state = chain_data.nb_state;
  nb_row = chain_data.nb_row;

  if (chain_data.initial) {
    initial = new Type[type == 'o' ? nb_state : nb_row];

    for (i = 0;i < (type == 'o' ? nb_state : nb_row);i++) {
      initial[i] = chain_data.initial[i];
    }
  }

  else {
    initial = NULL;
  }

  transition = new Type*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new Type[nb_state];

    for (j = 0;j < nb_state;j++) {
      transition[i][j] = chain_data.transition[i][j];
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Destruction des champs d'un objet ChainReestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
void ChainReestimation<Type>::remove()

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
 *  Destructeur de la classe ChainReestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ChainReestimation<Type>::~ChainReestimation()

{
  remove();
}


/*--------------------------------------------------------------*
 *
 *  Operateur d'assignement de la classe ChainReestimation.
 *
 *  argument : reference sur un objet ChainReestimation.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ChainReestimation<Type>& ChainReestimation<Type>::operator=(const ChainReestimation<Type> &chain_data)

{
  if (&chain_data != this) {
    remove();
    copy(chain_data);
  }

  return *this;
}


/*--------------------------------------------------------------*
 *
 *  Visualisation d'un objet ChainReestimation.
 *
 *  argument : stream.
 *
 *--------------------------------------------------------------*/

template <typename Type>
ostream& ChainReestimation<Type>::print(ostream &os) const

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


};  // namespace stat_tool



#endif
