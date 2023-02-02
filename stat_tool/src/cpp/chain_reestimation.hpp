/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
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



#ifndef CHAIN_REESTIMATION_C
#define CHAIN_REESTIMATION_C



namespace stat_tool {



/*--------------------------------------------------------------*/
/**
 *  \brief Default constructor of a ChainReestimation object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
ChainReestimation<Type>::ChainReestimation()

{
  type = ORDINARY;
  nb_state = 0;
  nb_row = 0;

  initial = NULL;
  transition = NULL;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Initialization of a ChainReestimation object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void ChainReestimation<Type>::init()

{
  int i , j;


  if (initial) {
    for (i = 0;i < (type == ORDINARY ? nb_state : nb_row);i++) {
      initial[i] = 0;
    }
  }

  for (i = 0;i < nb_row;i++) {
    for (j = 0;j < nb_state;j++) {
      transition[i][j] = 0;
    }
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the ChainReestimation class.
 *
 *  \param[in] itype     type,
 *  \param[in] inb_state number of states,
 *  \param[in] inb_row   number of rows of the transition probability matrix,
 *  \param[in] init_flag flag initialization.
 */
/*--------------------------------------------------------------*/

template <typename Type>
ChainReestimation<Type>::ChainReestimation(process_type itype , int inb_state ,
                                           int inb_row , bool init_flag)

{
  int i;


  type = itype;
  nb_state = inb_state;
  nb_row = inb_row;

  switch (type) {
  case ORDINARY :
    initial = new Type[nb_state];
    break;
  case EQUILIBRIUM :
    initial = new Type[nb_row];
    break;
  }

  transition = new Type*[nb_row];
  for (i = 0;i < nb_row;i++) {
    transition[i] = new Type[nb_state];
  }

  if (init_flag) {
    init();
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Copy of a ChainReestimation object.
 *
 *  \param[in] chain_data reference on a ChainReestimation object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void ChainReestimation<Type>::copy(const ChainReestimation<Type> &chain_data)

{
  int i , j;


  nb_state = chain_data.nb_state;
  nb_row = chain_data.nb_row;

  if (chain_data.initial) {
    switch (type) {

    case ORDINARY : {
      initial = new Type[nb_state];
      for (i = 0;i < nb_state;i++) {
        initial[i] = chain_data.initial[i];
      }
      break;
    }

    case EQUILIBRIUM : {
      initial = new Type[nb_row];
      for (i = 0;i < nb_row;i++) {
        initial[i] = chain_data.initial[i];
      }
      break;
    }
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


/*--------------------------------------------------------------*/
/**
 *  \brief Destruction of the data members of a ChainReestimation object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
void ChainReestimation<Type>::remove()

{
  delete [] initial;

  if (transition) {
    int i;

    for (i = 0;i < nb_row;i++) {
      delete [] transition[i];
    }
    delete [] transition;
  }
}


/*--------------------------------------------------------------*/
/**
 *  \brief Destructor of the ChainReestimation class.
 */
/*--------------------------------------------------------------*/

template <typename Type>
ChainReestimation<Type>::~ChainReestimation()

{
  remove();
}


/*--------------------------------------------------------------*/
/**
 *  \brief Assignment operator of the ChainReestimation class.
 *
 *  \param[in] chain_data reference on a ChainReestimation object.
 *
 *  \return               ChainReestimation object.
 */
/*--------------------------------------------------------------*/

template <typename Type>
ChainReestimation<Type>& ChainReestimation<Type>::operator=(const ChainReestimation<Type> &chain_data)

{
  if (&chain_data != this) {
    remove();
    copy(chain_data);
  }

  return *this;
}


/*--------------------------------------------------------------*/
/**
 *  \brief Writing of a ChainReestimation object.
 *
 *  \param[in,out] os stream.
 */
/*--------------------------------------------------------------*/

template <typename Type>
ostream& ChainReestimation<Type>::print(ostream &os) const

{
  int i , j;


  os << nb_state << " " << STAT_label[STATL_STATES] << endl;

  if (initial) {
    if (type == ORDINARY) {
      os << "\ninitial states" << endl;
    }
    else {
      os << "\nmemories" << endl;
    }

    for (i = 0;i < (type == ORDINARY ? nb_state : nb_row);i++) {
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
