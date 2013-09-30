/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2013 CIRAD/INRA/Inria Virtual Plants
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



#ifndef CHAIN_REESTIMATION_H
#define CHAIN_REESTIMATION_H



/****************************************************************
 *
 *  Definition de la classe :
 */


template <typename Type>
class ChainReestimation {  // structure de donnees correspondant a une chaine de Markov

    friend std::ostream& operator<<(std::ostream &os , const ChainReestimation<Type> &chain_data)
    { return chain_data.print(os); }

public :

    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int nb_state;           // nombre d'etats
    int nb_row;             // nombre de lignes de la matrice des frequences
                            // de transition
    Type *initial;          // frequences des etats initiaux
    Type **transition;      // matrice des frequences des transition

    void init();
    void copy(const ChainReestimation<Type> &chain_data);
    void remove();

    ChainReestimation();
    ChainReestimation(char itype , int inb_state , int inb_row , bool init_flag = false);
    ChainReestimation(const ChainReestimation<Type> &chain_data)
    { copy(chain_data); }
    ~ChainReestimation();
    ChainReestimation<Type>& operator=(const ChainReestimation<Type> &chain_data);

    std::ostream& print(std::ostream &os) const;
};



#include "chain_reestimation.hpp"



#endif
