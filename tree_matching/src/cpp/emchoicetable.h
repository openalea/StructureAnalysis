/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
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


#ifndef SB_EXTENDED_MATCHING_CHOICE_TABLE_HEADER
#define SB_EXTENDED_MATCHING_CHOICE_TABLE_HEADER

#include <iostream>
#include <vector>
#include <list>

typedef std::list<int> ChoiceList;
typedef std::vector<ChoiceList*> ChoiceListVector;
typedef std::vector<ChoiceListVector*> ChoiceListArray;

class ExtendedMatchingChoiceTable
{
        public :
                ExtendedMatchingChoiceTable() {};
                ExtendedMatchingChoiceTable(int ,int );
                ~ExtendedMatchingChoiceTable();
                void resize(int ,int );
                void putFirst(int ,int ,int );
                void putLast(int ,int ,int );
                void putIMFRMFLast(int ,int ,int );
                void putIPFRPFLast(int ,int ,int );
                void putIMFRFLast(int ,int ,int );
                void putIFRMFLast(int ,int ,int );
                int getFirst(int ,int) const ;
                void createList(int ,int );
                void destoyList(int ,int );
                ChoiceList* getList(int ,int ) const;
                ChoiceList* getIMFRMFList(int ,int ) const;
                ChoiceList* getIPFRPFList(int ,int ) const;
                ChoiceList* getIMFRFList(int ,int ) const;
                ChoiceList* getIFRMFList(int ,int ) const;

        private :
                int _i_size;
                int _r_size;
                ChoiceListArray _matchingChoice;
                ChoiceListArray _ipfrpfMatchingChoice;
                ChoiceListArray _imfrmfMatchingChoice;
                ChoiceListArray _imfrfMatchingChoice;
                ChoiceListArray _ifrmfMatchingChoice;
};

#endif






