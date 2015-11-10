/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *
 *       $Source$
 *       $Id: choicetable.h 3258 2007-06-06 13:18:26Z dufourko $
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


/**
 *\class ChoiceTable
 *\brief Choice table for a matching
 *\author Pascal ferraro
 *\date 2009
 */

#ifndef SB_CHOICE_TABLE_HEADER
#define SB_CHOICE_TABLE_HEADER

#include "treematch_config.h"

// #include "sequence.h"
#include "treegraph.h"

#include <iostream>
#include <list>
#include <vector>

typedef std::list<int> ChoiceList;
typedef std::vector<ChoiceList> ChoiceListVector;
typedef std::vector<ChoiceListVector> ChoiceListArray;

typedef std::pair<int,int> MatchRecord;
typedef std::list<MatchRecord> MatchRecordList;

class TREEMATCH_API  ChoiceTable
{
  public :

    /** Default ChoiceTable constructor. */
    ChoiceTable() {};

    /** Construct a ChoiceTable with {\e i_size, \e r_size}. */
    ChoiceTable(int ,int );

    /** ChoiceTable destructor. */
    ~ChoiceTable();

  void print();

    /** Change the size of a ChoiceTable with {\e i_size, \e r_size}. */
    void resize(int ,int );

    void putFirst(int ,int ,int );

    void putLast(int ,int ,int );

    int getFirst(int ,int) const ;

    void createList(int ,int );

    void destroyList(int ,int );

    ChoiceList* getList(int ,int ) ;



   void getList(int ,int , TreeGraphPtr T1, TreeGraphPtr T2, MatchRecordList*);
   void ForestList(int ,int , TreeGraphPtr T1, TreeGraphPtr T2,  MatchRecordList& );
   void TreeList(int ,int , TreeGraphPtr T1, TreeGraphPtr T2, MatchRecordList& );
   int Lat(ChoiceList* L, int vertex);

   void dump(const std::string& fname) const;
   static ChoiceTable load(const std::string& fname);

  private :
    int _i_size;
    int _r_size;
    ChoiceListArray _data;
};

#endif






