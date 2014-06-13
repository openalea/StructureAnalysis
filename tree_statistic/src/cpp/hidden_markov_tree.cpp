/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *
 *       $Source$
 *       $Id: hidden_markov_tree.cpp 3193 2007-05-29 10:03:19Z dufourko $
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
 *       MERCHANTABILITY or FITNESS for A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */

// #include <deque>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"   // definition of DiscreteParametricModel class
#include "stat_tool/vectors.h"
#include "stat_tool/discrete_mixture.h"
#include "sequence_analysis/sequences.h"

#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "typed_edge_trees.h"
#include "hidden_markov_tree.h"
#include "tree_labels.h"

#include <sstream>
#include <iomanip>

extern char* label(const char*);
extern void log_computation(int nb_value, const double *pmass, double *plog);
extern int column_width(int);
extern int column_width(int nb_value , const double *value , double scale = 1.);

using namespace Stat_trees;

/*****************************************************************
 *
 *  Default constructor of CategoricalTreeProcess class
 *  using the state number, the number of values and a flag
 *  on emission distributions (indicating whether the obsevation
 *  distributions are allocated)
 *
 **/

CategoricalTreeProcess::CategoricalTreeProcess(int inb_state,
                                               int inb_value,
                                               int observation_flag)
 : CategoricalProcess(inb_state, inb_value, observation_flag)

{
   size= NULL;
   no_occurrence= NULL;
   leave= NULL;
   absorption= NULL;
   first_occurrence_root= NULL;
   first_occurrence_leaves= NULL;
   sojourn_size= NULL;
   nb_zones= NULL;
   nb_occurrences= NULL;
}

/*****************************************************************
 *
 *  Constructor of CategoricalTreeProcess class
 *  usgin the state number and the state occupancy distributions
 *
 **/

/*
CategoricalTreeProcess::CategoricalTreeProcess(int inb_state,
                                               Distribution** occupancy)
 :  nb_state(inb_state)
 ,  nb_value(inb_state)
 ,  size(NULL)
 ,  no_occurrence(NULL)
 ,  leave(NULL)
 ,  first_occurrence_root(NULL)
 ,  first_occurrence_leaves(NULL)
 ,  nb_zones(NULL)
 ,  nb_occurrences(NULL)
{
   register int i;

   observation= NULL;

   absorption= new double[nb_state];
   sojourn_size= new (Distribution*)[nb_state];
   for (i= 0; i < nb_state; i++) {
     if (occupancy[i]) {
       sojourn_size[i] = new Distribution(*occupancy[i]);
     }
     else {
       sojourn_size[i]= NULL;
     }
   }
}
*/

/*****************************************************************
 *
 *  Constructor of CategoricalTreeProcess class
 *  from a CategoricalProcess object
 *
 **/

CategoricalTreeProcess::CategoricalTreeProcess(const CategoricalProcess& process)
 : CategoricalProcess(process)
 , size(NULL)
 , no_occurrence(NULL)
 , leave(NULL)
 , absorption(NULL)
 , first_occurrence_root(NULL)
 , first_occurrence_leaves(NULL)
 , sojourn_size(NULL)
 , nb_zones(NULL)
 , nb_occurrences(NULL)
{}
// { observation= NULL; } // weird !

/*****************************************************************
 *
 *  Copy constructor CategoricalTreeProcess class
 *  using a reference on a CategoricalTreeProcess object,
 *  the kind of copy ('c' : copy, 'o' :occupancy)
 *  and a parameter (flag on the caracteristic distribution computation /
 *  number of allocated values for state occupancy distributions /
 *  reference state).
 *
 **/

CategoricalTreeProcess::CategoricalTreeProcess(const CategoricalTreeProcess& process,
                                               char manip, int param)

 : CategoricalProcess(process)
 , size(NULL)
 , no_occurrence(NULL)
 , leave(NULL)
 , absorption(NULL)
 , first_occurrence_root(NULL)
 , first_occurrence_leaves(NULL)
 , sojourn_size(NULL)
 , nb_zones(NULL)
 , nb_occurrences(NULL)
{
   switch (manip)
   {

      case 'c':
      {
         // CategoricalProcess::copy(process);
         // should already be done by constructor
         copy(process, param);
         break;
      }

      case 'o':
      {
         init_occupancy(process, param);
         break;
      }
   }
}

/*****************************************************************
 *
 *  Destructor for CategoricalTreeProcess class
 *
 **/

CategoricalTreeProcess::~CategoricalTreeProcess()
{ remove(); }

/*****************************************************************
 *
 *  Assignement operator for CategoricalTreeProcess class
 *
 **/

Stat_trees::CategoricalTreeProcess& CategoricalTreeProcess::operator=(const CategoricalTreeProcess &process)

{
  if (&process != this)
  {
     remove();
     CategoricalProcess::remove();

     CategoricalProcess::copy(process);
     copy(process);
  }

  return *this;
}

/*****************************************************************
 *
 *  Permutation of the states for CategoricalTreeProcess
 *  using a given permutation perm (validity checked before call)
 *
 **/

void CategoricalTreeProcess::state_permutation(int* perm) const
{
   if (observation != NULL)
      CategoricalProcess::state_permutation(perm);

   // if the characteristic distributions were actually computed,
   // they should be permuted as well
}


/*****************************************************************
 *
 *  Access to members of CategoricalTreeProcess class
 *
 **/

double CategoricalTreeProcess::get_double_characteristic(double* d, int value) const
{
   double res= -1.;

   if (d != NULL)
   {
      assert(value < nb_value);
      res= d[value];
   }
   return res;
}

Distribution** CategoricalTreeProcess::get_ptptDistribution_characteristic(Distribution** d) const
{
   Distribution **res= NULL;
   int v;

   if (d != NULL)
   {
      res= new Distribution*[nb_value];
      for(v= 0; v < nb_value; v++)
      {
         if (d[v] != NULL)
            res[v]= new Distribution(*(d[v]));
         else
            res[v]= NULL;
      }
   }
   return res;
}

Distribution* CategoricalTreeProcess::get_ptDistribution_characteristic(Distribution** d,
                                                                        int value) const
{
   Distribution *res= NULL;

   if (d != NULL)
   {
      assert(value < nb_value);
      if (d[value] != NULL)
         res= new Distribution(*(d[value]));
   }
   return res;
}

Distribution* CategoricalTreeProcess::get_size() const
{
   Distribution *res= NULL;

   if (size != NULL)
      res= new Distribution(*size);

   return res;
}

double CategoricalTreeProcess::get_no_occurrence(int value) const
{ return get_double_characteristic(no_occurrence, value); }

double CategoricalTreeProcess::get_leave(int value) const
{ return get_double_characteristic(leave, value); }

double CategoricalTreeProcess::get_absorption(int value) const
{ return get_double_characteristic(absorption, value); }

Distribution** CategoricalTreeProcess::get_first_occurrence_root() const
{ return get_ptptDistribution_characteristic(first_occurrence_root); }

Distribution* CategoricalTreeProcess::get_first_occurrence_root(int value) const
{ return get_ptDistribution_characteristic(first_occurrence_root, value); }

Distribution** CategoricalTreeProcess::get_first_occurrence_leaves() const
{ return get_ptptDistribution_characteristic(first_occurrence_leaves); }

Distribution* CategoricalTreeProcess::get_first_occurrence_leaves(int value) const
{ return get_ptDistribution_characteristic(first_occurrence_leaves, value); }

Distribution** CategoricalTreeProcess::get_sojourn_size() const
{ return get_ptptDistribution_characteristic(sojourn_size); }

Distribution* CategoricalTreeProcess::get_sojourn_size(int value) const
{ return get_ptDistribution_characteristic(sojourn_size, value); }

Distribution** CategoricalTreeProcess::get_nb_zones() const
{ return get_ptptDistribution_characteristic(nb_zones); }

Distribution* CategoricalTreeProcess::get_nb_zones(int value) const
{ return get_ptDistribution_characteristic(nb_zones, value); }

Distribution** CategoricalTreeProcess::get_nb_occurrences() const
{ return get_ptptDistribution_characteristic(nb_occurrences); }

Distribution* CategoricalTreeProcess::get_nb_occurrences(int value) const
{ return get_ptDistribution_characteristic(nb_occurrences, value); }


/*****************************************************************
 *
 *  Copy operator for CategoricalTreeProcess class
 *  using a flag on the copy of caracteristic distributions
 *
 **/

void CategoricalTreeProcess::copy(const CategoricalTreeProcess& process,
                                  bool characteristic_flag)

{
   // register int i;

   if (characteristic_flag)
   {

     if (process.size != NULL)
        size= new Distribution(*(process.size));
     else
       size= NULL;

      copy_double_array(no_occurrence, process.no_occurrence, nb_value);

      copy_double_array(leave, process.leave, nb_value);

      copy_double_array(absorption, process.absorption, nb_value);

      copy_Distribution_array(first_occurrence_root,
                              process.first_occurrence_root,
                              nb_value);

      copy_Distribution_array(first_occurrence_leaves,
                              process.first_occurrence_leaves,
                              nb_value);

      copy_Distribution_array(sojourn_size,
                              process.sojourn_size,
                              nb_value);

      copy_Distribution_array(nb_zones,
                              process.nb_zones,
                              nb_value);

      copy_Distribution_array(nb_occurrences,
                              process.nb_occurrences,
                              nb_value);
   }

   else
   {
      size= NULL;
      no_occurrence= NULL;
      leave= NULL;
      absorption= NULL;
      first_occurrence_root= NULL;
      first_occurrence_leaves= NULL;
      sojourn_size= NULL;
      nb_zones= NULL;
      nb_occurrences= NULL;
   }
}

/*****************************************************************
 *
 *  Deallocation of the double and Distribution* arrays
 *  for CategoricalTreeProcess' characteristic distributions
 *
 **/

void CategoricalTreeProcess::remove_double_array(double*& d)
{
   if (d != NULL)
      delete [] d;

   d= NULL;
}

/*****************************************************************
 *
 *  Deallocation of arrays of pointers on Distributions
 *
 **/

void CategoricalTreeProcess::remove_Distribution_array(Distribution**& d,
                                                       int inb_value)
{
   int i;

   if (d != NULL)
   {
      for(i= 0; i < inb_value; i++)
          if (d[i] != NULL)
          {
             delete d[i];
             d[i]= NULL;
          }
      delete [] d;
      d= NULL;
   }
}

/*****************************************************************
 *
 *  Deallocation of the pointer for CategoricalTreeProcess
 *
 **/

void CategoricalTreeProcess::remove()

{
   // register int i;


   if (size != NULL)
      delete size;

   size= NULL;

   remove_double_array(no_occurrence);

   remove_double_array(leave);

   remove_double_array(absorption);

   remove_Distribution_array(first_occurrence_root, nb_value);

   remove_Distribution_array(first_occurrence_leaves, nb_value);

   remove_Distribution_array(sojourn_size, nb_value);

   remove_Distribution_array(nb_zones, nb_value);

   remove_Distribution_array(nb_occurrences, nb_value);
}

/*****************************************************************
 *
 *  Allocation of Distribution* arrays for
 *  CategoricalTreeProcess class
 *
 **/

void CategoricalTreeProcess::init_Distribution_array(Distribution**& d,
                                                     int inb_value)
{
   int i;

   if (d != NULL)
      remove_Distribution_array(d, nb_value);

   d= new Distribution*[inb_value];

   for(i= 0; i < nb_value; i++)
       d[i] = new Distribution(NB_VALUE);
}

void CategoricalTreeProcess::init_occupancy(const CategoricalTreeProcess& process,
                                            int occupancy_nb_value)
{
   register int i;

   nb_state= process.nb_state;
   nb_value= process.nb_value;

   observation= NULL;

   size= NULL;
   no_occurrence= NULL;
   leave= NULL;
   first_occurrence_root= NULL;
   first_occurrence_leaves= NULL;
   sojourn_size= new Distribution*[nb_value];


   absorption= new double[nb_value];
   for(i= 0; i < nb_value; i++)
   {
      absorption[i]= process.absorption[i];
      if (process.sojourn_size[i] != NULL)
         sojourn_size[i]= new Distribution(*(process.sojourn_size[i]),
                                           'c', occupancy_nb_value);
      else
         sojourn_size[i]= NULL;
   }

   nb_zones= NULL;
   nb_occurrences= NULL;
}

/*****************************************************************
 *
 *  Allocation of the distributions for CategoricalTreeProcess
 *  class, except the observation distributions, using the distribution
 *  of the tree size (and maybe several other things)
 *  and flags on couting distributions and on sojourn_size
 *
 **/

void CategoricalTreeProcess::create_characteristic(const Distribution& isize,
                                                   bool sojourn_size_flag,
                                                   bool counting_flag)
{
   register int i;
   int max_size= isize.nb_value - 1;


   if (size != NULL)
      delete size;

   size= new Distribution(isize);

   if (no_occurrence == NULL)
      no_occurrence= new double[nb_value];

   for(i= 0; i < nb_value; i++)
      no_occurrence[i]= 0.;

   if (first_occurrence_root != NULL)
   {
      for(i= 0; i < nb_value; i++)
         if (first_occurrence_root[i] != NULL)
            delete first_occurrence_root[i];
   }
   else
      first_occurrence_root= new Distribution*[nb_value];

   for(i= 0; i < nb_value; i++)
      first_occurrence_root[i]=new Distribution(max_size+1); //(NB_VALUE);

   if (first_occurrence_leaves != NULL)
   {
      for(i= 0; i < nb_value; i++)
         if (first_occurrence_leaves[i] != NULL)
            delete first_occurrence_leaves[i];
   }
   else
      first_occurrence_leaves= new Distribution*[nb_value];

   for(i= 0; i < nb_value; i++)
      first_occurrence_leaves[i]=new Distribution(max_size+1); //(NB_VALUE);

   if (leave == NULL)
      leave= new double[nb_value];

   for(i= 0; i < nb_value; i++)
     leave[i]= 0.;

   if (sojourn_size_flag)
   {
      if (absorption == NULL)
         absorption= new double[nb_value];

     for(i= 0; i < nb_value; i++)
        absorption[i]= 0.;

     if (sojourn_size != NULL)
     {
        for(i= 0; i < nb_value; i++)
           if (sojourn_size[i] != NULL)
              delete sojourn_size[i];
     }
     else
        sojourn_size= new Distribution*[nb_value];

     for(i= 0; i < nb_value; i++)
       sojourn_size[i]=new Distribution(max_size+1); //(NB_VALUE);
   }

   if (counting_flag)
   {
      if (nb_zones != NULL)
      {
         for(i= 0; i < nb_value; i++)
            if (nb_zones[i] != NULL)
               delete nb_zones[i];
      }
      else
         nb_zones= new Distribution*[nb_value];

      for(i= 0; i < nb_value; i++)
      {
         nb_zones[i]=new Distribution((max_size % 2 == 0 ?
                                       max_size / 2 : max_size / 2 + 1) + 1);
      }

     if (nb_occurrences != NULL)
     {
        for(i = 0;i < nb_value;i++)
           if (nb_occurrences[i] != NULL)
              delete nb_occurrences[i];
     }
     else
        nb_occurrences= new Distribution*[nb_value];

     for(i= 0; i < nb_value; i++)
        nb_occurrences[i]=new Distribution(max_size+1);
   }
}


/*****************************************************************
 *
 *  Copy of double or Distribution* arrays for
 *  CategoricalTreeProcess class
 *
 **/

void CategoricalTreeProcess::copy_double_array(double*& dest,
                                               const double* source,
                                               int inb_value)
{
   int i;

   if (source != NULL)
   {
      dest= new double[inb_value];

      for(i= 0; i < inb_value; i++)
         dest[i] = source[i];
   }
   else
      dest= NULL;
}

void CategoricalTreeProcess::copy_Distribution_array(Distribution**& dest,
                                                     Distribution** const source,
                                                     int inb_value)
{
   int i;

   if (source != NULL)
   {
      dest= new Distribution*[inb_value];

      for(i= 0; i < inb_value; i++)
         if (source[i] != NULL)
            dest[i]= new Distribution(*(source[i]));
         else
            dest[i]= NULL;
   }
   else
      dest= NULL;
}

/*****************************************************************
 *
 *  Prints a CategoricalTreeProcess
 *  using an output stream, the process identifier,
 *  empirical observation distributions, tree characteristics
 *  and flags on the level of detail and on the file use
 *
 **/

ostream& CategoricalTreeProcess::ascii_print(ostream &os, int process,
                                             FrequencyDistribution** empirical_observation,
                                             const TreeCharacteristics* characteristics,
                                             bool exhaustive, bool file_flag) const
{
   bool no_characteristic_print = false;
   register int i , j;
   int buff , width[2];
   double *pmass;

   if (observation != NULL)
   {
      for(i= 0; i < nb_state; i++)
      {
         os << "\n" << STAT_word[STATW_STATE] << " " << i << " "
            << STAT_word[STATW_OBSERVATION_DISTRIBUTION] << endl;
         pmass= observation[i]->mass+observation[i]->offset;

         for(j= observation[i]->offset; j < observation[i]->nb_value; j++)
         {
            if (*pmass > 0.)
               os << STAT_word[STATW_OUTPUT] << " " << j << " : " << *pmass << endl;
            pmass++;
         }

         if ((empirical_observation != NULL) && (exhaustive))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << "   | " << STAT_label[STATL_STATE] << " " << i << " "
               << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << " | " << STAT_label[STATL_STATE] << " " << i << " "
               << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION] << endl;

            observation[i]->ascii_print(os, file_flag, false, false, empirical_observation[i]);
         }
      }

      // column width computation

      width[0]= column_width(nb_state-1);

      width[1]= 0;
      for(i= 0; i < nb_state; i++)
      {
         buff= column_width(observation[i]->nb_value, observation[i]->mass);
         if (buff > width[1])
            width[1] = buff;
      }
      width[1]+= ASCII_SPACE;

      // prints the output probability matrix
      os << "\n";
      if (file_flag)
         os << "# ";

      os << STAT_label[STATL_OBSERVATION_PROBABILITIY_MATRIX] << endl;

      os << "\n";
      if (file_flag)
         os << "# ";

      os << setw(width[0]+width[1]) << 0;
      for(i= 1; i < nb_value; i++)
         os << setw(width[1]) << i;

      for(i= 0; i < nb_state; i++)
      {
         os << "\n";
         if (file_flag)
            os << "# ";

         os << setw(width[0]) << i ;
         for(j= 0; j < nb_value; j++)
            os << setw(width[1]) << observation[i]->mass[j];
      }
      os << endl;
   }

   if ((characteristics != NULL) && (exhaustive))
   {
      os << "\n";
      if (file_flag)
         os << "# ";

      os << "  ";

      for(i= 0; i < nb_value; i++)
      {
         os << " | " << STAT_TREES_label[TREESTATL_OBSERVED] << " "
            << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i;
      }

      os << " | " << STAT_label[STATL_FREQUENCY];
      os << endl;

   }

   if ((first_occurrence_root != NULL) || (characteristics != NULL))
   {
      for(i= 0; i < nb_value; i++)
      {
         // if the distribution of first occurrence of a value starting
         // from root is known, print it
         if (first_occurrence_root != NULL)
         {
            // if the probability of non-occurrence for value i is positive,
            // print it
            if (no_occurrence[i] > 0.)
            {
               os << "\n";
               if (file_flag)
                  os << "# ";

               os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NO_OCCURRENCE : TREESTATL_OUTPUT_NO_OCCURRENCE]
                  << " " << i <<  ": " << no_occurrence[i] << endl;
            }

            if (first_occurrence_root[i] != NULL)
            {
               os << "\n";
               if (file_flag)
                  os << "# ";

               os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                  << " " << i << " " << STAT_label[STATL_DISTRIBUTION] << endl;
               first_occurrence_root[i]->ascii_characteristic_print(os, false, file_flag);
            }
         }

         // if the empirical distribution is available, print it
         if ((characteristics != NULL) && (characteristics->first_occurrence_root != NULL))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
               << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
            characteristics->first_occurrence_root[i]->ascii_characteristic_print(os, false, file_flag);
         }

         // print the details if required
         if ((((first_occurrence_root != NULL) && (first_occurrence_root[i] != NULL))
              || ((characteristics != NULL) && (characteristics->first_occurrence_root != NULL)
               && (characteristics->first_occurrence_root[i]->nb_element > 0))) && (exhaustive))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << "  ";
            if ((characteristics != NULL) && (characteristics->first_occurrence_root != NULL))
            {
               no_characteristic_print = false;
               if (characteristics->first_occurrence_root[i]->nb_element > 0)
                  os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                     << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
            }
            else
               no_characteristic_print = true;

            if ((first_occurrence_root != NULL) && (first_occurrence_root[i] != NULL))
            {
               os << " | " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                  << " " << i << " " << STAT_label[STATL_DISTRIBUTION];

               if ((characteristics != NULL) && (characteristics->first_occurrence_root != NULL)
                    && (characteristics->first_occurrence_root[i]->nb_element > 0))
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION];

               os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                  << STAT_label[STATL_FUNCTION] << endl;

               first_occurrence_root[i]->ascii_print(os, file_flag, true, true,
                                                     ((!no_characteristic_print) ? characteristics->first_occurrence_root[i] : NULL));
            }
            else
            // first_occurrence_root == NULL
            {
              os << endl;
              characteristics->first_occurrence_root[i]->ascii_print(os, file_flag);
            }
         }
      }
   }

   if ((first_occurrence_leaves != NULL) || (characteristics != NULL))
   {
      for(i= 0; i < nb_value; i++)
      {
         // if the distribution of first occurrence of a value starting
         // from root is known, print it
         if (first_occurrence_leaves != NULL)
         {
            // if the probability of non-occurrence for value i is positive,
            // print it
            if (no_occurrence[i] > 0.)
            {
               os << "\n";
               if (file_flag)
                  os << "# ";

               os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NO_OCCURRENCE : TREESTATL_OUTPUT_NO_OCCURRENCE]
                  << " " << i <<  ": " << no_occurrence[i] << endl;
            }

            if (first_occurrence_leaves[i] != NULL)
            {
               os << "\n";
               if (file_flag)
                  os << "# ";

               os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES]
                  << " " << i << " " << STAT_label[STATL_DISTRIBUTION] << endl;
               first_occurrence_leaves[i]->ascii_characteristic_print(os, false, file_flag);
            }
         }

         // if the empirical distribution is available, print it
         if ((characteristics != NULL) && (characteristics->first_occurrence_leaves != NULL))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES]
               << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
            characteristics->first_occurrence_leaves[i]->ascii_characteristic_print(os, false, file_flag);
         }

         // print the details if required
         if ((((first_occurrence_leaves != NULL) && (first_occurrence_leaves[i] != NULL)) ||
              ((characteristics != NULL) && ((characteristics->first_occurrence_leaves != NULL))
              && (characteristics->first_occurrence_leaves[i]->nb_element > 0))) && (exhaustive))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << "  ";
            if ((characteristics != NULL) && (characteristics->first_occurrence_leaves != NULL))
            {
               no_characteristic_print = false;
               if (characteristics->first_occurrence_leaves[i]->nb_element > 0)
                  os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT]
                     << " " << i << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
            }
            else
               no_characteristic_print = true;

            if ((first_occurrence_leaves != NULL) && (first_occurrence_leaves[i] != NULL))
            {
               os << " | " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES]
                  << " " << i << " " << STAT_label[STATL_DISTRIBUTION];

               if ((characteristics != NULL) && (characteristics->first_occurrence_leaves != NULL)
                   && (characteristics->first_occurrence_leaves[i]->nb_element > 0))
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION];

               os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                  << STAT_label[STATL_FUNCTION] << endl;

               first_occurrence_leaves[i]->ascii_print(os, file_flag, true, true,
                                                     ((!no_characteristic_print) ? characteristics->first_occurrence_leaves[i] : NULL));
            }
            else
            // first_occurrence_leaves == NULL
            {
              os << endl;
              characteristics->first_occurrence_leaves[i]->ascii_print(os, file_flag);
            }
         }
      }
   }

   // print sojourn size distribution and frequency distribution
   if ((sojourn_size != NULL) || (characteristics != NULL))
   {
      for(i= 0;i < nb_value;i++)
      {
         // if the distribution of sojourn size in a value is known, print it
         if (sojourn_size != NULL)
         {
            // if the probability of non-occurrence for value i is positive,
            // print it
            if (absorption[i] > 0.)
            {
               os << "\n";
               if (file_flag)
                  os << "# ";

               os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_ABSORPTION : TREESTATL_OUTPUT_ABSORPTION]
                  << " " << i <<  ": " << absorption[i] << endl;
            }

            if (sojourn_size[i] != NULL)
            {
               os << "\n";
               if (file_flag)
                  os << "# ";

               os << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
                  << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
               sojourn_size[i]->ascii_characteristic_print(os, false, file_flag);
            }
         }

         // if the empirical sojourn size is available, print it
         if ((characteristics != NULL) && (characteristics->sojourn_size != NULL))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
               << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
            characteristics->sojourn_size[i]->ascii_characteristic_print(os, false, file_flag);
         }

         // print the details if required
         if ((((sojourn_size != NULL) && (sojourn_size[i] != NULL)) ||
              ((characteristics != NULL) && (characteristics->sojourn_size != NULL)
              && (characteristics->sojourn_size[i]->nb_element > 0))) && (exhaustive))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << "  ";
            if ((characteristics != NULL) && (characteristics->sojourn_size != NULL)
                && (characteristics->sojourn_size[i]->nb_element > 0))

            {
               no_characteristic_print = false;
               os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
                  << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
            }
            else
               no_characteristic_print = true;

            if ((sojourn_size != NULL) && (sojourn_size[i] != NULL))
            {
               os << " | " << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << i << " "
                  << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_DISTRIBUTION];
               if (!no_characteristic_print)
               {
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION];
               }
               os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                  << STAT_label[STATL_FUNCTION] << endl;

               sojourn_size[i]->Distribution::ascii_print(os, file_flag, true, false,
                                                          ((!no_characteristic_print) ?
                                                           characteristics->sojourn_size[i] : NULL));
            }
            else
               if (!no_characteristic_print)
               {
                  os << endl;
                  characteristics->sojourn_size[i]->ascii_print(os, file_flag);
               }
         }
      } // end for
   } // end if

   if ((nb_zones != NULL) || (characteristics != NULL))
   {
      for(i= 0; i < nb_value; i++)
      {
         if (nb_zones != NULL)
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            if (size->variance == 0.)
            {
               os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                  << " " << i << " " << STAT_TREES_label[TREESTATL_PER_SIZE] << " " << size->offset << " "
                  << STAT_TREES_label[TREESTATL_TREE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
            }
            else
            {
               os << STAT_TREES_label[TREESTATL_MIXTURE_OF]
                  << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES] << " " << i << " "
                  << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
            }
            nb_zones[i]->ascii_characteristic_print(os, (((size != NULL) && (size->variance > 0.)) ? false : true), file_flag);
         }

         if ((characteristics != NULL) && (characteristics->nb_zones != NULL))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES] << " " << i << " "
               << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
            characteristics->nb_zones[i]->ascii_characteristic_print(os, (((size != NULL) && (size->variance > 0.)) ? false : true), file_flag);
         }

         if (exhaustive)
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << "  ";
            if (characteristics != NULL)
            {
               os << " | " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                  << " " << i << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
            }

            if (nb_zones != NULL)
            {
               if ((size != NULL) && (size->variance == 0.))
               {
                  os << " | " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                     << " " << i << " " << STAT_TREES_label[TREESTATL_PER_SIZE] << " " << size->offset << " "
                     << STAT_TREES_label[TREESTATL_TREE] << " " << STAT_label[STATL_DISTRIBUTION];
               }
               else
               {
                  os << " | " << STAT_TREES_label[TREESTATL_MIXTURE_OF]
                     << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES] << " " << i << " "
                     << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_DISTRIBUTIONS];
               }
               if (characteristics != NULL)
               {
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION];
                  if (characteristics->nb_zones != NULL)
                     no_characteristic_print = false;
                  else
                     no_characteristic_print = true;
               }
               else
                  no_characteristic_print = true;
               if ((size != NULL) && (size->variance == 0.))
               {
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION] << endl;
               }
               else
               {
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
                     << STAT_label[STATL_FUNCTION] << endl;
               }
               nb_zones[i]->ascii_print(os, file_flag, true, false,
                                        ((!no_characteristic_print) ? characteristics->nb_zones[i] : NULL));
            }
            else
            {
               if ((!no_characteristic_print) && (characteristics->nb_zones[i] != NULL))
               {
                  os << endl;
                  characteristics->nb_zones[i]->ascii_print(os, file_flag);
               }
            }
         }
      }
   }

   if ((nb_occurrences != NULL) || (characteristics != NULL))
   {
      for(i= 0; i < nb_value; i++)
      {
         if (nb_occurrences != NULL)
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            if ((size != NULL) && (size->variance == 0.))
            {
               os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                  << " " << i << " " << STAT_TREES_label[TREESTATL_PER_SIZE] << " " << size->offset << " "
                  << STAT_TREES_label[TREESTATL_TREE] << " " << STAT_label[STATL_DISTRIBUTION] << endl;
            }
            else
            {
               os << STAT_TREES_label[TREESTATL_MIXTURE_OF]
                  << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                  << " " << i << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_DISTRIBUTIONS] << endl;
            }
            nb_occurrences[i]->ascii_characteristic_print(os, (((size != NULL) && (size->variance > 0.)) ? false : true), file_flag);
         }

         if ((characteristics != NULL) && (characteristics->nb_occurrences != NULL))
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
               << " " << i << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
            characteristics->nb_occurrences[i]->ascii_characteristic_print(os , (((size != NULL) && (size->variance > 0.)) ? false : true) , file_flag);
         }

         if (exhaustive)
         {
            os << "\n";
            if (file_flag)
               os << "# ";

            os << "  ";
            if (characteristics != NULL)
            {
               os << " | " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                  << " " << i << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
            }

            if (nb_occurrences != NULL)
            {
               if ((size != NULL) && (size->variance == 0.))
               {
                  os << " | " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                     << " " << i << " " << STAT_TREES_label[TREESTATL_PER_SIZE] << " " << size->offset << " "
                     << STAT_TREES_label[TREESTATL_TREE] << " " << STAT_label[STATL_DISTRIBUTION];
               }
               else
               {
                  os << " | " << STAT_TREES_label[TREESTATL_MIXTURE_OF]
                     << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                     << " " << i << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_DISTRIBUTIONS];
               }
               if (characteristics)
               {
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION];
                  if (characteristics->nb_zones != NULL)
                     no_characteristic_print = false;
                  else
                     no_characteristic_print = true;
               }
               else
                  no_characteristic_print = true;
               if ((size != NULL) && (size->variance == 0.))
               {
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_DISTRIBUTION] << " "
                     << STAT_label[STATL_FUNCTION] << endl;
               }
               else
               {
                  os << " | " << STAT_label[STATL_CUMULATIVE] << " " << STAT_label[STATL_MIXTURE] << " "
                     << STAT_label[STATL_FUNCTION] << endl;
               }
               nb_occurrences[i]->ascii_print(os, file_flag, true, false ,
                                             ((!no_characteristic_print) ? characteristics->nb_occurrences[i] : NULL));
            }
            else
            {
               if ((!no_characteristic_print) && (characteristics->nb_occurrences[i] != NULL))
               {
                  os << endl;
                  characteristics->nb_occurrences[i]->ascii_print(os, file_flag);
               }
            }
         }
      }
   }

   return os;
}

ostream& CategoricalTreeProcess::spreadsheet_print(ostream &os, int process,
                                                   FrequencyDistribution** empirical_observation,
                                                   const TreeCharacteristics* characteristics) const
{
#ifdef __GNUC__
   #warning CategoricalTreeProcess::spreadsheet_print not implemented
#endif
   return os;
}

/*****************************************************************
 *
 *  Gnuplot output of a CategoricalTreeProcess, using
 *  a file prefix, the title of the output figures,
 *  the identifier of the considered process
 *  the empirical observation distributions,
 *  the empirical tree characteristics and
 *  the tree size frequency distribution
 *
 **/

bool CategoricalTreeProcess::plot_print(const char * prefix, const char * title,
                                        int process,
                                        FrequencyDistribution** empirical_observation,
                                        const TreeCharacteristics * characteristics,
                                        const FrequencyDistribution * hsize) const
{
   bool status= false, start;
   register int val, i, j= 0, k= 0;
   int nb_histo, nb_dist, histo_index,
       dist_index, plot_index, *dist_nb_value, *index_dist; //index_length ,
   double *scale;
   // Curves *smoothed_curves =NULL;
   const Distribution **pdist;
   const FrequencyDistribution **phisto;
   ostringstream data_file_name[2];

   // data file printing

   // used only for smoothed probabilities ?
   data_file_name[0] << prefix << process << 0 << ".dat";


   // name of the file containing the data for all the graphs
   data_file_name[1] << prefix << process << 1 << ".dat";

   // the 6, 7, etc. would have to be changed
   // list of the successive distributions to be printed
   pdist= new const Distribution*[6 * nb_value + nb_state];
   dist_nb_value= new int[6 * nb_value + nb_state];
   // scales for printing distributions and frequency distributions together
   scale= new double[6 * nb_value + nb_state];
   // list of the successive frequency distributions to be printed
   phisto= new const FrequencyDistribution*[1 + 7 * nb_value + nb_state];
   // correspondance between the frequency distributions and the distributions ?
   index_dist= new int[1 + 7 * nb_value + nb_state];

   nb_histo= 0;
   nb_dist= 0;

   if (hsize != NULL)
   {
      phisto[nb_histo]= hsize;
      index_dist[nb_histo++]= I_DEFAULT;
   }

   // create a list of the successive distributions to be printed
   if ((first_occurrence_root != NULL) || (characteristics != NULL))
   {
      for(val= 0; val < nb_value; val++)
      {
         if ((first_occurrence_root != NULL) && (first_occurrence_root[val] != NULL))
         {
            pdist[nb_dist]= first_occurrence_root[val];

            if ((characteristics!= NULL) && (characteristics->first_occurrence_root != NULL))
            {
               phisto[nb_histo]= characteristics->first_occurrence_root[val];
               index_dist[nb_histo]= nb_dist;
               dist_nb_value[nb_dist] = MIN(first_occurrence_root[val]->nb_value,
                                            phisto[nb_histo]->nb_value*3);
               if (first_occurrence_root[val]->cumul[first_occurrence_root[val]->nb_value-1] > 0)
                  scale[nb_dist++]= phisto[nb_histo++]->nb_element /
                                    first_occurrence_root[val]->cumul[first_occurrence_root[val]->nb_value-1];
               else
                  scale[nb_dist++]= phisto[nb_histo++]->nb_element;
            }
            else
            {
               dist_nb_value[nb_dist]= first_occurrence_root[val]->nb_value;
               scale[nb_dist++]= 1.;
            }
         }
         else
         {
            if ((characteristics != NULL) && (characteristics->first_occurrence_root != NULL))
            {
               phisto[nb_histo]= characteristics->first_occurrence_root[val];
               index_dist[nb_histo++]= I_DEFAULT;
            }
         }
      }
   }

   if ((first_occurrence_leaves != NULL) || (characteristics != NULL))
   {
      for(val= 0; val < nb_value; val++)
      {
         if ((first_occurrence_leaves != NULL) && (first_occurrence_leaves[val] != NULL))
         {
            pdist[nb_dist]= first_occurrence_leaves[val];

            if ((characteristics!= NULL) && (characteristics->first_occurrence_leaves != NULL))
            {
               phisto[nb_histo]= characteristics->first_occurrence_leaves[val];
               index_dist[nb_histo]= nb_dist;
               dist_nb_value[nb_dist] = MIN(first_occurrence_leaves[val]->nb_value,
                                            phisto[nb_histo]->nb_value*3);
               if (first_occurrence_leaves[val]->cumul[first_occurrence_leaves[val]->nb_value-1] > 0)
                  scale[nb_dist++]= phisto[nb_histo++]->nb_element /
                                    first_occurrence_leaves[val]->cumul[first_occurrence_leaves[val]->nb_value-1];
               else
                  scale[nb_dist++]= phisto[nb_histo++]->nb_element;
            }
            else
            {
               dist_nb_value[nb_dist]= first_occurrence_leaves[val]->nb_value;
               scale[nb_dist++]= 1.;
            }
         }

         else
         {
            if ((characteristics != NULL) && (characteristics->first_occurrence_leaves != NULL))
            {
               phisto[nb_histo]= characteristics->first_occurrence_leaves[val];
               index_dist[nb_histo++]= I_DEFAULT;
            }
         }
      }
   }

   if ((sojourn_size != NULL) || (characteristics != NULL))
   {
      for(val= 0; val < nb_value; val++)
      {
         if ((sojourn_size != NULL) && (sojourn_size[val] != NULL))
         {
            pdist[nb_dist]= sojourn_size[val];
            dist_nb_value[nb_dist]= sojourn_size[val]->nb_value;

            if ((characteristics != NULL) && (characteristics->sojourn_size != NULL) &&
                (characteristics->sojourn_size[val]->nb_element > 0))
            {
               phisto[nb_histo]= characteristics->sojourn_size[val];
               index_dist[nb_histo]= nb_dist;
               if ((sojourn_size[val]->cumul[sojourn_size[val]->nb_value-1] < CUMUL_THRESHOLD) &&
                   (sojourn_size[val]->cumul[sojourn_size[val]->nb_value-1] > 0))
                  scale[nb_dist++]= phisto[nb_histo++]->nb_element /
                                    sojourn_size[val]->cumul[sojourn_size[val]->nb_value-1];
               else
                  scale[nb_dist++] = phisto[nb_histo++]->nb_element;
           }
           else
             scale[nb_dist++] = 1.;
         }
         else
         {
            if ((characteristics != NULL) && (characteristics->sojourn_size != NULL) &&
                (characteristics->sojourn_size[val]->nb_element > 0))
            {
               phisto[nb_histo]= characteristics->sojourn_size[val];
               index_dist[nb_histo++]= I_DEFAULT;
            }
         }
      }
   }

   if ((nb_zones != NULL) || (nb_occurrences != NULL) ||
       (characteristics != NULL))
   {
      for(val= 0; val < nb_value; val++)
      {
         if ((nb_zones != NULL) && (nb_zones[val] != NULL))
         {
            pdist[nb_dist] = nb_zones[val];

            if ((characteristics != NULL) && (characteristics->nb_zones != NULL))
            {
               phisto[nb_histo]= characteristics->nb_zones[val];
               index_dist[nb_histo]= nb_dist;
               dist_nb_value[nb_dist]= nb_zones[val]->plot_nb_value_computation(phisto[nb_histo]);
               scale[nb_dist++]= phisto[nb_histo++]->nb_element;
            }
            else
            {
               dist_nb_value[nb_dist] = nb_zones[val]->plot_nb_value_computation();
               scale[nb_dist++] = 1.;
            }
         }
         else
            if ((characteristics != NULL) && (characteristics->nb_zones != NULL))
            {
               phisto[nb_histo]= characteristics->nb_zones[val];
               index_dist[nb_histo++]= I_DEFAULT;
            }

        if ((nb_occurrences != NULL) && (nb_occurrences[val] != NULL))
        {
           pdist[nb_dist] = nb_occurrences[val];

           if ((characteristics != NULL) && (characteristics->nb_occurrences != NULL))
           {
              phisto[nb_histo]= characteristics->nb_occurrences[val];
              index_dist[nb_histo]= nb_dist;
              dist_nb_value[nb_dist]= nb_occurrences[val]->plot_nb_value_computation(phisto[nb_histo]);
              scale[nb_dist++]= phisto[nb_histo++]->nb_element;
           }
           else
           {
              dist_nb_value[nb_dist] = nb_occurrences[val]->plot_nb_value_computation();
              scale[nb_dist++] = 1.;
           }
        }

        else
           if ((characteristics != NULL) && (characteristics->nb_occurrences != NULL))
           {
              phisto[nb_histo]= characteristics->nb_occurrences[val];
              index_dist[nb_histo++]= I_DEFAULT;
           }
      }
   }

   if (observation != NULL)
   {
      for(val= 0; val < nb_state; val++)
      {
         pdist[nb_dist]= observation[val];
         dist_nb_value[nb_dist]= observation[val]->nb_value;

         if (empirical_observation != NULL)
         {
           phisto[nb_histo]= empirical_observation[val];
           index_dist[nb_histo]= nb_dist;
           scale[nb_dist++]= phisto[nb_histo++]->nb_element;
         }
         else
            scale[nb_dist++]= 1.;
      }
   }

   // create the file, named prefix(variable+1)1.dat, containing
   // the values for all the gnuplot graphs

   status= ::plot_print((data_file_name[1].str()).c_str(), nb_dist, pdist,
                         scale, dist_nb_value, nb_histo, phisto);

   // write the command and printing files
   if (status)
   {
      // i == 0 -> output on terminal (.plot file)
      // i == 1 -> output into a postscript file  (.plot file)

      // what follows supposedly writes the smoothed probabilities
      /*
      for(i= 0; i < 2; i++)
      {
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << process << 1 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << process << 1 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << process << 1 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n";

         if (characteristics != NULL)
         {
            if (characteristics->index_value->frequency[index_length-1] < MAX_FREQUENCY)
            {
                out_file << "set title \"";
                if (title != NULL)
                   out_file << title << " - ";

                if (process > 0)
                   out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process << " - ";

                out_file << SEQ_label[SEQL_SMOOTHED_OBSERVED_PROBABILITIES] << "\"\n\n";

                if (index_length - 1 < TIC_THRESHOLD) {
                  out_file << "set xtics 0,1" << endl;
                }

                out_file << "plot [0:" << index_length - 1 << "] [0:1] ";

                j= 0;
                for (k = 0;k < nb_value;k++) {
                  j++;
                  out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                           << (index_value ? 2 * nb_value : nb_value) + k + 1
                           << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                           << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                           << k << "\" with linespoints";
                  if (index_value) {
                    j++;
                    out_file << ",\\" << endl;
                    out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                             << k + 1 << " title \"" << SEQ_label[SEQL_THEORETICAL] << " "
                             << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                             << k << "\" with linespoints";
                  }

                  if ((j == PLOT_NB_CURVE) && (k < nb_value - 1)) {
                    out_file << endl;
                    if (i == 0) {
                      out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    }
                    out_file << "\nplot [0:" << index_length - 1 << "] [0:1] ";
                  }

                  else {
                    if (k < nb_value - 1) {
                      out_file << ",\\";
                    }
                    out_file << endl;
                  }
                }

                if (index_length - 1 < TIC_THRESHOLD) {
                  out_file << "set xtics autofreq" << endl;
                }

                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << endl;
            }

            out_file << "set title \"";
            if (title != NULL) {
              out_file << title;
              if (process > 0) {
                out_file << " - ";
              }
            }
            if (process > 0) {
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            }
            out_file << "\"\n\n";

            if (index_length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << index_length - 1 << "] [0:1] ";

            j= 0;
            for (k = 0;k < nb_value;k++) {
              j++;
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << (index_value ? nb_value : 0) + k + 1 << " title \"" << SEQ_label[SEQL_OBSERVED] << " "
                       << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                       << k << "\" with linespoints";
              if (index_value) {
                j++;
                out_file << ",\\" << endl;
                out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                         << k + 1 << " title \"" << SEQ_label[SEQL_THEORETICAL] << " "
                         << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                         << k << "\" with linespoints";
              }

              if ((j == PLOT_NB_CURVE) && (k < nb_value - 1)) {
                out_file << endl;
                if (i == 0) {
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                }
                out_file << "\nplot [0:" << index_length - 1 << "] [0:1] ";
              }

              else {
                if (k < nb_value - 1) {
                  out_file << ",\\";
                }
                out_file << endl;
              }
            }

            if (index_length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }

            if (i == 0) {
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
            }
            out_file << endl;

            if (hlength->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }
            if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics 0,1" << endl;
            }

            out_file << "plot [0:" << hlength->nb_value - 1 << "] [0:"
                     << (int)(hlength->max * YSCALE) + 1 << "] \""
                     << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                     << SEQ_label[SEQL_SEQUENCE_LENGTH] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                     << "\" with impulses" << endl;

            if (hlength->nb_value - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
            if ((int)(hlength->max * YSCALE) + 1 < TIC_THRESHOLD) {
              out_file << "set ytics autofreq" << endl;
            }
          }

          else {
            out_file << "set title";
            if ((title != NULL) || (process > 0)) {
              out_file << " \"";
              if (title != NULL) {
                out_file << title;
                if (process > 0) {
                  out_file << " - ";
                }
              }
              if (process > 0) {
                out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
              }
              out_file << "\"";
            }
            out_file << "\n\n";

            if (index_value->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics 0,1" << endl;
            }

            out_file << "plot [0:" << index_value->length - 1 << "] [0:1] ";

            for (j = 0;j < nb_value;j++) {
              out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                       << j + 1 << " title \"" << SEQ_label[SEQL_THEORETICAL] << " "
                       << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                       << j << "\" with linespoints";
              if (j < nb_value - 1) {
                out_file << ",\\";
              }
              out_file << endl;
            }

            if (index_value->length - 1 < TIC_THRESHOLD) {
              out_file << "set xtics autofreq" << endl;
            }
          }

          if (i == 1) {
            out_file << "\nset terminal x11" << endl;
         }

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      } // end for i
      */

      histo_index= 1;
      dist_index= 0;

      // first occurrence_root
      for(i= 0; i < 2; i++)
      {
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << process << 1 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << process << 1 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << process << 1 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";
         if ((title != NULL) || (process > 0))
         {
            out_file << " \"";
            if (title != NULL)
            {
               out_file << title;
               if (process > 0)
                  out_file << " - ";
            }
            if (process > 0)
              out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;

            out_file << "\"";
         }
         out_file << "\n\n";

         j= histo_index; // j refers to the current frequency distribution
         k= dist_index;  // k refers to the current distribution

         start= true;
         for(val= 0; val < nb_value; val++)
         {
            if ((first_occurrence_root != NULL) && (first_occurrence_root[val] != NULL))
            {
               if (!start)
               {
                  if (i == 0)
                    out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;

                  out_file << endl;
               }
               else
                  start= false;

               if (MAX(dist_nb_value[k] , 2) - 1 < TIC_THRESHOLD)
                  out_file << "set xtics 0,1" << endl;

               if ((characteristics != NULL) && (characteristics->first_occurrence_root != NULL))
               {
                  out_file << "plot [0:" << MAX(dist_nb_value[k], 2)-1 << "] [0:"
                           << (int)(MAX(phisto[j]->max, pdist[k]->max*scale[k]) * YSCALE)+1
                           << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                           << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses,\\" << endl;
                  out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                           << " " << val << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;
                  j++;
               }

               else
               {
                  out_file << "plot [0:" << MAX(dist_nb_value[k], 2)-1 << "] [0:"
                           << MIN(pdist[k]->max*YSCALE , 1.) << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                           << " " << val << " " << STAT_label[STATL_DISTRIBUTION] << "\" with linespoints" << endl;
               }

               if (MAX(dist_nb_value[k], 2)-1 < TIC_THRESHOLD)
                 out_file << "set xtics autofreq" << endl;

               k++;
            }
            else // first_occurrence_root == NULL
               if (characteristics != NULL)
               {
                  if (!start)
                  {
                    if (i == 0)
                       out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                  }
               else
                 start = false;

               if (MAX(phisto[j]->nb_value, 2)-1 < TIC_THRESHOLD)
                  out_file << "set xtics 0,1" << endl;

               if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                  out_file << "set ytics 0,1" << endl;

               out_file << "plot [0:" << MAX(phisto[j]->nb_value, 2)-1 << "] [0:"
                        << (int)(phisto[j]->max * YSCALE)+1 << "] \""
                        << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                        << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                        << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

               if (MAX(phisto[j]->nb_value, 2)-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;

               if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                  out_file << "set ytics autofreq" << endl;

               j++;
            }
         }

         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }

      histo_index = j;
      dist_index = k;

      plot_index= 1;

      // first occurrence_leaves
      if ((first_occurrence_leaves != NULL) || (characteristics != NULL))
      {
         plot_index++;

         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
               case 0 :
                  file_name[0] << prefix << process << plot_index << ".plot";
                  break;
               case 1 :
                  file_name[0] << prefix << process << plot_index << ".print";
                  break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << process << plot_index << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";

            if ((title != NULL) || (process > 0))
            {
               out_file << " \"";
               if (title != NULL)
               {
                  out_file << title;
                  if (process > 0)
                     out_file << " - ";
               }
               if (process > 0)
                  out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
              out_file << "\"";
            }
            out_file << "\n\n";

            j= histo_index;
            k= dist_index;

            start= true;
            for(val= 0; val < nb_value; val++)
            {
               if ((first_occurrence_leaves != NULL) &&
                   (first_occurrence_leaves[val] != NULL))
               {
                 if (!start)
                 {
                    if (i == 0)
                       out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                 }
                 else
                   start= false;

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                    out_file << "set xtics 0,1" << endl;

               if ((characteristics != NULL) &&
                   (characteristics->first_occurrence_leaves != NULL) &&
                   (characteristics->first_occurrence_leaves[val]->nb_element > 0))
               {
                  out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                           << (int)(MAX(phisto[j]->max, pdist[k]->max*scale[k])*YSCALE)+1
                           << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses,\\" << endl;
                  out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES] << " " << STAT_label[STATL_DISTRIBUTION];
                  first_occurrence_leaves[val]->plot_title_print(out_file);
                  out_file << "\" with linespoints" << endl;
                  j++;
               }

               else
               {
                  out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                           << MIN(pdist[k]->max*YSCALE, 1.) << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES] << " " << STAT_label[STATL_DISTRIBUTION];
                  first_occurrence_leaves[val]->plot_title_print(out_file);
                  out_file << "\" with linespoints" << endl;
               }

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;
               k++;
            }
            else // first_occurrence_leaves != NULL
               if ((characteristics != NULL) &&
                   (characteristics->first_occurrence_leaves[val]->nb_element > 0))
               {
                  if (!start)
                  {
                    if (i == 0)
                      out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                  }
                  else
                     start = false;

                  if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics 0,1" << endl;

                  out_file << "plot [0:" << phisto[j]->nb_value-1 << "] [0:"
                           << (int)(phisto[j]->max*YSCALE)+1 << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses" << endl;

                  if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;

                  if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics autofreq" << endl;
                  j++;
               }
            } // end for val
            if (i == 1)
               out_file << "\nset terminal x11" << endl;

            out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
         } // end for i
         histo_index= j;
         dist_index= k;
      }


      // sojourn_size
      if ((sojourn_size != NULL) || (characteristics != NULL))
      {
         plot_index++;

         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
               case 0 :
                  file_name[0] << prefix << process << plot_index << ".plot";
                  break;
               case 1 :
                  file_name[0] << prefix << process << plot_index << ".print";
                  break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << process << plot_index << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";

            if ((title != NULL) || (process > 0))
            {
               out_file << " \"";
               if (title != NULL)
               {
                  out_file << title;
                  if (process > 0)
                     out_file << " - ";
               }
               if (process > 0)
                  out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
              out_file << "\"";
            }
            out_file << "\n\n";

            j= histo_index;
            k= dist_index;

            start= true;
            for(val= 0; val < nb_value; val++)
            {
               if ((sojourn_size != NULL) && (sojourn_size[val] != NULL))
               {
                 if (!start)
                 {
                    if (i == 0)
                       out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                 }
                 else
                   start= false;

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                    out_file << "set xtics 0,1" << endl;

               if ((characteristics != NULL) && (characteristics->sojourn_size != NULL) &&
                   (characteristics->sojourn_size[val]->nb_element > 0))
               {
                  out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                           << (int)(MAX(phisto[j]->max, pdist[k]->max*scale[k])*YSCALE)+1
                           << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses,\\" << endl;
                  out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_DISTRIBUTION];
                  sojourn_size[val]->plot_title_print(out_file);
                  out_file << "\" with linespoints" << endl;
                  j++;
               }

               else
               {
                  out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                           << MIN(pdist[k]->max*YSCALE, 1.) << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_DISTRIBUTION];
                  sojourn_size[val]->plot_title_print(out_file);
                  out_file << "\" with linespoints" << endl;
               }

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;
               k++;
            }
            else // sojourn_size != NULL
               if ((characteristics != NULL) && (characteristics->sojourn_size[val]->nb_element > 0))
               {
                  if (!start)
                  {
                    if (i == 0)
                      out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                  }
                  else
                     start = false;

                  if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics 0,1" << endl;

                  out_file << "plot [0:" << phisto[j]->nb_value-1 << "] [0:"
                           << (int)(phisto[j]->max*YSCALE)+1 << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses" << endl;

                  if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;

                  if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics autofreq" << endl;
                  j++;
               }
            } // end for val
            if (i == 1)
               out_file << "\nset terminal x11" << endl;

            out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
         } // end for i
         histo_index= j;
         dist_index= k;
      }

      // counting
      if ((nb_zones != NULL) || (nb_occurrences != NULL) ||
          (characteristics != NULL))
      {
         plot_index++;

         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
               case 0 :
                  file_name[0] << prefix << process << plot_index << ".plot";
                  break;
               case 1 :
                  file_name[0] << prefix << process << plot_index << ".print";
                  break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << process << plot_index << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";
            if ((title != NULL) || (process > 0))
            {
               out_file << " \"";
               if (title != NULL)
               {
                  out_file << title;
                  if (process > 0)
                     out_file << " - ";
               }
               if (process > 0)
                  out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
               out_file << "\"";
            }
            out_file << "\n\n";

            j= histo_index;
            k= dist_index;

            start= true;
            for(val= 0; val < nb_value; val++)
            {
               if (nb_zones != NULL)
               {
                  if (!start)
                  {
                     if (i == 0)
                        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                     out_file << endl;
                  }
                  else
                     start= false;

                  if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if ((characteristics != NULL) && (characteristics->nb_zones != NULL))
                  {
                     out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                              << (int)(MAX(phisto[j]->max, pdist[k]->max * scale[k])*YSCALE)+1
                              << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                              << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                              << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                              << "\" with impulses,\\" << endl;
                     out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1;
                     if ((size != NULL) && (size->variance == 0.))
                     {
                        out_file << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                                 << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << size->offset << " "
                                 << STAT_TREES_label[TREESTATL_TREE] << " " << STAT_label[STATL_DISTRIBUTION];
                     }
                     else
                     {
                        out_file << " title \"" << STAT_TREES_label[TREESTATL_MIXTURE_OF]
                                 << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                                 << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_DISTRIBUTION];
                        out_file << "\" with linespoints" << endl;
                        j++;
                     }
                  }
                  // uncomment if a distribution for the size is added
                  /*
                  else // characteristics == NULL
                  {
                     out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                              << MIN(pdist[k]->max*YSCALE, 1.) << "] \""
                              << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                              << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                              << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << length->offset << " "
                              << STAT_TREES_label[TREESTATL_TREE] << " " << STAT_label[STATL_DISTRIBUTION]
                              << "\" with linespoints" << endl;
                  }*/
                  if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;
                  k++;
               }
               else // nb_zones == NULL
                  if (characteristics != NULL)
                  {
                     if (!start)
                     {
                        if (i == 0)
                           out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                        out_file << endl;
                     }
                     else
                        start = false;

                     if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                        out_file << "set xtics 0,1" << endl;

                     if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                        out_file << "set ytics 0,1" << endl;

                     out_file << "plot [0:" << phisto[j]->nb_value-1 << "] [0:"
                              << (int)(phisto[j]->max*YSCALE)+1 << "] \""
                              << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                              << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                              << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                              << "\" with impulses" << endl;

                     if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                        out_file << "set xtics autofreq" << endl;

                     if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                        out_file << "set ytics autofreq" << endl;

                     j++;
                  }

                 if (nb_occurrences != NULL)
                 {
                    if (!start)
                    {
                       if (i == 0)
                          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                       out_file << endl;
                    }
                    else
                       start= false;

                    if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                       out_file << "set xtics 0,1" << endl;

                    if ((characteristics != NULL) && (characteristics->nb_occurrences != NULL))
                    {
                       out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                                << (int)(MAX(phisto[j]->max, pdist[k]->max * scale[k])*YSCALE)+1
                                << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                                << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                                << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                                << "\" with impulses,\\" << endl;
                       out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1;
                       if ((size != NULL) && (size->variance == 0.))
                       {
                         out_file << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                                  << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << size->offset << " "
                                  << STAT_TREES_label[TREESTATL_TREE] << " " << STAT_label[STATL_DISTRIBUTION];
                       }
                       else
                       {
                          out_file << " title \"" << STAT_TREES_label[TREESTATL_MIXTURE_OF]
                                   << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                                   << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_DISTRIBUTION];
                          out_file << "\" with linespoints" << endl;
                          j++;
                       }
                    }
                    else
                    // uncomment if a distribution of the size is added
                    /*
                    {
                       out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                                << MIN(pdist[k]->max*YSCALE, 1.) << "] \""
                                << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                                << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                                << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << length->offset << " "
                                << STAT_TREES_label[TREESTATL_TREE] << " " << STAT_label[STATL_DISTRIBUTION]
                                << "\" with linespoints" << endl;
                    }*/
                    if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                       out_file << "set xtics autofreq" << endl;
                    k++;
                 }
                 else // nb_occurrences == NULL
                    if (characteristics != NULL)
                    {
                       if (i == 0)
                          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                       out_file << endl;

                       if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                          out_file << "set xtics 0,1" << endl;

                       if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                          out_file << "set ytics 0,1" << endl;

                       out_file << "plot [0:" << phisto[j]->nb_value-1 << "] [0:"
                                << (int)(phisto[j]->max*YSCALE)+1 << "] \""
                                << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                                << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                                << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                                << "\" with impulses" << endl;

                       if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                          out_file << "set xtics autofreq" << endl;

                       if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                          out_file << "set ytics autofreq" << endl;
                       j++;
                    }
               } // end nb_zones == NULL

               if (characteristics != NULL)
               {
                  if (i == 0)
                     out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                  out_file << endl;

                  if (hsize->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if ((int)(hsize->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics 0,1" << endl;

                  out_file << "plot [0:" << hsize->nb_value - 1 << "] [0:"
                           << (int)(hsize->max*YSCALE)+1 << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                           << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses" << endl;

                  if (hsize->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;

                  if ((int)(hsize->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics autofreq" << endl;
               }

               if (i == 1)
                  out_file << "\nset terminal x11" << endl;

               out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
            } // end for val

         histo_index= j;
         dist_index= k;
      }

      // observation
      if (observation != NULL)
      {
         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
              case 0 :
                 file_name[0] << prefix << process << 0 << ".plot";
                 break;
              case 1 :
                 file_name[0] << prefix << process << 0 << ".print";
                 break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << process << 0 << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";

            if (title != NULL)
               out_file << " \"" << title << " - " << STAT_label[STATL_OUTPUT_PROCESS]
                        << " " << process << "\"";
            out_file << "\n\n";

            j= histo_index;
            k= dist_index;

            for(val= 0; val < nb_state; val++)
            {
               if (dist_nb_value[k] - 1 < TIC_THRESHOLD)
                  out_file << "set xtics 0,1" << endl;

               if (empirical_observation != NULL)
               {
                  out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                           << (int)(MAX(phisto[j]->max, pdist[k]->max * scale[k])*YSCALE)+1
                           << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[STATL_STATE] << " " << val << " "
                           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses,\\" << endl;
                  out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[STATL_STATE] << " " << val << " "
                           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION]
                           << "\" with linespoints" << endl;
                  j++;
               }
               else
               {
                  out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                           << MIN(pdist[k]->max*YSCALE, 1.) << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[STATL_STATE] << " " << val << " "
                           << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION]
                           << "\" with linespoints" << endl;
               }

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;
               k++;

               if ((i == 0) && (val < nb_state - 1))
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
               out_file << endl;
            }

            if (i == 1)
               out_file << "\nset terminal x11" << endl;

            out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
         }
      }

      delete [] pdist;
      delete [] dist_nb_value;
      delete [] scale;
      delete [] phisto;
      delete [] index_dist;
   }

   return status;
}

/*****************************************************************
 *
 *  Matplotlib output of CategoricalTreeProcess, using
 *  a MultiPlotSet instance, the current Plot index,
 *  the identifier of the considered process
 *  the empirical observation distributions,
 *  the empirical tree characteristics and
 *  the tree size frequency distribution
 *
 **/

MultiPlotSet* CategoricalTreeProcess::plotable_write(MultiPlotSet &plot, int &index,
                                                     int process, FrequencyDistribution * const * empirical_observation,
                                                     const TreeCharacteristics * characteristics,
                                                     const FrequencyDistribution * hsize) const
{
   register int val, i, j, var;
   int dist_nb_value;
   double scale, max;
   ostringstream title, legend;


   // create a list of the successive distributions to be printed
   if ((first_occurrence_root != NULL) || (characteristics != NULL))
   {
      plot.variable_nb_viewpoint[process]++; // add a view for first_occurrence_root
      for(val = 0; val < nb_value; val++)
      {
         if ((first_occurrence_root != NULL) && (first_occurrence_root[val] != NULL))
         {
            plot.variable[index] = process;
            plot.viewpoint[index] = FIRST_OCCURRENCE_ROOT;

            if (process > 0)
            {
               title.str("");
               title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
               plot[index].title = title.str();
            }

            plot[index].xrange = Range(0 , MAX(first_occurrence_root[val]->nb_value, 2) - 1);
            if (MAX(first_occurrence_root[val]->nb_value, 2) - 1 < TIC_THRESHOLD)
               plot[index].xtics = 1;

            if ((characteristics != NULL) && (val < characteristics->get_nb_values())
                && (characteristics->first_occurrence_root[val]->nb_element > 0))
            {
               scale = characteristics->first_occurrence_root[val]->nb_element /
                          first_occurrence_root[val]->cumul[first_occurrence_root[val]->nb_value-1];

               plot[index].yrange = Range(0, ceil(MAX(characteristics->first_occurrence_root[val]->max,
                                                      first_occurrence_root[val]->max * scale) * YSCALE));

               plot[index].resize(2);
               legend.str("");
               legend << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                      << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
               plot[index][0].legend = legend.str();

               plot[index][0].style = "impulses";

               characteristics->first_occurrence_root[val]->plotable_frequency_write(plot[index][0]);
               j = 1;
            }
            else
            {
               scale = 1.;
               plot[index].yrange = Range(0, ceil(MAX(first_occurrence_root[val]->max * YSCALE,
                                                      1.)));
               plot[index].resize(1);
               j = 0;
            }
            legend.str("");
            legend << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                   << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
            plot[index][j].legend = legend.str();
            plot[index][j].style = "linespoints";

            first_occurrence_root[val]->plotable_mass_write(plot[index][j], scale);
            index++;
         } // ((first_occurrence_root == NULL) || (first_occurrence_root[val] == NULL))
         else
            if ((characteristics != NULL) && (val < characteristics->get_nb_values())
                && (characteristics->first_occurrence_root[val]->nb_element > 0))
            {
               plot.variable[index] = process;
               plot.viewpoint[index] = FIRST_OCCURRENCE_ROOT;

               if (process > 0)
               {
                  title.str("");
                  title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
                  plot[index].title = title.str();
               }

               plot[index].xrange = Range(0 , MAX(characteristics->first_occurrence_root[val]->nb_value, 2) - 1);
               plot[index].yrange = Range(0, ceil(MAX(characteristics->first_occurrence_root[val]->max,
                                                      characteristics->first_occurrence_root[val]->max * scale) * YSCALE));
               if (MAX(characteristics->first_occurrence_root[val]->nb_value, 2) - 1 < TIC_THRESHOLD)
                  plot[index].xtics = 1;

               if (ceil(characteristics->first_occurrence_root[val]->max * YSCALE) < TIC_THRESHOLD)
                  plot[index].ytics = 1;

               plot[index].resize(1);
               legend.str("");
               legend << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                      << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
               plot[index][0].legend = legend.str();

               plot[index][0].style = "impulses";

               characteristics->first_occurrence_root[val]->plotable_frequency_write(plot[index][0]);
               index++;
            }

      } // end for val
   }

   /*
   if ((first_occurrence_leaves != NULL) || (characteristics != NULL))
   {
      plot.variable_nb_viewpoint[process]++; // add a view for first_occurrence_leaves
      for(val= 0; val < nb_value; val++)
      {
         if ((first_occurrence_leaves != NULL) && (first_occurrence_leaves[val] != NULL))
         {
            pdist[nb_dist]= first_occurrence_leaves[val];

            if ((characteristics!= NULL) && (characteristics->first_occurrence_leaves != NULL))
            {
               phisto[nb_histo]= characteristics->first_occurrence_leaves[val];
               index_dist[nb_histo]= nb_dist;
               dist_nb_value[nb_dist] = MIN(first_occurrence_leaves[val]->nb_value,
                                            phisto[nb_histo]->nb_value*3);
               if (first_occurrence_leaves[val]->cumul[first_occurrence_leaves[val]->nb_value-1] > 0)
                  scale[nb_dist++]= phisto[nb_histo++]->nb_element /
                                    first_occurrence_leaves[val]->cumul[first_occurrence_leaves[val]->nb_value-1];
               else
                  scale[nb_dist++]= phisto[nb_histo++]->nb_element;
            }
            else
            {
               dist_nb_value[nb_dist]= first_occurrence_leaves[val]->nb_value;
               scale[nb_dist++]= 1.;
            }
         }

         else
         {
            if ((characteristics != NULL) && (characteristics->first_occurrence_leaves != NULL))
            {
               phisto[nb_histo]= characteristics->first_occurrence_leaves[val];
               index_dist[nb_histo++]= I_DEFAULT;
            }
         }
      }
   }

   if ((sojourn_size != NULL) || (characteristics != NULL))
   {
      plot.variable_nb_viewpoint[process]++; // add a view for sojourn_size
      for(val= 0; val < nb_value; val++)
      {
         if ((sojourn_size != NULL) && (sojourn_size[val] != NULL))
         {
            pdist[nb_dist]= sojourn_size[val];
            dist_nb_value[nb_dist]= sojourn_size[val]->nb_value;

            if ((characteristics != NULL) && (characteristics->sojourn_size != NULL) &&
                (characteristics->sojourn_size[val]->nb_element > 0))
            {
               phisto[nb_histo]= characteristics->sojourn_size[val];
               index_dist[nb_histo]= nb_dist;
               if ((sojourn_size[val]->cumul[sojourn_size[val]->nb_value-1] < CUMUL_THRESHOLD) &&
                   (sojourn_size[val]->cumul[sojourn_size[val]->nb_value-1] > 0))
                  scale[nb_dist++]= phisto[nb_histo++]->nb_element /
                                    sojourn_size[val]->cumul[sojourn_size[val]->nb_value-1];
               else
                  scale[nb_dist++] = phisto[nb_histo++]->nb_element;
           }
           else
             scale[nb_dist++] = 1.;
         }
         else
         {
            if ((characteristics != NULL) && (characteristics->sojourn_size != NULL) &&
                (characteristics->sojourn_size[val]->nb_element > 0))
            {
               phisto[nb_histo]= characteristics->sojourn_size[val];
               index_dist[nb_histo++]= I_DEFAULT;
            }
         }
      }
   }

   if ((nb_zones != NULL) || (nb_occurrences != NULL) ||
       (characteristics != NULL))
   {
      plot.variable_nb_viewpoint[process]++; // add a view for nb_zones or nb_occurrences
      for(val= 0; val < nb_value; val++)
      {
         if ((nb_zones != NULL) && (nb_zones[val] != NULL))
         {
            pdist[nb_dist] = nb_zones[val];

            if ((characteristics != NULL) && (characteristics->nb_zones != NULL))
            {
               phisto[nb_histo]= characteristics->nb_zones[val];
               index_dist[nb_histo]= nb_dist;
               dist_nb_value[nb_dist]= nb_zones[val]->plot_nb_value_computation(phisto[nb_histo]);
               scale[nb_dist++]= phisto[nb_histo++]->nb_element;
            }
            else
            {
               dist_nb_value[nb_dist] = nb_zones[val]->plot_nb_value_computation();
               scale[nb_dist++] = 1.;
            }
         }
         else
            if ((characteristics != NULL) && (characteristics->nb_zones != NULL))
            {
               phisto[nb_histo]= characteristics->nb_zones[val];
               index_dist[nb_histo++]= I_DEFAULT;
            }

        if ((nb_occurrences != NULL) && (nb_occurrences[val] != NULL))
        {
           pdist[nb_dist] = nb_occurrences[val];

           if ((characteristics != NULL) && (characteristics->nb_occurrences != NULL))
           {
              phisto[nb_histo]= characteristics->nb_occurrences[val];
              index_dist[nb_histo]= nb_dist;
              dist_nb_value[nb_dist]= nb_occurrences[val]->plot_nb_value_computation(phisto[nb_histo]);
              scale[nb_dist++]= phisto[nb_histo++]->nb_element;
           }
           else
           {
              dist_nb_value[nb_dist] = nb_occurrences[val]->plot_nb_value_computation();
              scale[nb_dist++] = 1.;
           }
        }

        else
           if ((characteristics != NULL) && (characteristics->nb_occurrences != NULL))
           {
              phisto[nb_histo]= characteristics->nb_occurrences[val];
              index_dist[nb_histo++]= I_DEFAULT;
           }
      }
   }
   */


   if (observation != NULL)
   {
      if (empirical_observation != NULL)
      {
         for(val = 0; val < nb_state; val++)
         {
            plot.variable[index] = process;
            plot.viewpoint[index] = OBSERVATION;

            title.str("");
            title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
            plot[index].title = title.str();

            plot[index].xrange = Range(0 , observation[val]->nb_value - 1);
            if (observation[val]->nb_value - 1 < TIC_THRESHOLD)
               plot[index].xtics = 1;

            if (empirical_observation[val]->nb_element > 0)
            {
               scale = empirical_observation[val]->nb_element;
               plot[index].yrange = Range(0, ceil(MAX(empirical_observation[val]->max,
                                                      observation[val]->max * scale) * YSCALE));
               plot[index].resize(2);
               legend.str("");
               legend <<STAT_label[STATL_STATE] << " " << val << " "
                        << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
               plot[index][j].legend = legend.str();
               plot[index][j].style = "impulses";

               empirical_observation[val]->plotable_frequency_write(plot[index][0]);
               j = 1;
            }
            else
            {
               scale = 1.;
               plot[index].yrange = Range(0, MIN(observation[val]->max * YSCALE,
                                                 1.));
               plot[index].resize(1);
               j = 0;
            }

            legend.str("");
            legend <<STAT_label[STATL_STATE] << " " << val << " "
                     << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
            plot[index][j].legend = legend.str();
            plot[index][j].style = "linepoints";

            observation[val]->plotable_mass_write(plot[index][j], scale);
            index++;
         }
      }
      else // (empirical_observation == NULL)
      {
         plot.variable[index] = process;
         plot.viewpoint[index] = OBSERVATION;

         title.str("");
         title << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
         plot[index].title = title.str();

         plot[index].xrange = Range(0, nb_value - 1);
         if (nb_value - 1 < TIC_THRESHOLD)
            plot[index].xtics = 1;

         max = observation[0]->max;
         for(val = 1; val < nb_state; val++)
            if (observation[val]->max > max)
               max = observation[val]->max;

         plot[index].yrange = Range(0, MIN(max * YSCALE, 1.));
         plot[index].resize(nb_state);

         for(val = 0; val < nb_state; val++)
         {
            legend.str("");
            legend <<STAT_label[STATL_STATE] << " " << val << " "
                     << STAT_label[STATL_OBSERVATION] << " " << STAT_label[STATL_DISTRIBUTION];
            plot[index][val].legend = legend.str();
            plot[index][val].style = "linepoints";

            observation[val]->plotable_mass_write(plot[index][val]);
         }
         index++;
      }
   }

   /*

      // first occurrence_leaves
      if ((first_occurrence_leaves != NULL) || (characteristics != NULL))
      {
         plot_index++;

         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
               case 0 :
                  file_name[0] << prefix << process << plot_index << ".plot";
                  break;
               case 1 :
                  file_name[0] << prefix << process << plot_index << ".print";
                  break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << process << plot_index << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";

            if ((title != NULL) || (process > 0))
            {
               out_file << " \"";
               if (title != NULL)
               {
                  out_file << title;
                  if (process > 0)
                     out_file << " - ";
               }
               if (process > 0)
                  out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
              out_file << "\"";
            }
            out_file << "\n\n";

            j= histo_index;
            k= dist_index;

            start= true;
            for(val= 0; val < nb_value; val++)
            {
               if ((first_occurrence_leaves != NULL) &&
                   (first_occurrence_leaves[val] != NULL))
               {
                 if (!start)
                 {
                    if (i == 0)
                       out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                 }
                 else
                   start= false;

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                    out_file << "set xtics 0,1" << endl;

               if ((characteristics != NULL) &&
                   (characteristics->first_occurrence_leaves != NULL) &&
                   (characteristics->first_occurrence_leaves[val]->nb_element > 0))
               {
                  out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                           << (int)(MAX(phisto[j]->max, pdist[k]->max*scale[k])*YSCALE)+1
                           << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses,\\" << endl;
                  out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES] << " " << STAT_label[STATL_DISTRIBUTION];
                  first_occurrence_leaves[val]->plot_title_print(out_file);
                  out_file << "\" with linespoints" << endl;
                  j++;
               }

               else
               {
                  out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                           << MIN(pdist[k]->max*YSCALE, 1.) << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES] << " " << STAT_label[STATL_DISTRIBUTION];
                  first_occurrence_leaves[val]->plot_title_print(out_file);
                  out_file << "\" with linespoints" << endl;
               }

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;
               k++;
            }
            else // first_occurrence_leaves != NULL
               if ((characteristics != NULL) &&
                   (characteristics->first_occurrence_leaves[val]->nb_element > 0))
               {
                  if (!start)
                  {
                    if (i == 0)
                      out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                  }
                  else
                     start = false;

                  if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics 0,1" << endl;

                  out_file << "plot [0:" << phisto[j]->nb_value-1 << "] [0:"
                           << (int)(phisto[j]->max*YSCALE)+1 << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[process == 0 ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses" << endl;

                  if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;

                  if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics autofreq" << endl;
                  j++;
               }
            } // end for val
            if (i == 1)
               out_file << "\nset terminal x11" << endl;

            out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
         } // end for i
         histo_index= j;
         dist_index= k;
      }


      // sojourn_size
      if ((sojourn_size != NULL) || (characteristics != NULL))
      {
         plot_index++;

         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
               case 0 :
                  file_name[0] << prefix << process << plot_index << ".plot";
                  break;
               case 1 :
                  file_name[0] << prefix << process << plot_index << ".print";
                  break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << process << plot_index << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";

            if ((title != NULL) || (process > 0))
            {
               out_file << " \"";
               if (title != NULL)
               {
                  out_file << title;
                  if (process > 0)
                     out_file << " - ";
               }
               if (process > 0)
                  out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
              out_file << "\"";
            }
            out_file << "\n\n";

            j= histo_index;
            k= dist_index;

            start= true;
            for(val= 0; val < nb_value; val++)
            {
               if ((sojourn_size != NULL) && (sojourn_size[val] != NULL))
               {
                 if (!start)
                 {
                    if (i == 0)
                       out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                 }
                 else
                   start= false;

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                    out_file << "set xtics 0,1" << endl;

               if ((characteristics != NULL) && (characteristics->sojourn_size != NULL) &&
                   (characteristics->sojourn_size[val]->nb_element > 0))
               {
                  out_file << "plot [0:" << dist_nb_value[k] - 1 << "] [0:"
                           << (int)(MAX(phisto[j]->max, pdist[k]->max*scale[k])*YSCALE)+1
                           << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses,\\" << endl;
                  out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_DISTRIBUTION];
                  sojourn_size[val]->plot_title_print(out_file);
                  out_file << "\" with linespoints" << endl;
                  j++;
               }

               else
               {
                  out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                           << MIN(pdist[k]->max*YSCALE, 1.) << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_DISTRIBUTION];
                  sojourn_size[val]->plot_title_print(out_file);
                  out_file << "\" with linespoints" << endl;
               }

               if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;
               k++;
            }
            else // sojourn_size != NULL
               if ((characteristics != NULL) && (characteristics->sojourn_size[val]->nb_element > 0))
               {
                  if (!start)
                  {
                    if (i == 0)
                      out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                    out_file << endl;
                  }
                  else
                     start = false;

                  if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics 0,1" << endl;

                  out_file << "plot [0:" << phisto[j]->nb_value-1 << "] [0:"
                           << (int)(phisto[j]->max*YSCALE)+1 << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                           << " title \"" << STAT_label[process == 0 ? STATL_STATE : STATL_OUTPUT] << " " << val
                           << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses" << endl;

                  if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;

                  if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics autofreq" << endl;
                  j++;
               }
            } // end for val
            if (i == 1)
               out_file << "\nset terminal x11" << endl;

            out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
         } // end for i
         histo_index= j;
         dist_index= k;
      }

      // counting
      if ((nb_zones != NULL) || (nb_occurrences != NULL) ||
          (characteristics != NULL))
      {
         plot_index++;

         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
               case 0 :
                  file_name[0] << prefix << process << plot_index << ".plot";
                  break;
               case 1 :
                  file_name[0] << prefix << process << plot_index << ".print";
                  break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << process << plot_index << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";
            if ((title != NULL) || (process > 0))
            {
               out_file << " \"";
               if (title != NULL)
               {
                  out_file << title;
                  if (process > 0)
                     out_file << " - ";
               }
               if (process > 0)
                  out_file << STAT_label[STATL_OUTPUT_PROCESS] << " " << process;
               out_file << "\"";
            }
            out_file << "\n\n";

            j= histo_index;
            k= dist_index;

            start= true;
            for(val= 0; val < nb_value; val++)
            {
               if (nb_zones != NULL)
               {
                  if (!start)
                  {
                     if (i == 0)
                        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                     out_file << endl;
                  }
                  else
                     start= false;

                  if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if ((characteristics != NULL) && (characteristics->nb_zones != NULL))
                  {
                     out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                              << (int)(MAX(phisto[j]->max, pdist[k]->max * scale[k])*YSCALE)+1
                              << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                              << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                              << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                              << "\" with impulses,\\" << endl;
                     out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1;
                     // uncomment if a distribution for the size is added
                     out_file << " title \"" << STAT_TREES_label[TREESTATL_MIXTURE_OF]
                              << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                              << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_DISTRIBUTION];
                     out_file << "\" with linespoints" << endl;
                     j++;
                  }
                  if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;
                  k++;
               }
               else // nb_zones == NULL
                  if (characteristics != NULL)
                  {
                     if (!start)
                     {
                        if (i == 0)
                           out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                        out_file << endl;
                     }
                     else
                        start = false;

                     if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                        out_file << "set xtics 0,1" << endl;

                     if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                        out_file << "set ytics 0,1" << endl;

                     out_file << "plot [0:" << phisto[j]->nb_value-1 << "] [0:"
                              << (int)(phisto[j]->max*YSCALE)+1 << "] \""
                              << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                              << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                              << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                              << "\" with impulses" << endl;

                     if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                        out_file << "set xtics autofreq" << endl;

                     if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                        out_file << "set ytics autofreq" << endl;

                     j++;
                  }

                 if (nb_occurrences != NULL)
                 {
                    if (!start)
                    {
                       if (i == 0)
                          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                       out_file << endl;
                    }
                    else
                       start= false;

                    if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                       out_file << "set xtics 0,1" << endl;

                    if ((characteristics != NULL) && (characteristics->nb_occurrences != NULL))
                    {
                       out_file << "plot [0:" << dist_nb_value[k]-1 << "] [0:"
                                << (int)(MAX(phisto[j]->max, pdist[k]->max * scale[k])*YSCALE)+1
                                << "] \"" << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                                << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                                << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                                << "\" with impulses,\\" << endl;
                       out_file << "\"" << label((data_file_name[1].str()).c_str()) << "\" using " << nb_histo+k+1;
                       out_file << " title \"" << STAT_TREES_label[TREESTATL_MIXTURE_OF]
                                << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                                << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_DISTRIBUTION];
                       out_file << "\" with linespoints" << endl;
                       j++;
                    }
                    // uncomment if a distribution for the size is added
                    if (dist_nb_value[k]-1 < TIC_THRESHOLD)
                       out_file << "set xtics autofreq" << endl;
                    k++;
                 }
                 else // nb_occurrences == NULL
                    if (characteristics != NULL)
                    {
                       if (i == 0)
                          out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                       out_file << endl;

                       if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                          out_file << "set xtics 0,1" << endl;

                       if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                          out_file << "set ytics 0,1" << endl;

                       out_file << "plot [0:" << phisto[j]->nb_value-1 << "] [0:"
                                << (int)(phisto[j]->max*YSCALE)+1 << "] \""
                                << label((data_file_name[1].str()).c_str()) << "\" using " << j+1
                                << " title \"" << STAT_TREES_label[process == 0 ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                                << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                                << "\" with impulses" << endl;

                       if (phisto[j]->nb_value-1 < TIC_THRESHOLD)
                          out_file << "set xtics autofreq" << endl;

                       if ((int)(phisto[j]->max*YSCALE)+1 < TIC_THRESHOLD)
                          out_file << "set ytics autofreq" << endl;
                       j++;
                    }
               } // end nb_zones == NULL

               if (characteristics != NULL)
               {
                  if (i == 0)
                     out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                  out_file << endl;

                  if (hsize->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if ((int)(hsize->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics 0,1" << endl;

                  out_file << "plot [0:" << hsize->nb_value - 1 << "] [0:"
                           << (int)(hsize->max*YSCALE)+1 << "] \""
                           << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                           << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                           << "\" with impulses" << endl;

                  if (hsize->nb_value-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;

                  if ((int)(hsize->max*YSCALE)+1 < TIC_THRESHOLD)
                     out_file << "set ytics autofreq" << endl;
               }

               if (i == 1)
                  out_file << "\nset terminal x11" << endl;

               out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
            } // end for val

         histo_index= j;
         dist_index= k;
      } */

}

/*****************************************************************
 *
 *  Return the number of views (i.e. the size) in Matplotlib output
 *  of CategoricalTreeProcess, using
 *  the identifier of the considered process
 *  the empirical observation distributions,
 *  the empirical tree characteristics and
 *  the tree size frequency distribution
 *
 **/

unsigned int CategoricalTreeProcess::nb_plot_set_computation(int process, FrequencyDistribution * const * empirical_observation,
                                                             const TreeCharacteristics * characteristics,
                                                             const FrequencyDistribution * hsize) const
{
   unsigned int nb_plot_set = 0;
   register int val;

   // create a list of the successive distributions to be printed
   if ((first_occurrence_root != NULL) || (characteristics != NULL))
      for(val = 0; val < nb_value; val++)
      {
         if ((first_occurrence_root != NULL) || (first_occurrence_root[val] != NULL)
             || ((characteristics != NULL) && (val < characteristics->get_nb_values())
                && (characteristics->first_occurrence_root[val]->nb_element > 0)))
               nb_plot_set++;
      }

   // add other viewpoints, etc.

   if (observation != NULL)
      nb_plot_set += nb_state;

   return nb_plot_set;

}

/*****************************************************************
 *
 *  Test on the hidden character of the model
 *  using the number of output processes and pointeurs on these
 *
 **/

bool Stat_trees::test_hidden(int nb_output_process, CategoricalTreeProcess** process)

{
   bool hidden= false;
   register int i;


   for(i= 1; i <= nb_output_process; i++)
   {
       if (process[i] != NULL)
          hidden= process[i]->test_hidden();
       if (hidden)
         break;
   }

   return hidden;
}

/*****************************************************************
 *
 *  Default constructor of HiddenMarkovTree class
 *
 **/

HiddenMarkovTree::HiddenMarkovTree()
 : markov_data(NULL)
 , _ch_order(1)
 , self_row(NULL)
 , _nb_ioutput_process(0)
 , _nb_doutput_process(0)
 , npprocess(NULL)
 , piprocess(NULL)
 , pdprocess(NULL)
{}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovTree class
 *  using the type of markov tree ('o'rdinary or 'e'quilibrium),
 *  the number of states, the number of observed integral
 *  and floating processes, the (children) order
 *  and the number of observed values for each variable
 *  (ending with the floating variables)
 *
 **/

HiddenMarkovTree::HiddenMarkovTree(char itype, int inb_state,
                                   int ich_order,
                                   int inb_ioutput_process,
                                   int inb_doutput_process,
                                   int* nb_value,
                                   bool* force_param)
 : Chain(itype, inb_state, (int)pow((double)inb_state, (double)ich_order), true)
 , markov_data(NULL)
 , _ch_order(ich_order)
 , self_row(NULL)
 , _nb_ioutput_process(inb_ioutput_process)
 , _nb_doutput_process(inb_doutput_process)
 , npprocess(NULL)
 , piprocess(NULL)
 , pdprocess(NULL)
{
   register int i, s, v;
   double sum= 0.0, inc;
   bool *fparam= NULL;

   self_row= new int[nb_state];
   self_row_computation();

   fparam= new bool[_nb_ioutput_process];
   if (force_param == NULL)
   {
      for (i= 0; i < _nb_ioutput_process; i++)
         fparam[i]= false;
   }
   else
   {
      for (i= 0; i < _nb_ioutput_process; i++)
         fparam[i]= force_param[i];
   }

   npprocess= new CategoricalTreeProcess*[_nb_ioutput_process+1];
   npprocess[0]= new CategoricalTreeProcess(nb_state, nb_state, false);
   // since the number of values is equal to the number of states,
   // npprocess[0] represents the hidden state process

   piprocess= new DiscreteParametricProcess*[_nb_ioutput_process+1];
   piprocess[0]= NULL;

   for(i= 0; i < _nb_ioutput_process; i++)
      if ((*nb_value <= NB_OUTPUT) && !(fparam[i]))
      {
         npprocess[i+1]= new CategoricalTreeProcess(nb_state, *nb_value++, true);
         piprocess[i+1]= NULL;
      }
      else
      {
         npprocess[i+1]= NULL;
         piprocess[i+1]= new DiscreteParametricProcess(nb_state, (int)(*nb_value++ * SAMPLE_NB_VALUE_COEFF));
      }

   pdprocess= new DiscreteParametricProcess*[_nb_doutput_process+1];
   pdprocess[0]= NULL;
   for(i= 0; i < _nb_doutput_process; i++)
      pdprocess[i+1]= new DiscreteParametricProcess(nb_state, (int)(*nb_value++ * SAMPLE_NB_VALUE_COEFF));

   delete [] fparam;
   fparam= NULL;
}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovTree class
 *  using a Chain object, the "children order",
 *  the number of observed integral processes, the processes themselves,
 *  the tree size and a flag on the counting distribution computation
 *
 **/

HiddenMarkovTree::HiddenMarkovTree(const Chain * pchain, int ich_order,
                                   int inb_ioutput_process,
                                   CategoricalProcess** pobservation,
                                   int size, bool counting_flag)
 : Chain(*pchain)
 , markov_data(NULL)
 , _ch_order(ich_order)
 , self_row(NULL)
 , _nb_ioutput_process(inb_ioutput_process)
 , _nb_doutput_process(0)
 , npprocess(NULL)
 , piprocess(NULL)
 , pdprocess(NULL)
{
   register int i;

   self_row= new int[nb_state];
   self_row_computation();

   npprocess= new CategoricalTreeProcess*[_nb_ioutput_process+1];
   npprocess[0]= new CategoricalTreeProcess(nb_state, nb_state, false);
   // index 0 corresponds to the hidden process for npprocess

   for(i= 0; i < _nb_ioutput_process; i++)
      npprocess[i+1]= new CategoricalTreeProcess(*(pobservation[i]));

   // index 0 corresponds to the first observed process for pobservation

   /* for(i= 0; i < nb_state; i++)
   {
      if (transition[i][i] < 1.)
         npprocess[0]->absorption[i]= 0.;
      else
         npprocess[0]->absorption[i]= 1.;
   }*/
   // This may be replaced by something more relevant

   if (size > COUNTING_MAX_SIZE)
      counting_flag= false;

   characteristic_computation(size, counting_flag);

}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovTree class
 *  using a Chain object, the "children order",
 *  the number of observed integral processes, the processes themselves,
 *  the tree size and a flag on the counting distribution computation
 *  and a flag on the counting distribution computation
 *
 **/

HiddenMarkovTree::HiddenMarkovTree(const Chain * pchain, int ich_order,
                                   int inb_ioutput_process,
                                   int inb_doutput_process,
                                   CategoricalProcess** categorical_observation,
                                   DiscreteParametricProcess** iparametric_observation,
                                   DiscreteParametricProcess** dparametric_observation,
                                   int size, bool counting_flag)
 : Chain(*pchain)
 , markov_data(NULL)
 , _ch_order(ich_order)
 , self_row(NULL)
 , _nb_ioutput_process(inb_ioutput_process)
 , _nb_doutput_process(inb_doutput_process)
 , npprocess(NULL)
 , piprocess(NULL)
 , pdprocess(NULL)
{
   register int i;

   self_row= new int[nb_state];
   self_row_computation();

   npprocess= new CategoricalTreeProcess*[_nb_ioutput_process+1];
   npprocess[0]= new CategoricalTreeProcess(nb_state, nb_state, false);

   piprocess= new DiscreteParametricProcess*[_nb_ioutput_process+1];
   piprocess[0]= NULL;

   for(i= 0; i < _nb_ioutput_process; i++)
   {
      if (categorical_observation[i] != NULL)
      {
         npprocess[i+1]= new CategoricalTreeProcess(*(categorical_observation[i]));
         piprocess[i+1]= NULL;
      }
      else
      {
         npprocess[i+1]= NULL;
         piprocess[i+1]= new DiscreteParametricProcess(*(iparametric_observation[i]));
      }
   }

   /* for(i= 0; i < nb_state; i++)
   {
      if (transition[i][i] < 1.)
         npprocess[0]->absorption[i]= 0.;
      else
         npprocess[0]->absorption[i]= 1.;
   }*/
   // This might be replaced by something more relevant

   if (size > COUNTING_MAX_SIZE)
      counting_flag= false;

   characteristic_computation(size, counting_flag);
}

/*****************************************************************
 *
 *  Copy constructor of HiddenMarkovTree class
 *  using a flag on the HiddenMarkovTreeData copy
 *  and one on the characteristic distribution copy
 *
 **/

HiddenMarkovTree::HiddenMarkovTree(const HiddenMarkovTree& markov, bool data_flag,
                                   bool characteristic_flag)
 : StatInterface()
 , Chain(markov)
 , markov_data(NULL)
 , _ch_order(0)
 , self_row(NULL)
 , _nb_ioutput_process(0)
 , _nb_doutput_process(0)
 , npprocess(NULL)
 , piprocess(NULL)
 , pdprocess(NULL)

{ copy(markov, data_flag, characteristic_flag); }

/*****************************************************************
 *
 *  Destructor for HiddenMarkovTree class
 *
 **/

HiddenMarkovTree::~HiddenMarkovTree()
{ remove(); }

/*****************************************************************
 *
 *  Assignement operator for HiddenMarkovTree class
 *
 **/

HiddenMarkovTree& HiddenMarkovTree::operator=(const HiddenMarkovTree& markov)
{
  if (&markov != this)
  {
     remove();
     Chain::remove();

     Chain::copy(markov);
     copy(markov);
  }

  return *this;
}

/*****************************************************************
 *
 *  Return a copy of an instance of HiddenMarkovTree with the same
 *  dynamic class
 *
 **/

HiddenMarkovTree* HiddenMarkovTree::HiddenMarkovTreeCopy(bool data_flag,
                                                         bool characteristic_flag) const
{
   HiddenMarkovTree *res = new HiddenMarkovTree(*this, data_flag, characteristic_flag) ;
   return res;
}


/*****************************************************************
 *
 *  Extraction of a distribution for HiddenMarkovTree class
 *  using a StatError object, the type of distribution,
 *  the considered variable and value (or state)
 *
 **/

DiscreteParametricModel* HiddenMarkovTree::extract(StatError& error,
                                                   int type,
                                                   int ivariable,
                                                   int value) const
{
   bool status= true;
   int hvariable= 0;
   int variable= ivariable;
   Distribution *pdist;
   // DiscreteParametric *pparam;
   DiscreteParametricModel *dist;
   FrequencyDistribution *phisto;


   dist= NULL;
   error.init();

   pdist= NULL;

   if (type == OBSERVATION)
   {
      if ((ivariable < 1) || (ivariable > _nb_ioutput_process))
      {
         status= false;
         error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
      }
      else
      {
         if ((value < 0) || (value >= nb_state))
         {
            status= false;
            std::ostringstream error_message;
            error_message << STAT_label[STATL_STATE] << " " << value << " "
                          << STAT_error[STATR_NOT_PRESENT];
            error.update((error_message.str()).c_str());
         }
         else
            if (npprocess[variable] == NULL)
               pdist= piprocess[variable]->observation[value];
            else
               pdist= npprocess[variable]->observation[value];
      }
   }
   else // type != OBSERVATION
   {
      if ((ivariable < 0) || (ivariable > _nb_ioutput_process))
      {
         status= false;
         error.update(STAT_error[STATR_OUTPUT_PROCESS_INDEX]);
      }

      else
      {
         if ((value < 0) || (value >= npprocess[variable]->nb_value))
         {
            status= false;
            ostringstream error_message;
            error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                          << value << " " << STAT_error[STATR_NOT_PRESENT];
            error.update((error_message.str()).c_str());
         }
         else
         {
            switch (type)
            {
               case FIRST_OCCURRENCE_ROOT :
                  pdist= npprocess[variable]->first_occurrence_root[value];
                  break;
               case FIRST_OCCURRENCE_LEAVES :
                  pdist= npprocess[variable]->first_occurrence_leaves[value];
                  break;
               case SOJOURN_SIZE :
                  pdist= npprocess[variable]->sojourn_size[value];
                  break;
               case NB_ZONES :
                  pdist= npprocess[variable]->nb_zones[value];
                  break;
               case NB_OCCURRENCES :
                  pdist= npprocess[variable]->nb_occurrences[value];
                  break;
            }

            if (pdist == NULL)
            {
               status= false;
               error.update(STAT_TREES_error[TREESTATR_NON_EXISTING_CHARACTERISTIC_DISTRIBUTION]);
            }
         }
      }
   }

   if (status)
   {
      phisto= NULL;

      if (markov_data != NULL)
      {
         switch (markov_data->_type[0])
         {
            case INT_VALUE :
               hvariable= variable - 1;
               break;
            case STATE :
               hvariable= variable;
               break;
         }

         if (hvariable >= 0)
         {
            switch (type)
            {
               case OBSERVATION :
               {
                 if ((markov_data->observation_distribution != NULL)
                    && (markov_data->observation_distribution[hvariable] != NULL))
                   phisto= markov_data->observation_distribution[hvariable][value];
                 break;
               }

               case FIRST_OCCURRENCE_ROOT :
               {
                  phisto= markov_data->characteristics[hvariable]->first_occurrence_root[value];
                  break;
               }

               case FIRST_OCCURRENCE_LEAVES :
               {
                  phisto= markov_data->characteristics[hvariable]->first_occurrence_leaves[value];
                  break;
               }

               case SOJOURN_SIZE :
               {
                  phisto= markov_data->characteristics[hvariable]->sojourn_size[value];
                  break;
               }

               case NB_ZONES :
               {
                  phisto= markov_data->characteristics[hvariable]->nb_zones[value];
                  break;
               }

               case NB_OCCURRENCES :
               {
                  phisto= markov_data->characteristics[hvariable]->nb_occurrences[value];
                  break;
               }
            }
         }
      }

      if (pdist != NULL)
        dist= new DiscreteParametricModel(*pdist, phisto);
   }

   return dist;
}

/*****************************************************************
 *
 *  Return the data part of a HiddenMarkovTree
 *  using a StatError object, keeping a reference on self
 *
 **/

HiddenMarkovTreeData* HiddenMarkovTree::extract_data(StatError& error) const
{
   bool status= true;
   HiddenMarkovTreeData *tree= NULL;

   error.init();

   if (markov_data == NULL)
   {
      status= false;
      error.update(STAT_error[STATR_NO_DATA]);
   }

   if ((_nb_ioutput_process != markov_data->_nb_integral)
       || (_nb_doutput_process != markov_data->_nb_float))
   {
      status= false;
      error.update(STAT_TREES_error[TREESTATR_STATE_TREES]);
   }

   if (status)
   {
      tree= new HiddenMarkovTreeData(*markov_data);
      tree->markov= this->HiddenMarkovTreeCopy(false);
   }

   return tree;
}

/*****************************************************************
 *
 *  Parameter thresholding for HiddenMarkovTree class
 *  using the min probability
 *
 **/

HiddenMarkovTree* HiddenMarkovTree::thresholding(double min_probability) const
{
   register int i;
   HiddenMarkovTree *markov;

   markov= this->HiddenMarkovTreeCopy(false, false);
   markov->Chain::thresholding(min_probability);

   for(i= 0; i < _nb_ioutput_process; i++)
      npprocess[i+1]->thresholding(min_probability);

   markov->component_computation();
   // check the meaning of this

   return markov;
}

/*****************************************************************
 *
 *  Permutation of the states of a HiddenMarkovTree based on
 *  a given permutation perm (which must be the same length as
 *  the number of states)
 *
 **/

void HiddenMarkovTree::state_permutation(StatError& error,
                                         int* perm) const
{
   register int i, j;
   bool status= true;
   // each element of check_perm must be used exactly once
   bool * check_perm= new bool[nb_state];
   double* pinitial= NULL;
   double** ptransition= NULL;

   // check permutation
   error.init();

   for (i= 0; i < nb_state; i++)
      check_perm[i]= false; // indicates if ith element was used in perm

   for (i= 0; i < nb_state; i++)
   {
      if (check_perm[perm[i]])
         status= false;
      else
         check_perm[perm[i]]= true;
   }
   delete [] check_perm;
   check_perm= NULL;

   if (!status)
      error.update(STAT_TREES_error[TREESTATR_NO_PERMUTATION]);
   else
   {
      // permutation of initial probabilities
      pinitial= new double[nb_state];
      for (i= 0; i < nb_state; i++)
         pinitial[perm[i]]= initial[i];
      for (i= 0; i < nb_state; i++)
         initial[i]= pinitial[i];
      delete [] pinitial;
      pinitial= NULL;

      // permutation of transition probabilities
      ptransition= new double*[nb_state];
      for (i= 0; i < nb_state; i++)
         ptransition[i]= new double[nb_state];
      for (i= 0; i < nb_state; i++)
         for (j= 0; j < nb_state; j++)
            ptransition[perm[i]][perm[j]]= transition[i][j];
      for (i= 0; i < nb_state; i++)
         for (j= 0; j < nb_state; j++)
            transition[i][j]= ptransition[i][j];
      for (i= 0; i < nb_state; i++)
         delete [] ptransition[i];
      delete [] ptransition;
      ptransition= NULL;

      if (self_row != NULL)
      {
         // permutation of auto-transition probabilities
         int* pself_row= new int[nb_state];
         for (i= 0; i < nb_state; i++)
            pself_row[perm[i]]= self_row[i];
         for (i= 0; i < nb_state; i++)
            self_row[i]= pself_row[i];
         delete [] pself_row;
         pself_row= NULL;
      }
      // permutation of observation distributions
      if (npprocess[0] != NULL)
         npprocess[0]->state_permutation(perm);
      for (i= 1; i <= _nb_ioutput_process; i++)
      {
         if (npprocess[i] != NULL)
            npprocess[i]->state_permutation(perm);
         if (piprocess[i] != NULL)
            piprocess[i]->state_permutation(perm);
      }
      for (i= 1; i <= _nb_doutput_process; i++)
         if (pdprocess[i] != NULL)
            pdprocess[i]->state_permutation(perm);

      if (markov_data != NULL)
         markov_data->state_permutation(perm);
   }
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovTree on a single line
 *  using an output stream
 *
 **/

ostream& HiddenMarkovTree::line_write(ostream& os) const
{
   os << nb_state << " " << STAT_word[STATW_STATES];

//    os << "   " << SEQ_label[SEQL_MAX_ORDER] << " " << _ch_order;

   return os;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovTree
 *  using an output stream and a flag on the level of detail
 *
 **/

ostream& HiddenMarkovTree::ascii_write(ostream& os, bool exhaustive) const
{ return ascii_write(os, markov_data, exhaustive, false); }


/*****************************************************************
 *
 *  Prints a HiddenMarkovTree into a file
 *  using a StatError object, the path
 *  and a flag on the level of detail
 *
 **/

bool HiddenMarkovTree::ascii_write(StatError& error, const char * path,
                                   bool exhaustive) const
{
   bool status;
   ofstream out_file(path);

   error.init();

   if (!out_file)
   {
      status= false;
      error.update(STAT_error[STATR_FILE_NAME]);
   }

   else
   {
      status= true;
      ascii_write(out_file, markov_data, exhaustive, true);
   }

   return status;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovTree in a spreadsheet fashion
 *  using a StatError object and the path
 *
 **/

bool HiddenMarkovTree::spreadsheet_write(StatError& error,
                                         const char * path) const
{
   bool status;
   ofstream out_file(path);


   error.init();

   if (!out_file)
   {
      status= false;
      error.update(STAT_error[STATR_FILE_NAME]);
   }

   else
   {
      status= true;
      spreadsheet_write(out_file, markov_data);
   }

   return status;
}

/*****************************************************************
 *
 *  Gnuplot output for HiddenMarkovTree class
 *  using a StatError object, a prefix for the files
 *  and the title of output figures
 *
 **/

bool HiddenMarkovTree::plot_write(StatError& error,
                                  const char * prefix,
                                  const char * title) const
{
   bool status= plot_write(prefix, title, markov_data);

   error.init();

   if (!status)
     error.update(STAT_error[STATR_FILE_PREFIX]);

   return status;
}

/*****************************************************************
 *
 *  Computation of the characteristic distributions
 *  for HiddenMarkovTree class
 *  using a HiddenMarkovTreeData object,
 *  a flag on the computation of counting distributions,
 *  the considered variable
 *  and a flag indicating whether the size of the trees (or other
 *  quantities for which the characteristic distributions are invariant
 *  - if any) are taken into account
 *
 **/

void HiddenMarkovTree::characteristic_computation(const HiddenMarkovTreeData& tree,
                                                  bool counting_flag,
                                                  int variable,
                                                  bool size_flag)
{
   if (nb_component > 0)
   {
      register int i, j, k;
      int nb_value;
      double *memory= NULL;
      Distribution dsize;

      if (tree.hsize != NULL)
         dsize= *(tree.hsize);

      // computation of the characteristic distributions at state level

      if (((variable == I_DEFAULT) || (variable == 0)) &&
          // the 4 conditions beloww are nested
          ((!size_flag) || ((size_flag) && ((npprocess[0]->size == NULL)
           || (dsize != *(npprocess[0]->size))))))
      {
         npprocess[0]->create_characteristic(dsize, true, counting_flag);

         memory= memory_computation();
         //index_state_distribution();

         for (i= 0; i < nb_state; i++)
         {
            state_no_occurrence_probability(i);
            if (tree.state_characteristics != NULL)
               nb_value= tree.state_characteristics->_max_value;
            else
               nb_value= 1;

            state_first_occurrence_root_distribution(i, nb_value);
            state_first_occurrence_leaves_distribution(i, nb_value);

            // state_first_occurrence_root_distribution(i);
            // state_first_occurrence_leaves_distribution(i);

            state_leave_probability(memory, i);

            if (state_type[i] != 'a')
            {
               if (tree.state_characteristics != NULL)
                  nb_value= tree.state_characteristics->_max_value;
               else
                  nb_value= 1;

               state_sojourn_size_distribution(memory, i, nb_value);

               // state_sojourn_size_distribution(memory, i);
            }
            else
            {
               npprocess[0]->absorption[i]= 1.;
               delete npprocess[0]->sojourn_size[i];
               npprocess[0]->sojourn_size[i]= NULL;
            }

            if (counting_flag)
            {
               state_nb_pattern_mixture(i, 'r');
               state_nb_pattern_mixture(i, 'o');
            }
         }
      }

      // computation of the characteristic distributions at observation level

      for(i= 1; i <= _nb_ioutput_process; i++)
      {
         if ((npprocess[i] != NULL) && (((variable == I_DEFAULT) || (i == variable) &&
             ((!size_flag) || ((size_flag) && ((npprocess[i]->size == NULL) ||
             ((npprocess[i]->size != NULL) && (dsize != *(npprocess[i]->size)))))))))
         {
            npprocess[i]->create_characteristic(dsize, true, counting_flag);

            if (memory == NULL)
               memory = memory_computation();
            // index_output_distribution(i);

            for (j= 0; j < npprocess[i]->nb_value; j++)
            {
               output_no_occurrence_probability(i, j);
               if (tree.characteristics[i-1] != NULL)
                  nb_value= tree.characteristics[i-1]->_max_value;
               else
                  nb_value= 1;

               output_first_occurrence_root_distribution(i, j, nb_value);
               output_first_occurrence_leaves_distribution(i, j, nb_value);

               output_leave_probability(memory, i, j);

               for (k= 0; k < nb_state; k++)
                  if ((npprocess[i]->observation[k]->mass[j] > 0.) &&
                     ((state_type[k] != 'a') || (npprocess[i]->observation[k]->mass[j] < 1.)))
                     break;

               if (k < nb_state)
               {
                  if (tree.characteristics[i-1] != NULL)
                     nb_value= tree.characteristics[i-1]->_max_value;
                  else
                     nb_value= 1;

                  output_sojourn_size_distribution(memory, i, j, nb_value);
               }
               else
               {
                  npprocess[i]->absorption[j]= 1.;
                  delete npprocess[i]->sojourn_size[j];
                  npprocess[i]->sojourn_size[j]= NULL;
               }

               if (counting_flag)
               {
                  output_nb_zones_mixture(i, j);
                  output_nb_occurrences_mixture(i, j);
               }
            }
         }
      }
      delete [] memory;
   }
}

// virtual functions common to all hidden_markov_trees

double HiddenMarkovTree::likelihood_computation(const Trees& otrees,
                                                int index) const
{ return D_INF; }

double HiddenMarkovTree::likelihood_computation(HiddenMarkovTreeData& otrees) const
{ return D_INF; }

HiddenMarkovTreeData*
HiddenMarkovTree::state_tree_computation(StatError& error,
                                         const Trees& trees,
                                         int algorithm,
                                         bool characteristic_flag,
                                         int index) const
{ return NULL; }

// access to class members

/*****************************************************************
 *
 *  Access to the HiddenMarkovTreeData object
 *  for HiddenMarkovTree class
 *
 **/

HiddenMarkovTreeData* HiddenMarkovTree::get_markov_data() const
{
   HiddenMarkovTreeData *res= NULL;

   if (markov_data != NULL)
      res= new HiddenMarkovTreeData(*markov_data);
      // check the flag argument of the constructor
   return res;
}

/*****************************************************************
 *
 *  Access to children order of a HiddenMarkovTree
 *
 **/

int HiddenMarkovTree::get_ch_order() const
{ return _ch_order; }

/*****************************************************************
 *
 *  Access to the number of states of a HiddenMarkovTree
 *
 **/

int HiddenMarkovTree::get_nb_state() const
{ return nb_state; }

/*****************************************************************
 *
 *  Access to the probabilities to stay in a given state
 *  for HiddenMarkovTree class
 *
 **/

int HiddenMarkovTree::get_self_row(int state) const
{ return self_row[state]; }

/*****************************************************************
 *
 *  Access to the number of observed processes
 *  for HiddenMarkovTree class
 *
 **/

int HiddenMarkovTree::get_nb_ioutput_process() const
{ return _nb_ioutput_process; }

int HiddenMarkovTree::get_nb_doutput_process() const
{ return _nb_doutput_process; }

/*****************************************************************
 *
 * Return the number of values of a given variable
 * with finite values for a HiddenMarkovTree
 *
 **/

int HiddenMarkovTree::get_nb_values(int variable) const
{
  int nb_values;

  assert((variable >= 0) && (variable < _nb_ioutput_process));
  // the number of values for a float random variable does not make much sense

  if (npprocess[variable+1] != NULL)
     nb_values= npprocess[variable+1]->nb_value;
  else
     nb_values= piprocess[variable+1]->nb_value;

  return nb_values;
}

/*****************************************************************
 *
 * Return "true" if and only if process ivariable is parametric
 *
 **/

bool HiddenMarkovTree::is_parametric(int ivariable) const
{
  assert((ivariable >= 0) && (ivariable < _nb_ioutput_process));
  // all float random processes are parametric

  return (piprocess[ivariable+1] != NULL);
}


/*****************************************************************
 *
 *  Access to the non-parametric observation processes
 *  for HiddenMarkovTree class (copy of the process)
 *
 **/

CategoricalTreeProcess* HiddenMarkovTree::get_categorical_process(int variable) const
{
   CategoricalTreeProcess *res= NULL;

   if (npprocess != NULL)
   {
      assert((variable > 0) && (variable <= _nb_ioutput_process));
      // npprocess[0] points to the (hidden) state process
      if (npprocess[variable] != NULL)
         res= new CategoricalTreeProcess(*(npprocess[variable]));
   }
   return res;
}

/*****************************************************************
 *
 *  Access to the discrete parametric observation processes
 *  for HiddenMarkovTree class (copy of the process)
 *
 **/

DiscreteParametricProcess* HiddenMarkovTree::get_iparametric_process(int variable) const
{
   DiscreteParametricProcess *res= NULL;

   if (piprocess != NULL)
   {
      assert((variable > 0) && (variable <= _nb_ioutput_process));
      // piprocess[0] == NULL
      if (piprocess[variable] != NULL)
         res= new DiscreteParametricProcess(*(piprocess[variable]));
   }
   return res;
}

/*****************************************************************
 *
 *  Access to the continuous parametric observation processes
 *  for HiddenMarkovTree class (copy of the process)
 *
 **/

DiscreteParametricProcess* HiddenMarkovTree::get_dparametric_process(int variable) const
{
   DiscreteParametricProcess *res= NULL;

   if (pdprocess != NULL)
   {
      assert((variable > 0) && (variable <= _nb_doutput_process));
      // pdprocess[0] == NULL
      if (pdprocess[variable] != NULL)
         res= new DiscreteParametricProcess(*(pdprocess[variable]));
   }
   return res;
}

/*****************************************************************
 *
 *  Access to the non-parametric observation processes
 *  for HiddenMarkovTree class (actual pointer)
 *
 **/

CategoricalTreeProcess** HiddenMarkovTree::get_categorical_process() const
{ return npprocess; }

/*****************************************************************
 *
 *  Access to the discrete parametric observation processes
 *  for HiddenMarkovTree class (actual pointer)
 *
 **/

DiscreteParametricProcess** HiddenMarkovTree::get_iparametric_process() const
{ return piprocess; }

/*****************************************************************
 *
 *  Access to the continuous parametric observation processes
 *  for HiddenMarkovTree class (actual pointer)
 *
 **/

DiscreteParametricProcess** HiddenMarkovTree::get_dparametric_process() const
{ return pdprocess; }


/*****************************************************************
 *
 *  Computation of the characteristic distributions
 *  for HiddenMarkovTree class
 *  using the size of the trees (or other quantities for which
 *  the characteristic distributions are invariant - if any),
 *  a flag on the computation of counting distributions
 *  and the considered variable
 *
 **/

void HiddenMarkovTree::characteristic_computation(int size,
                                                    bool counting_flag,
                                                    int variable)
{
   bool computation[NB_OUTPUT_PROCESS+1];
   register int i, j, k;
   double sum, *memory= NULL;
   DiscreteParametric dsize(UNIFORM, size, size, D_DEFAULT, D_DEFAULT);

   if (nb_component > 0)
   {

     // computation of the intensity and interval distributions
     // at state level

     if (((variable == I_DEFAULT) || (variable == 0)) &&
         ((!(npprocess[0]->size != NULL)) || (dsize != *(npprocess[0]->size))))
     {
        computation[0]= true;
        npprocess[0]->create_characteristic((Distribution)dsize, false, counting_flag);

        if (memory == NULL)
        {
           switch (type)
           {
              case 'o' :
              {
                 memory= memory_computation();
                 break;
              }

              case 'e' :
              {
                 memory= new double[nb_state];
                 for(j= 0; j < nb_state; j++)
                    memory[j]= initial[j];
                 break;
              }
           }
        }

        for(i= 0; i < nb_state; i++)
           state_leave_probability(memory, i);

       // computation of the stationnary distribution in the case of
       // an equilibrium state process

       if (type == 'e')
       {
          sum= 0.;
          for(i= 0; i < nb_state; i++)
          {
             initial[i]= 1. / nb_state;
             sum += initial[i];
          }
          // the true formula seems to imply the computation of the
          // transposed transition matrix eigenvectors

          for(i= 0; i < nb_state; i++)
             initial[i] /= sum;
       }

       // index_state_distribution();

       for(i= 0; i < nb_state; i++)
       {
          state_no_occurrence_probability(i);
          state_first_occurrence_root_distribution(i);
          state_first_occurrence_leaves_distribution(i);

          /* if ((transition[i][i] > 0.) && (transition[i][i] < 1.))
          {
             if ((npprocess[0]->sojourn_size[i] != NULL) &&
                 (npprocess[0]->sojourn_size[i]->parameter != 1.-transition[i][i]))
             {
                delete npprocess[0]->sojourn_size[i];
                npprocess[0]->sojourn_size[i]= NULL;
             }

             if (npprocess[0]->sojourn_size[i] == NULL)
                npprocess[0]->sojourn_size[i]= new DiscreteParametric(NEGATIVE_BINOMIAL , 1 , I_DEFAULT , 1. ,
                                                                      1.-transition[i][i] , OCCUPANCY_THRESHOLD);
             }
             npprocess[0]->sojourn_size[i]->ident= CATEGORICAL;
          }*/

          // valid for semi-markov processes
          state_sojourn_size_distribution(memory, i);
       }
   }
   else // nb_component < 1
      computation[0]= false;

     // computation of the intensity and interval distributions
     // at observation level

   for(i= 1; i <= _nb_ioutput_process; i++)
   {
      if ((npprocess[i] != NULL) && ((variable == I_DEFAULT) || (i == variable)) &&
          ((npprocess[i]->size == NULL) || ((npprocess[i]->size != NULL) &&
            (dsize != *(npprocess[i]->size)))))
      {
         computation[i]= true;
         npprocess[i]->create_characteristic(dsize, true, counting_flag);

         // index_output_distribution(i);

         for(j= 0; j < npprocess[i]->nb_value; j++)
         {
            output_no_occurrence_probability(i, j);
            output_first_occurrence_root_distribution(i, j);
            output_first_occurrence_leaves_distribution(i, j);

            output_leave_probability(memory, i, j);

            for(k= 0; k < nb_state; k++)
            {
               if ((npprocess[i]->observation[k]->mass[j] > 0.) &&
                  ((state_type[k] != 'a') || (npprocess[i]->observation[k]->mass[j] < 1.)))
                break;
            }

            if (k < nb_state)
               output_sojourn_size_distribution(memory, i, j);
            else
            {
               npprocess[i]->absorption[j]= 1.;
               delete npprocess[i]->sojourn_size[j];
               npprocess[i]->sojourn_size[j]= NULL;
            }
         }
      }
      else
        computation[i]= false;
   }

      delete [] memory;

      if (counting_flag)
      {
        // computation of the counting distributions at state level

        if (computation[0])
           for (i = 0;i < nb_state;i++)
           {
              state_nb_pattern_mixture(i, 'r');
              state_nb_pattern_mixture(i, 'o');
           }

        // computation of the counting distributions at observation level

        for(i= 1; i <= _nb_ioutput_process; i++)
           if (computation[i])
              for (j= 0; j < npprocess[i]->nb_value; j++)
              {
                output_nb_zones_mixture(i, j);
                output_nb_occurrences_mixture(i, j);
              }
      }
   }
}

/*****************************************************************
 *
 * Compute the state profiles for a given tree,
 * including the smoothed probabilities,
 * the upward-downward and the generalized Viterbi algorithms
 *
 **/

bool HiddenMarkovTree::state_profile(StatError& error,
                                     const HiddenMarkovTreeData& trees,
                                     int identifier,
                                     HiddenMarkovTreeData*& smoothed_state_tree,
                                     HiddenMarkovTreeData*& nstate_trees,
                                     HiddenMarkovTreeData*& viterbi_upward_downward,
                                     HiddenMarkovTreeData*& generalized_restoration,
                                     std::vector<ostringstream*>& messages,
                                     int state_tree,
                                     unsigned int nb_state_trees,
                                     int entropy_algo,
                                     int root) const
{ return false; }

/*****************************************************************
 *
 *  Write Gnuplot output of state and Viterbi profiles
 *
 **/

bool HiddenMarkovTree::tree_state_profile_plot_write(StatError &error,
                                                     const char *prefix,
                                                     const HiddenMarkovTreeData& trees,
                                                     int identifier,
                                                     int vertex,
                                                     const char *title,
                                                     int entropy_algo) const
{ return false; }

/*****************************************************************
 *
 *  Write Gnuplot output of state and entropy profiles
 *
 **/

bool HiddenMarkovTree::state_profile_plot_write(StatError &error,
                                                const char *prefix,
                                                int identifier,
                                                int vertex,
                                                const char *title,
                                                int entropy_algo) const
{
   bool status;

   error.init();

   if (markov_data == NULL)
   {
      status= false;
      error.update(STAT_error[STATR_NO_DATA]);
   }
   else
      status= tree_state_profile_plot_write(error , prefix , *markov_data ,
                                            identifier, vertex, title ,
                                            entropy_algo);
   return status;
}

/*****************************************************************
 *
 *  Copy operator of HiddenMarkovTree class
 *  with a flag on the HiddenMarkovTreeData copy
 *  and one on the characteristic distribution copy
 *
 **/

void HiddenMarkovTree::copy(const HiddenMarkovTree& markov, bool data_flag,
                            bool characteristic_flag)
{
   register int i;

   // Chain::copy(markov) must be used before (or Chain::Chain)
   // as well as remove()

   if ((data_flag) && (markov.markov_data != NULL))
      markov_data= new HiddenMarkovTreeData(*(markov.markov_data), false);
   else
      markov_data= NULL;

   _ch_order= markov._ch_order;

   if (markov.self_row != NULL)
   {
      self_row= new int[nb_state];
      for(i = 0;i < nb_state;i++)
         self_row[i] = markov.self_row[i];
   }

   _nb_ioutput_process= markov._nb_ioutput_process;
   _nb_doutput_process= markov._nb_doutput_process;

   if ((markov.npprocess != NULL) && (markov.piprocess != NULL))
   {
      npprocess= new CategoricalTreeProcess*[_nb_ioutput_process+1];
      piprocess= new DiscreteParametricProcess*[_nb_ioutput_process+1];
      piprocess[0]= NULL;

      if (markov.npprocess[0] != NULL)
         npprocess[0]= new CategoricalTreeProcess(*(markov.npprocess[0]),
                                                  'c',
                                                   characteristic_flag);
      for(i= 0; i < _nb_ioutput_process; i++)
      {
         if (markov.npprocess[i+1] != NULL)
         {
            npprocess[i+1]= new CategoricalTreeProcess(*(markov.npprocess[i+1]),
                                                       'c',
                                                        characteristic_flag);
            piprocess[i+1]= NULL;
         }
         else
         {
            npprocess[i+1]= NULL;
            piprocess[i+1]= new DiscreteParametricProcess(*(markov.piprocess[i+1]));
         }
      }
   }
   else
   {
      npprocess= NULL;
      piprocess= NULL;
   }

   if (markov.pdprocess != NULL)
   {
      pdprocess= new DiscreteParametricProcess*[_nb_doutput_process+1];
      pdprocess[0]= NULL;
      for(i= 0; i < _nb_doutput_process; i++)
         pdprocess[i+1]= new DiscreteParametricProcess(*(markov.pdprocess[i+1]));
   }
   else
      pdprocess= NULL;
}

/*****************************************************************
 *
 *  Deallocation of the pointers for HiddenMarkovTree class
 *
 **/

void HiddenMarkovTree::remove()
{
   register int i;

   if (markov_data != NULL)
      delete markov_data;
   markov_data= NULL;

   if (self_row != NULL)
      delete [] self_row;
   self_row= NULL;

   if (npprocess != NULL)
   {
      for(i= 0; i <= _nb_ioutput_process; i++)
      {
         if (npprocess[i] != NULL)
            delete npprocess[i];
         npprocess[i]= NULL;

         if (piprocess[i] != NULL)
            delete piprocess[i];
         piprocess[i]= NULL;
      }
      delete [] npprocess;
      delete [] piprocess;

      npprocess= NULL;
      piprocess= NULL;
   }

   if (pdprocess != NULL)
   {
      for(i= 0; i <= _nb_doutput_process; i++)
      {
         if (pdprocess[i] != NULL)
            delete pdprocess[i];
         pdprocess[i]= NULL;
      }
      delete [] pdprocess;
      pdprocess= NULL;
   }
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovTree and the corresponding data structure
 *  using an output stream, a HiddenMarkovTreeData object,
 *  a flag on the level of detail, a flag on the file use,
 *  a Test object and a flag on the (children) order printing
 *
 **/

ostream& HiddenMarkovTree::ascii_write(ostream& os,
                                       const HiddenMarkovTreeData * otrees,
                                       bool exhaustive,
                                       bool file_flag,
                                       const Test* test,
                                       bool ch_order_flag) const
{
   register int i, j, k;
   int variable, cumul_size, nb_variable, buff;
   int width[2];
   int nb_output_process = _nb_ioutput_process+_nb_doutput_process;
   bool **logic_transition = NULL;
   double **distance = NULL;
   FrequencyDistribution **observation = NULL, *marginal = NULL;
   TreeCharacteristics *characteristics = NULL;

   // printing of the Markov tree parameters

   Chain::ascii_print(os, file_flag);

   if ((otrees != NULL) && (otrees->_type[0] == STATE))
   {
      variable= 0;
      characteristics= otrees->characteristics[variable];
   }

   npprocess[0]->ascii_print(os, 0, NULL, characteristics, exhaustive, file_flag);

   // printing of the (characteristic ?) distributions of each
   // observed process

   if (nb_output_process > 0)
   {
      os << "\n" << nb_output_process << " "
         << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;

      for(i = 1; i <= _nb_ioutput_process; i++)
      {
         os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
         os << " " << i;

         if (npprocess[i] != NULL)
            os << " : " << STAT_word[STATW_CATEGORICAL];
         else
            os << " : " << STAT_word[STATW_DISCRETE_PARAMETRIC];

         os << endl;

         if (otrees != NULL)
         {
            nb_variable = otrees->get_nb_int() + otrees->get_nb_float();
            switch (otrees->_type[0])
            {
               case STATE:
               // first variable of otrees is not an output process
               // but the state process
                  if ((nb_output_process > 1) && (nb_variable > nb_output_process))
                     variable= i;
                  else
                     variable= i-1;
                  break;
               /* case INT_VALUE:
                  variable= i - 1;
                  break; */
               // _type is generally not defined for the moment
               default :
                  variable= i - 1;
            }

           if (otrees->observation_distribution != NULL)
              observation = otrees->observation_distribution[variable];
           else
              observation = NULL;

           if (otrees->characteristics[variable] != NULL)
           {
              characteristics = otrees->characteristics[variable];
              marginal = otrees->characteristics[variable]->marginal_distribution;
           }
           else
           {
              characteristics = NULL;
              marginal = NULL;
           }

         }

         logic_transition = logic_transition_computation();

         distance = new double*[nb_state];
         for(j = 0; j < nb_state; j++)
            distance[j] = new double[nb_state];

         if (npprocess[i] != NULL)
         {
            npprocess[i]->ascii_print(os, i, observation, characteristics,
                                      exhaustive, file_flag);

            for(j = 0; j < nb_state; j++)
            {
               distance[j][j] = 0.;

               for(k = j + 1;k < nb_state;k++)
               {
                  if ((logic_transition[j][k]) || (logic_transition[k][j]))
                     distance[j][k] =
                         npprocess[i]->observation[j]->overlap_distance_computation(*(npprocess[i]->observation[k]));
                  else
                     distance[j][k] = 1.;

                  distance[k][j] = distance[j][k];
               }
            }
         }
         else
         {
            // does not seem to work when observation == NULL and exhaustiv == true
            piprocess[i]->ascii_print(os, observation, marginal,
                                      exhaustive, file_flag);

            for(j = 0; j < nb_state; j++)
            {
               distance[j][j] = 0.;

               for(k = j + 1; k < nb_state; k++)
               {
                  if ((logic_transition[j][k]) || (logic_transition[k][j]))
                     distance[j][k] =
                         piprocess[i]->observation[j]->sup_norm_distance_computation(*(piprocess[i]->observation[k]));
                  else
                     distance[j][k] = 1.;

                  distance[k][j] = distance[j][k];
               }
            }
         }

         width[0] = column_width(nb_state , distance[0]);
         for(j = 1; j < nb_state; j++)
         {
            buff = column_width(nb_state , distance[j]);
            if (buff > width[0])
                      width[0] = buff;
         }
         width[0] += ASCII_SPACE;

         os.setf(ios::left , ios::adjustfield);

         os << "\n";
         if (file_flag)
            os << "# ";

         os << STAT_TREES_label[TREESTATL_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

         for(j = 0; j < nb_state; j++)
         {
            if (file_flag)
               os << "# ";

            for(k = 0; k < nb_state; k++)
            {
               if ((k != j) && (logic_transition[j][k]))
                  os << setw(width[0]) << distance[j][k];
               else
                  os << setw(width[0]) << "_";
            }
            os << endl;
         }

         for(j = 0; j < nb_state; j++)
         {
            delete [] logic_transition[j];
            logic_transition[j] = NULL;
         }

         delete [] logic_transition;
         logic_transition = NULL;

         for(j = 0; j < nb_state; j++)
         {
            delete [] distance[j];
            distance[j] = NULL;
         }
         delete [] distance;
      } // end for i <= _nb_ioutput_process

      for(i = 0; i < _nb_doutput_process; i++)
      {
         os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
         os << " " << i;

         os << " : " << STAT_word[STATW_DISCRETE_PARAMETRIC] << endl;

         if (otrees != NULL)
         {
           if (otrees->observation_distribution != NULL)
              observation= otrees->observation_distribution[i];
           // This seems weird : how is observation supposed to be defined
           // for parametric processes ? What happens if otrees->observation_distribution == NULL ?
         }
         if (otrees->characteristics[variable] != NULL)
            pdprocess[i]->ascii_print(os, observation, otrees->characteristics[variable]->marginal_distribution,
                                      exhaustive, file_flag);
         else
            pdprocess[i]->ascii_print(os, observation, NULL, exhaustive, file_flag);

         for(j = 0; j < nb_state; j++)
         {
            distance[j][j] = 0.;

            for(k = j + 1; k < nb_state; k++)
            {
               if ((logic_transition[j][k]) || (logic_transition[k][j]))
                  distance[j][k] =
                      pdprocess[i]->observation[j]->sup_norm_distance_computation(*(pdprocess[i]->observation[k]));
               else
                  distance[j][k] = 1.;

               distance[k][j] = distance[j][k];
            }
         }
         if (pdprocess[i] != NULL)
         {
            width[0] = column_width(nb_state , distance[0]);
            for(j = 1; j < nb_state; j++)
            {
               buff = column_width(nb_state , distance[j]);
               if (buff > width[0])
                         width[0] = buff;
            }
            width[0] += ASCII_SPACE;

            os.setf(ios::left , ios::adjustfield);

            os << "\n";
            if (file_flag)
               os << "# ";

            os << STAT_TREES_label[TREESTATL_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

            for(j = 0; j < nb_state; j++)
            {
               if (file_flag)
                  os << "# ";

               for(k = 0; k < nb_state; k++)
               {
                  if ((k != j) && (logic_transition[j][k]))
                     os << setw(width[0]) << distance[j][k];
                  else
                     os << setw(width[0]) << "_";
               }
               os << endl;
            }

            for(j = 0; j < nb_state; j++)
            {
               delete [] logic_transition[j];
               logic_transition[j] = NULL;
            }

            delete [] logic_transition;
            logic_transition = NULL;

            for(j = 0; j < nb_state; j++)
            {
               delete [] distance[j];
               distance[j] = NULL;
            }
            delete [] distance;
         }
      }
   }

   if (otrees != NULL)
   {
      int nb_parameter= nb_parameter_computation(MIN_PROBABILITY);
      double information, likelihood, hidden_likelihood;

      hidden_likelihood= otrees->hidden_likelihood;
      likelihood= otrees->likelihood;

      // print the quantities for which the characteristic distributions
      // are invariant - if any

      os << "\n";
      if (file_flag)
         os << "# ";
      os << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
      otrees->hsize->ascii_characteristic_print(os, false, file_flag);

      if (exhaustive)
      {
         os << "\n";
         if (file_flag)
            os << "# ";
         os << "   | " << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
         otrees->hsize->ascii_print(os, file_flag);
      }

      os << "\n";
      if (file_flag)
        os << "# ";
      os << STAT_TREES_label[TREESTATL_CUMULATIVE_SIZE] << ": " << otrees->cumul_size_computation() << endl;

      // print the tree information quantity in the iid case

      information= otrees->iid_information_computation();

      os << "\n";
      if (file_flag)
         os << "# ";
      cumul_size= otrees->cumul_size_computation();
      os << STAT_TREES_label[TREESTATL_TREES_IID_INFORMATION] << ": " << information << " ("
         << information / cumul_size << ")" << endl;

      // print the likelihood

      if (likelihood != D_INF)
      {
         os << "\n";
         if (file_flag)
            os << "# ";
         os << STAT_TREES_label[TREESTATL_STATE_TREES_LIKELIHOOD] << ": " << hidden_likelihood << "   ("
            << STAT_label[STATL_NORMALIZED] << ": " << hidden_likelihood / cumul_size << ")" << endl;
      }

      if (otrees->sample_entropy != D_DEFAULT) {
        os << "\n";
        if (file_flag) {
          os << "# ";
        }
        os << STAT_TREES_label[TREESTATL_STATE_TREE_ENTROPY] << ": " << otrees->sample_entropy << "   ("
           << STAT_label[STATL_NORMALIZED] << ": " << otrees->sample_entropy / cumul_size << ")" << endl;
      }

      if (likelihood != D_INF)
      {
         os << "\n";
         if (file_flag)
           os << "# ";
         os << STAT_TREES_label[TREESTATL_OBSERVED_TREES_LIKELIHOOD] << ": " << likelihood << "   ("
            << STAT_label[STATL_NORMALIZED] << ": " << likelihood / cumul_size << ")" << endl;
      }

      // print AIC, AICc and BIC.
      if (likelihood != D_INF)
      {
         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AIC] << "): "
            << 2 * (likelihood - nb_parameter) << endl;

         if (nb_parameter < cumul_size-1)
         {
            os << "\n";
            if (file_flag)
               os << "# ";
            os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
               << 2*(likelihood-(double)(nb_parameter*cumul_size) /
                 (double)(cumul_size-nb_parameter-1)) << endl;
         }

         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
            << 2 * likelihood - nb_parameter * log((double)cumul_size) << endl;

         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BICc] << "): "
            << 2 * likelihood-penalty_computation(MIN_PROBABILITY) << endl;
      }

      // print BIC-ICL
      if (otrees->hidden_likelihood != D_INF)
      {
         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICL] << "): "
            // << 2*hidden_likelihood-nb_parameter*log((double)cumul_size) << endl;
            << ((otrees->sample_entropy != D_INF) ? (2 * (likelihood - otrees->sample_entropy) - nb_parameter
                  * log((double)cumul_size)) : D_INF) << endl;
         // otrees->likelihood == completed likelihood

         os << "\n";
         if (file_flag)
            os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[ICLc] << "): "
            // << 2*hidden_likelihood-penalty_computation(MIN_PROBABILITY) << endl;
            << ((otrees->sample_entropy != D_INF) ? (2 * (hidden_likelihood - otrees->sample_entropy)
               - penalty_computation(MIN_PROBABILITY)) : D_INF) << endl;
      }


      /* if ((likelihood != D_INF) && (nb_component == 1))
      {
         if (nb_parameter < cumul_size-1)
         {
            os << "\n";
            if (file_flag)
              os << "# ";
            os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
               << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << "): "
               << 2 * (likelihood - (double)(nb_parameter * cumul_size) /
                 (double)(cumul_size-nb_parameter-1)) << endl;
         }

         os << "\n";
         if (file_flag)
           os << "# ";
         os << nb_parameter << " " << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS]
            << "   2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << "): "
            << 2 * likelihood - nb_parameter * log((double)cumul_size) << endl;
      } */
   }

   if (test)
   {
      os << "\n";
      test->ascii_print(os, file_flag);
   }

   return os;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovTree and the corresponding data structure
 *  in a spreadsheet fashion
 *  using an output stream, a HiddenMarkovTreeData object
 *  and a Test object
 *
 **/

ostream& HiddenMarkovTree::spreadsheet_write(ostream& os, const HiddenMarkovTreeData * otrees,
                                             const Test * test) const
{
   register int i, j, k;
   int variable = 0, cumul_size,
       nb_output_process = _nb_ioutput_process+_nb_doutput_process;
   bool **logic_transition = NULL;
   double **distance = NULL;
   FrequencyDistribution **observation = NULL;
   TreeCharacteristics *characteristics = NULL;

   switch (type)
   {
      case 'o' :
         os << STAT_TREES_word[TREESTATW_HIDDEN_MARKOV_TREE] << endl;
         break;
      case 'e' :
         os << STAT_TREES_word[TREESTATW_EQUILIBRIUM_HIDDEN_MARKOV_TREE] << endl;
         break;
   }

   // printing of the Markov tree parameters

   spreadsheet_print(os);

   if ((otrees != NULL) && (otrees->_type[0] == STATE))
   {
      variable= 0;
      characteristics= otrees->characteristics[variable];
   }

   npprocess[0]->spreadsheet_print(os, 0, 0, characteristics);

   // printing of the (characteristic ?) distributions of each
   // observed process

   if (nb_output_process > 0)
   {
      os << "\n" << nb_output_process << "\t"
         << STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES] << endl;

      for(i= 1; i <= _nb_ioutput_process; i++)
      {
         os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
         os << "\t" << i;

         if (npprocess[i] != NULL)
            os << "\t" << STAT_word[STATW_CATEGORICAL];
         else
            os << "\t" << STAT_word[STATW_DISCRETE_PARAMETRIC];

         os << endl;

         if (otrees != NULL)
         {
            switch (otrees->_type[0])
            {
               case INT_VALUE :
                  variable= i - 1;
                  break;
               case STATE :
                  variable= i;
                  break;
            }

            if (otrees->observation_distribution != NULL)
               observation = otrees->observation_distribution[variable];
            if (otrees->characteristics[variable] != NULL)
               characteristics = otrees->characteristics[variable];
         }

         logic_transition = logic_transition_computation();

         distance = new double*[nb_state];
         for(j = 0; j < nb_state; j++)
            distance[j] = new double[nb_state];

         if (npprocess[i] != NULL)
         {
            npprocess[i]->spreadsheet_print(os, i, observation, characteristics);
            for(j = 0; j < nb_state; j++)
            {
               distance[j][j] = 0.;

               for(k = j + 1;k < nb_state;k++)
               {
                  if ((logic_transition[j][k]) || (logic_transition[k][j]))
                     distance[j][k] =
                         npprocess[i]->observation[j]->overlap_distance_computation(*(npprocess[i]->observation[k]));
                  else
                     distance[j][k] = 1.;

                  distance[k][j] = distance[j][k];
               }
            }
         }
         else
         {
            piprocess[i]->spreadsheet_print(os, observation);

            for(j = 0; j < nb_state; j++)
            {
               distance[j][j] = 0.;

               for(k = j + 1; k < nb_state; k++)
               {
                  if ((logic_transition[j][k]) || (logic_transition[k][j]))
                     distance[j][k] =
                         piprocess[i]->observation[j]->sup_norm_distance_computation(*(piprocess[i]->observation[k]));
                  else
                     distance[j][k] = 1.;

                  distance[k][j] = distance[j][k];
               }
            }
         }
         os << "\n";
         os << STAT_TREES_label[TREESTATL_OBSERVATION_DISTRIBUTION_DISTANCE] << endl;

         for(j = 0; j < nb_state; j++)
         {
            for(k = 0; k < nb_state; k++)
            {
               if ((k != j) && (logic_transition[j][k]))
                  os << distance[j][k];
               else
                  os << "\t";
            }
            os << endl;
         }

         for(j = 0; j < nb_state; j++)
         {
            delete [] logic_transition[j];
            logic_transition[j] = NULL;
         }

         delete [] logic_transition;
         logic_transition = NULL;

         for(j = 0; j < nb_state; j++)
         {
            delete [] distance[j];
            distance[j] = NULL;
         }
         delete [] distance;


      }

      for(i= 0; i < _nb_ioutput_process; i++)
      {
         os << "\n" << STAT_word[STATW_OUTPUT_PROCESS];
         os << "\t" << i;

         os << "\t" << STAT_word[STATW_DISCRETE_PARAMETRIC] << endl;

         if (otrees != NULL)
         {
            if (otrees->observation_distribution != NULL)
               observation= otrees->observation_distribution[i];
            // warning : akward !
         }
         piprocess[i]->spreadsheet_print(os, observation);
      }
   }

   if (otrees != NULL)
   {
      int nb_parameter= nb_parameter_computation(MIN_PROBABILITY);
      double information, likelihood;


      // printing of the quantities for which the characteristic distributions
      // are invariant - if any

      os << "\n" << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\t";
      otrees->hsize->spreadsheet_characteristic_print(os);

      os << "\n\t" << STAT_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
      otrees->hsize->spreadsheet_print(os);

      cumul_size= otrees->cumul_size_computation();
      os << "\n" << STAT_TREES_label[TREESTATL_CUMULATIVE_SIZE] << "\t" << cumul_size << endl;

      // printing of the tree information quantity in the iid case

      information= otrees->iid_information_computation();

      os << "\n" << STAT_TREES_label[TREESTATL_TREES_IID_INFORMATION] << "\t" << information << "\t"
         << information / cumul_size << endl;

      // printing of the likelihood

      if (otrees->likelihood != D_INF)
      {
         os << "\n" << STAT_TREES_label[TREESTATL_STATE_TREES_LIKELIHOOD] << "\t" << otrees->likelihood << "\t"
            << STAT_label[STATL_NORMALIZED] << "\t" << otrees->likelihood / cumul_size << endl;
      }

      likelihood= otrees->hidden_likelihood;

      if (likelihood != D_INF)
      {
         os << "\n" << STAT_TREES_label[TREESTATL_OBSERVED_TREES_LIKELIHOOD] << "\t" << likelihood << "\t"
            << STAT_label[STATL_NORMALIZED] << "\t" << likelihood / cumul_size << endl;
      }

      if ((likelihood != D_INF) && (nb_component == 1))
      {
         if (nb_parameter < cumul_size-1)
         {
            os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
               << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[AICc] << ")\t"
               << 2 * (likelihood - (double)(nb_parameter * cumul_size) /
                  (double)(cumul_size-nb_parameter-1)) << endl;
         }

         os << "\n" << nb_parameter << "\t" << STAT_label[nb_parameter == 1 ? STATL_FREE_PARAMETER : STATL_FREE_PARAMETERS] << "\t"
            << "2 * " << STAT_label[STATL_PENALIZED_LIKELIHOOD] << " ("  << STAT_criterion_word[BIC] << ")\t"
            << 2 * likelihood - nb_parameter * log((double)cumul_size) << endl;
      }
   }

   if (test)
   {
      os << "\n";
      test->spreadsheet_print(os);
   }

   return os;
}

/*****************************************************************
 *
 *  Gnuplot output for HiddenMarkovTree class
 *  using a prefix for the files, the title of figures
 *  and the observed trees
 *
 **/

bool HiddenMarkovTree::plot_write(const char * prefix, const char * title,
                                  const HiddenMarkovTreeData * otrees) const
{
   bool status;
   register int i;
   int variable= 0; //, cumul_size, nb_output_process= _nb_ioutput_process+_nb_doutput_process;
   FrequencyDistribution *hsize= NULL, **observation= NULL;
   TreeCharacteristics *characteristics= NULL;

   if ((otrees != NULL) && (otrees->_type[0] == STATE))
   {
      variable= 0;
      characteristics= otrees->state_characteristics;
      // characteristics= otrees->characteristics[variable];
      hsize= otrees->hsize;
   }

   status= npprocess[0]->plot_print(prefix, title, 0, NULL, characteristics, hsize);

   // print the (characteristic ?) distributions of each
   // observed process

   if (status)
   {
      if (otrees != NULL)
         hsize= otrees->hsize;

      for(i= 1; i <= _nb_ioutput_process; i++)
      {
         if (otrees != NULL)
         {
            switch (otrees->_type[0])
            {
               case INT_VALUE :
                  variable= i - 1;
                  break;
               case STATE :
                  variable= i;
                  break;
            }

            if (otrees->observation_distribution != NULL)
               observation= otrees->observation_distribution[variable];

            if (otrees->characteristics[variable] != NULL)
               characteristics= otrees->characteristics[variable];
         }

         if (npprocess[i] != NULL)
            npprocess[i]->plot_print(prefix, title, i, observation,
                                     characteristics, hsize);
         else
            piprocess[i]->plot_print(prefix, title, i, observation);
      }

      for(i= 0; i < _nb_doutput_process; i++)
      {
         if (otrees->observation_distribution != NULL)
            observation= otrees->observation_distribution[i];

         piprocess[i]->plot_print(prefix, title, i, observation);
      }
   }
   return status;
}


/*****************************************************************
 *
 *  Computation of the indices of the transition probability matrix
 *  corresponding to the probabilities of staying in each state
 *  for HiddenMarkovTree class
 *
 **/

void HiddenMarkovTree::self_row_computation()
{
   register int i, j;
   int state_index[ORDER], *pself_row;

   for (i= 0; i < _ch_order; i++)
      state_index[i]= 0;

   pself_row= self_row;

   for(i= 0; i < nb_row; i++)
   {
      for(j= 1; j < _ch_order; j++)
         if (state_index[j] != state_index[j-1])
            break;

      if (j == _ch_order)
         *pself_row++ = i;

      for (j = 0; j < _ch_order; j++)
      {
         if (state_index[j] < nb_state-1)
         {
            state_index[j]++;
            break;
         }
         else
           state_index[j]= 0;
      }
   }
}

/*****************************************************************
 *
 *  Extraction of the HiddenMarkovTree classes
 *  from state accessibility
 *
 **/

void HiddenMarkovTree::component_computation()

{
   bool **logic_transition;
   register int i, j, k;
   int power= nb_row / nb_state;
   double sum, **ptransition;

   // computation of the matrix of the allowed state transitions

   logic_transition= new bool*[nb_state];

   for (i= 0; i < nb_state; i++)
   {
      logic_transition[i] = new bool[nb_state];

      for (j= 0; j < nb_state; j++)
      {
         if (j == i)
            logic_transition[i][j]= false;
         else
         {
            ptransition= transition + i * power;
            sum= 0.;
            for(k = 0;k < power;k++)
            {
               sum+= *(*ptransition+j);
               ptransition++;
            }
            logic_transition[i][j]= (sum == 0. ? false : true);
         }
      }
   }

   Chain::component_computation(logic_transition);

   for (i= 0; i < nb_state; i++)
     delete [] logic_transition[i];

   delete [] logic_transition;
}

/*****************************************************************
 *
 *  Computation of the distribution of each memory for a HiddenMarkovTree
 *  given the tree size distribution
 *
 **/

double* HiddenMarkovTree::memory_computation() const
{
   register int i , j , k;
   int power[ORDER], state_index[ORDER];
   double *memory, *state_tree, *pstate_tree, *states, *pstates, **ptransition;


   memory= new double[nb_row];
   for(i= 0; i < nb_row; i++)
      memory[i] = 0.;

   i= 1;
   for(j= 0; j < _ch_order; j++)
   {
      power[j]= i;
      i*= nb_state;
   }

   // initialization of the tree state probabilities

   state_tree= new double[nb_row];
   pstate_tree= new double[nb_row];

   pstates= pstate_tree;
   i= 0;

   for(j = 0;j < nb_row;j++)
   {
      if (j == self_row[i])
         *pstates++ = initial[i++];
      else
        *pstates++ = 0.;
   }

   // computation of the probability of each memory, depending on the index
   for(i= 1; i < npprocess[0]->size->nb_value-2; i++)
   {
     // computation of the state tree probabilities having a size
     // (or tree depth ?) = _ch_order

     for(j= 0; j < _ch_order; j++)
        state_index[j]= 0;

     states= state_tree;

     for(j = 0;j < nb_row;j++)
     {
        ptransition= transition;
        pstates= pstate_tree;
        for (k= 0; k < _ch_order-1; k++)
        {
           ptransition+= state_index[k] * power[k+1];
           pstates+= state_index[k] * power[k+1];
        }

        *states= 0.;
        for(k= 0; k < nb_state; k++)
        {
           *states+= *(*ptransition + state_index[_ch_order-1]) * *pstates++;
           ptransition++;
        }
        states++;

        // update of the state indices

        for(k= 0; k < _ch_order; k++)
        {
           if (state_index[k] < nb_state - 1)
           {
              state_index[k]++;
              break;
           }
           else
             state_index[k]= 0;
        }
     }

     // update of the tree state probabilities and
     // of the cumulative memory probabilities

     states= state_tree;
     pstates= pstate_tree;

     for(j= 0; j < nb_row; j++)
     {
        memory[j]+= *states * (1. - npprocess[0]->size->cumul[i]);
        *pstates++ = *states++;
     }
   }

   delete [] state_tree;
   delete [] pstate_tree;

   return memory;
}

/*****************************************************************
 *
 *  Initialization of the HiddenMarkovTree parameters
 *  using a flag on the HiddenMarkovTree nature
 *  and a fixed probability of staying in any state
 *
 **/

void HiddenMarkovTree::init(bool left_right, double self_transition)
{
   register int i, j;
   int power= nb_row / nb_state, state, state_index[ORDER];

   accessibility= new bool*[nb_state];
   for(i= 0; i < nb_state; i++)
      accessibility[i]= new bool[nb_state];

   state_type= new char[nb_state];

   switch (left_right)
   {
      // case where each transition is allowed

      case false :
      {
         nb_component= 1;
         component_nb_state= new int[nb_component];
         component_nb_state[0]= nb_state;
         component= new int*[nb_component];
         component[0]= new int[component_nb_state[0]];

         for(i= 0; i < nb_state; i++)
         {
            for(j= 0; j < nb_state; j++)
               accessibility[i][j] = true;

            component[0][i]= i;
            state_type[i]= 'r';
         }

         for(i= 0; i < nb_state; i++)
            initial[i] = 1. / (double)nb_state;

         for(i= 0; i < nb_row; i++)
         {
            state= i / power;
            for(j= 0; j < state; j++)
               transition[i][j]= (1. - self_transition) / (nb_state - 1);

            transition[i][state]= self_transition;
            for(j= state+1; j < nb_state; j++)
               transition[i][j] = (1. - self_transition) / (nb_state - 1);
         }
         break;
      }
      // case of left-right models
      case true :
      {
         nb_component= nb_state;
         component_nb_state= new int[nb_component];
         component= new int*[nb_component];

         for(i= 0; i < nb_state; i++)
         {
            for(j= 0; j <= i; j++)
               accessibility[i][j]= false;

            for(j=i+1; j < nb_state; j++)
               accessibility[i][j] = true;

            component_nb_state[i]= 1;
            component[i]= new int[component_nb_state[i]];
            component[i][0]= i;

            if (i < nb_state - 1)
               state_type[i]= 't';
            else
               state_type[i] = 'a';
         }

         for(i= 0; i < nb_state-1; i++)
            initial[i]= 1. / (double)(nb_state - 1);

         initial[nb_state-1]= 0.;

         for(i= 0; i < _ch_order; i++)
            state_index[i]= 0;

         for(i= 0; i < nb_row; i++)
         {
            for(j= 0;j < state_index[_ch_order-1]; j++)
               transition[i][j]= 0.;

            for(j= 1;j < _ch_order;j++)
               if (state_index[j] < state_index[j-1])
                  break;

            if ((j == _ch_order) && (state_index[_ch_order-1] < nb_state-1))
            {
               transition[i][state_index[_ch_order - 1]] = self_transition;
               for(j= state_index[_ch_order-1]+1; j < nb_state; j++)
                  transition[i][j]
                     = (1. - self_transition) / (nb_state - (state_index[_ch_order-1] + 1));
            }
            else
            {
               transition[i][state_index[_ch_order-1]]= 1.;
               for(j= state_index[_ch_order-1]+1; j < nb_state; j++)
                  transition[i][j] = 0.;
            }

            for(j= 0; j < _ch_order; j++)
            {
               if (state_index[j] < nb_state-1)
               {
                  state_index[j]++;
                  break;
               }
               else
                  state_index[j] = 0;
            }
         }
      break;
      }
   }
}

/*****************************************************************
 *
 *  Initialization of the HiddenMarkovTree parameters
 *  (compliant with null probabilities)
 *  using a fixed probability of staying in any state
 *
 **/

void HiddenMarkovTree::init(double self_transition)
{
   register int i, j;
   int nb_positive, power= nb_row / nb_state, state;
   double self;


   nb_positive= 0;
   // number of states with positive initial probability
   for(i= 0; i < nb_state; i++)
      if (initial[i] > 0.)
         nb_positive++;

   for(i= 0; i < nb_state; i++)
      if (initial[i] > 0.)
         initial[i] = 1. / (double)nb_positive;

   for(i= 0; i < nb_row; i++)
   {
      state= i / power;
      if (transition[i][state] > 0.)
         self= self_transition;
      else
         self= 0.;

      nb_positive= 0;

      // number of possible transitions from current state to another state

      for(j= 0; j < state; j++)
         if (transition[i][j] > 0.)
            nb_positive++;

      for(j= state+1; j < nb_state; j++)
         if (transition[i][j] > 0.)
            nb_positive++;

      for(j= 0; j < state; j++)
         if (transition[i][j] > 0.)
            transition[i][j] = (1. - self) / nb_positive;

      if (transition[i][j] < 1.)
         transition[i][j]= self;

      for(j= state+1; j < nb_state; j++)
         if (transition[i][j] > 0.)
            transition[i][j] = (1. - self) / nb_positive;
   }
}

/*****************************************************************
 *
 *  Computation of the log parameters of a HiddenMarkovTree
 *
 **/

void HiddenMarkovTree::log_computation()
{
   register int i, j= 0;

   Chain::log_computation();

   for(i= 0; i < _nb_ioutput_process; i++)
   {
      if (npprocess[i+1] != NULL)
         for(j= 0; j < nb_state; j++)
            npprocess[i+1]->observation[j]->log_computation();
      else
         for(j= 0; j < nb_state; j++)
            ::log_computation(piprocess[i+1]->nb_value,
                              piprocess[i+1]->observation[j]->mass,
                              piprocess[i+1]->observation[j]->cumul);
   }

   for(i= 0; i < _nb_doutput_process; i++)
   {
      ::log_computation(pdprocess[i+1]->nb_value,
                        pdprocess[i+1]->observation[j]->mass,
                        pdprocess[i+1]->observation[j]->cumul);
   }


}

/*****************************************************************
 *
 *  Virtual functions of a HiddenMarkovTree which are only
 *  implemented in derived classes
 *
 **/

void HiddenMarkovTree::state_no_occurrence_probability(int state, double increment)
{}

void HiddenMarkovTree::state_first_occurrence_root_distribution(int state,
                                                                int min_nb_value,
                                                                double cumul_threshold)
{}

void HiddenMarkovTree::state_first_occurrence_leaves_distribution(int state, int min_nb_value,
                                                                  double cumul_threshold)
{}

void HiddenMarkovTree::state_leave_probability(const double * memory, int state,
                                               double increment)
{}

void HiddenMarkovTree::state_sojourn_size_distribution(const double * memory, int state,
                                                       int min_nb_value,
                                                       double cumul_threshold)
{}

void HiddenMarkovTree::state_nb_pattern_mixture(int state, char pattern)
{}

void HiddenMarkovTree::output_no_occurrence_probability(int variable, int output,
                                                        double increment)
{}

void HiddenMarkovTree::output_first_occurrence_root_distribution(int variable, int output,
                                                                 int min_nb_value,
                                                                 double cumul_threshold)
{}

void HiddenMarkovTree::output_first_occurrence_leaves_distribution(int variable, int output,
                                                                   int min_nb_value,
                                                                   double cumul_threshold)
{}

void HiddenMarkovTree::output_leave_probability(const double * memory,
                                                int variable, int output,
                                                double increment)
{}

void HiddenMarkovTree::output_sojourn_size_distribution(const double * memory, int variable,
                                                        int output, int min_nb_value,
                                                        double cumul_threshold)
{}

void HiddenMarkovTree::output_nb_zones_mixture(int variable, int output)
{}

void HiddenMarkovTree::output_nb_occurrences_mixture(int variable, int output)
{}

int HiddenMarkovTree::nb_parameter_computation(double min_probability) const
{ return I_DEFAULT; }

/****************************************************************
 *
 *  Compute an adaptative penalty used for model selection
 *  of HiddenMarkovTrees
 *
 **/

double HiddenMarkovTree::penalty_computation(double min_probability) const
{
   register int i, val, j, var;
   int nb_parameter, cumul_size; // sample_size,
   double sum, *memory= NULL , *state_marginal= NULL;
   double penalty= 0.;

   if (markov_data != NULL)
   {
      memory= memory_computation();

      state_marginal= new double[nb_state];

      switch (type)
      {
         case 'o' :
         {
            sum= 0.;
            for(val= 0; val < npprocess[0]->size->nb_value-2; val++)
               sum+= (1.-npprocess[0]->size->cumul[val+1]);

            for(j= 0; j < nb_state; j++)
               memory[j] /= sum;

            for(j= 0; j < nb_state; j++)
               state_marginal[j] = 0.;

            /*
            for(val= 0; val < npprocess[0]->size->nb_value-1; val++)
               for(j = 0;j < nb_state;j++)
                  state_marginal[j] += npprocess[0]->index_value->point[j][val] *
                                       (1.-categorical_process[0]->size->cumul[val]);
            */
            sum= 0.;

            for(j= 0; j < nb_state; j++)
               sum += state_marginal[j];

            for(j= 0; j < nb_state; j++)
               state_marginal[j] /= sum;

            break;
         }

         case 'e' :
         {
            for(j= 0; j < nb_state; j++)
               state_marginal[j] = initial[j];
            break;
         }
      } // end switch

      cumul_size= markov_data->cumul_size_computation();
      for(i= 0; i < nb_state; i++)
      {
         nb_parameter= 0;
         for(j= 0; j < nb_state; j++)
            if (transition[i][j] > min_probability)
               nb_parameter++;

         nb_parameter--;

         if ((nb_parameter > 0) && (memory[i] > 0.))
            penalty += nb_parameter * log(memory[i] * cumul_size);
      }

      for(j= 0; j < nb_state; j++)
      {
         // for parametric sojourn size distributions
         /*
         if (transition[j][j] == 0.)
         {
            nb_parameter= npprocess[0]->sojourn_size[j]->nb_parameter_computation();
            if (npprocess[0]->sojourn_size[j]->inf_bound == 1)
               nb_parameter--;

            penalty += nb_parameter * log(state_marginal[j] * cumul_size);
         }
         */
      }

      if (_nb_ioutput_process > 0)
      {
         for(var= 1; var <= _nb_ioutput_process; var++)
         {
            if (npprocess[var] != NULL)
            {
               for(j= 0; j < nb_state; j++)
               {
                  nb_parameter= 0;
                  for(val= 0; val < npprocess[var]->nb_value; val++)
                     if (npprocess[var]->observation[j]->mass[val] > min_probability)
                        nb_parameter++;

                  nb_parameter--;

                  if (nb_parameter > 0)
                     penalty+= nb_parameter * log(state_marginal[j] * cumul_size);
               }
            }
            else
            {
               for(j= 0; j < nb_state; j++)
               {
                  nb_parameter= piprocess[var]->observation[j]->nb_parameter_computation();

                  penalty+= nb_parameter * log(state_marginal[j] * cumul_size);
               }
            }
         }
      }

      delete [] memory;
      memory= NULL;
      delete [] state_marginal;
      state_marginal= NULL;
   }
   return penalty;
}


double*** HiddenMarkovTree::state_marginal_distribution(const HiddenMarkovTreeData& trees) const
{ return NULL; }

double** HiddenMarkovTree::state_marginal_distribution(const Trees& trees,
                                                       int index) const
{ return NULL; }

/*****************************************************************
 *
 *  Computation of the observed_data conditional distributions
 *  for HiddenMarkovTree class
 *  using a HiddenMarkovTreeData object,
 *  the stored conditional probabilities,
 *  a flag on the computation of the log probabilities
 *  and the index of considered tree
 *
 **/
void HiddenMarkovTree::output_conditional_distribution(const HiddenMarkovTreeData& trees,
                                                       double_array_3d& output_cond,
                                                       bool log_computation,
                                                       int index) const
{
   typedef HiddenMarkovTreeData::tree_type tree_type;
   typedef tree_type::vertex_iterator vertex_iterator;
   typedef tree_type::value value;

   int nb_trees= trees._nb_trees, t, current_size;
   register int j, var;
   Typed_edge_int_fl_tree<Int_fl_container> *current_tree;
   vertex_iterator it, end;
   value val;

   // output_cond[t][j][u] corresponds to conditional distribution, given the state variable
   // is equal to j, taken at the value of node u of tree t

   assert((_nb_ioutput_process == trees._nb_integral)
          && (_nb_doutput_process == trees._nb_float));
   if (output_cond == NULL)
   {
      output_cond= new double_array_2d[nb_trees];
      for(t= 0; t < nb_trees; t++)
      output_cond[t]= NULL;
   }

   for(t= 0; t < nb_trees; t++)
      if ((index == I_DEFAULT) || (index == t))
      {
         if (output_cond[t] == NULL)
         {
            output_cond[t]= new double*[nb_state];
            for(j= 0; j < nb_state; j++)
            output_cond[t][j]= NULL;
         }

         current_tree= trees.trees[t];
         current_size= current_tree->get_size();
         for(j= 0; j < nb_state; j++)
            if (output_cond[t][j] == NULL)
               output_cond[t][j]= new double[current_size];

         Tree_tie::tie(it, end)= current_tree->vertices();
         while (it < end)
         {
            val= current_tree->get(*it);
            for(j= 0; j < nb_state; j++)
            {
               if (log_computation)
                  output_cond[t][j][*it]= 0.;
               else
                  output_cond[t][j][*it]= 1.;
               for(var= 0; var < _nb_ioutput_process; var++)
                  if (npprocess[var+1] != NULL)
                  {
                     if (log_computation)
                     {
#                       ifdef DEBUG
                        if (((npprocess[var+1]->observation[j]->mass[val.Int(var)] > 0)
                             && (abs(log(npprocess[var+1]->observation[j]->mass[val.Int(var)])-npprocess[var+1]->observation[j]->cumul[val.Int(var)]) > DOUBLE_ERROR))
                           || ((npprocess[var+1]->observation[j]->mass[val.Int(var)] == 0) && (abs(npprocess[var+1]->observation[j]->cumul[val.Int(var)]-D_INF) > DOUBLE_ERROR)))
                           cout << "Warning: computation error at observation[" << var+1 << "]["
                                << j << "]." << endl;
#                       endif
                        // output_cond[t][j][*it]+= npprocess[var+1]->observation[j]->cumul[val.Int(var)];
                        if (npprocess[var+1]->observation[j]->mass[val.Int(var)] != 0)
                           output_cond[t][j][*it]+= log(npprocess[var+1]->observation[j]->mass[val.Int(var)]);
                        else
                           output_cond[t][j][*it]= D_INF;
                     }
                     else
                        output_cond[t][j][*it]*= npprocess[var+1]->observation[j]->mass[val.Int(var)];
                  }
                  else
                  {
                     if (log_computation)
                     {
#                       ifdef DEBUG
                        if (((piprocess[var+1]->observation[j]->mass[val.Int(var)] > 0)
                             && (abs(log(piprocess[var+1]->observation[j]->mass[val.Int(var)])-piprocess[var+1]->observation[j]->cumul[val.Int(var)]) > DOUBLE_ERROR))
                           || ((piprocess[var+1]->observation[j]->mass[val.Int(var)] == 0) && (abs(piprocess[var+1]->observation[j]->cumul[val.Int(var)]-D_INF) > DOUBLE_ERROR)))
                           cout << "Warning: computation error at transition[" << var+1 << "]["
                                << j << "]." << endl;
#                       endif
                           // output_cond[t][j][*it]+= piprocess[var+1]->observation[j]->cumul[val.Int(var)];
                           if (piprocess[var+1]->observation[j]->mass[val.Int(var)] != 0)
                              output_cond[t][j][*it]+= log(piprocess[var+1]->observation[j]->mass[val.Int(var)]);
                           else
                              output_cond[t][j][*it]= D_INF;
                     }
                     else
                        output_cond[t][j][*it]*= piprocess[var+1]->observation[j]->mass[val.Int(var)];
                  }

               /* for(var= 0; var < _nb_doutput_process; var++)
                  output_cond[t][j][*it]*= pdprocess[var]->observation[j]->mass[val.Double(var)];*/
               // case of floating observed processes, not implemented for the moment
            }
            it++;
         }
   }
}

double HiddenMarkovTree::likelihood_correction(const HiddenMarkovTreeData& trees) const
{ return D_INF; }

double HiddenMarkovTree::upward_downward(const HiddenMarkovTreeData& trees,
                                         double& max_marginal_entropy,
                                         double& entropy1,
                                         double& likelihood,
                                         std::deque<int>*& vd,
                                         int index, ostream* os,
                                         char format, int vertex,
                                         int entropy_algo) const
{ return D_INF; }

/*****************************************************************
 *
 *  Compute the entropy of partial state processes
 *
 **/

double HiddenMarkovTree::partial_entropy_computation(const HiddenMarkovTreeData& trees,
                                                     int t,
                                                     double_array_3d output_cond_prob,
                                                     double_array_3d marginal_prob,
                                                     double_array_3d upward_parent_prob,
                                                     double_array_3d downward_prob,
                                                     double_array_4d downward_pair_prob,
                                                     double_array_3d state_entropy,
                                                     double_array_3d conditional_entropy,
                                                     double_array_4d conditional_prob,
                                                     double_array_2d& partial_entropy,
                                                     int entropy_algo) const
{
   double res;
   int tr;

   switch (entropy_algo)
   {
      case UPWARD:
      {
         if (partial_entropy == NULL)
         {
            partial_entropy = new double*[trees.get_nb_trees()];
            for(tr = 0; tr < trees.get_nb_trees(); tr++)
               partial_entropy[tr] = NULL;
         }

         res= upward_partial_entropy_computation(trees, t, downward_prob,
                                                 state_entropy, partial_entropy[t]);
         break;
      }

      case DOWNWARD:
      {
         res= downward_partial_entropy_computation(trees, t,
                                                   output_cond_prob,
                                                   marginal_prob,
                                                   upward_parent_prob,
                                                   downward_prob,
                                                   downward_pair_prob,
                                                   state_entropy,
                                                   conditional_entropy,
                                                   conditional_prob,
                                                   partial_entropy);
         break;
      }
   }
   return res;
}

/*************************************************************************
 *
 *  Compute the entropy of partial state processes by a downward algorithm
 *
 **/

double HiddenMarkovTree::downward_partial_entropy_computation(const HiddenMarkovTreeData& trees,
                                                              int t, double_array_3d output_cond_prob,
                                                              double_array_3d marginal_prob,
                                                              double_array_3d upward_parent_prob,
                                                              double_array_3d downward_prob,
                                                              double_array_4d downward_pair_prob,
                                                              double_array_3d state_entropy,
                                                              double_array_3d conditional_entropy,
                                                              double_array_4d conditional_prob,
                                                              double_array_2d& partial_entropy) const
{ return D_INF; }

/*************************************************************************
 *
 *  Compute the entropy of partial state processes by an upward algorithm
 *
 **/

double HiddenMarkovTree::upward_partial_entropy_computation(const HiddenMarkovTreeData& trees,
                                                            int t,
                                                            double_array_3d downward_prob,
                                                            double_array_3d state_entropy,
                                                            double*& partial_entropy) const
{ return D_INF; }

/*************************************************************************
 *
 *  Compute the state process' conditional entropy
 *
 **/

void HiddenMarkovTree::conditional_entropy_computation(const HiddenMarkovTreeData& trees,
                                                       double_array_3d marginal_prob,
                                                       double_array_3d upward_prob,
                                                       double_array_3d upward_parent_prob,
                                                       double_array_3d downward_prob,
                                                       double_array_2d& expected_conditional_entropy,
                                                       double_array_3d& conditional_entropy,
                                                       double_array_4d& conditional_prob,
                                                       double_array_3d& state_entropy,
                                                       int index,
                                                       int entropy_algo) const
{
   switch (entropy_algo)
   {

      case UPWARD:
      {
         upward_conditional_entropy_computation(trees, marginal_prob, upward_prob,
                                                upward_parent_prob,
                                                downward_prob,
                                                expected_conditional_entropy,
                                                index);
         break;
      }

      case DOWNWARD:
      {
         downward_conditional_entropy_computation(trees, marginal_prob,
                                                  downward_prob,
                                                  upward_prob,
                                                  upward_parent_prob,
                                                  expected_conditional_entropy,
                                                  conditional_entropy,
                                                  conditional_prob,
                                                  state_entropy,
                                                  index);
         break;
      }
   }
}

/*************************************************************************
 *
 *  Compute the state process' conditional entropy by a downward algorithm
 *
 **/

double HiddenMarkovTree::downward_conditional_entropy_computation(const HiddenMarkovTreeData& trees,
                                                                  double_array_3d marginal_prob,
                                                                  double_array_3d downward_prob,
                                                                  double_array_3d upward_prob,
                                                                  double_array_3d upward_parent_prob,
                                                                  double_array_2d& expected_conditional_entropy,
                                                                  double_array_3d& conditional_entropy,
                                                                  double_array_4d& conditional_prob,
                                                                  double_array_3d& state_entropy,
                                                                  int index) const
{ return -D_INF; }

/*************************************************************************
 *
 *  Compute the state process' conditional entropy by an upward algorithm
 *
 **/

double HiddenMarkovTree::upward_conditional_entropy_computation(const HiddenMarkovTreeData& trees,
                                                                double_array_3d marginal_prob,
                                                                double_array_3d upward_prob,
                                                                double_array_3d upward_parent_prob,
                                                                double_array_3d downward_prob,
                                                                double_array_2d& conditional_entropy,
                                                                int index) const
{ return -D_INF; }

double HiddenMarkovTree::smoothed_probabilities(const HiddenMarkovTreeData& trees,
                                                double_array_3d& smoothed_prob,
                                                double_array_2d& marginal_entropy,
                                                double_array_2d& conditional_entropy,
                                                double_array_2d& partial_entropy,
                                                int index,
                                                int entropy_algo) const
{ return D_INF; }

double HiddenMarkovTree::viterbi(const HiddenMarkovTreeData& trees,
                                 int index) const
{ return D_INF; }

long double HiddenMarkovTree::nb_state_trees(const HiddenMarkovTreeData& trees,
                                             int index) const
{ return D_DEFAULT; }

HiddenMarkovTreeData* HiddenMarkovTree::viterbi_upward_downward(const HiddenMarkovTreeData& trees,
                                                                std::vector<ostringstream*>& messages,
                                                                double likelihood,
                                                                double& state_likelihood,
                                                                std::deque<int>*& vd,
                                                                int index,
                                                                std::ostream* os,
                                                                char format,
                                                                int vertex) const
{ return NULL; }

HiddenMarkovTreeData* HiddenMarkovTree::generalized_viterbi(const HiddenMarkovTreeData& trees,
                                                            std::vector<ostringstream*>& messages,
                                                            int nb_state_trees,
                                                            double likelihood,
                                                            int index,
                                                            int sroot) const
{ return NULL; }

double HiddenMarkovTree::state_likelihood_computation(const HiddenMarkovTreeData& trees) const
{ return D_DEFAULT; }

double HiddenMarkovTree::state_likelihood_computation(const HiddenMarkovTreeData& trees,
                                                      int index) const
{ return D_DEFAULT; }

/*****************************************************************
 *
 *  Creates a HiddenMarkovTree from a file
 *  using a StatError object, the path,
 *  and the tree sizes (or other quantities for which
 *  the characteristic distributions are invariant - if any),
 *
 **/

HiddenMarkovTree* Stat_trees::hidden_markov_tree_ascii_read(StatError& error,
                                                            const char * path,
                                                            int size, bool counting_flag,
                                                            double cumul_threshold)
{
   RWLocaleSnapshot locale("en");
   RWCString buffer, token, tree_ident;
   size_t position;
   char type= 'v';
   bool status, lstatus, categorical= false;
   register int i;
   int line, ch_order, nb_state, index, // j, n,
       nb_ioutput_process= 0, nb_doutput_process, nb_output_process= 0;
   long value;
   const Chain *chain;
   CategoricalProcess **np_observation;
   DiscreteParametricProcess **ip_observation;
   DiscreteParametricProcess **dp_observation= NULL;
   HiddenMarkovTree *markov= NULL;
   ifstream in_file(path);

   error.init();

   if (!in_file)
      error.update(STAT_error[STATR_FILE_NAME]);
   else
   {
      status= true;
      line= 0;

      if (size < 2)
      {
         status= false;
         error.update(STAT_TREES_error[TREESTATR_SMALL_TREE_SIZE]);
      }

      if (size > MAX_SIZE)
      {
         status= false;
         error.update(STAT_TREES_error[TREESTATR_BIG_TREE_SIZE]);
      }

      while (buffer.readLine(in_file, false))
      {
         line++;

#        ifdef DEBUG
         cout << line << "  " << buffer << endl;
#        endif

         position= buffer.first('#');
         if (position != RW_NPOS)
            buffer.remove(position);

         i= 0;

         RWCTokenizer next(buffer);

         while (!((token = next()).isNull()))
         {
           // test keyword (EQUILIBRIUM) HIDDEN_MARKOV_TREE

            if (i == 0)
            {
               tree_ident= token;
               if (token == STAT_TREES_word[TREESTATW_HIDDEN_MARKOV_IND_OUT_TREE])
                   // || (token == STAT_TREES_word[TREESTATW_HIDDEN_MARKOV_IN_TREE]))
                  type= 'o';
               else
                  if (token == STAT_TREES_word[TREESTATW_EQUILIBRIUM_HIDDEN_MARKOV_IND_OUT_TREE])
                      // || (token == STAT_TREES_word[TREESTATW_EQUILIBRIUM_HIDDEN_MARKOV_IN_TREE]))
                     type= 'e';
                  else
                  {
                     status= false;
                     ostringstream correction_message;
                     correction_message << STAT_TREES_word[TREESTATW_HIDDEN_MARKOV_TREE];
                     error.correction_update (STAT_parsing[STATP_KEY_WORD],
                                             (correction_message.str()).c_str(),
                                             line);
                  }
            }
           i++;
         }

         if (i > 0)
         {
            if (i != 1)
            {
               status= false;
               error.update(STAT_parsing[STATP_FORMAT], line);
            }
            break;
         }
      }

      // analysis of the Markov tree format and processing

      if (type != 'v')
      {
         if ((tree_ident == STAT_TREES_word[TREESTATW_HIDDEN_MARKOV_IND_OUT_TREE]) ||
             (tree_ident == STAT_TREES_word[TREESTATW_EQUILIBRIUM_HIDDEN_MARKOV_IND_OUT_TREE]))
            {
               ch_order= 1;
               chain= chain_parsing(error, in_file, line, type);
            }
         else
            chain= chain_parsing(error, in_file, line, type);

         if (chain != NULL)
         {
            nb_ioutput_process= I_DEFAULT;
            nb_doutput_process= I_DEFAULT;

            nb_state= chain->nb_state;

            np_observation= NULL;
            ip_observation= NULL;
            dp_observation= NULL;

            while (buffer.readLine(in_file , false))
            {
               line++;

#              ifdef DEBUG
               cout << line << "  " << buffer << endl;
#              endif

               position = buffer.first('#');
               if (position != RW_NPOS)
                  buffer.remove(position);

               i= 0;

               RWCTokenizer next(buffer);

               while (!((token = next()).isNull()))
               {
                  switch (i)
                  {
                     // test on the number of observed processes

                     case 0 :
                     {
                        lstatus= locale.stringToNum(token, &value);
                        if (lstatus)
                        {
                           if ((value < 1) || (value > NB_OUTPUT_PROCESS))
                              lstatus = false;
                           else
                              nb_output_process= value;
                        }

                        if (!lstatus)
                        {
                           status= false;
                           error.update(STAT_parsing[STATP_NB_OUTPUT_PROCESS],
                                        line, i+1);
                        }
                        break;
                     }

                     // test on keyword OUTPUT_PROCESS(ES)

                     case 1 :
                     {
                        if (token != STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES])
                        {
                           status = false;
                           error.correction_update(STAT_parsing[STATP_KEY_WORD] ,
                                                   STAT_word[nb_output_process == 1 ? STATW_OUTPUT_PROCESS : STATW_OUTPUT_PROCESSES],
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

            if (nb_output_process == I_DEFAULT)
            {
               status = false;
               error.update(STAT_parsing[STATP_FORMAT] , line);
            }
            else
            {
               np_observation= new CategoricalProcess*[nb_output_process];
               ip_observation = new DiscreteParametricProcess*[nb_output_process];
               for(i= 0; i < nb_output_process; i++)
               {
                  np_observation[i]= NULL;
                  ip_observation[i]= NULL;
               }

               index= 0;

               while (buffer.readLine(in_file , false))
               {
                  line++;

#                 ifdef DEBUG
                  cout << line << "  " << buffer << endl;
#                 endif

                  position = buffer.first('#');
                  if (position != RW_NPOS)
                     buffer.remove(position);

                  i= 0;

                  RWCTokenizer next(buffer);

                  while (!((token = next()).isNull()))
                  {
                     switch (i)
                     {

                     // test on keyword OUTPUT_PROCESS

                        case 0 :
                        {
                           categorical= true;

                           if (token == STAT_word[STATW_OUTPUT_PROCESS])
                              index++;
                           else
                           {
                              status= false;
                              error.correction_update(STAT_parsing[STATP_KEY_WORD], STAT_word[STATW_OUTPUT_PROCESS],
                                                      line, i+1);
                           }
                           break;
                        }

                        // test on the observation process index

                        case 1 :
                        {
                           lstatus= locale.stringToNum(token, &value);
                           if ((lstatus) && ((value != index) || (value > nb_output_process)))
                              lstatus= false;

                           if (!lstatus)
                           {
                              status= false;
                              error.update(STAT_parsing[STATP_OUTPUT_PROCESS_INDEX],
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

                        // test on keyword CATEGORICAL / DISCRETE_PARAMETRIC

                        case 3 :
                        {
              if ((token == STAT_word[STATW_CATEGORICAL]) || (token == STAT_word[STATW_NONPARAMETRIC]))
                              categorical= true;
                           else
                           {
                             if ((token == STAT_word[STATW_DISCRETE_PARAMETRIC]) || (token == STAT_word[STATW_PARAMETRIC]))
                                 categorical= false;
                              else
                              {
                                 status= false;
                                 ostringstream correction_message;
                                 correction_message << STAT_word[STATW_CATEGORICAL] << " or "
                                                    << STAT_word[STATW_DISCRETE_PARAMETRIC] << " or "
                                                    << STAT_word[STATW_NONPARAMETRIC] << " or "  // pour compatibilite ascendante
                                                    << STAT_word[STATW_PARAMETRIC];
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

                     switch (categorical)
                     {

                        case true :
                        {
                           np_observation[index-1]= categorical_observation_parsing(error, in_file, line,
                                                                                    chain->nb_state,
                                                                                    HIDDEN_MARKOV, true);
                           // ip_observation[index-1]= NULL;
                           if (np_observation[index-1] == NULL)
                              status= false;
                           break;
                        }

                        case false :
                        {
                           ip_observation[index-1]= discrete_observation_parsing(error, in_file, line,
                                                                                 chain->nb_state, HIDDEN_MARKOV,
                                                                                 cumul_threshold);
                           // np_observation[index-1]= NULL;
                           if (ip_observation[index-1] == NULL)
                              status = false;
                           break;
                        }
                     }
                  }
               } // end : while (buffer.readLine(in_file , false))

               if (index != nb_output_process)
               {
                  status= false;
                  error.update(STAT_parsing[STATP_FORMAT] , line);
               }

               if (status)
                  markov= new HiddenMarkovTree(chain, ch_order, nb_output_process, 0,
                                               np_observation, ip_observation,
                                               dp_observation, size, counting_flag);

               for(i= 0; i < nb_output_process; i++)
               {
                  if (np_observation[i] != NULL)
                     delete np_observation[i];
                  if (ip_observation[i] != NULL)
                     delete ip_observation[i];
               }

               if (nb_output_process > 0)
               {
                  delete [] np_observation;
                  delete [] ip_observation;
               }

               if (chain != NULL)
                  delete chain;

            } // end if (nb_output_process == I_DEFAULT)
         }
      }
   }
   return markov;
}

/*****************************************************************
 *
 * Left (bit) shift operator of a HiddenMarkovTree
 *
 **/

std::ostream& Stat_trees::operator<<(std::ostream& os, const HiddenMarkovTree& markov)
{ return markov.ascii_write(os, markov.markov_data); }


/*****************************************************************
 *
 *  Default constructor of HiddenMarkovTreeData class
 *
 **/

HiddenMarkovTreeData::HiddenMarkovTreeData()
 : Trees() // unsure
 , markov(NULL)
 , chain_data(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(I_DEFAULT)
 , sample_entropy(D_INF)
 , entropy(NULL)
 , state_trees(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
{}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovTreeData class
 *  using the number of integral and float variables,
 *  frequency distributions for the size and the number of children,
 *  a flag on the possibility for a node to have no child due to
 *  random and an initialization flag
 *
 **/

HiddenMarkovTreeData::HiddenMarkovTreeData(int inb_integral,
                                           int inb_float,
                                           const FrequencyDistribution& ihsize,
                                           const FrequencyDistribution& ihnb_children,
                                           bool no_child_flag,
                                           bool init_flag)
 : Trees(inb_integral, inb_float,
         ihsize, ihnb_children,
         no_child_flag, init_flag)
 , markov(NULL)
 , chain_data(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(I_DEFAULT)
 , sample_entropy(D_INF)
 , entropy(NULL)
 , state_trees(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
{
   typedef Typed_edge_one_int_tree::value value;
   value default_value;
   int t;
   Unlabelled_typed_edge_tree *utree;

   if (_nb_integral+_nb_float > 0)
   {
      this->state_trees= new Typed_edge_one_int_tree*[_nb_trees];

      default_value.Int()= I_DEFAULT;
      // with appropriate namespace

      if (init_flag)
      {
         for(t= 0; t < _nb_trees; t++)
         {
            this->state_trees[t]= new Typed_edge_one_int_tree;
            utree= trees[t]->get_structure();
            this->state_trees[t]->set_structure(*utree, default_value);
            delete utree;
         }
      }
      // it might be more relevant to simulate the data using markov, if available
      else
         for(t= 0; t < _nb_trees; t++)
            this->state_trees[t]= NULL;
   }
}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovTreeData class
 *  using the number of integral and float variables
 *  and the number of trees
 *
 **/

HiddenMarkovTreeData::HiddenMarkovTreeData(int inb_integral,
                                           int inb_float,
                                           int inb_trees)
 : Trees(inb_integral, inb_float, inb_trees)
 , markov(NULL)
 , chain_data(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(I_DEFAULT)
 , sample_entropy(D_INF)
 , entropy(NULL)
 , state_trees(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
{}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovTreeData class
 *  using the number of trees, the type of each variable
 *  and the multidimensional observed trees
 *
 **/

HiddenMarkovTreeData::HiddenMarkovTreeData(int inb_trees,
                                           int* itype,
                                           Default_tree** otrees)
 : Trees(inb_trees, itype, otrees)
 , markov(NULL)
 , chain_data(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(I_DEFAULT)
 , sample_entropy(D_INF)
 , entropy(NULL)
 , state_trees(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
{}

/*****************************************************************
 *
 *  Constructor of HiddenMarkovTreeData class
 *  from a Trees object
 *
 **/

HiddenMarkovTreeData::HiddenMarkovTreeData(const Trees& otrees)
 : Trees(otrees)
 , markov(NULL)
 , chain_data(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(I_DEFAULT)
 , sample_entropy(D_INF)
 , entropy(NULL)
 , state_trees(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
{}


/*****************************************************************
 *
 *  Copy constructor of HiddenMarkovTreeData class
 *  using a flag on the HiddenMarkovTree copy
 *
 **/

HiddenMarkovTreeData::HiddenMarkovTreeData(const HiddenMarkovTreeData& trees,
                                           bool model_flag,
                                           bool characteristic_flag)

 : Trees(trees)
 , markov(NULL)
 , chain_data(NULL)
 , likelihood(D_INF)
 , hidden_likelihood(D_INF)
 , _nb_states(I_DEFAULT)
 , sample_entropy(D_INF)
 , entropy(NULL)
 , state_trees(NULL)
 , observation_distribution(NULL)
 , observation_histogram(NULL)
 , state_characteristics(NULL)
{ copy(trees , model_flag, characteristic_flag); }

/*****************************************************************
 *
 *  Destructor of HiddenMarkovTreeData class
 *
 **/

HiddenMarkovTreeData::~HiddenMarkovTreeData()
{
   remove();
   Trees::remove(); // unsure
}

/*****************************************************************
 *
 *  Assignement operator of HiddenMarkovTreeData class
 *
 **/

HiddenMarkovTreeData& HiddenMarkovTreeData::operator=(const HiddenMarkovTreeData& trees)
{
   if (&trees != this)
   {
      remove();
      Trees::remove();

      Trees::copy(trees);
      copy(trees);
   }
   return *this;
}

/*****************************************************************
 *
 *  FrequencyDistribution extraction for Hidden_markov_data class,
 *  using a StatError object, the kind of frequency distribution,
 *  the considered variable or state and value
 *  (case of observation frequency distributions)
 *
 **/

DiscreteDistributionData* HiddenMarkovTreeData::extract(StatError& error, int type,
                                                        int variable, int value) const

{
   bool status= true;
   Distribution *pdist= NULL;
   FrequencyDistribution *phisto= NULL;
   DiscreteDistributionData *histo= NULL;

   error.init();

   // not to be used for HiddenMarkovTreeData objects with embedded
   // state tree (i.e. the state is not considered as a variable)

   // consequently, the number of observed processes for markov and
   // for Trees differ by one

   if (_type[0] != STATE)
   {
      status= false;
      error.update(STAT_TREES_error[TREESTATR_STATE_TREES]);
   }
   else
   {
      // frequency distribution part
      if (type == OBSERVATION)
      {
         if ((variable < 2) || (variable > _nb_integral))
         // I guess this is nonsense for a floating variable ...
         // Furthermore, variable 1 supposedly refers to the state variable
         {
            status= false;
            error.update(STAT_error[STATR_VARIABLE_INDEX]);
         }

         else
         {
             variable--;

             if ((value < 0) || (value >= state_characteristics->get_nb_values()))
             // value represents a (hidden) state
             {
                status= false;
                ostringstream error_message;
                error_message << STAT_label[STATL_STATE] << " " << value << " "
                              << STAT_error[STATR_NOT_PRESENT];
                error.update((error_message.str()).c_str());
             }

             else
             {
                phisto= observation_distribution[variable][value];

                if (phisto->nb_element == 0)
                {
                   status= false;
                   error.update(STAT_error[STATR_EMPTY_SAMPLE]);
                }
             }
         }

         if (status)
         {
            pdist= NULL;

            if (markov->npprocess[variable] != NULL) // ->observation
               pdist= markov->npprocess[variable]->observation[value];
            else
               pdist= markov->piprocess[variable]->observation[value];

            if (pdist != NULL)
               histo= new DiscreteDistributionData(*phisto, pdist);
            // else
            //     there might be something with a Param*

         }
      }
      else // type != OBSERVATION
      {
         phisto= Trees::extract(error, type, variable, value);
         if (phisto == NULL)
            status= false;
      }

      if (status)
      // distribution part
      {
         pdist= NULL;
         // variable-1 corresponds to variable for markov
         switch (type)
         {
            case FIRST_OCCURRENCE_ROOT :
               pdist= markov->npprocess[variable-1]->first_occurrence_root[value];
               break;
            case FIRST_OCCURRENCE_LEAVES :
               pdist= markov->npprocess[variable-1]->first_occurrence_leaves[value];
               break;
            case SOJOURN_SIZE :
               pdist= markov->npprocess[variable-1]->sojourn_size[value];
               break;
            case NB_ZONES :
               pdist= markov->npprocess[variable-1]->nb_zones[value];
               break;
            case NB_OCCURRENCES :
               pdist= markov->npprocess[variable-1]->nb_occurrences[value];
               break;
            case OBSERVATION :
               if (markov->npprocess[variable] != NULL) // ->observation
                  pdist= markov->npprocess[variable]->observation[value];
               else
                  pdist= markov->piprocess[variable]->observation[value];
         }

         if (pdist != NULL)
            histo= new DiscreteDistributionData(*phisto, pdist);
         // else
         //     there might be something with a Param*

      }

   } // type != STATE
   return histo;
}

/*****************************************************************
 *
 *  Return the model part of a HiddenMarkovTreeData
 *  using a StatError object, keeping a reference on self
 *
 **/

HiddenMarkovTree* HiddenMarkovTreeData::extract_model(StatError& error) const
{
   bool status = true, ignore_state_variable = false;
   HiddenMarkovTree *model = NULL;
   HiddenMarkovTreeData *hmt_data = NULL;

   error.init();

   if (markov == NULL)
   {
      status = false;
      error.update(STAT_TREES_error[TREESTATR_NO_MODEL]);
   }

   // if the first variable corresponds to the state variable,
   // ignore this variable
   if ((_nb_integral > 0)
       && (_nb_integral-1 == markov->_nb_ioutput_process)
       && (this->_type[0] == STATE))
       ignore_state_variable = true;
   if (!(ignore_state_variable)
         && ((_nb_integral != markov->_nb_ioutput_process))
      || (_nb_float != markov->_nb_doutput_process))
   {
      status = false;
      error.update(STAT_error[STATR_NB_OUTPUT_PROCESS]);
   }

   if (status)
   {
      model = markov->HiddenMarkovTreeCopy(true, true);
      if (ignore_state_variable)
      {
         hmt_data = this->remove_state_variable();
         model->markov_data = new HiddenMarkovTreeData(*hmt_data, false);
         delete hmt_data;
         hmt_data = NULL;
      }
      else
         model->markov_data = new HiddenMarkovTreeData(*this, false);
     model->markov_data->likelihood = model->likelihood_computation(*model->markov_data);

   }

   return model;
}

/*****************************************************************
 *
 *  Extract marginal Histogram with observation distributions
 *  for Hidden_markov_data class, using a StatError object
 *  and the considered variable
 *
 **/

DiscreteMixtureData* HiddenMarkovTreeData::extract_marginal(StatError& error,
                                                            int variable) const

{
   register int j;
   int ivariable = I_DEFAULT;
   bool status = true, mixture = false;
   double *pweight = NULL;
   FrequencyDistribution *pstate_marginal = NULL; // state_marginal
   DiscreteDistributionData *phisto = NULL, // marginal histogram
                            *pstate_marginal_ddd = NULL; // state_marginal
   DiscreteMixtureData *histo = NULL;
   DiscreteMixture *pmixt = NULL;
   // DiscreteParametric **pcomp = NULL;
   const DiscreteParametric **pcomp = new const DiscreteParametric*[_nb_states];
   // DiscreteParametricModel *pm = NULL;

   error.init();

   // not to be used for HiddenMarkovTreeData objects with embedded
   // state tree (i.e. the state is not considered as a variable)

   // consequently, the number of observed processes for markov and
   // for Trees differ by one

   if ((_type[0] != STATE) && (state_trees == NULL))
   {
      status = false;
      error.update(STAT_TREES_error[TREESTATR_STATE_TREES]);
   }
   else
   {
      if ((variable < 1) || (variable > _nb_integral))
      {
         status = false;
         error.update(STAT_error[STATR_VARIABLE_INDEX]);
      }
      else
      {
         if ((status) && (_type[0] == STATE) && (variable == 1))
         {
            // cannot plot state variable versus itself
            status = false;
            error.update(STAT_error[STATR_VARIABLE_INDEX]);
         }
         else
         {
            // marginal histogram of target variable
            phisto = Trees::extract(error, variable);

            if (phisto == NULL)
            {
               status = false;
               error.update(STAT_error[STATR_EMPTY_SAMPLE]);
            }

            if (_type[0] == STATE)
            {
               pstate_marginal_ddd = Trees::extract(error, 1);
               if (pstate_marginal_ddd != NULL)
               {
                  pstate_marginal = new FrequencyDistribution(*pstate_marginal_ddd);
                  delete pstate_marginal_ddd;
                  pstate_marginal_ddd = NULL;
               }
               ivariable = variable-1; // HMT variable != trees variable
            }
            else
            {
               ivariable = variable; // HMT variable == trees variable
               if (state_characteristics)
                  pstate_marginal = state_characteristics->marginal_distribution;
            }

            if (pstate_marginal == NULL)
            {
               status = false;
               error.update(STAT_error[STATR_EMPTY_SAMPLE]);
            }
         }
      }
      if (status)
      {
         // mixture weights
         pweight = new double[_nb_states];

         for(j = 0; j < _nb_states; j++)
            pweight[j] = (double)(pstate_marginal->frequency[j]) /
                            pstate_marginal->nb_element;

         // add observation distributions iif markov part is present
         if (this->markov != NULL)
            mixture = true;

         if (mixture)
         {
            // pcomp = new DiscreteParametric*[_nb_states];
            if (markov->npprocess[ivariable] != NULL)
            {
               for(j = 0; j < _nb_states; j++)
                  pcomp[j] =
                     new DiscreteParametric(*markov->npprocess[ivariable]->observation[j]);
            }
            else
            {
               for(j = 0; j < _nb_states; j++)
                  pcomp[j] =
                     new DiscreteParametric(*markov->piprocess[ivariable]->observation[j]);
            }
         }

         if (mixture)
            pmixt = new DiscreteMixture(_nb_states, pweight, pcomp);

         if (pmixt != NULL)
            histo = new DiscreteMixtureData(*phisto, pmixt);
         else
            histo = new DiscreteMixtureData(*phisto, _nb_states);

         if (pmixt != NULL)
         {
            delete pmixt;
            pmixt = NULL;
         }

         if (mixture)
         {
            for(j = 0; j < _nb_states; j++)
            {
               delete pcomp[j];
               pcomp[j]= NULL;
            }
            delete [] pcomp;
            pcomp = NULL;
         }
         delete [] pweight;
         pweight = NULL;
         delete pmixt;
         pmixt = NULL;
         // delete pm;
         // pm = NULL;
      }
      if (phisto != NULL)
      {
         delete phisto;
         phisto = NULL;
      }
      if ((pstate_marginal != NULL) && (_type[0] == STATE))
      {
         delete pstate_marginal;
         pstate_marginal = NULL;
      }
   }
   return histo;
}

/*****************************************************************
 *
 *  Merge HiddenMarkovTreeData, using a StatError object
 *  and a collection of HiddenMarkovTreeData objects.
 *  The collection otrees of HiddenMarkovTreeData objects is
 *  merged with *this to create a new object
 *  beginning with *this and adding the other trees by increasing index.
 *  State_trees are also merged if possible.
 *  Model part is not taken into account.
 *
 **/

HiddenMarkovTreeData*
HiddenMarkovTreeData::merge(StatError& error,
                            const pt_hmtd_vector& otrees) const
{
   bool status= true;
   unsigned int t, i, // index of trees object in otrees
                tree_index; // index of current tree object in result trees
   const unsigned int nb_sample= otrees.size();
   ostringstream error_message;
   state_tree_type state_res_t;
   Trees *res_t; // result as tree_type
   HiddenMarkovTreeData *res= NULL;
   pt_observed_trees_array hmtd_array= NULL;

   hmtd_array= new observed_trees*[nb_sample];

   for(t= 0; t < nb_sample; t++)
      hmtd_array[t]= otrees[t];

   error.init();

   res_t= this->Trees::merge(error, nb_sample, hmtd_array);

   if (error.get_nb_error() == 0)
   {
      res= new HiddenMarkovTreeData(*res_t);
      if (this->state_trees != NULL)
      {
         for(t= 0; t < nb_sample; t++)
         {
            if (otrees[t]->state_trees == NULL)
            {
               status= false;
               error_message << STAT_label[STATL_SAMPLE] << " " << t+2 << ": "
                             << STAT_TREES_error[TREESTATR_STATE_TREES];
            }
         }
         if (status)
         {
            res->state_trees= new Typed_edge_one_int_tree*[res->_nb_trees];
            // copy state trees for this
            tree_index= 0;
            for(t= 0; t < (unsigned int)this->_nb_trees; t++)
            {
               assert (res->get_size(tree_index) == this->get_size(t));
               res->state_trees[tree_index]=
                  new Typed_edge_one_int_tree(*(this->state_trees[t]));
               tree_index++;
            }
            for(i= 0; i < nb_sample; i++)
               for(t= 0; t < (unsigned int)otrees[i]->_nb_trees; t++)
               {
                  assert (res->get_size(tree_index) == otrees[i]->get_size(t));
                  res->state_trees[tree_index]=
                     new Typed_edge_one_int_tree(*(otrees[i]->state_trees[t]));
                  tree_index++;
               }
         }
         else
         {
            delete res;
            res= NULL;
         }
      }
   }
   delete [] hmtd_array;
   hmtd_array= NULL;

   if (res_t != NULL)
   {
      delete res_t;
      res_t= NULL;
   }

   return res;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovTreeData
 *  using an output stream and a flag on the level of detail
 *
 **/

ostream& HiddenMarkovTreeData::ascii_write(ostream &os, bool exhaustive) const

{
   Test *test= new Test(CHI2);

   if (markov != NULL)
   {
      markov->chi2_fit(*chain_data, *test);
      markov->ascii_write(os, this, exhaustive, false, test);
   }
   delete test;
   return os;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovTreeData into a file
 *  using a StatError object, the path
 *  and a flag on the level of detail
 *
 **/

bool HiddenMarkovTreeData::ascii_write(StatError& error,
                                       const char * path,
                                       bool exhaustive) const

{
   Test *test;
   ofstream out_file(path);

   bool status= false;

   if (markov != NULL)
   {
      error.init();

      if (!out_file)
      {
         status= false;
         error.update(STAT_error[STATR_FILE_NAME]);
      }

      else
      {
         status= true;

         test = new Test(CHI2);
         markov->chi2_fit(*chain_data , *test);

         markov->ascii_write(out_file, this, exhaustive, true, test);
         delete test;
      }
   }
   return status;
}

/*****************************************************************
 *
 *  Prints a HiddenMarkovTreeData in a spreadsheet fashion
 *  using a StatError object and the path
 *
 **/

bool HiddenMarkovTreeData::spreadsheet_write(StatError& error, const char * path) const
{
   bool status= false;
   Test *test;
   ofstream out_file(path);

   if (markov != NULL)
   {
      error.init();

      if (!out_file)
      {
         status= false;
         error.update(STAT_error[STATR_FILE_NAME]);
      }

      else
      {
         status= true;

         test= new Test(CHI2);
         markov->chi2_fit(*chain_data, *test);

         markov->spreadsheet_write(out_file, this, test);
         delete test;
      }
   }
   return status;
}

/*****************************************************************
 *
 *  Gnuplot output for HiddenMarkovTreeData class
 *  using a StatError object, a prefix for the files
 *  and the title of figures
 *
 **/

bool HiddenMarkovTreeData::plot_write(StatError& error,
                                      const char * prefix,
                                      const char * title) const
{
   bool status= false;

   if (markov != NULL)
   {
      status= markov->plot_write(prefix, title, this);

      error.init();

      if (!status)
         error.update(STAT_error[STATR_FILE_PREFIX]);
   }
   return status;
}

/*****************************************************************
 *
 *  Access to the class members for HiddenMarkovTreeData class
 *
 **/

HiddenMarkovTree* HiddenMarkovTreeData::get_markov() const
{
   HiddenMarkovTree *res= NULL;

   if (markov != NULL)
      res= markov->HiddenMarkovTreeCopy();

   return res;
}

ChainData* HiddenMarkovTreeData::get_chain_data() const
{
   ChainData *res= NULL;

   if (chain_data != NULL)
      res= new ChainData(*chain_data);

   return res;
}

double HiddenMarkovTreeData::get_likelihood() const
{ return likelihood; }

double HiddenMarkovTreeData::get_hidden_likelihood() const
{ return hidden_likelihood; }

int HiddenMarkovTreeData::get_nb_states() const
{ return _nb_states; }

TreeCharacteristics* HiddenMarkovTreeData::get_state_characteristics() const
{
   TreeCharacteristics *res;

   if (state_characteristics != NULL)
      res= new TreeCharacteristics(*state_characteristics);
   else
      res= NULL;

   return res;
}

HiddenMarkovTreeData::ptFrequencyDistribution_array_2d HiddenMarkovTreeData::get_observation() const
{  // frequency distributions corresponding to the conditional observation distributions
   // for each variable and state
   register int var, j;
   ptFrequencyDistribution_array_2d res= NULL;

   if (observation_distribution != NULL)
   {
      res= new ptFrequencyDistribution_array[_nb_integral];
      for(var= 0; var < _nb_integral; var++)
         if (observation_distribution[var] != NULL)
         {
            res[var]= new FrequencyDistribution*[state_characteristics->marginal_distribution->nb_value];
            for(j= 0; j < state_characteristics->marginal_distribution->nb_value; j++)
            {
               if (observation_distribution[var][j] != NULL)
                  res[var][j]= new FrequencyDistribution(*(observation_distribution[var][j]));
               else
                  res[var][j]= NULL;
            }
         }
         else
            res[var]= NULL;
   }
   else
      res= NULL;

   return res;
}

HiddenMarkovTreeData::ptFrequencyDistribution_array HiddenMarkovTreeData::get_observation(int variable) const
{  // frequency distributions corresponding to the conditional observation distributions
   // for each state and one given variable
   register int j;
   ptFrequencyDistribution_array res= NULL;

   if (observation_distribution != NULL)
   {
      assert(variable < _nb_integral);
      if (observation_distribution[variable] != NULL)
      {
         res= new FrequencyDistribution*[state_characteristics->marginal_distribution->nb_value];
         for(j= 0; j < state_characteristics->marginal_distribution->nb_value; j++)
         {
            if (observation_distribution[variable][j] != NULL)
               res[j]= new FrequencyDistribution(*(observation_distribution[variable][j]));
            else
               res[j]= NULL;
         }
      }
      else
         res= NULL;
   }
   else
      res= NULL;

   return res;
}

FrequencyDistribution* HiddenMarkovTreeData::get_observation(int variable,
                                                             int state) const
{  // frequency distributions corresponding to the conditional observation distributions
   // for one given state and one given variable
   // register int j;
   FrequencyDistribution *res= NULL;

   if (observation_distribution != NULL)
   {
      assert(variable < _nb_integral);
      if (observation_distribution[variable] != NULL)
      {
         assert(state < state_characteristics->marginal_distribution->nb_value);
         if (observation_distribution[variable][state] != NULL)
            res= new FrequencyDistribution(*(observation_distribution[variable][state]));
         else
            res= NULL;
      }
      else
         res= NULL;
   }
   else
      res= NULL;

   return res;
}

HiddenMarkovTreeData::ptOne_int_tree_set HiddenMarkovTreeData::get_state_trees() const
{
   int t;
   ptOne_int_tree_set res_trees= NULL;

   if (this->state_trees != NULL)
   {
      res_trees= new Typed_edge_one_int_tree*[_nb_trees];
      for(t= 0; t < _nb_trees; t++)
      {
         if (this->state_trees[t] != NULL)
            res_trees[t]= new Typed_edge_one_int_tree(*(this->state_trees[t]));
         else
            res_trees[t]= NULL;
      }
   }

   return res_trees;
}

HiddenMarkovTreeData::ptOne_int_tree_set HiddenMarkovTreeData::get_state_trees_ptr() const
{ return state_trees; }

Typed_edge_one_int_tree* HiddenMarkovTreeData::get_state_tree(int itree) const
{
   Typed_edge_one_int_tree *res_tree;

   assert(itree < _nb_trees);

   if (this->state_trees != NULL)
   {
      if (this->state_trees[itree] != NULL)
         res_tree= new Typed_edge_one_int_tree(*(this->state_trees[itree]));
      else
         res_tree= NULL;
   }
   else
      res_tree= NULL;

   return res_tree;
}

Typed_edge_one_int_tree* HiddenMarkovTreeData::get_state_tree_ptr(int itree) const
{
   assert((itree < _nb_trees) && (this->state_trees != NULL));
   return state_trees[itree];
}

/*****************************************************************
 *
 *  Return a HiddenMarkovTreeData containing the states
 *  as a variable
 *
 **/

HiddenMarkovTreeData*
HiddenMarkovTreeData::get_state_hidden_markov_tree_data() const
{
   register int offset;
   const short int added_int_variables= 1;
   HiddenMarkovTreeData *res= NULL;
   Unlabelled_typed_edge_tree *tmp_utree= NULL;
   bool return_this= false;
   int inb_variables, inb_trees, var, t;
   int *itype= NULL;
   // register int st;
   // Int_fl_container i(_nb_integral+added_int_variables, 0), s;
   Int_fl_container i(_nb_integral+added_int_variables, _nb_float), s;
   vertex_iterator it, end;
   Default_tree** otrees= NULL;

   if ((this->state_trees == NULL) || (markov == NULL))
   {
      return_this= true;
      res= new HiddenMarkovTreeData(*this, true);
   }
   else
   {
      inb_trees= _nb_trees;
      inb_variables= _nb_integral + 0 + added_int_variables;
      itype= new int[inb_variables];

      for(var=0; var < added_int_variables; var++)
         itype[var]= STATE;

      offset= added_int_variables;
      for(var=0; var < _nb_integral; var++)
         itype[var+offset]= _type[var];

      offset+= _nb_integral;
      for(var=0; var < _nb_float; var++)
         itype[var+offset]= _type[var+_nb_integral];

      otrees= new Default_tree*[inb_trees];
      for(t= 0; t < inb_trees; t++)
      {
         tmp_utree= trees[t]->get_structure();
         otrees[t]= new Default_tree(_nb_integral+added_int_variables,
                                     _nb_float+0,
                                     trees[t]->root(), 1);
         otrees[t]->set_structure(*tmp_utree, i);
         Tree_tie::tie(it, end)= trees[t]->vertices();
         while (it < end)
         {
            s= trees[t]->get(*it);
            // add integer variables
            i.Int(0)= (this->state_trees[t]->get(*it)).Int();
            // copy existing integer variables
            for(var= 0; var < _nb_integral; var++)
               i.Int(var+added_int_variables)=s.Int(var);
            // copy existing floating variables
            for(var= 0; var < _nb_float; var++)
               i.Double(var+0)=s.Double(var);
            otrees[t]->put(*it++, i);
         }
         delete tmp_utree;
         tmp_utree= NULL;
      }

      res= new HiddenMarkovTreeData(inb_trees, itype, otrees);
      res->markov= markov->HiddenMarkovTreeCopy(false, false);

      res->state_trees= new Typed_edge_one_int_tree*[inb_trees];
      for(t= 0; t < inb_trees; t++)
         res->state_trees[t]= new Typed_edge_one_int_tree(*(this->state_trees[t]));
      res->likelihood= likelihood;
      res->hidden_likelihood= hidden_likelihood;
      res->_nb_states= _nb_states;

      // res->chain_data= new ChainData(*res, 0, 1, markov);
      res->chain_data= new ChainData(markov->type, markov->nb_state,
                                     markov->nb_state);
      res->build_characteristics();
      res->build_size_frequency_distribution();
      res->build_nb_children_frequency_distribution();
      res->build_observation_frequency_distribution();

      res->markov->characteristic_computation(*res, true);

      delete [] itype;
      itype= NULL;
      for(t= 0; t < inb_trees; t++)
      {
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete [] otrees;
      otrees= NULL;
   }
   return res;
}

/*****************************************************************
 *
 *  Return a HiddenMarkovTreeData containing the states
 *  (restored by the specified algorithm),
 *  the smoothed probabilities and entropy as variables,
 *  given a tree index and a type of entropy profile
 *
 **/

HiddenMarkovTreeData*
HiddenMarkovTreeData::get_state_smoothed_hidden_markov_tree_data(int index,
                                                                 int algorithm,
                                                                 int entropy_algo) const
{
   typedef HiddenMarkovTree::double_array_3d double_array_3d;
   typedef HiddenMarkovTree::double_array_2d double_array_2d;

   const int added_int_variables= 1;
   register int offset, added_float_variables= 0, st, j, cstate,
                added_variables= added_int_variables+added_float_variables;
   bool return_this= false;
   int inb_variables, inb_trees, var, t, ti; // index of current tree
   double pmax;
   Int_fl_container i(_nb_integral+added_int_variables,
                      _nb_float+added_float_variables),
                    s;
   state_value val;
   HiddenMarkovTreeData *res= NULL, *this_cp= NULL; // copy of this
   const HiddenMarkovTree *hmarkovt= NULL;
   Unlabelled_typed_edge_tree *tmp_utree= NULL;
   // key v;
   int *itype= NULL;
   vertex_iterator it, end;
   Default_tree** otrees= NULL;
   double computed_likelihood= D_INF;
   double_array_3d smoothed_prob= NULL;
   double_array_2d marginal_entropy= NULL, conditional_entropy= NULL,
                   partial_entropy= NULL;

   assert(index < get_nb_trees());

   hmarkovt= this->markov;
   if ((this->state_trees == NULL) || (hmarkovt == NULL))
   {
      return_this= true;
      res= new HiddenMarkovTreeData(*this, true);
   }
   else
   {
      added_float_variables= hmarkovt->nb_state + 3;
      added_variables= added_int_variables+added_float_variables;
      // determine dimension of result
      i.reset(_nb_integral+added_int_variables,
              _nb_float+added_float_variables);
      if (index == I_DEFAULT)
         inb_trees= _nb_trees;
      else
         inb_trees= 1;
      inb_variables= _nb_integral + _nb_float + added_variables;
      itype= new int[inb_variables];

      // optimal state
      for(var=0; var < added_int_variables; var++)
         itype[var]= STATE;

      // copy existing integer variables
      offset= added_int_variables;
      for(var=0; var < _nb_integral; var++)
         itype[var+offset]= _type[var];

      // smoothed probabilities and entropy
      offset+= _nb_integral;
      for(var=0; var < added_float_variables; var++)
         itype[var+offset]= REAL_VALUE;

      // copy existing floating variables
      offset= added_float_variables;
      for(var=0; var < _nb_float; var++)
         itype[var+offset]= _type[var+_nb_integral];

      // compute the smoothed probabilities
      computed_likelihood= hmarkovt->smoothed_probabilities(*this,
                                                            smoothed_prob,
                                                            marginal_entropy,
                                                            conditional_entropy,
                                                            partial_entropy,
                                                            index,
                                                            entropy_algo);

      if ((smoothed_prob != NULL) && (computed_likelihood != D_INF))
      {
         otrees= new Default_tree*[inb_trees];
         for(t= 0; t < _nb_trees; t++)
            if ((index == I_DEFAULT) || (t == index))
            {
               // if only one tree is involved, build only otrees[0]
               if (index == I_DEFAULT)
                  ti= t;
               else
                  ti= 0;
               tmp_utree= trees[t]->get_structure();
               otrees[ti]= new Default_tree(_nb_integral+added_int_variables,
                                            _nb_float+added_float_variables,
                                            trees[t]->root(), 1);
               // set topology and dimension of the result
               otrees[ti]->set_structure(*tmp_utree, i);
               Tree_tie::tie(it, end)= trees[t]->vertices();
               while (it < end)
               {
                  s= trees[t]->get(*it);
                  // add integer variables
                  // add restored state (smoothing)
                  i.Int(0)= (this->state_trees[t]->get(*it)).Int();
                  // copy existing integer variables
                  for(var= 0; var < _nb_integral; var++)
                     i.Int(var+added_int_variables)=s.Int(var);
                  // add floating variables
                  for(st= 0; st < added_float_variables-3; st++)
                     i.Double(st)= smoothed_prob[t][st][*it];
                  i.Double(st)= conditional_entropy[t][*it]; st++;
                  i.Double(st)= marginal_entropy[t][*it]; st++;
                  i.Double(st)= partial_entropy[t][*it];
                  // copy existing floating variables
                  for(var= 0; var < _nb_float; var++)
                     i.Double(var+added_float_variables)=s.Double(var);
                  otrees[ti]->put(*it++, i);
               }
               delete tmp_utree;
               tmp_utree= NULL;
            }

         res = new HiddenMarkovTreeData(inb_trees, itype, otrees);
         res->markov = hmarkovt->HiddenMarkovTreeCopy(false, false);
         // new HiddenMarkovTree(*markov, false, false);

         res->state_trees= new Typed_edge_one_int_tree*[inb_trees];
         for(t= 0; t < _nb_trees; t++)
            if ((index == I_DEFAULT) || (t == index))
            {
               if (index == I_DEFAULT)
                  ti= t;
               else
                  ti= 0;
               res->state_trees[ti]=
                  new Typed_edge_one_int_tree(*(this->state_trees[t]));
               // default state tree restoration : viterbi algorithm
               // restored state tree must be updated if algorithm == FORWARD_BACKWARD
               if (algorithm == FORWARD_BACKWARD)
               {
                  Tree_tie::tie(it, end)= res->state_trees[ti]->vertices();
                  while (it < end)
                  {
                     pmax= 0.;
                     val= res->state_trees[ti]->get(*it);
                     for(j= 0; j < hmarkovt->nb_state; j++)
                        if (smoothed_prob[t][j][*it] > pmax)
                        {
                           pmax= smoothed_prob[t][j][*it];
                           cstate= j;
                        }
                     val.Int()= cstate;
                     res->state_trees[ti]->put(*it, val);
                     i= res->trees[ti]->get(*it);
                     i.Int(0)= cstate;
                     res->trees[ti]->put(*it++, i);
                  }
               }
            }
         // res->markov->likelihood_computation(*res) impossible
         // since the number of variables of res has been increased
         res->likelihood = likelihood;
         res->sample_entropy = sample_entropy;
         if (entropy != NULL)
         {
            res->entropy = new double[_nb_trees];
            for(t = 0; t < _nb_trees; t++)
               res->entropy[t] = entropy[t];
         }


         if (algorithm == FORWARD_BACKWARD)
         {
            this_cp= new HiddenMarkovTreeData(*this, false);
            for(t= 0; t < _nb_trees; t++)
               if (this_cp->state_trees[t] != NULL)
               {
                  delete this_cp->state_trees[t];
                  this_cp->state_trees[t]= NULL;
               }
               if (this_cp->state_trees != NULL)
            delete [] this_cp->state_trees;
            this_cp->state_trees= new Typed_edge_one_int_tree*[_nb_trees];
            for(t= 0; t < _nb_trees; t++)
               if ((index == I_DEFAULT) || (t == index))
               {
                  if (index == I_DEFAULT)
                     ti= t;
                  else
                     ti= 0;
                  this_cp->state_trees[t]=
                     new Typed_edge_one_int_tree(*res->state_trees[ti]);
               }
               else
                  this_cp->state_trees[t]= NULL;
            res->hidden_likelihood= hmarkovt->state_likelihood_computation(*this_cp,
                                                                           index);
            delete this_cp;
            this_cp= NULL;
         }
         else
            res->hidden_likelihood= hidden_likelihood;
         res->_nb_states= _nb_states;

         // res->chain_data= new ChainData(*res, 0, 1, markov);
         res->chain_data= new ChainData(markov->type, markov->nb_state,
                                        markov->nb_state);
         res->build_characteristics();
         res->build_size_frequency_distribution();
         res->build_nb_children_frequency_distribution();
         res->build_observation_frequency_distribution();

         res->markov->characteristic_computation(*res, true);

         delete [] itype;
         itype= NULL;
         for(t= 0; t < _nb_trees; t++)
            if ((index == I_DEFAULT) || (t == index))
            {
               if (index == I_DEFAULT)
                  ti= t;
               else
                  ti= 0;
               delete otrees[ti];
               otrees[ti]= NULL;
            }
         delete [] otrees;
         otrees= NULL;

         for(t= 0; t < _nb_trees; t++)
            if ((index == I_DEFAULT) || (t == index))
            {
               for(st= 0; st < hmarkovt->nb_state; st++)
               {
                  delete [] smoothed_prob[t][st];
                  smoothed_prob[t][st]= NULL;
               }
               delete [] smoothed_prob[t];
               smoothed_prob[t]= NULL;
               delete [] conditional_entropy[t];
               conditional_entropy[t]= NULL;
               delete [] marginal_entropy[t];
               marginal_entropy[t]= NULL;
               delete [] partial_entropy[t];
               partial_entropy[t]= NULL;
            }
         if (smoothed_prob != NULL)
         {
            delete [] smoothed_prob;
            smoothed_prob= NULL;
         }
         if (conditional_entropy != NULL)
         {
            delete [] conditional_entropy;
            conditional_entropy= NULL;
         }
         if (marginal_entropy != NULL)
         {
            delete [] marginal_entropy;
            marginal_entropy= NULL;
         }
         if (partial_entropy != NULL)
         {
            delete [] partial_entropy;
            partial_entropy= NULL;
         }
      }
   }

   if ((computed_likelihood == D_INF) && !(return_this))
   {
      delete res;
      res= new HiddenMarkovTreeData(*this, true);
   }
   return res;
}

std::ostream& HiddenMarkovTreeData::ascii_write_observation(std::ostream &os,
                                                            bool exhaustive,
                                                            bool file_flag) const
{
   int value, var;

   if (observation_distribution != NULL)
      for(var= 0; var < _nb_integral; var++)
         if (observation_distribution[var] != NULL)
            for(value= 0; value < _nb_states; value++)
               if (observation_distribution[var][value] != NULL)
               {
                  os << "Observation frequency distribution for state " << value
                     << " and variable " << var << " : " << endl;
                  observation_distribution[var][value]->ascii_write(os, exhaustive, file_flag);
                  os << endl;
               }
   return os;
}

void HiddenMarkovTreeData::nb_state_computation()
{  // computation of the number of states
   int t, res= 0;
   Typed_edge_one_int_tree::vertex_iterator it, end;
   value v;

   for(t= 0; t < _nb_trees; t++)
   {
        Tree_tie::tie(it, end)= this->state_trees[t]->vertices();
        while (it < end)
           res= max(res, (this->state_trees[t]->get(*it++)).Int());
   }
   _nb_states= res+1;
}

/*****************************************************************
 *
 *  Computation of the frequency distributions corresponding to the conditional
 *  distribution of one observed (given) variable
 *  of HiddenMarkovTreeData
 *
 **/

void HiddenMarkovTreeData::observation_frequency_distribution_computation(int ivariable)
{
   register int t, s;
   int val;
   Typed_edge_one_int_tree::vertex_iterator it, end;

   assert(ivariable < _nb_integral);

   for(t= 0; t < _nb_trees; t++)
      if (this->state_trees[t] != NULL)
      {
         Tree_tie::tie(it, end)= this->state_trees[t]->vertices();
         while (it < end)
         {
            s= (trees[t]->get(*it)).Int(ivariable);
            s= (this->state_trees[t]->get(*it)).Int();
            val= (trees[t]->get(*it)).Int(ivariable);
            observation_distribution[ivariable][s]->frequency[val]++;
            it++;
         }
      }

   for(s= 0; s < state_characteristics->marginal_distribution->nb_value; s++)
   {
      observation_distribution[ivariable][s]->nb_value_computation();
      observation_distribution[ivariable][s]->offset_computation();
      observation_distribution[ivariable][s]->nb_element_computation();
      observation_distribution[ivariable][s]->max_computation();
      observation_distribution[ivariable][s]->mean_computation();
      observation_distribution[ivariable][s]->variance_computation();
   }
}

/*****************************************************************
 *
 *  Allocation of the frequency distributions corresponding to the observation
 *  (conditional) distributions for HiddenMarkovTreeData class
 *
 **/

void HiddenMarkovTreeData::create_observation_frequency_distribution(int nb_state)
{
   register int var, j;

   if (_nb_integral > 0)
   {
      if (observation_distribution == NULL)
      {
         observation_distribution= new ptFrequencyDistribution_array[_nb_integral];
         for(var= 0; var < _nb_integral; var++)
         {
            observation_distribution[var]= new FrequencyDistribution*[nb_state];
            for(j= 0; j < nb_state; j++)
               observation_distribution[var][j]= new FrequencyDistribution(get_max_int_value(var)+1);
         }
      }
      else // the number of values may have changed
         for(var= 0; var < _nb_integral; var++)
         {
            if (observation_distribution[var] != NULL)
               for(j= 0; j < nb_state; j++)
               {
                  if (observation_distribution[var][j] != NULL)
                  {
                     delete observation_distribution[var][j];
                     observation_distribution[var][j]= NULL;
                  }
                  observation_distribution[var][j]= new FrequencyDistribution(get_max_int_value(var)+1);
               }
            else
            {
               observation_distribution[var]= new FrequencyDistribution*[nb_state];
               for(j= 0; j < nb_state; j++)
                  observation_distribution[var][j]= new FrequencyDistribution(get_max_int_value(var)+1);
            }
         }
   }
   else // delete the potential frequency distributions
   {
      if (observation_distribution != NULL)
         for(var= 0; var < _nb_integral; var++)
         {
            if (observation_distribution[var] != NULL)
            {
               for(j= 0; j < nb_state; j++)
               {
                  if (observation_distribution[var][j] != NULL)
                  {
                     delete observation_distribution[var][j];
                     observation_distribution[var][j]= NULL;
                  }
               }
               delete [] observation_distribution[var];
               observation_distribution[var]= NULL;
            }
         }
      observation_distribution= NULL;
   }
}

/*****************************************************************
 *
 *  Computation of the frequency distributions corresponding to the observation
 *  (conditional) distributions for Hidden_markov_data_tree class
 *  - without prior allocation
 *
 **/

void HiddenMarkovTreeData::observation_frequency_distribution_computation()
{
   register int var;

   for(var= 0; var < _nb_integral; var++)
      observation_frequency_distribution_computation(var);
}

/*****************************************************************
 *
 *  Allocation and computation of the frequency distributions corresponding
 *  to the observation (conditional) distributions for
 *  HiddenMarkovTree data class
 *
 **/


void HiddenMarkovTreeData::build_observation_frequency_distribution()
{
   build_state_characteristics();
   create_observation_frequency_distribution(state_characteristics->marginal_distribution->nb_value);
   observation_frequency_distribution_computation();
}

/*****************************************************************
 *
 *  Allocation of the hidden state trees
 *  for HiddenMarkovTreeData class
 *
 **/

void HiddenMarkovTreeData::build_state_trees()
{
   typedef Typed_edge_one_int_tree::value value;

   register int t;
   Unlabelled_typed_edge_tree *utree;
   value default_value;

   default_value.Int()= 0; // I_DEFAULT;

   if (this->state_trees == NULL)
   {
      this->state_trees= new Typed_edge_one_int_tree*[_nb_trees];
      for(t= 0; t < _nb_trees; t++)
      {
         this->state_trees[t]= new Typed_edge_one_int_tree;
         utree= trees[t]->get_structure();
         this->state_trees[t]->set_structure(*utree, default_value);
         delete utree;
      }
   }
}

/*****************************************************************
 *
 * Computation of the characteristic quantity frequency distributions
 * of HiddenMarkovTreeData class for the (hidden) state variable
 *
 **/

void HiddenMarkovTreeData::build_state_characteristics()
{

   Typed_edge_one_int_tree **otrees1= new Typed_edge_one_int_tree*[_nb_trees];
   int t; // i

   if (_nb_states == I_DEFAULT)
      nb_state_computation();

   if (state_characteristics != NULL)
   {
      delete state_characteristics;
      state_characteristics= NULL;
   }
   for(t= 0; t < _nb_trees; t++)
      otrees1[t]= this->state_trees[t];

   state_characteristics= new TreeCharacteristics(0,
                                                  _nb_states-1,
                                                  _max_size,
                                                  _max_depth,
                                                  _nb_trees,
                                                  otrees1,
                                                  0);

   delete [] otrees1;
   otrees1= NULL;
}

/*****************************************************************
 *
 * Return a HiddenMarkovTreeData without the state variable
 * If the first variable corresponds to the state variable,
 * remove it. Otherwise, return a copy of this.
 *
 **/

HiddenMarkovTreeData* HiddenMarkovTreeData::remove_state_variable() const
{
   register int offset;
   const short int remove_int_variable = 1;
   int state;
   HiddenMarkovTreeData *res = NULL, *pre_res = NULL;
   Unlabelled_typed_edge_tree *tmp_utree = NULL;
   bool return_this = false;
   int inb_variables, inb_trees, var, t;
   int *itype = NULL;
   // register int st;
   Int_fl_container i(MAX(1,
                          _nb_integral-remove_int_variable),
                      _nb_float);
   Int_fl_container s;
   vertex_iterator it, end;
   Default_tree** otrees = NULL;

   if ((this->state_trees == NULL) || (markov == NULL))
   {
      return_this = true;
      // res = new HiddenMarkovTreeData(*this, true);
   }

   // check whether the first variable is the state variable
   if ((!return_this) && (this->_type[0] != STATE))
      return_this = true;

   if (!return_this)
   {
      inb_trees= _nb_trees;
      inb_variables = _nb_integral - remove_int_variable;
      itype = new int[inb_variables];

      offset = remove_int_variable;
      for(var=offset ; var < _nb_integral; var++)
         itype[var-offset] = _type[var];

      for(var=0; var < _nb_float; var++)
         itype[var+_nb_integral-offset] = _type[var+_nb_integral];

      otrees = new Default_tree*[inb_trees];
      for(t = 0; t < inb_trees; t++)
      {
         tmp_utree = trees[t]->get_structure();
         otrees[t] = new Default_tree(_nb_integral-remove_int_variable,
                                      _nb_float+0,
                                      trees[t]->root(), 1);
         otrees[t]->set_structure(*tmp_utree, i);
         Tree_tie::tie(it, end) = trees[t]->vertices();
         while (it < end)
         {
            s = trees[t]->get(*it);
            // remove state variable
            state = (this->state_trees[t]->get(*it)).Int();
            if (s.Int(0) != state)
               return_this = true;
            // copy existing integer variables
            for(var = offset; var < _nb_integral; var++)
               i.Int(var-offset) = s.Int(var);
            // copy existing floating variables
            for(var=0; var < _nb_float; var++)
               i.Double(var+0) = s.Double(var);
            otrees[t]->put(*it++, i);
         }
         delete tmp_utree;
         tmp_utree = NULL;
      }

      pre_res = new HiddenMarkovTreeData(inb_trees, itype, otrees);
      pre_res->markov = markov->HiddenMarkovTreeCopy(false, false);

      pre_res->state_trees = new Typed_edge_one_int_tree*[inb_trees];
      for(t=0; t < inb_trees; t++)
         pre_res->state_trees[t] =
            new Typed_edge_one_int_tree(*(this->state_trees[t]));
      pre_res->likelihood = likelihood;
      pre_res->hidden_likelihood = hidden_likelihood;
      pre_res->_nb_states = _nb_states;

      pre_res->chain_data= new ChainData(markov->type, markov->nb_state,
                                         markov->nb_state);
      pre_res->build_characteristics();
      pre_res->build_size_frequency_distribution();
      pre_res->build_nb_children_frequency_distribution();
      pre_res->build_observation_frequency_distribution();

      pre_res->markov->characteristic_computation(*pre_res, true);

      if (itype != NULL)
      {
         delete [] itype;
         itype = NULL;
         for(t=0; t < inb_trees; t++)
         {
            delete otrees[t];
            otrees[t] = NULL;
         }
         delete [] otrees;
         otrees = NULL;
      }
   }

   if (return_this)
      res = new HiddenMarkovTreeData(*this, true);
   else
   {
      res = new HiddenMarkovTreeData(*pre_res, true);
      delete pre_res;
      pre_res = NULL;
   }
   return res;
}

/*****************************************************************
 *
 *  Permutation of the states of a HiddenMarkovTreeData
 *  based on a given permutation perm.
 *  Validity of permutation must be ensured before call
 *
 **/

void HiddenMarkovTreeData::state_permutation(int* perm)
{
   register int i, var, t;
   TreeCharacteristics* pstate_char= new TreeCharacteristics[_nb_states];
   ptFrequencyDistribution_array_2d pobservation= new FrequencyDistribution**[_nb_integral];
   Typed_edge_one_int_tree::vertex_iterator it, end;
   One_int_container c;

   for(var= 0; var < _nb_integral; var++)
   {
      pobservation[var]= new FrequencyDistribution*[_nb_states];
      for (i= 0; i < _nb_states; i++)
         pobservation[var][perm[i]]= observation_distribution[var][i];
      for (i= 0; i < _nb_states; i++)
         observation_distribution[var][i]= pobservation[var][i];
      delete [] pobservation[var];
      pobservation[var]= NULL;
   }

   delete [] pstate_char;
   pstate_char= NULL;
   delete [] pobservation;
   pobservation= NULL;

   if (state_trees != NULL)
   {
      for(t= 0; t < _nb_trees; t++)
      {
         Tree_tie::tie(it, end)= state_trees[t]->vertices();
         while (it < end)
         {
            c= state_trees[t]->get(*it);
            c.Int()= perm[c.Int()];
            state_trees[t]->put(*it++, c);
         }
      }
   }
   build_state_characteristics();
}

std::ostream& HiddenMarkovTreeData::state_profile_ascii_print(std::ostream& os,
                                                              int index,
                                                              int nb_state,
                                                              double_array_2d smoothed) const
{
   return os;
}

std::ostream&
HiddenMarkovTreeData::state_profile_spreadsheet_print(std::ostream& os,
                                                      int index,
                                                      int nb_state,
                                                      double_array_2d smoothed) const
{
   return os;
}

std::ostream&
HiddenMarkovTreeData::state_profile_plot_print(std::ostream& os, int index,
                                               int nb_state,
                                               double_array_2d smoothed) const
{
   return os;
}

/*****************************************************************
 *
 *  Write state and entroy profiles (Gnuplot output)
 *  using a given output stream, the index of considered tree,
 *  the number of states of the model, the profiles
 *  the vertex defining the tree path to consider and the path itself
 *
 **/

std::ostream&
HiddenMarkovTreeData::profile_plot_print(std::ostream& os,
                                         int index,
                                         int nb_state,
                                         double_array_2d smoothed,
                                         double_array conditional_entropy,
                                         double_array marginal_entropy,
                                         double_array partial_entropy,
                                         key vertex,
                                         generic_visitor<tree_type>::vertex_deque*& vd) const
{  // validity of index and vertex must be checked before
   register unsigned int i, j;
   tree_type current_tree= *trees[index];
   generic_visitor<tree_type> *visitor= NULL;


   if (vd == NULL)
   {
      visitor= new generic_visitor<tree_type>;
      traverse_tree(current_tree.root(), current_tree, *visitor);
      vd= new generic_visitor<tree_type>::vertex_deque();
      *vd= visitor->get_vertex_ancestors(current_tree, vertex);
      delete visitor;
      visitor= NULL;
   }

   for(i= 0; i < vd->size(); i++)
   {
      for(j= 0; j < (unsigned int) nb_state; j++)
         os << smoothed[j][(*vd)[i]] << " ";
      os << conditional_entropy[(*vd)[i]] << " " << marginal_entropy[(*vd)[i]]
         << " " << partial_entropy[(*vd)[i]] << endl;
   }

   return os;
}

/*****************************************************************
 *
 *  Write Viterbi state profile (Gnuplot output)
 *  using a given output stream, the index of considered tree,
 *  the number of states of the model, the state profile
 *  the vertex defining the tree path to consider and the path itself
 *
 **/

std::ostream&
HiddenMarkovTreeData::profile_plot_print(std::ostream& os,
                                         int index,
                                         int nb_state,
                                         double_array_2d ratio,
                                         key vertex,
                                         generic_visitor<tree_type>::vertex_deque*& vd) const
{  // validity of index and vertex must be checked before
   typedef generic_visitor<tree_type>::vertex_deque vertex_deque;

   register unsigned int i, j;
   tree_type current_tree= *trees[index];
   generic_visitor<tree_type> *visitor= NULL;

   if (vd == NULL)
   {
      visitor= new generic_visitor<tree_type>;
      traverse_tree(current_tree.root(), current_tree, *visitor);
      vd= new generic_visitor<tree_type>::vertex_deque();
      *vd= visitor->get_vertex_ancestors(current_tree, vertex);
      delete visitor;
      visitor= NULL;
   }

   for(i= 0; i < vd->size(); i++)
   {
      for(j= 0; j < (unsigned int) nb_state; j++)
         os << ratio[j][(*vd)[i]] << " ";
      os << endl;
   }

   return os;
}

/*****************************************************************
 *
 *  Deallocation of the pointers for HiddenMarkovTreeData class
 *  The markov part is deleted if any
 *
 **/

void HiddenMarkovTreeData::remove()
{
   register int t, var, j;

   if (markov != NULL)
   {
      delete markov;
      markov = NULL;
   }

   if (chain_data != NULL)
   {
      delete chain_data;
      chain_data = NULL;
   }

   if (this->entropy != NULL)
   {
      delete [] this->entropy;
      this->entropy = NULL;
   }

   if (this->state_trees != NULL)
   {
       for(t = 0; t < _nb_trees; t++)
       {
          if (this->state_trees[t] != NULL)
          {
             delete this->state_trees[t];
             this->state_trees[t] = NULL;
          }
       }
       delete [] this->state_trees;
       this->state_trees = NULL;
   }

   if (observation_distribution != NULL)
   {
      for(var = 0; var < _nb_integral; var++)
      {
         if (observation_distribution[var] != NULL)
         {
            for(j = 0; j < state_characteristics->marginal_distribution->nb_value; j++)
            {
               if (observation_distribution[var][j] != NULL)
               {
                  delete observation_distribution[var][j];
                  observation_distribution[var][j] = NULL;
               }
            }
            delete [] observation_distribution[var];
         }
         observation_distribution[var] = NULL;
      }
      delete [] observation_distribution;
      observation_distribution = NULL;
   }

   if (observation_histogram != NULL)
   {
      for(var = 0; var < _nb_float; var++)
      {
         if (observation_histogram[var] != NULL)
         {
            for(j = 0; j < state_characteristics->marginal_distribution->nb_value; j++)
            {
               if (observation_histogram[var][j] != NULL)
               {
                  delete observation_histogram[var][j];
                  observation_histogram[var][j] = NULL;
               }
            }
            delete [] observation_histogram[var];
         }
         observation_histogram[var] = NULL;
      }
      delete [] observation_histogram;
      observation_histogram = NULL;
   }

   if (state_characteristics != NULL)
   {
      delete state_characteristics;
      state_characteristics = NULL;
   }

   // delete observation_histogram?

   // A posterior call to Trees::remove() is required...
}

/*****************************************************************
 *
 *  Copy operator of HiddenMarkovTreeData class
 *  using a flag on the HiddenMarkovTree copy
 *
 **/

void HiddenMarkovTreeData::copy(const HiddenMarkovTreeData& trees,
                                bool model_flag,
                                bool characteristic_flag)
{
   register int t, var, j;

   if ((model_flag) && (trees.markov != NULL))
   {
      if (markov != NULL)
         delete markov;
      markov = trees.markov->HiddenMarkovTreeCopy(false);
   }
   else
      markov = NULL;

   if (trees.chain_data != NULL)
   {
      if (chain_data != NULL)
         delete chain_data;
      chain_data = new ChainData(*(trees.chain_data));
   }
   else
      chain_data = NULL;

   likelihood = trees.likelihood;
   hidden_likelihood = trees.hidden_likelihood;
   sample_entropy = trees.sample_entropy;

   if (trees.entropy != NULL)
   {
      this->entropy = new double[_nb_trees];
      for(t = 0; t < _nb_trees; t++)
      {
         this->entropy[t] = trees.entropy[t];
      }
   }
   else
      this->entropy = NULL;


   if (trees.state_trees != NULL)
   {
      this->state_trees = new Typed_edge_one_int_tree*[_nb_trees];
      for(t = 0; t < _nb_trees; t++)
      {
          if (trees.state_trees[t] != NULL)
             this->state_trees[t] =
                new Typed_edge_one_int_tree(*(trees.state_trees[t]));
          else
             this->state_trees[t] = NULL;
      }
   }
   else
      this->state_trees= NULL;

   if (characteristic_flag)
   {
      if (trees.state_characteristics != NULL)
         state_characteristics = new TreeCharacteristics(*(trees.state_characteristics));
      else
         state_characteristics = NULL;

      if (trees.observation_distribution != NULL)
      {
         observation_distribution= new ptFrequencyDistribution_array[_nb_integral];
         for(var = 0; var < _nb_integral; var++)
            if (trees.observation_distribution[var] != NULL)
            {
               observation_distribution[var]= new FrequencyDistribution*[state_characteristics->marginal_distribution->nb_value];
               for(j = 0; j < state_characteristics->marginal_distribution->nb_value; j++)
               {
                  if (trees.observation_distribution[var][j] != NULL)
                     observation_distribution[var][j] = new FrequencyDistribution(*(trees.observation_distribution[var][j]));
                  else
                     observation_distribution[var][j] = NULL;
               }
            }
            else
               observation_distribution[var] = NULL;
      }
      else
         observation_distribution = NULL;

      if (trees.observation_histogram != NULL)
      {
         observation_histogram = new ptHistogram_array[_nb_integral];
         for(var = 0; var < _nb_integral; var++)
            if (trees.observation_histogram[var] != NULL)
            {
               observation_histogram[var]= new Histogram*[state_characteristics->marginal_distribution->nb_value];
               for(j = 0; j < state_characteristics->marginal_distribution->nb_value; j++)
               {
                  if (trees.observation_histogram[var][j] != NULL)
                     observation_histogram[var][j] = new Histogram(*(trees.observation_distribution[var][j]));
                  else
                     observation_histogram[var][j] = NULL;
               }
            }
            else
               observation_histogram[var] = NULL;
      }
      else
         observation_histogram = NULL;
   }
   else
   {
      if (state_characteristics != NULL)
         delete state_characteristics;
      if (observation_distribution != NULL)
         delete observation_distribution;
      state_characteristics = NULL;
      observation_distribution = NULL;
   }

}

/*****************************************************************
 *
 * Left (bit) shift operator for HiddenMarkovTreeData class
 *
 **/

ostream& Stat_trees::operator<<(ostream &os, const HiddenMarkovTreeData& trees)
{ return trees.ascii_write(os, false); }
