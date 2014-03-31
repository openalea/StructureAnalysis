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
 *       $Id: typed_edge_trees.cpp 3193 2007-05-29 10:03:19Z dufourko $
 *
 *       Forum for OpenAlea developers: Openalea-devlp@lists.gforge.inria.fr
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

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"    // definition of Curves class for Sequences
#include "stat_tool/markovian.h" // definition of constants
#include "stat_tool/distribution.h"
#include "stat_tool/vectors.h"
#include "sequence_analysis/sequences.h" // definition of constants
#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "typed_edge_trees.h"
#include <iostream>
#include <cstdlib>

using namespace Stat_trees;

/*****************************************************************
 *
 *  Default constructor of TreeCharacteristics class
 **/

TreeCharacteristics::TreeCharacteristics()
  :  _variable(I_DEFAULT)
  ,  _min_value(0)
  ,  _max_value(0)
  ,  marginal_distribution(NULL)
  ,  marginal_histogram(NULL)
  ,  first_occurrence_root(NULL)
  ,  first_occurrence_leaves(NULL)
  ,  sojourn_size(NULL)
  ,  nb_zones(NULL)
  ,  nb_occurrences(NULL)
  // , children_pairs(0)
{}

/*****************************************************************
 *
 *  Constructor of TreeCharacteristics class using observations
 *  using a flag on considering value I_DEFAULT as a missing value
 *  and ignore this value
 *
 **/

TreeCharacteristics::TreeCharacteristics(int imin_value,
                                         int imax_value,
                                         int imax_size,
                                         int imax_depth,
                                         int inb_trees,
                                         Typed_edge_one_int_tree** otrees1,
                                         int ivariable,
                                         bool final_value)
  :  _variable(ivariable)
  ,  _min_value(imin_value)
  ,  _max_value(imax_value)
  ,  marginal_distribution(NULL)
  ,  marginal_histogram(NULL)
  ,  first_occurrence_root(NULL)
  ,  first_occurrence_leaves(NULL)
  ,  sojourn_size(NULL)
  ,  nb_zones(NULL)
  ,  nb_occurrences(NULL)
  // , children_pairs(0)
{  build_characteristics(_min_value,
                         _max_value,
                         imax_size,
                         imax_depth,
                         inb_trees,
                         otrees1,
                         _variable,
                         final_value); }

/*****************************************************************
 *
 *  Copy constructor of TreeCharacteristics class
 *
 **/

TreeCharacteristics::TreeCharacteristics(const TreeCharacteristics& t_char)
 :  _variable(t_char._variable)
 ,  _min_value(t_char._min_value)
 ,  _max_value(t_char._max_value)
 ,  marginal_distribution(NULL)
 ,  marginal_histogram(NULL)
 ,  first_occurrence_root(NULL)
 ,  first_occurrence_leaves(NULL)
 ,  sojourn_size(NULL)
 ,  nb_zones(NULL)
 ,  nb_occurrences(NULL)

{ copy(t_char); }

/*****************************************************************
 *
 *  Destructor of TreeCharacteristics class
 *
 **/

TreeCharacteristics::~TreeCharacteristics()
{ remove(); }

/*****************************************************************
 *
 *  Copy operator of TreeCharacteristics class
 *
 **/

void TreeCharacteristics::copy(const TreeCharacteristics& t_char)
{
   int value, nb_values;

   _variable= t_char._variable;
   _min_value= t_char._min_value;
   _max_value= t_char._max_value;

   nb_values= _max_value - _min_value + 1;

   if (t_char.marginal_distribution != NULL)
      marginal_distribution= new FrequencyDistribution(*(t_char.marginal_distribution));
   else
      marginal_distribution= NULL;

   if (t_char.marginal_histogram != NULL)
      marginal_histogram= new Histogram(*(t_char.marginal_histogram));
   else
      marginal_histogram= NULL;

   if (t_char.first_occurrence_root != NULL)
   {
      first_occurrence_root= new FrequencyDistribution*[nb_values];
      for(value= 0; value < nb_values; value++)
      {
          if (t_char.first_occurrence_root[value] != NULL)
             first_occurrence_root[value]= new FrequencyDistribution(*(t_char.first_occurrence_root[value]));
          else
             first_occurrence_root[value]= NULL;
      }
   }
   else
      first_occurrence_root= NULL;

   if (t_char.first_occurrence_leaves != NULL)
   {
       first_occurrence_leaves= new FrequencyDistribution*[nb_values];
       for(value= 0; value < nb_values; value++)
       {
           if (t_char.first_occurrence_leaves[value] != NULL)
              first_occurrence_leaves[value]= new FrequencyDistribution(*(t_char.first_occurrence_leaves[value]));
           else
              first_occurrence_leaves[value]= NULL;
       }
   }
   else
      first_occurrence_leaves= NULL;

   if (t_char.sojourn_size != NULL)
   {
       sojourn_size= new FrequencyDistribution*[nb_values];
       for(value= 0; value < nb_values; value++)
       {
           if (t_char.sojourn_size[value] != NULL)
              sojourn_size[value]= new FrequencyDistribution(*(t_char.sojourn_size[value]));
           else
              sojourn_size[value]= NULL;
       }
   }
   else
      sojourn_size= NULL;

   if (t_char.nb_zones != NULL)
   {
       nb_zones= new FrequencyDistribution*[nb_values];
       for(value= 0; value < nb_values; value++)
       {
           if (t_char.nb_zones[value] != NULL)
              nb_zones[value]= new FrequencyDistribution(*(t_char.nb_zones[value]));
           else
              nb_zones[value]= NULL;
       }
   }
   else
      nb_zones= NULL;

   if (t_char.nb_occurrences != NULL)
   {
       nb_occurrences= new FrequencyDistribution*[nb_values];
       for(value= 0; value < nb_values; value++)
       {
           if (t_char.nb_occurrences[value] != NULL)
              nb_occurrences[value]= new FrequencyDistribution(*(t_char.nb_occurrences[value]));
           else
              nb_occurrences[value]= NULL;
       }
   }
   else
      nb_occurrences= NULL;
}

/*****************************************************************
 *
 *  Assignement operator of TreeCharacteristics class
 *
 **/

TreeCharacteristics&
TreeCharacteristics::operator=(const TreeCharacteristics &t_char)
{
   if (&t_char != this)
   {
      remove();
      copy(t_char);
   }
   return *this;
}


/*****************************************************************
 *
 *  Extract the characteristic quantities from a sample
 *  of trees and build frequency distributions:
 *     path length before the first occurrence of a given value:
 *        - starting from the root
 *        - starting from all the terminal vertices
 *     size of the homogeneous zones for a given value
 *     number of homogeneous zones in the tree for a given value
 *     occurrence number for a given value in the whole tree
 *
 **/

void TreeCharacteristics::build_characteristics(int imin_value,
                                                int imax_value,
                                                int imax_size,
                                                int imax_depth,
                                                int inb_trees,
                                                Typed_edge_one_int_tree** otrees1,
                                                int ivariable,
                                                bool final_value)

{
   bool build= true;
   int val, nb_values= _max_value - _min_value + 1,
       cp_nb_values= nb_values;

   _variable= ivariable;
   _min_value= imin_value;
   _max_value= imax_value;
   _nb_trees= inb_trees;
   remove();

   if (nb_values <= NB_OUTPUT)
      init(imax_size, imax_depth);
   else
      marginal_distribution= new FrequencyDistribution(_max_value+1);

   build_marginal_frequency_distribution(otrees1, final_value);

   nb_values= marginal_distribution->nb_value;

   if (nb_values > NB_OUTPUT)
      // the characteristics are not built if the number of values is too high...
      build=false;

   for (val= 0; val < nb_values; val++)
      if (marginal_distribution->frequency[val] == 0)
      {
         build= false;
         break;
      }

   if (build)
   {
      build_first_occurrence_root_frequency_distribution(otrees1, imax_depth, final_value);
      build_first_occurrence_leaves_frequency_distribution(otrees1, imax_depth, final_value);
      build_zone_frequency_distributions(otrees1, imax_size, final_value);
      build_nb_occurrences_frequency_distribution(otrees1, imax_size, final_value);
      // build_children_pairs_frequency_distribution(otrees1, final_value)
   }
   else
   {
      remove_characteristic(first_occurrence_root, cp_nb_values);
      remove_characteristic(first_occurrence_leaves, cp_nb_values);
      remove_characteristic(sojourn_size, cp_nb_values);
      remove_characteristic(nb_zones, cp_nb_values);
      remove_characteristic(nb_occurrences, cp_nb_values);
   }
}

/*****************************************************************
 *
 *  Access to class members of TreeCharacteristics
 *
 **/

int TreeCharacteristics::get_variable() const
{ return _variable; }

int TreeCharacteristics::get_nb_values() const
{ return(_max_value-_min_value+1); }

/*****************************************************************
 *
 *  Access to the characteristic quantity frequency distributions
 *
 **/

FrequencyDistribution* TreeCharacteristics::get_marginal_distribution() const
{
   FrequencyDistribution *res;

   if (marginal_distribution != NULL)
      res= new FrequencyDistribution(*marginal_distribution);
   else
      res= NULL;
   return res;
}
  // ? Tree_curves* get_index_value() const
  //  { return index_value; }
  // ? Curves* get_depth_value() const
  //  { return depth_value; }


FrequencyDistribution* TreeCharacteristics::get_characteristic(int value,
                                                               FrequencyDistribution** h) const
{
   FrequencyDistribution *res;

   assert((_min_value <= value) && (value <= _max_value));
   if (h != NULL)
   {
      if (h[value] != NULL)
         res= new FrequencyDistribution(*(h[value-_min_value]));
      else
         res= NULL;
   }
   else
      res= NULL;

   return res;
}

/*****************************************************************
 *
 *  Check whether given characteristic, refered to by charac, is present
 *
 **/

bool TreeCharacteristics::is_characteristic(int charac) const
{
   bool present= true;
   const int nb_values= _max_value - _min_value + 1;
   int value;

   switch (charac)
   {
      case FIRST_OCCURRENCE_ROOT :
         if (first_occurrence_root == NULL)
        present= false;
     else
        for(value= 0; value < nb_values; value++)
               if (first_occurrence_root[value] == NULL)
                  present= false;
         break;
      case FIRST_OCCURRENCE_LEAVES :
         if (first_occurrence_leaves == NULL)
        present= false;
     else
        for(value= 0; value < nb_values; value++)
               if (first_occurrence_leaves[value] == NULL)
                  present= false;
         break;
      case SOJOURN_SIZE :
         if (sojourn_size == NULL)
        present= false;
     else
        for(value= 0; value < nb_values; value++)
               if (sojourn_size[value] == NULL)
                  present= false;
         break;
      case NB_ZONES :
         if (nb_zones == NULL)
        present= false;
     else
        for(value= 0; value < nb_values; value++)
               if (nb_zones[value] == NULL)
                  present= false;
         break;
      case NB_OCCURRENCES :
         if (nb_occurrences == NULL)
        present= false;
     else
        for(value= 0; value < nb_values; value++)
               if (nb_occurrences[value] == NULL)
                  present= false;
         break;
  }
  return present;
}


FrequencyDistribution* TreeCharacteristics::get_first_occurrence_root(int value) const
{ return get_characteristic(value, first_occurrence_root); }

FrequencyDistribution* TreeCharacteristics::get_first_occurrence_leaves(int value) const
{ return get_characteristic(value, first_occurrence_leaves); }

FrequencyDistribution* TreeCharacteristics::get_sojourn_size(int value) const
{ return get_characteristic(value, sojourn_size); }

FrequencyDistribution* TreeCharacteristics::get_nb_zones(int value) const
{ return get_characteristic(value, nb_zones); }

FrequencyDistribution* TreeCharacteristics::get_nb_occurrences(int value) const
{ return get_characteristic(value, nb_occurrences); }

int TreeCharacteristics::get_nb_value_first_occurrence_root(int value) const
{ return first_occurrence_root[value]->nb_value; }

int TreeCharacteristics::get_nb_value_sojourn_size(int value) const
{ return sojourn_size[value]->nb_value; }

std::ostream& TreeCharacteristics::ascii_write_marginal_distribution(std::ostream &os,
                                                                     bool exhaustive,
                                                                     bool file_flag) const
{
   if (marginal_distribution != NULL)
      marginal_distribution->ascii_write(os, exhaustive, file_flag);
   return os;
}

void TreeCharacteristics::ascii_write_characteristic(FrequencyDistribution** c,
                                                     int inb_values,
                                                     std::ostream& os,
                                                     std::string c1,
                                                     std::string c2,
                                                     bool exhaustive,
                                                     bool file_flag) const
{
   int value;

   if (c != NULL)
      for(value= 0; value < inb_values; value++)
         if (c[value] != NULL)
         {
            os << c1 << value+_min_value << c2 << endl;
            c[value]->ascii_write(os, exhaustive, file_flag);
            os << endl;
         }
}

/*****************************************************************
 *
 *  Print the TreeCharacteristics using an output stream,
 *  the variable type (STATE / VALUE), a frequency distribution for the tree size,
 *  a flag on the level of detail and one on printing into a file
 *
 **/

std::ostream& TreeCharacteristics::ascii_print(std::ostream& os,
                                               int type,
                                               const FrequencyDistribution& hsize,
                                               bool exhaustive,
                                               bool comment_flag) const

{
   register int val;
   int nb_value= _max_value - _min_value + 1;

   if ((first_occurrence_root != NULL) || (first_occurrence_leaves != NULL)
       || (sojourn_size != NULL) || (nb_zones != NULL) || (nb_occurrences != NULL))
   {
      if (exhaustive)
      {
         os << "\n";
         if (comment_flag)
            os << "# ";

         os << "  ";

         for (val= 0; val < nb_value; val++)
            os << " | " << STAT_TREES_label[TREESTATL_OBSERVED] << " "
               << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << val;
         os << " | " << STAT_label[STATL_FREQUENCY] << endl;

         // index_value->ascii_print(os , comment_flag);
      }
   }

   // frequency distributions of first occurrence
   if (first_occurrence_root != NULL)
   {
      for (val= 0; val < nb_value; val++)
      {
         os << "\n";
         if (comment_flag)
            os << "# ";

         os << STAT_TREES_label[type == STATE ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
            << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
         first_occurrence_root[val]->ascii_characteristic_print(os, false, comment_flag);

         if (exhaustive)
         {
            os << "\n";
            if (comment_flag)
               os << "# ";

            os << "   | " << STAT_TREES_label[type == STATE ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
               << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
            first_occurrence_root[val]->ascii_print(os, comment_flag);
         }
      }
   }

   if (first_occurrence_leaves != NULL)
   {
      for (val= 0; val < nb_value; val++)
      {
         os << "\n";
         if (comment_flag)
            os << "# ";

         os << STAT_TREES_label[type == STATE ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES]
            << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
         first_occurrence_leaves[val]->ascii_characteristic_print(os, false, comment_flag);

         if (exhaustive)
         {
            os << "\n";
            if (comment_flag)
               os << "# ";

            os << "   | " << STAT_TREES_label[type == STATE ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES]
               << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
            first_occurrence_leaves[val]->ascii_print(os, comment_flag);
         }
      }
   }

   // sojourn size
   if (sojourn_size != NULL)
   {
      for (val= 0; val < nb_value; val++)
      {
         os << "\n";
         if (comment_flag)
            os << "# ";

         os << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << val << " "
            << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
         sojourn_size[val]->ascii_characteristic_print(os, false, comment_flag);

         if ((sojourn_size[val]->nb_element > 0) && (exhaustive))
         {
            os << "\n";
            if (comment_flag)
               os << "# ";

            os << "   | " << STAT_label[type == STATE ? STATL_STATE : STATL_VALUE] << " " << val << " "
               << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
            sojourn_size[val]->ascii_print(os, comment_flag);
          }
      }
   }

   // number of zones
   if (nb_zones != NULL)
   {
      for (val= 0; val < nb_value; val++)
      {
         os << "\n";
         if (comment_flag)
            os << "# ";

         os << STAT_TREES_label[type == STATE ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
            << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
         nb_zones[val]->ascii_characteristic_print(os , (hsize.variance > 0. ? false : true) , comment_flag);

         if (exhaustive)
         {
            os << "\n";
            if (comment_flag)
               os << "# ";

            os << "   | " << STAT_TREES_label[type == STATE ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
               << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
            nb_zones[val]->ascii_print(os, comment_flag);
         }
      }
   }

   // number of occurrences
   if (nb_occurrences != NULL)
   {
      for (val= 0; val < nb_value; val++)
      {
        os << "\n";
        if (comment_flag)
           os << "# ";

        os << STAT_TREES_label[type == STATE ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
           << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << " - ";
        nb_occurrences[val]->ascii_characteristic_print(os, (hsize.variance > 0. ? false : true) , comment_flag);

        if (exhaustive)
        {
           os << "\n";
           if (comment_flag)
              os << "# ";

           os << "   | " << STAT_TREES_label[type == STATE ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
              << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << endl;
           nb_occurrences[val]->ascii_print(os , comment_flag);
        }
      }
   }
   return os;
}

std::ostream& TreeCharacteristics::spreadsheet_print(std::ostream& os,
                                                     int type,
                                                     const FrequencyDistribution& hsize) const
{ return os; }

/*****************************************************************
 *  Gnuplot output of TreeCharacteristics using
 *  a prefix for the files, the title of the output figures,
 *  the considered variable, the total number of variables,
 *  the variable type (STATE / INT_VALUE)
 *  and the frequency distribution of the tree sizes.
 **/

bool TreeCharacteristics::plot_print(const char * prefix,
                                     const char * title,
                                     int variable,
                                     int nb_variables,
                                     int type,
                                     const FrequencyDistribution& hsize) const
{
   bool status= true, start;
   register int val, i, j, nb_values;
   int nb_histo, histo_index;
   const FrequencyDistribution **phisto;
   ostringstream data_file_name[2];

   nb_values= _max_value - _min_value + 1;


   // used only for smoothed probabilities
   // data_file_name[0] << prefix << variable+1 << 0 << ".dat";

   phisto= new const FrequencyDistribution*[NB_OUTPUT*5];

   // name of the file containing the data for all the graphs
   data_file_name[1] << prefix << variable+1 << 1 << ".dat";

   nb_histo= 0;
   // create a list of the successive frequency distributions to be printed
   if (first_occurrence_root != NULL)
      for (val= 0; val < nb_values; val++)
         phisto[nb_histo++]= first_occurrence_root[val];

   if (first_occurrence_leaves != NULL)
      for (val= 0; val < nb_values; val++)
         phisto[nb_histo++]= first_occurrence_leaves[val];

   if (sojourn_size != NULL)
      for (val= 0; val < nb_values; val++)
        if (sojourn_size[val]->nb_element > 0)
           phisto[nb_histo++] = sojourn_size[val];

   // note that the number of zones and the sojourn size
   // are alternate
   if ((nb_zones != NULL) && (nb_occurrences!= NULL))
      for (val= 0; val < nb_values; val++)
      {
         phisto[nb_histo++] = nb_zones[val];
         phisto[nb_histo++] = nb_occurrences[val];
      }

   hsize.plot_print((data_file_name[1].str()).c_str(), nb_histo, phisto);
   // create the file, named prefix(variable+1)1.dat, containing
   // the values for all the gnuplot graphs

   // what follows is useful for potential smoothed probabilities
   // and uses the file named prefix(variable+1)0.dat
   /*
   for (i= 0; i < 2; i++)
   {

      ostringstream file_name[2];

      switch (i)
      {
         case 0 :
            file_name[0] << prefix << variable + 1 << 0 << ".plot";
            break;
         case 1 :
            file_name[0] << prefix << variable + 1 << 0 << ".print";
            break;
      }

      ofstream out_file((file_name[0].str()).c_str());

      if (i == 1)
      {
         out_file << "set terminal postscript" << endl;
         file_name[1] << label(prefix) << variable + 1 << 1 << ".ps";
         out_file << "set output \"" << file_name[1].str() << "\"\n\n";
      }

      out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n";

      out_file << "set title \"";
      if (title) {
        out_file << title;
        if (nb_variables > 1) {
          out_file << " - ";
        }
      }
      if (nb_variables > 1) {
        out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
      }
      out_file << "\"\n\n";

      for (j = 0;j < nb_values;j++) {
        out_file << "\"" << label((data_file_name[0].str()).c_str()) << "\" using "
                 << j + 1 << " title \"" << STAT_TREES_label[TREESTATL_OBSERVED] << " "
                 << STAT_label[type == STATE ? STATL_STATE : STATL_OUTPUT] << " "
                 << j << "\" with linespoints";
        if (j < nb_values - 1) {
          out_file << ",\\";
        }
        out_file << endl;
      }

      if (i == 0) {
        out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
      }
      out_file << endl;

      if (hsize.nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics 0,1" << endl;
      }
      if ((int)(hsize.max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics 0,1" << endl;
      }

      out_file << "plot [0:" << hsize.nb_value - 1 << "] [0:"
               << (int)(hsize.max * YSCALE) + 1 << "] \""
               << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
               << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
               << "\" with impulses" << endl;

      if (hsize.nb_value - 1 < TIC_THRESHOLD) {
        out_file << "set xtics autofreq" << endl;
      }
      if ((int)(hsize.max * YSCALE) + 1 < TIC_THRESHOLD) {
        out_file << "set ytics autofreq" << endl;
      }

      if (i == 1) {
        out_file << "\nset terminal x11" << endl;
      }

      out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
   }
   */

   // first_occurrence_root
   histo_index= 2;
   j= histo_index;

   if (first_occurrence_root != NULL)
   {
      for (i= 0; i < 2; i++)
      {
         // i == 0 -> output on terminal (.plot file)
         // i == 1 -> output into a postscript file  (.plot file)
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << variable+1 << 1 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << variable+1 << 1 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << variable+1 << 1 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";
         if ((title != NULL) || (nb_variables > 1))
         {
            out_file << " \"";
            if (title)
            {
               out_file << title;
               if (nb_variables > 1)
                  out_file << " - ";
            }
            if (nb_variables > 1)
               out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;

            out_file << "\"";
         }
         out_file << "\n\n";

         j= histo_index;

         start= true;

         for (val= 0; val < nb_values; val++)
         {
            if (!start)
            {
              if (i == 0)
                 out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              out_file << endl;
            }
            else
               start= false;

            if (MAX(1, first_occurrence_root[val]->nb_value-1) < TIC_THRESHOLD)
               out_file << "set xtics 0,1" << endl;

            if ((int)(first_occurrence_root[val]->max * YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics 0,1" << endl;

            out_file << "plot [0:" << MAX(first_occurrence_root[val]->nb_value-1 , 1) << "] [0:"
                       << (int)(first_occurrence_root[val]->max*YSCALE)+1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                       << " title \"" << STAT_TREES_label[type == STATE ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                       << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

            if (MAX(1 , first_occurrence_root[val]->nb_value-1) < TIC_THRESHOLD)
               out_file << "set xtics autofreq" << endl;

            if ((int)(first_occurrence_root[val]->max*YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics autofreq" << endl;
         }

         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
   }

   // first_occurrence_leaves
   histo_index = j;

   if (first_occurrence_leaves != NULL)
   {
      for (i= 0; i < 2; i++)
      {
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << variable+1 << 2 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << variable+1 << 2 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << variable+1 << 2 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";
         if ((title != NULL) || (nb_variables > 1))
         {
            out_file << " \"";
            if (title)
            {
               out_file << title;
               if (nb_variables > 1)
                  out_file << " - ";
            }
            if (nb_variables > 1)
               out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;

            out_file << "\"";
         }
         out_file << "\n\n";

         j= histo_index;

         start= true;

         for (val= 0; val < nb_values; val++)
         {
            if (!start)
            {
              if (i == 0)
                 out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              out_file << endl;
            }
            else
               start= false;

            if (MAX(1, first_occurrence_leaves[val]->nb_value-1) < TIC_THRESHOLD)
               out_file << "set xtics 0,1" << endl;

            if ((int)(first_occurrence_leaves[val]->max * YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics 0,1" << endl;

            out_file << "plot [0:" << MAX(first_occurrence_leaves[val]->nb_value-1 , 1) << "] [0:"
                       << (int)(first_occurrence_leaves[val]->max*YSCALE)+1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                       << " title \"" << STAT_TREES_label[type == STATE ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES]
                       << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

            if (MAX(1 , first_occurrence_leaves[val]->nb_value-1) < TIC_THRESHOLD)
               out_file << "set xtics autofreq" << endl;

            if ((int)(first_occurrence_leaves[val]->max*YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics autofreq" << endl;
         }

         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
   }

   // sojourn_size
   histo_index= j;

   if (sojourn_size != NULL)
   {
      for (i= 0; i < 2; i++)
      {
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << variable+1 << 3 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << variable+1 << 3 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << variable+1 << 3 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";

         if ((title != NULL) || (nb_variables > 1))
         {
            out_file << " \"";
            if (title != NULL)
            {
               out_file << title;
               if (nb_variables > 1)
                  out_file << " - ";
            }
            if (nb_variables > 1)
               out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
            out_file << "\"";
         }
         out_file << "\n\n";

         j= histo_index;

         start= true;
         for (val= 0; val < nb_values; val++)
         {
            if (sojourn_size[val]->nb_element > 0)
            {
               if (!start)
               {
                  if (i == 0)
                     out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                  out_file << endl;
               }
               else
                  start= false;

               if (sojourn_size[val]->nb_value-1 < TIC_THRESHOLD)
                  out_file << "set xtics 0,1" << endl;

               if ((int)(sojourn_size[val]->max*YSCALE)+1 < TIC_THRESHOLD)
                  out_file << "set ytics 0,1" << endl;

               out_file << "plot [0:" << sojourn_size[val]->nb_value-1 << "] [0:"
                        << (int)(sojourn_size[val]->max*YSCALE)+1 << "] \""
                        << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                        << " title \"" << STAT_label[type == STATE ? STATL_STATE : STATL_OUTPUT] << " " << val
                        << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                        << "\" with impulses" << endl;

               if (sojourn_size[val]->nb_value-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;

               if ((int)(sojourn_size[val]->max*YSCALE)+1 < TIC_THRESHOLD)
                  out_file << "set ytics autofreq" << endl;
            }
         }
         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
   }

   // nb_zones and nb_occurrences
   histo_index= j;

   if ((nb_zones != NULL) && (nb_occurrences != NULL))
   {
      for (i= 0; i < 2; i++)
      {
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << variable+1 << 4 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << variable+1 << 4 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << variable + 1 << 4 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";

         if ((title != NULL) || (nb_variables > 1))
         {
            out_file << " \"";
            if (title != NULL)
            {
               out_file << title;
               if (nb_variables > 1)
                  out_file << " - ";
            }
            if (nb_variables > 1)
               out_file << STAT_label[STATL_VARIABLE] << " " << variable+1;

            out_file << "\"";
         }
         out_file << "\n\n";

         j= histo_index;

         for (val= 0; val < nb_values; val++)
         {
           if (nb_zones[val]->nb_value-1 < TIC_THRESHOLD)
              out_file << "set xtics 0,1" << endl;

           if ((int)(nb_zones[val]->max*YSCALE)+1 < TIC_THRESHOLD)
              out_file << "set ytics 0,1" << endl;

           out_file << "plot [0:" << nb_zones[val]->nb_value-1 << "] [0:"
                    << (int)(nb_zones[val]->max * YSCALE) + 1 << "] \""
                    << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                    << " title \"" << STAT_TREES_label[type == STATE ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                    << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                    << "\" with impulses" << endl;

           if (nb_zones[val]->nb_value-1 < TIC_THRESHOLD)
              out_file << "set xtics autofreq" << endl;

           if ((int)(nb_zones[val]->max*YSCALE)+1 < TIC_THRESHOLD)
              out_file << "set ytics autofreq" << endl;

           if (i == 0)
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;

           out_file << endl;

           if (nb_occurrences[val]->nb_value-1 < TIC_THRESHOLD)
              out_file << "set xtics 0,1" << endl;

           if ((int)(nb_occurrences[val]->max*YSCALE)+1 < TIC_THRESHOLD)
              out_file << "set ytics 0,1" << endl;

           out_file << "plot [0:" << nb_occurrences[val]->nb_value-1 << "] [0:"
                    << (int)(nb_occurrences[val]->max*YSCALE)+1 << "] \""
                    << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                    << " title \"" << STAT_TREES_label[type == STATE ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                    << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                    << "\" with impulses" << endl;

           if (nb_occurrences[val]->nb_value-1 < TIC_THRESHOLD)
              out_file << "set xtics autofreq" << endl;

           if ((int)(nb_occurrences[val]->max*YSCALE)+1 < TIC_THRESHOLD)
              out_file << "set ytics autofreq" << endl;

           if (i == 0)
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;

           out_file << endl;
         }

         if (hsize.nb_value-1 < TIC_THRESHOLD)
            out_file << "set xtics 0,1" << endl;

         if ((int)(hsize.max*YSCALE)+1 < TIC_THRESHOLD)
            out_file << "set ytics 0,1" << endl;

         out_file << "plot [0:" << hsize.nb_value - 1 << "] [0:"
                  << (int)(hsize.max*YSCALE)+1 << "] \""
                  << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                  << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                  << "\" with impulses" << endl;

         if (hsize.nb_value-1 < TIC_THRESHOLD)
            out_file << "set xtics autofreq" << endl;

         if ((int)(hsize.max*YSCALE)+1 < TIC_THRESHOLD)
            out_file << "set ytics autofreq" << endl;

         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
   }

   delete [] phisto;
   phisto= NULL;

   return status;
}

/*****************************************************************
 *
 *  Graphical output of TreeCharacteristics
 *  using a plotable data structure and the type of characteristic
 *
 **/

MultiPlotSet* TreeCharacteristics::get_plotable(int plot_type,
                                                int variable,
                                                int nb_variables,
                                                int type) const
{
   bool status= true, start;
   const int nb_values= _max_value - _min_value + 1;;
   int val;
   MultiPlotSet *plotset= new MultiPlotSet(nb_values);
   MultiPlotSet &plot= *plotset;
   std::ostringstream vstring;

   switch (plot_type)
   {
      case FIRST_OCCURRENCE_ROOT :
      {
         plot.border= "15 lw 0";
         // plot.bidule= "set tics out";
         // plot.bidule= "nomirror";
         vstring.str("");
         if (nb_variables > 1)
            vstring << STAT_label[STATL_VARIABLE] << " " << variable + 1;
         plot.title= vstring.str();

         for (val= 0; val < nb_values; val++)
         {
            plot[val].resize(1);
            plot[val][0].style= "impulses";
            vstring.str("");
            if (MAX(1, first_occurrence_root[val]->nb_value-1) < TIC_THRESHOLD)
               plot[val].xtics= 1;

            if ((int)(first_occurrence_root[val]->max * YSCALE)+1 < TIC_THRESHOLD)
               plot[val].ytics= 1;

            plot[val].xrange= Range(0, MAX(first_occurrence_root[val]->nb_value-1 , 1));
            plot[val].yrange= Range(0, (int)(first_occurrence_root[val]->max*YSCALE)+1);
            vstring << STAT_TREES_label[type == STATE ? TREESTATL_STATE_FIRST_OCCURRENCE_ROOT : TREESTATL_OUTPUT_FIRST_OCCURRENCE_ROOT]
                    << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION];
            plot[val].title= vstring.str();
            first_occurrence_root[val]->plotable_frequency_write(plot[val][0]);

            // plot[val].bidule ="autofreq";
         }
         break;
      }
      case FIRST_OCCURRENCE_LEAVES :
      {
         break;
      }
      case SOJOURN_SIZE :
      {
         break;
      }
      case NB_ZONES :
      {
         break;
      }
      case NB_OCCURRENCES :
      {
         break;
      }
   }
   /*
   // first_occurrence_leaves
   histo_index = j;

   if (first_occurrence_leaves != NULL)
   {
      for (i= 0; i < 2; i++)
      {
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << variable+1 << 2 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << variable+1 << 2 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << variable+1 << 2 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";
         if ((title != NULL) || (nb_variables > 1))
         {
            out_file << " \"";
            if (title)
            {
               out_file << title;
               if (nb_variables > 1)
                  out_file << " - ";
            }
            if (nb_variables > 1)
               out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;

            out_file << "\"";
         }
         out_file << "\n\n";

         j= histo_index;

         start= true;

         for (val= 0; val < nb_values; val++)
         {
            if (!start)
            {
              if (i == 0)
                 out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
              out_file << endl;
            }
            else
               start= false;

            if (MAX(1, first_occurrence_leaves[val]->nb_value-1) < TIC_THRESHOLD)
               out_file << "set xtics 0,1" << endl;

            if ((int)(first_occurrence_leaves[val]->max * YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics 0,1" << endl;

            out_file << "plot [0:" << MAX(first_occurrence_leaves[val]->nb_value-1 , 1) << "] [0:"
                       << (int)(first_occurrence_leaves[val]->max*YSCALE)+1 << "] \""
                       << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                       << " title \"" << STAT_TREES_label[type == STATE ? TREESTATL_STATE_FIRST_OCCURRENCE_LEAVES : TREESTATL_OUTPUT_FIRST_OCCURRENCE_LEAVES]
                       << " " << val << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION] << "\" with impulses" << endl;

            if (MAX(1 , first_occurrence_leaves[val]->nb_value-1) < TIC_THRESHOLD)
               out_file << "set xtics autofreq" << endl;

            if ((int)(first_occurrence_leaves[val]->max*YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics autofreq" << endl;
         }

         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
   }

   // sojourn_size
   histo_index= j;

   if (sojourn_size != NULL)
   {
      for (i= 0; i < 2; i++)
      {
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << variable+1 << 3 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << variable+1 << 3 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << variable+1 << 3 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";

         if ((title != NULL) || (nb_variables > 1))
         {
            out_file << " \"";
            if (title != NULL)
            {
               out_file << title;
               if (nb_variables > 1)
                  out_file << " - ";
            }
            if (nb_variables > 1)
               out_file << STAT_label[STATL_VARIABLE] << " " << variable + 1;
            out_file << "\"";
         }
         out_file << "\n\n";

         j= histo_index;

         start= true;
         for (val= 0; val < nb_values; val++)
         {
            if (sojourn_size[val]->nb_element > 0)
            {
               if (!start)
               {
                  if (i == 0)
                     out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;
                  out_file << endl;
               }
               else
                  start= false;

               if (sojourn_size[val]->nb_value-1 < TIC_THRESHOLD)
                  out_file << "set xtics 0,1" << endl;

               if ((int)(sojourn_size[val]->max*YSCALE)+1 < TIC_THRESHOLD)
                  out_file << "set ytics 0,1" << endl;

               out_file << "plot [0:" << sojourn_size[val]->nb_value-1 << "] [0:"
                        << (int)(sojourn_size[val]->max*YSCALE)+1 << "] \""
                        << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                        << " title \"" << STAT_label[type == STATE ? STATL_STATE : STATL_OUTPUT] << " " << val
                        << " " << STAT_TREES_label[TREESTATL_SOJOURN_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                        << "\" with impulses" << endl;

               if (sojourn_size[val]->nb_value-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;

               if ((int)(sojourn_size[val]->max*YSCALE)+1 < TIC_THRESHOLD)
                  out_file << "set ytics autofreq" << endl;
            }
         }
         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
   }

   // nb_zones and nb_occurrences
   histo_index= j;

   if ((nb_zones != NULL) && (nb_occurrences != NULL))
   {
      for (i= 0; i < 2; i++)
      {
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
               file_name[0] << prefix << variable+1 << 4 << ".plot";
               break;
            case 1 :
               file_name[0] << prefix << variable+1 << 4 << ".print";
               break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << variable + 1 << 4 << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";

         if ((title != NULL) || (nb_variables > 1))
         {
            out_file << " \"";
            if (title != NULL)
            {
               out_file << title;
               if (nb_variables > 1)
                  out_file << " - ";
            }
            if (nb_variables > 1)
               out_file << STAT_label[STATL_VARIABLE] << " " << variable+1;

            out_file << "\"";
         }
         out_file << "\n\n";

         j= histo_index;

         for (val= 0; val < nb_values; val++)
         {
           if (nb_zones[val]->nb_value-1 < TIC_THRESHOLD)
              out_file << "set xtics 0,1" << endl;

           if ((int)(nb_zones[val]->max*YSCALE)+1 < TIC_THRESHOLD)
              out_file << "set ytics 0,1" << endl;

           out_file << "plot [0:" << nb_zones[val]->nb_value-1 << "] [0:"
                    << (int)(nb_zones[val]->max * YSCALE) + 1 << "] \""
                    << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                    << " title \"" << STAT_TREES_label[type == STATE ? TREESTATL_STATE_NB_ZONES : TREESTATL_OUTPUT_NB_ZONES]
                    << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                    << "\" with impulses" << endl;

           if (nb_zones[val]->nb_value-1 < TIC_THRESHOLD)
              out_file << "set xtics autofreq" << endl;

           if ((int)(nb_zones[val]->max*YSCALE)+1 < TIC_THRESHOLD)
              out_file << "set ytics autofreq" << endl;

           if (i == 0)
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;

           out_file << endl;

           if (nb_occurrences[val]->nb_value-1 < TIC_THRESHOLD)
              out_file << "set xtics 0,1" << endl;

           if ((int)(nb_occurrences[val]->max*YSCALE)+1 < TIC_THRESHOLD)
              out_file << "set ytics 0,1" << endl;

           out_file << "plot [0:" << nb_occurrences[val]->nb_value-1 << "] [0:"
                    << (int)(nb_occurrences[val]->max*YSCALE)+1 << "] \""
                    << label((data_file_name[1].str()).c_str()) << "\" using " << j++
                    << " title \"" << STAT_TREES_label[type == STATE ? TREESTATL_STATE_NB_OCCURRENCES : TREESTATL_OUTPUT_NB_OCCURRENCES]
                    << " " << val << " " << STAT_TREES_label[TREESTATL_PER_TREE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                    << "\" with impulses" << endl;

           if (nb_occurrences[val]->nb_value-1 < TIC_THRESHOLD)
              out_file << "set xtics autofreq" << endl;

           if ((int)(nb_occurrences[val]->max*YSCALE)+1 < TIC_THRESHOLD)
              out_file << "set ytics autofreq" << endl;

           if (i == 0)
              out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;

           out_file << endl;
         }

         if (hsize.nb_value-1 < TIC_THRESHOLD)
            out_file << "set xtics 0,1" << endl;

         if ((int)(hsize.max*YSCALE)+1 < TIC_THRESHOLD)
            out_file << "set ytics 0,1" << endl;

         out_file << "plot [0:" << hsize.nb_value - 1 << "] [0:"
                  << (int)(hsize.max*YSCALE)+1 << "] \""
                  << label((data_file_name[1].str()).c_str()) << "\" using 1 title \""
                  << STAT_TREES_label[TREESTATL_TREE_SIZE] << " " << STAT_label[STATL_FREQUENCY_DISTRIBUTION]
                  << "\" with impulses" << endl;

         if (hsize.nb_value-1 < TIC_THRESHOLD)
            out_file << "set xtics autofreq" << endl;

         if ((int)(hsize.max*YSCALE)+1 < TIC_THRESHOLD)
            out_file << "set ytics autofreq" << endl;

         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
   }

   */
   return plotset;
}


std::ostream& TreeCharacteristics::ascii_write_first_occurrence_root(std::ostream &os,
                                                                     bool exhaustive,
                                                                     bool file_flag) const
{
   ascii_write_characteristic(first_occurrence_root,
                              _max_value-_min_value+1,
                              os,
                              "First occurrence of value ",
                              " (from root) : ",
                              exhaustive,
                              file_flag);
   return os;
}

std::ostream& TreeCharacteristics::ascii_write_first_occurrence_leaves(std::ostream &os,
                                                                       bool exhaustive,
                                                                       bool file_flag) const
{
   ascii_write_characteristic(first_occurrence_leaves,
                              _max_value-_min_value+1,
                              os,
                              "First occurrence of value ",
                              " (from the leaves) : ",
                              exhaustive,
                              file_flag);
   return os;
}

std::ostream& TreeCharacteristics::ascii_write_sojourn_size(std::ostream &os,
                                                            bool exhaustive,
                                                            bool file_flag) const
{
   ascii_write_characteristic(sojourn_size,
                              _max_value-_min_value+1,
                              os,
                              "Size of homogeneous zones for value ",
                              " : ",
                              exhaustive,
                              file_flag);
   return os;
}

std::ostream& TreeCharacteristics::ascii_write_nb_zones(std::ostream &os,
                                                        bool exhaustive,
                                                        bool file_flag) const
{
   ascii_write_characteristic(nb_zones,
                              _max_value-_min_value+1,
                              os,
                              "Number of homogeneous zones for value ",
                              " : ",
                              exhaustive,
                              file_flag);
   return os;
}

std::ostream& TreeCharacteristics::ascii_write_nb_occurrences(std::ostream &os,
                                                              bool exhaustive,
                                                              bool file_flag) const
{
   ascii_write_characteristic(nb_occurrences,
                              _max_value-_min_value+1,
                              os,
                              "Number of occurrences for value ",
                              " : ",
                              exhaustive,
                              file_flag);
   return os;
}

/*****************************************************************
 *
 *  Deallocation of the pointers for TreeCharacteristics
 *
 **/

void TreeCharacteristics::remove_characteristic(FrequencyDistribution**& c,
                                                int inb_values)
{
   int value;

   if (c != NULL)
   {
      for(value= 0; value < inb_values; value++)
         if (c[value] != NULL)
         {
            delete c[value];
            c[value]= NULL;
         }
      delete [] c;
      c= NULL;
   }
}

void TreeCharacteristics::remove()
{ // deallocation
  int nb_values= _max_value - _min_value + 1; //value,

  if (marginal_distribution != NULL)
  {
     delete marginal_distribution;
     marginal_distribution= NULL;
  }

  if (marginal_histogram != NULL)
  {
     delete marginal_histogram;
     marginal_histogram= NULL;
  }


  remove_characteristic(first_occurrence_root, nb_values);

  remove_characteristic(first_occurrence_leaves, nb_values);

  remove_characteristic(sojourn_size, nb_values);

  remove_characteristic(nb_zones, nb_values);

  remove_characteristic(nb_occurrences, nb_values);
}

/*****************************************************************
 *
 *  Allocation of the frequency distributions given the number of values
 *  of the variable
 *
 **/

void TreeCharacteristics::init_characteristic(FrequencyDistribution**& c,
                                              int inb_values,
                                              int imax_val)
{
   int value;

   c= new FrequencyDistribution*[inb_values];

   for(value= 0; value < inb_values; value++)
      c[value]= new FrequencyDistribution(imax_val+1);
};


void TreeCharacteristics::init(int imax_size, int imax_depth)
{

   // int value;
   int nb_values= _max_value - _min_value + 1;

   assert(nb_values > 0);
   // allocating non existing frequency distributions is pointless

   marginal_distribution= new FrequencyDistribution(_max_value+1);

   init_characteristic(first_occurrence_root, nb_values, imax_depth);

   init_characteristic(first_occurrence_leaves, nb_values, imax_depth);

   // init_characteristic(sojourn_size, nb_values, imax_size);

   init_characteristic(nb_zones, nb_values, imax_size);

   init_characteristic(nb_occurrences, nb_values, imax_size);
}


/*****************************************************************
 *
 *  Construction of the frequency distribution for the marginal distribution
 *
 **/

void TreeCharacteristics::build_marginal_frequency_distribution(Typed_edge_one_int_tree** otrees1,
                                                                bool final_value)
{ // marginal frequency distribution

   typedef Typed_edge_one_int_tree::vertex_iterator vertex_iterator;

   int t, v;
   vertex_iterator it, end;

   if ((otrees1 != NULL) &&(_min_value >= 0) && (_max_value <= MARGINAL_DISTRIBUTION_MAX_VALUE))
   {
      if (marginal_distribution == NULL)
         marginal_distribution = new FrequencyDistribution(_max_value+1);

      if (final_value)
      {
         for(t = 0; t < _nb_trees; t++)
         {
            if (otrees1[t] != NULL)
            {
               Tree_tie::tie(it, end)= otrees1[t]->vertices();
               while (it < end)
               {
                  v= (otrees1[t]->get(*it++)).Int();
                  if (v != I_DEFAULT)
                     (marginal_distribution->frequency[v])++;
               }
            }
         }
      }
      else
      {
         for(t = 0; t < _nb_trees; t++)
         {
            if (otrees1[t] != NULL)
            {
               Tree_tie::tie(it, end)= otrees1[t]->vertices();
               while (it < end)
               {
                  v= (otrees1[t]->get(*it++)).Int();
                  (marginal_distribution->frequency[v])++;
               }
            }
         }
      }


      marginal_distribution->offset = _min_value;
      marginal_distribution->nb_element_computation();
      marginal_distribution->max_computation();
      marginal_distribution->mean_computation();
      marginal_distribution->variance_computation();
   }
}

/*****************************************************************
 *
 *  Construction of the frequency distribution for the distribution of the
 *  1st occurrence depth for a given value, starting from root
 *
 **/

void TreeCharacteristics::build_first_occurrence_root_frequency_distribution(Typed_edge_one_int_tree** otrees1,
                                                                             int imax_depth,
                                                                             bool final_value)
{ //  frequency distribution of first occurrence (root)

   typedef Typed_edge_one_int_tree::vertex_iterator vertex_iterator;
   typedef Typed_edge_one_int_tree::vertex_array vertex_array;

   int t, v, curr_val, nb_values= _max_value - _min_value + 1;
   int size, index;
   vertex_array va;
   generic_visitor<Typed_edge_one_int_tree> visitor;
   bool *occurrence;

//      first_occurrence_root= new  FrequencyDistribution*[_nb_values];
//     for(value= 0; value < nb_values; value++)
//       first_occurrence_root[value]= new FrequencyDistribution(imax_depth);

   if (otrees1 != NULL)
   {
      if (first_occurrence_root == NULL)
      {
         first_occurrence_root= new FrequencyDistribution*[nb_values];
         for(v=0; v < nb_values; v++)
            first_occurrence_root[v]= NULL;
      }

      for(v= 0; v < nb_values; v++)
         if (first_occurrence_root[v] == NULL)
            first_occurrence_root[v]= new FrequencyDistribution(imax_depth+1);

      occurrence= new bool[nb_values];

      for(t= 0; t < _nb_trees; t++)
      {
         if (otrees1[t] != NULL)
         {
            for(v= 0; v < nb_values; v++)
               occurrence[v]= false;

            size= otrees1[t]->get_size();
            va= visitor.get_breadthorder(*otrees1[t]);
            index= 0;

            curr_val= _min_value; // current value

            if (!final_value)
            {
               while ((index < size) && (curr_val <= _max_value))
               {
                  v= (otrees1[t]->get(va[index])).Int() - _min_value;
                  if (!occurrence[v])
                  {
                     occurrence[v]= true;
                     (first_occurrence_root[v]->frequency[otrees1[t]->get_depth(va[index])])++;
                     curr_val++;
                  }
                  index++;
               }
            }
            else
            {
               while ((index < size) && (curr_val <= _max_value))
               {
                  v= (otrees1[t]->get(va[index])).Int();
                  if (v != I_DEFAULT)
                  {
                     v-= _min_value;
                     if (!occurrence[v])
                     {
                        occurrence[v]= true;
                        (first_occurrence_root[v]->frequency[otrees1[t]->get_depth(va[index])])++;
                        curr_val++;
                     }
                  }
                  index++;
               }
            }
         }
      }

      for(v= 0; v < nb_values; v++)
      {
         first_occurrence_root[v]->nb_value_computation();
         first_occurrence_root[v]->offset_computation();
         first_occurrence_root[v]->nb_element_computation();
         first_occurrence_root[v]->max_computation();
         first_occurrence_root[v]->mean_computation();
         first_occurrence_root[v]->variance_computation();
      }
      delete [] occurrence;
      occurrence= NULL;
   }
}

/*****************************************************************
 *
 *  Construction of the frequency distribution for the distribution of the
 *  1st occurrence depth for a given value, starting from the leaves
 *
 **/

void TreeCharacteristics::build_first_occurrence_leaves_frequency_distribution(Typed_edge_one_int_tree** otrees1,
                                                                               int imax_depth,
                                                                               bool final_value)
{ //  frequency distribution of first occurrence (leaves)

   typedef Typed_edge_one_int_tree::vertex_iterator vertex_iterator;
   typedef Typed_edge_one_int_tree::vertex_array vertex_array;

   int t, v, curr_val, nb_values= _max_value - _min_value + 1;
   int size, index;
   vertex_array va;
   std::vector<int> depth;
   generic_visitor<Typed_edge_one_int_tree> visitor;
   bool *occurrence;

   if (otrees1 != NULL)
   {
      if (first_occurrence_leaves == NULL)
      {
         first_occurrence_leaves= new FrequencyDistribution*[nb_values];
         for(v= 0; v < nb_values; v++)
            first_occurrence_leaves[v]= NULL;
      }

      for(v= 0; v < nb_values; v++)
         if (first_occurrence_leaves[v] == NULL)
            first_occurrence_leaves[v]= new FrequencyDistribution(imax_depth+1);

      occurrence= new bool[nb_values];

      for(t= 0; t < _nb_trees; t++)
      {
         if (otrees1[t] != NULL)
         {
            for(v= 0; v < nb_values; v++)
               occurrence[v]= false;

            size= otrees1[t]->get_size();
            va= visitor.get_leavesfirstorder(*otrees1[t], depth);
            index= 0;

            curr_val= _min_value; // current value

            if (!final_value)
            {
               while ((index < size) && (curr_val <= _max_value))
               {
                  v= (otrees1[t]->get(va[index])).Int() - _min_value;
                  if (!occurrence[v])
                  {
                     occurrence[v]= true;
                     (first_occurrence_leaves[v]->frequency[depth[index]])++;
                     curr_val++;
                  }
                  index++;
               }
            }
            else
            {
               while ((index < size) && (curr_val <= _max_value))
               {

                  v= (otrees1[t]->get(va[index])).Int();
                  if (v != I_DEFAULT)
                  {
                     v -= _min_value;
                     if (!occurrence[v])
                     {
                        occurrence[v]= true;
                        (first_occurrence_leaves[v]->frequency[depth[index]])++;
                        curr_val++;
                     }
                  }
                  index++;
               }
            }
         }
      }

      for(v= 0; v < nb_values; v++)
      {
         first_occurrence_leaves[v]->nb_value_computation();
         first_occurrence_leaves[v]->offset_computation();
         first_occurrence_leaves[v]->nb_element_computation();
         first_occurrence_leaves[v]->max_computation();
         first_occurrence_leaves[v]->mean_computation();
         first_occurrence_leaves[v]->variance_computation();
      }
      delete [] occurrence;
      occurrence= NULL;
   }
}

/*****************************************************************
 *
 *  Construction of the frequency distributions for the distribution of the
 *  numbers and sizes of the homogeneous zones
 *
 **/

void TreeCharacteristics::build_zone_frequency_distributions(Typed_edge_one_int_tree** otrees1,
                                                             int imax_size,
                                                             bool final_value)
{ //  frequency distribution of homogeneous zones

   typedef Typed_edge_one_int_tree::children_iterator children_iterator;
   typedef Typed_edge_one_int_tree::key key;

   int val, size, nb_values= _max_value - _min_value + 1;
   int t, current_zone, nb_nodes;
   int *nb_zones_t= new int[nb_values];
   key v;
   std::vector<std::deque<key> > zones;
   std::deque<key> node_list, zone_id;
   // set of the homogeneous zones
   std::deque<key> tmp_zone;
   children_iterator it, end;

   if (otrees1 != NULL)
   {
      if (sojourn_size == NULL)
      {
         sojourn_size= new FrequencyDistribution*[nb_values];
         for(val= 0; val < nb_values; val++)
            sojourn_size[val]= NULL;
      }

      for(val= 0; val < nb_values; val++)
         if (sojourn_size[val] == NULL)
            sojourn_size[val]= new FrequencyDistribution(imax_size+1);

      if (nb_zones == NULL)
      {
         nb_zones= new FrequencyDistribution*[nb_values];
         for(val= 0; val < nb_values; val++)
            nb_zones[val]= NULL;
      }

      for(val= 0; val < nb_values; val++)
         if (nb_zones[val] == NULL)
            nb_zones[val]= new FrequencyDistribution(imax_size+1);

      for(t= 0; t < _nb_trees; t++)
      {
         if (otrees1[t] != NULL)
         {
            nb_nodes= 1;

            v= otrees1[t]->root();
            size= otrees1[t]->get_size();
            val= (otrees1[t]->get(v)).Int();
            current_zone= 0;
            tmp_zone.push_back(v);
            zones.push_back(tmp_zone);
            node_list.push_back(v);
            zone_id.push_back(current_zone);
            tmp_zone.pop_back();

            while (!node_list.empty())
            {
               v= node_list.front();
               current_zone= zone_id.front();
               Tree_tie::tie(it, end)= otrees1[t]->children(v);
               //cout << "Traitement du sommet " << v << " de la zone "
               //     << zones[current_zone]
               if (!final_value)
               {
                  while (it < end)
                  {
                     val= (otrees1[t]->get(*it)).Int();
                     if (val == (otrees1[t]->get(v)).Int())
                     // *it belongs to the current homogeneous zone
                     {
                         (zones[current_zone]).push_back(*it);
                         node_list.push_back(*it);
                         zone_id.push_back(current_zone);
                     }
                     else
                     // beginning of a new zone
                     {
                        tmp_zone.push_back(*it);
                        zones.push_back(tmp_zone);
                        tmp_zone.pop_back();
                        node_list.push_back(*it);
                        zone_id.push_back(zones.size()-1);
                     }
                     it++;
                  } // each child handled
               }
               else
               {
                  while (it < end)
                  {
                     val= (otrees1[t]->get(*it)).Int();
                     if (val != I_DEFAULT)
                     {
                        if (val == (otrees1[t]->get(v)).Int())
                        // *it belongs to the current homogeneous zone
                        {
                            (zones[current_zone]).push_back(*it);
                            node_list.push_back(*it);
                            zone_id.push_back(current_zone);
                        }
                        else
                        // beginning of a new zone
                        {
                           tmp_zone.push_back(*it);
                           zones.push_back(tmp_zone);
                           tmp_zone.pop_back();
                           node_list.push_back(*it);
                           zone_id.push_back(zones.size()-1);
                        }
                     }
                     it++;
                  } // each child handled
               }
               node_list.pop_front();
               zone_id.pop_front();
            }

            for(val= 0; val < nb_values; val++)
               nb_zones_t[val]= 0;

            if (!final_value)
            {
               while (!zones.empty())
               {
                  v= (zones.back()).front();
                  val= (otrees1[t]->get(v)).Int();
                  nb_zones_t[val-_min_value]++;
                  // current tree has one more zone of value "val"

                  size= (zones.back()).size();
                  // cout << "Zone de valeur " << val << " enracinee au sommet " << v
                  //      << " de taille " << size << endl;
                  (sojourn_size[val-_min_value]->frequency[size])++;
                  // ... and of size "size'

                  zones.pop_back();
               }
            }
            else
            {
               while (!zones.empty())
               {
                  v= (zones.back()).front();
                  val= (otrees1[t]->get(v)).Int();
                  if (val != I_DEFAULT)
                  {
                     nb_zones_t[val-_min_value]++;
                     // current tree has one more zone of value "val"

                     size= (zones.back()).size();
                     // cout << "Zone de valeur " << val << " enracinee au sommet " << v
                     //      << " de taille " << size << endl;
                     (sojourn_size[val-_min_value]->frequency[size])++;
                     // ... and of size "size'
                  }
                  zones.pop_back();
               }
            }

            for(val= 0; val < nb_values; val++)
               (nb_zones[val]->frequency[nb_zones_t[val]])++;

         }
      } // each tree handled

      for(val= 0; val < nb_values; val++)
      {
         nb_zones[val]->nb_value_computation();
         nb_zones[val]->offset_computation();
         nb_zones[val]->nb_element_computation();
         nb_zones[val]->max_computation();
         nb_zones[val]->mean_computation();
         nb_zones[val]->variance_computation();

         sojourn_size[val]->nb_value_computation();
         sojourn_size[val]->offset_computation();
         sojourn_size[val]->nb_element_computation();
         sojourn_size[val]->max_computation();
         sojourn_size[val]->mean_computation();
         sojourn_size[val]->variance_computation();
      }
   }
   delete [] nb_zones_t;
   nb_zones_t= NULL;
}

/*****************************************************************
 *
 *  Construction of the frequency distributions for the number of occurrences
 *  of a given value
 *
 **/

void TreeCharacteristics::build_nb_occurrences_frequency_distribution(Typed_edge_one_int_tree** otrees1,
                                                                      int imax_size,
                                                                      bool final_value)
{
   typedef Typed_edge_one_int_tree::vertex_iterator vertex_iterator;

   int t, val, nb_values= _max_value - _min_value + 1;
   int *nb_occurrences_t= new int[nb_values];
   // number of occurrences of each value for current tree
   vertex_iterator it, end;

   if (nb_occurrences == NULL)
   {
      nb_occurrences= new FrequencyDistribution*[nb_values];
      for(val= 0; val < nb_values; val++)
         nb_occurrences[val]= NULL;
   }

   for(val= 0; val < nb_values; val++)
      if (nb_occurrences[val] == NULL)
         nb_occurrences[val]= new FrequencyDistribution(imax_size+1);

   for(t = 0; t < _nb_trees; t++)
   {
      if (otrees1[t] != NULL)
      {
         for(val= 0; val < nb_values; val++)
            nb_occurrences_t[val]= 0;

         Tree_tie::tie(it, end)= otrees1[t]->vertices();
         if (!final_value)
            while (it < end)
            {
               val= (otrees1[t]->get(*it++)).Int() - _min_value;
               nb_occurrences_t[val]++;
            }
         else
            while (it < end)
            {
               val= (otrees1[t]->get(*it++)).Int();
               if (val != I_DEFAULT)
               {
                  val -= _min_value;
                  nb_occurrences_t[val]++;
               }
            }

         for(val= 0; val < nb_values; val++)
            (nb_occurrences[val]->frequency[nb_occurrences_t[val]])++;
      }
   }

   for(val= 0; val < nb_values; val++)
   {
      nb_occurrences[val]->nb_value_computation();
      nb_occurrences[val]->offset_computation();
      nb_occurrences[val]->nb_element_computation();
      nb_occurrences[val]->max_computation();
      nb_occurrences[val]->mean_computation();
      nb_occurrences[val]->variance_computation();
   }

   delete [] nb_occurrences_t;
   nb_occurrences_t= NULL;
}

