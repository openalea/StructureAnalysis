/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@cirad.fr)
 *
 *       $Source$
 *       $Id: observed_trees.hpp 3186 2007-05-25 15:10:30Z dufourko $
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
#ifndef OBSERVED_TREES_TCC
#define OBSERVED_TREES_TCC


/*****************************************************************
 *
 *  Default constructor of Int_fl_tree class
 *
 **/

template<class Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>::Int_fl_tree(key root,
                                                   int n,
                                                   const value& default_value)
  : Generic_tree<Generic_Int_fl_container>(root, n, default_value)
{
   _nb_integral= default_value.nb_int();
   _nb_float= default_value.nb_float();
}

/*****************************************************************
 *
 *  Constructor of Int_fl_tree class specifying the number of
 *  float and integral variables
 *
 **/

template<class Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>::Int_fl_tree(int inb_integral,
                                                   int inb_float,
                                                   key root,
                                                   int n)
  : _nb_integral(inb_integral)
  , _nb_float(inb_float)
{
   value v;

   v.reset(_nb_integral, _nb_float);
   Int_fl_tree(root, n, v);
}


/*****************************************************************
 *
 *  Copy constructor of Int_fl_tree class
 *
 **/

template<class Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>::Int_fl_tree(const Int_fl_tree& tree)
  : Generic_tree<Generic_Int_fl_container>(tree)
{
   _nb_integral= tree._nb_integral;
   _nb_float= tree._nb_float;
}

/*****************************************************************
 *
 *  Constructor of Int_fl_tree class using a structure and
 *  a default value
 *
 **/

template<class Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>::Int_fl_tree(Unlabelled_tree& utree, const value& default_value)
  :  Generic_tree<Generic_Int_fl_container>(utree, default_value)
{
   _nb_integral= default_value.nb_int();
   _nb_float= default_value.nb_float();
}

/*****************************************************************
 *
 *  Constructor of Int_fl_tree class using a structure and
 *  the number of float and integral variables
 *
 **/

template<class Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>::Int_fl_tree(const
Generic_tree<Generic_Int_fl_container>& tree)
{
   value v= tree.get(tree.root());

   _nb_integral= v.nb_int();
   _nb_float= v.nb_float();

   copy(tree);
}


/*****************************************************************
 *
 *  Constructor of Int_fl_tree class using a structure and
 *  the number of float and integral variables
 *
 **/

template<class Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>::Int_fl_tree(int inb_integral,
                                                   int inb_float,
                                                   Unlabelled_tree& utree)
{
   Int_fl_tree<Generic_Int_fl_container> *res;
   value v;

   v.reset(inb_integral, inb_float);
   res= new Int_fl_tree<Generic_Int_fl_container>(utree, v);

   Generic_tree<Generic_Int_fl_container>::copy(*res);
   copy(*res);

   delete res;
}

/*****************************************************************
 *
 *  Destructor of Int_fl_tree class
 *
 **/

template<class Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>::~Int_fl_tree()
{}

/*****************************************************************
 *
 *  Assignement operator of Int_fl_tree class
 *
 **/

template<class Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>&
Int_fl_tree<Generic_Int_fl_container>::operator=(const Int_fl_tree<Generic_Int_fl_container>& tree)
{
   if (&tree != this)
   {
      Generic_tree<Generic_Int_fl_container>::copy(tree);
      copy(tree);
   }
   return *this;
}

/*****************************************************************
 *
 *  Adding a vertex to a Int_fl_tree
 *
 **/

template<class Generic_Int_fl_container>
typename Int_fl_tree<Generic_Int_fl_container>::key
Int_fl_tree<Generic_Int_fl_container>::add_vertex(const value& data)
{
    if (this->_size > 0)
       assert((data.nb_int() == this->_nb_integral) && (data.nb_float() == this->_nb_float));
    else
    {
       this->_nb_integral= data.nb_int();
       this->_nb_float= data.nb_float();
    }
    return (Generic_tree<Generic_Int_fl_container>::add_vertex(data));
}
/*****************************************************************
 *
 *  Selection of one single variable for Int_fl_tree
 *
 **/

template<class Generic_Int_fl_container>
One_int_tree* Int_fl_tree<Generic_Int_fl_container>::select_int_variable(int variable)
{

   typedef typename  Int_fl_tree<Generic_Int_fl_container>::vertex_iterator gvertex_iterator;
   typedef One_int_tree::value value;

   Unlabelled_tree* utree;
   One_int_tree* res_tree;
   gvertex_iterator it, end;
   value v;

   // assert(blabla);
   // does Int(variable exist ?)
   utree= this->get_structure();
   res_tree= new One_int_tree(*utree);
   tie(it, end)= this->vertices();
   while ( it != end )
     {
       v.Int()= get(*it).Int(variable);
       res_tree->put(*it, v);
       it++;
     }
   return res_tree;
}

/*****************************************************************
 *
 *  Access to the minimal and maximal values
 *  of integral and float variables
 *
 **/

template<class Generic_Int_fl_container>
Generic_Int_fl_container* Int_fl_tree<Generic_Int_fl_container>::get_min_value()
{
    int var;
    vertex_iterator it, end;
    value v;
    Generic_Int_fl_container *_min_value;

    _min_value= new Generic_Int_fl_container(_nb_integral, _nb_float);
    for(var= 0; var < _nb_integral; var++)
        _min_value->Int(var)= INT_MAX;

    for(var= 0; var < _nb_float; var++)
        _min_value->Double(var)= FLT_MAX;

    tie(it, end)= this->vertices();
    while (it < end)
    {
       v= get(*it);
       for(var= 0; var < this->_nb_integral; var++)
           _min_value->Int(var)= min(_min_value->Int(var),
                                     v.Int(var));

       for(var= 0; var < this->_nb_float; var++)
            _min_value->Double(var)= min(_min_value->Double(var),
                                         v.Double(var));
       it++;
    }
    return _min_value;
}

template<class Generic_Int_fl_container>
Generic_Int_fl_container* Int_fl_tree<Generic_Int_fl_container>::get_max_value()
{
    int var;
    vertex_iterator it, end;
    value v;
    Generic_Int_fl_container *_max_value;

    _max_value= new Generic_Int_fl_container(_nb_integral, _nb_float);
    for(var= 0; var < _nb_integral; var++)
        _max_value->Int(var)= INT_MIN;

    for(var= 0; var < _nb_float; var++)
        _max_value->Double(var)= -FLT_MAX;

    tie(it, end)= this->vertices();
    while (it < end)
    {
       v= get(*it);
       for(var= 0; var < _nb_integral; var++)
           _max_value->Int(var)= max(_max_value->Int(var),
                                     v.Int(var));


       for(var= 0; var < _nb_float; var++)
            _max_value->Double(var)= max(_max_value->Double(var),
                                         v.Double(var));
       it++;
    }
    return _max_value;
}


/*****************************************************************
 *
 *  Random determination of the integral variables for
 *  Int_fl_tree class using a set of Distributions (one for
 *  each integral variable)
 *
 **/

template<class Generic_Int_fl_container>
void Int_fl_tree<Generic_Int_fl_container>::iid_simulation(Distribution ** const distrib)
{
   int var;
   vertex_iterator it, end;
   value v(_nb_integral, 0);

   tie(it, end)= this->vertices();
   while ( it != end )
   {
       v= get(*it);
       for(var= 0; var < _nb_integral; var++)
          v.Int(var)= distrib[var]->simulation();
       put(*it++, v);
   }
}

/*****************************************************************
 *
 *  Access to the number of integral and float variables
 *
 **/

template<class Generic_Int_fl_container>
int Int_fl_tree<Generic_Int_fl_container>::get_nb_int() const
{ return _nb_integral; }

template<class Generic_Int_fl_container>
int Int_fl_tree<Generic_Int_fl_container>::get_nb_float() const
{ return _nb_float; }

/*****************************************************************
 *
 *  Copy operator of Int_fl_tree class
 *
 **/

template<class Generic_Int_fl_container>
void Int_fl_tree<Generic_Int_fl_container>::copy(const Int_fl_tree<Generic_Int_fl_container>& tree)
{
   _nb_integral= tree._nb_integral;
   _nb_float= tree._nb_float;
}

/*****************************************************************
 *
 *  Default constructor of Observed_tree class
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>::Observed_trees(int inb_integral,
                                                         int inb_float,
                                                         int inb_trees)
 : _nb_integral(inb_integral)
 , _nb_float(inb_float)
 , _type(NULL)
 , _max_size(0)
 , _max_depth(0)
 , _nb_trees(inb_trees)
 , hsize(NULL)
 , hnb_children(NULL)
 , trees(NULL)
 , characteristics(NULL)
{
   if (_nb_trees)
      trees= new Int_fl_tree<Generic_Int_fl_container>*[_nb_trees];

   _min_value.reset(_nb_integral, _nb_float);
   _max_value.reset(_nb_integral, _nb_float);
}

/*****************************************************************
 *
 *  Copy constructor of Observed_trees class
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>::Observed_trees(const Observed_trees& otrees)
{ copy(otrees); }

/*****************************************************************
 *
 *  Constructor of Observed_trees class using the number of trees,
 *  the type of each variable, the (multidimensional) observed trees
 *  and a flag on the computation of the histograms
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>::Observed_trees(int inb_trees,
                                                         int const * itype,
                                                         Int_fl_tree<Generic_Int_fl_container>** otrees,
                                                         bool histogram_flag)
 : _nb_integral(0)
 , _nb_float(0)
 , _type(NULL)
 , _max_size(0)
 , _max_depth(0)
 , _nb_trees(inb_trees)
 , hsize(NULL)
 , hnb_children(NULL)
 , trees(NULL)
 , characteristics(NULL)
{
   typedef Int_fl_tree<Generic_Int_fl_container> treetype;
   typedef typename treetype::value value;

   treetype tparadigm;
   // Tree_characteristics tc;
   value vparadigm;
   int var, t, nb_variables;

   if ((otrees != NULL) && (inb_trees > 0))
   {
       tparadigm= *(otrees[0]);
       vparadigm= tparadigm.get(tparadigm.root());
       _nb_integral= vparadigm.nb_int();
       // number of integral variables
       _nb_float= vparadigm.nb_float();
       // number of float variables
       nb_variables= _nb_integral + _nb_float;
       trees= new Int_fl_tree<Generic_Int_fl_container>*[_nb_trees];
       _type= new int[nb_variables];

       for(var= 0; var < _nb_integral; var++)
          _type[var]= itype[var];

       for(var= _nb_integral; var < nb_variables; var++)
          _type[var]= REAL_VALUE;

       for(t= 0; t < _nb_trees; t++)
       {
          assert((otrees[t]->get_nb_int() == _nb_integral)
             && (otrees[t]->get_nb_float() == _nb_float));
          trees[t]= new Int_fl_tree<Generic_Int_fl_container>(*(otrees[t]));
       }

       if (histogram_flag)
       {
          build_characteristics();
          build_size_histogram();
          build_nb_children_histogram();
       }
   }
   else
     Observed_trees<Generic_Int_fl_container>(0);

}

/*****************************************************************
 *
 *  Constructor of Observed_trees class
 *  using the number of integral and float variables,
 *  histograms for the size and number of children of the trees,
 *  a flag on the possibility for a node to have no child due to random
 *  and a flag on the (random) tree initialization
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>::Observed_trees(int inb_integral,
                                                         int inb_float,
                                                         const Histogram& ihsize,
                                                         const Histogram& ihnb_children,
                                                         bool no_child_flag,
                                                         bool init_flag)
 : _nb_integral(inb_integral)
 , _nb_float(inb_float)
 , _type(NULL)
 , _max_size(0)
 , _max_depth(0)
 , _nb_trees(0)
 , hsize(NULL)
 , hnb_children(NULL)
 , trees(NULL)
 , characteristics(NULL)
{
   typedef Int_fl_tree<Generic_Int_fl_container> treetype;
   typedef typename treetype::value value;

   register int i, j, t;
   int *psize, *size;
   value default_value;
   Unlabelled_tree utree;

   if (_nb_integral+_nb_float > 0)
   {
      _min_value.reset(_nb_integral, _nb_float);
      _max_value.reset(_nb_integral, _nb_float);
      default_value.reset(_nb_integral, _nb_float);

      _type= new int[_nb_integral+_nb_float];

      for(i=0; i < _nb_integral; i++)
      {
         _type[i]= INT_VALUE;
         default_value.Int(i)= 0;
      }

      for(i=0; i < _nb_float; i++)
      {
         _type[_nb_integral+i]= REAL_VALUE;
         default_value.Double(i)= 0.;
      }
      _nb_trees= ihsize.nb_element;


      trees= new Int_fl_tree<Generic_Int_fl_container>*[_nb_trees];

      hnb_children= new Histogram(ihnb_children);

      if (init_flag)
      {
         size= new int[_nb_trees];
         psize= size;
         // array of the size of each tree
         for(i= ihsize.offset; i < ihsize.nb_value; i++)
            for(j= 0; j < ihsize.frequency[i]; j++)
               *psize++ = i;

         if (!no_child_flag)
            hnb_children->frequency[0]= 0;

         for(t= 0; t < _nb_trees; t++)
         {
            trees[t]= new Int_fl_tree<Generic_Int_fl_container>();
            utree.simulation(*hnb_children, size[t], size[t]);
            trees[t]->set_structure(utree, default_value);
         }
         build_characteristics();
         build_size_histogram();
         build_nb_children_histogram();
      }
      else
      {
         for(t= 0; t < _nb_trees; t++)
            trees[t]= NULL;
         hsize= new Histogram(ihsize);
      }
   }
}

/*****************************************************************
 *
 *  Constructor of Observed_trees class by selecting some trees,
 *  whose number and list are given in argument
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>::Observed_trees(const Observed_trees& otrees,
                                                         int inb_trees,
                                                         int_array index)
 : _nb_integral(0)
 , _nb_float(0)
 , _type(NULL)
 , _max_size(0)
 , _max_depth(0)
 , _nb_trees(inb_trees)
 , hsize(NULL)
 , hnb_children(NULL)
 , trees(NULL)
 , characteristics(NULL)
{
   typedef Int_fl_tree<Generic_Int_fl_container> treetype;
   typedef typename treetype::value value;

   treetype tparadigm;
   // Tree_characteristics tc;
   value vparadigm;
   int var, t, nb_variables;

   if ((index != NULL) && (inb_trees > 0))
   {
       _nb_integral= otrees._nb_integral;
       _nb_float= otrees._nb_float;

       nb_variables= _nb_integral + _nb_float;
       trees= new Int_fl_tree<Generic_Int_fl_container>*[_nb_trees];
       _type= new int[nb_variables];

       for(var= 0; var < _nb_integral; var++)
          _type[var]= otrees._type[var];

       for(var= _nb_integral; var < nb_variables; var++)
          _type[var]= REAL_VALUE;

       for(t= 0; t < _nb_trees; t++)
          trees[t]= new Int_fl_tree<Generic_Int_fl_container>(*(otrees.trees[index[t]]));

       build_characteristics();
       build_size_histogram();
       build_nb_children_histogram();
   }
   else
     Observed_trees<Generic_Int_fl_container>(0);
}

/*****************************************************************
 *
 *  Destructor of Observed_trees class
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>::~Observed_trees()
{ remove(); }

/*****************************************************************
 *
 *  Assignement operator of Observed_trees class
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>&
Observed_trees<Generic_Int_fl_container>::operator=(const Observed_trees<Generic_Int_fl_container> &otrees)
{ // assignement operator
  if (&otrees != this)
  {
     remove();
     copy(otrees);
  }
}

/*****************************************************************
 *
 *  Extract marginal distribution histogram from Observed_trees
 *  using a Format_error and the considered integral variable
 *
 **/

template<typename Generic_Int_fl_container>
Distribution_data* Observed_trees<Generic_Int_fl_container>::extract(Format_error &error,
                                                                     int variable) const
{
  bool status= true;
  Distribution_data *histo= NULL;

  error.init();

  if ((variable < 1) || (variable > _nb_integral))
  {
     status= false;
     error.update(STAT_error[STATR_VARIABLE_INDEX]);
  }
  else
  {
     variable--;

     if ((_type[variable] != INT_VALUE) && (_type[variable] != STATE))
     {
        status= false;
        ostringstream error_message , correction_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable+1 << ": "
                      << STAT_TREES_error[STATR_VARIABLE_TYPE];
        correction_message << STAT_TREES_word[INT_VALUE] << " or " << STAT_TREES_word[STATE];
        error.correction_update((error_message.str()).c_str() , (correction_message.str()).c_str());
     }
     else if (characteristics[variable] == NULL)
     {
        status= false;
        ostringstream error_message;
        error_message << STAT_label[STATL_VARIABLE] << " " << variable+1 << ": "
                      << STAT_error[STATR_MARGINAL_HISTOGRAM];
        error.update((error_message.str()).c_str());
     }
  }

  if (status)
     histo= new Distribution_data(*(characteristics[variable]->marginal));

  return histo;

}

/*****************************************************************
 *
 *  Extract a characteristic distribution histogram from Observed_trees
 *  using a Format_error, the characteristic type,
 *  the considered variable and value
 *
 **/

template<typename Generic_Int_fl_container> Distribution_data*
Observed_trees<Generic_Int_fl_container>::extract(Format_error &error, int type,
                                                  int variable, int value) const
{
   bool status= true;
   // Distribution *pdist= NULL;
   Histogram *phisto= NULL;
   Distribution_data *histo= NULL;

   error.init();

   if ((variable < 1) || (variable > _nb_integral))
   {
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }

   if ((status) &&
         ((characteristics != NULL) && (characteristics[variable-1] != NULL)))
   {
      variable--;

      if ((value < 0) || (value >= characteristics[variable]->get_nb_values()))
      {
         status= false;
         ostringstream error_message;
         error_message << STAT_label[variable == 0 ? STATL_STATE : STATL_OUTPUT] << " "
                       << value << " " << STAT_TREES_error[STATR_NOT_PRESENT];
         error.update((error_message.str()).c_str());
      }

      else
      {
         switch (type)
         {
            case FIRST_OCCURRENCE_ROOT :
            {
               if (characteristics[variable]->first_occurrence_root != NULL)
                  phisto= characteristics[variable]->first_occurrence_root[value];
               else
                  status= false;
               break;
            }
            case FIRST_OCCURRENCE_LEAVES :
            {
               if (characteristics[variable]->first_occurrence_leaves != NULL)
                  phisto= characteristics[variable]->first_occurrence_leaves[value];
               else
                  status= false;
               break;
            }
            case SOJOURN_SIZE :
            {
               if (characteristics[variable]->sojourn_size != NULL)
                  phisto= characteristics[variable]->sojourn_size[value];
               else
                  status= false;
               break;
            }
            case NB_ZONES :
            {
               if (characteristics[variable]->nb_zones != NULL)
                  phisto = characteristics[variable]->nb_zones[value];
               else
                  status= false;
               break;
            }
            case NB_OCCURRENCES :
            {
               if (characteristics[variable]->nb_occurrences != NULL)
                  phisto = characteristics[variable]->nb_occurrences[value];
               else
                  status= false;
               break;
            }
         }

         if (!status)
            error.update(STAT_TREES_error[STATR_CHARACTERISTICS_NOT_COMPUTED]);
         else
         {
            if (phisto->nb_element == 0)
            {
               status= false;
               error.update(STAT_error[STATR_EMPTY_HISTOGRAM]);
            }
            else
               histo= new Distribution_data(*phisto);
         }
      }
   }
   else // the caracteristics have not been computed
   {
      status= false;
      error.update(STAT_TREES_error[STATR_CHARACTERISTICS_NOT_COMPUTED]);
   }

   return histo;
}


/*****************************************************************
 *
 *  Extract sequences from Observed_trees, using a Format_error object
 *  and a flag on considering all possible sequences from root
 *  to leaf nodes or on avoiding redundancy
 *
 **/

template<typename Generic_Int_fl_container> Sequences*
Observed_trees<Generic_Int_fl_container>::extract_sequences(Format_error& error,
                                                            bool all_paths) const
{
   // only the case all_paths == True is implemented
   typedef typename generic_visitor<tree_type>::vertex_array vertex_array;

   Sequences* res= NULL;
   int *itype= NULL, *ilength= NULL;
   int ***isequence= NULL; // isequence[identifier][variable][index]
   int *iidentifier= NULL;
   int var, nb_sequences= 0, t, s, pos, sequence_id= 0;
   bool status= true;
   ostringstream error_message;
   std::vector<vertex_array> leaf_set;
   vertex_array leaves;
   generic_visitor<tree_type> visitor;
   key current_leaf, current_vertex;
   value val;

   error.init();

   if (_nb_integral == 0)
   {
      status= false;
      error_message << STAT_TREES_error[STATR_NB_INT_OUTPUT_PROCESS] << " ";

      error.correction_update((error_message.str()).c_str(), _nb_integral);
      error_message.str("");
   }
   if (status)
   {
      if (trees == NULL)
      {
         status= false;
         error.update(STAT_error[STATR_EMPTY_SAMPLE]);
      }
   }
   if (status)
   {
      itype= new int[_nb_integral];
      for (var= 0; var < _nb_integral; var++)
         itype[var]= _type[var];

      for (t= 0; t < _nb_trees; t++)
      {
         leaves.clear();
         // determine the set of leaves for each tree
         visitor.find_leaves(*(trees[t]), trees[t]->root(), leaves);
         leaf_set.push_back(leaves);
         nb_sequences+= leaves.size();
      }

      isequence= new int**[nb_sequences];
      ilength= new int[nb_sequences];
      iidentifier= new int[nb_sequences];
      for (s= 0; s < nb_sequences; s++)
         isequence[s]= new int*[_nb_integral];

      for (t= 0; t < _nb_trees; t++)
      {
         while(!leaf_set[t].empty())
         {
            current_leaf= leaf_set[t].back();
            leaf_set[t].pop_back();
            ilength[sequence_id]= trees[t]->get_depth(current_leaf)+1;
            iidentifier[sequence_id]= sequence_id;
            for (var= 0; var < _nb_integral; var++)
               isequence[sequence_id][var]= new int[ilength[sequence_id]];

            current_vertex= current_leaf;
            for (pos= ilength[sequence_id]-1; pos > -1; pos--)
            {
               val= trees[t]->get(current_vertex);
               for (var= 0; var < _nb_integral; var++)
                  isequence[sequence_id][var][pos]= val.Int(var);
               current_vertex= trees[t]->parent(current_vertex);
            }
            sequence_id++;
         }
      }
      res= new Sequences(_nb_integral, itype,  nb_sequences, ilength,
                         isequence, iidentifier);

      for (s= 0; s < nb_sequences; s++)
      {
         for (var= 0; var < _nb_integral; var++)
         {
            delete [] isequence[s][var];
            isequence[s][var]= NULL;
         }

         delete[] isequence[s];
         isequence[s]= NULL;
      }
      delete [] isequence;
      isequence= NULL;

      delete [] itype;
      itype= NULL;

      delete [] ilength;
      ilength =NULL;

      delete [] iidentifier;
      iidentifier =NULL;
   }


   return res;
}

/*****************************************************************
 *
 *  Merge Observed_trees, using a Format_error object,
 *  the number of Observed_trees objects and a collection
 *  of Observed_trees objects.
 *  The collection otrees of Observed_trees objects is
 *  merged with *this to create a new object
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::merge(Format_error& error,
                                                int nb_sample,
                                                const pt_observed_trees_array otrees) const
{
   bool status= true;
   register int o, var, t;
   int inb_trees, nb_variables= _nb_integral+_nb_float, tree_id;
   Observed_trees<Generic_Int_fl_container> *res= NULL;
   ostringstream error_message;
   pt_tree_type_array ptta;

   error.init();

   for(o= 0; o < nb_sample; o++)
   {
      if ((otrees[o]->_nb_integral+otrees[o]->_nb_float) !=
         (_nb_integral+_nb_float))
      {
         status= false;
         error_message << STAT_label[STATL_SAMPLE] << " " << o+2 << ": "
                       << STAT_error[STATR_NB_VARIABLE] << " ";

         error.correction_update((error_message.str()).c_str(), nb_variables);
         error_message.str("");
      }
      else
      {
         for(var= 0; var < nb_variables; var++)
            if (otrees[o]->_type[var] != _type[var])
            {
               status= false;
               error_message << STAT_label[STATL_SAMPLE] << " " << o+2 << ": "
                             << STAT_label[STATL_VARIABLE] << " " << var+1 << ": "
                             << STAT_TREES_error[STATR_VARIABLE_TYPE];

               error.correction_update((error_message.str()).c_str(), STAT_TREES_type[_type[var]]);
               error_message.str("");
           }
      }
   }

   if (status)
   {
      // computation of the number of trees
      inb_trees= _nb_trees;

      for(o= 0; o < nb_sample; o++)
         inb_trees+= otrees[o]->_nb_trees;


      ptta= new tree_type*[inb_trees];
      tree_id= 0;

      for(t= 0; t < _nb_trees; t++)
         ptta[tree_id++]= get_tree(t);

      for(o= 0; o < nb_sample; o++)
         for(t= 0; t < otrees[o]->_nb_trees; t++)
            ptta[tree_id++]= otrees[o]->get_tree(t);

      res= new Observed_trees(inb_trees, _type, ptta);

      for(t= 0; t < inb_trees; t++)
      {
         if (ptta[t] != NULL);
            delete ptta[t];
         ptta[t]= NULL;
      }
      delete [] ptta;
      ptta= NULL;
   }

   return res;
}

/*****************************************************************
 *
 *  Translation of the values of a given integral variable
 *  for Observed_trees, using a Format_error object,
 *  the variable index and the translation parameter
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::shift(Format_error& error,
                                                int variable,
                                                int shift_param) const
{
   const int nb_variables= _nb_integral + _nb_integral;
   bool status= true;
   register int t;
   Observed_trees<Generic_Int_fl_container> *otrees= NULL;
   vertex_iterator it, end;
   value v;

   error.init();

   if ((variable < 1) && (variable > nb_variables))
   {
      // variable index must be valid
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }
   else
   {
      variable--;

      if ((_type[variable] != INT_VALUE) && (_type[variable] != TIME)
          && (_type[variable] != STATE))
      {
         status= false;
         ostringstream correction_message;
         correction_message << STAT_TREES_type[INT_VALUE]
                            << " or " << STAT_TREES_type[TIME]
                            << " or " << STAT_TREES_type[STATE];
         error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
      }
      else
      {
         if (shift_param+_min_value.Int(variable) < INT_MIN)
         {
            status= false;
            ostringstream correction_message;
            correction_message << STAT_error[STATR_GREATER_THAN] << " "
                               << INT_MIN-_min_value.Int(variable);
            error.correction_update(STAT_error[STATR_SHIFT_VALUE], (correction_message.str()).c_str());
         }
         if (shift_param+_max_value.Int(variable) > INT_MAX)
         {
            status = false;
            ostringstream correction_message;
            correction_message << STAT_error[STATR_SMALLER_THAN] << " "
                               << INT_MAX-_max_value.Int(variable);
            error.correction_update(STAT_error[STATR_SHIFT_VALUE], (correction_message.str()).c_str());
         }
      }
   }

   if (status)
   {
      otrees= new Observed_trees(_nb_trees, _type, trees);

      for(t= 0; t < _nb_trees; t++)
      {
         tie(it, end)= (otrees->trees[t])->vertices();
         while (it < end)
         {
            // translation of the variable
            v= (otrees->trees[t])->get(*it);
            v.Int(variable)+= shift_param;
            (otrees->trees[t])->put(*it++, v);
         }
      }
      otrees->_max_value.Int(variable)+= shift_param;
      otrees->_min_value.Int(variable)+= shift_param;
      otrees->build_characteristics(variable);
   }
   return otrees;
}

/*****************************************************************
 *
 *  Translation of the values of a given floating variable
 *  for Observed_trees, using a Format_error object,
 *  the variable index and the translation parameter
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::shift(Format_error& error,
                                                int variable,
                                                double shift_param) const
{
   const int nb_variables= _nb_integral + _nb_float;
   bool status= true;
   register int t, cvariable= 0;
   Observed_trees<Generic_Int_fl_container> *otrees= NULL;
   vertex_iterator it, end;
   value v;

   error.init();

   if ((variable < 1) && (variable > nb_variables))
   {
      // variable index must be valid
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }
   else
   {
      variable--;
      cvariable= variable - _nb_integral;
      // position of the variable among the floating ones

      if (_type[variable] != REAL_VALUE)
      {
         status= false;
         ostringstream correction_message;
         correction_message << STAT_TREES_type[REAL_VALUE];
         error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
      }
      else
      {
         if (shift_param+_min_value.Double(cvariable) < -FLT_MAX)
         {
            status= false;
            ostringstream correction_message;
            correction_message << STAT_error[STATR_GREATER_THAN] << " "
                               << -FLT_MAX-_min_value.Double(cvariable);
            error.correction_update(STAT_error[STATR_SHIFT_VALUE], (correction_message.str()).c_str());
         }
         if (shift_param+_max_value.Double(cvariable) > FLT_MAX)
         {
            status= false;
            ostringstream correction_message;
            correction_message << STAT_error[STATR_SMALLER_THAN] << " "
                               << FLT_MAX-_max_value.Double(cvariable);
            error.correction_update(STAT_error[STATR_SHIFT_VALUE], (correction_message.str()).c_str());
         }
      }
   }

   if (status)
   {
      otrees= new Observed_trees(_nb_trees, _type, trees);

      for(t= 0; t < _nb_trees; t++)
      {
         tie(it, end)= (otrees->trees[t])->vertices();
         while (it < end)
         {
            // translation of the variable
            v= (otrees->trees[t])->get(*it);
            v.Double(cvariable)+= shift_param;
            (otrees->trees[t])->put(*it++, v);
         }
      }
      otrees->_max_value.Double(cvariable)+= shift_param;
      otrees->_min_value.Double(cvariable)+= shift_param;
   }
   return otrees;
}

/*****************************************************************
 *
 *  Cluster the values of a given integral variable
 *  for Observed_trees, using a Format_error object,
 *  the variable index and the clustering step
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::cluster(Format_error &error,
                                                  int variable,
                                                  int step) const
{
   // this function should be extended to floating values
   bool status= true;
   Observed_trees<Generic_Int_fl_container> *otrees= NULL;

   error.init();

   if ((variable < 1) || (variable > _nb_integral))
   {
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }
   else
   {
      variable--;

      if ((_type[variable] != INT_VALUE) && (_type[variable] != STATE))
      {
         status= false;
         ostringstream correction_message;
         correction_message << STAT_TREES_type[INT_VALUE] << " or " << STAT_TREES_type[STATE];
         error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE], (correction_message.str()).c_str());
      }
   }

   if (step < 1)
   {
      status= false;
      error.update(STAT_error[STATR_CLUSTERING_STEP]);
   }

   if (status)
   {
      otrees= new observed_trees(_nb_trees, _type, trees);
      otrees->cluster(*this, variable, step);
   }

   return otrees;
}

/*****************************************************************
 *
 *  Transcoding of the values of a given integral variable
 *  for Observed_trees, using a Format_error object,
 *  the variable index and a transcoding table
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::transcode(Format_error& error,
                                                    int variable,
                                                    int_array symbol) const
{
   bool status= true, *presence;
   register int i;
   int min_symbol, max_symbol;
   observed_trees *otrees= NULL;

   error.init();

   if ((variable < 1) || (variable > _nb_integral))
   {
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }
   else
   {
      variable--;

      if ((_type[variable] != INT_VALUE) && (_type[variable] != STATE))
      {
         status= false;
         ostringstream correction_message;
         correction_message << STAT_TREES_type[INT_VALUE] << " or " << STAT_TREES_type[STATE];
         error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
      }
      else
      {
         min_symbol= INT_MAX;
         max_symbol= INT_MIN;

         for(i= 0; i <= _max_value.Int(variable)-_min_value.Int(variable); i++)
         {
            if (symbol[i] < min_symbol)
               min_symbol= symbol[i];

            if (symbol[i] > max_symbol)
               max_symbol= symbol[i];
         }

         if (max_symbol-min_symbol == 0)
         {
            status= false;
            error.update(STAT_error[STATR_NB_SYMBOL]);
         }

         if (max_symbol-min_symbol > _max_value.Int(variable)-_min_value.Int(variable))
         {
            status = false;
            error.update(STAT_error[STATR_NON_CONSECUTIVE_SYMBOLS]);
         }
      }

      if (status)
      {
         presence= new bool[max_symbol-min_symbol+1];
         for(i= 0; i <= max_symbol-min_symbol; i++)
            presence[i]= false;

         for(i= 0; i <= _max_value.Int(variable)-_min_value.Int(variable); i++)
            presence[symbol[i]-min_symbol]= true;

         for (i= 0; i <= max_symbol-min_symbol; i++)
            if (!presence[i])
            {
               status= false;
               ostringstream error_message;
               error_message << STAT_error[STATR_MISSING_SYMBOL] << " " << i+min_symbol;
               error.update((error_message.str()).c_str());
            }
         delete [] presence;
      }

      if (status)
      {
         for (i= 0; i <= _max_value.Int(variable)-_min_value.Int(variable); i++)
            symbol[i]-= min_symbol;

         otrees= new observed_trees(_nb_trees, _type, trees);
         otrees->transcode(*this, variable, min_symbol, max_symbol, symbol);
      }
   }
   return otrees;
}

/*****************************************************************
 *
 *  Cluster the values (e.g. symbols) of a given integral
 *  variable for Observed_trees, using a Format_error object,
 *  the variable index, the number of clusters and the cluster limits.
 *  If the size of the array exceeds the number of clusters,
 *  the exceeding clusters are ignored (i.e. the array is truncated).
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::cluster(Format_error& error,
                                                  int variable,
                                                  int nb_class,
                                                  int_array ilimit) const
{
   bool status= true;
   register int i, j, k;
   int *limit, *symbol;
   observed_trees *otrees= NULL;

   error.init();

   if ((variable < 1) || (variable > _nb_integral))
   {
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }
   else
   {
      variable--;

      if ((_type[variable] != INT_VALUE) && (_type[variable] != STATE))
      {
         status= false;
         ostringstream correction_message;
         correction_message << STAT_TREES_type[INT_VALUE] << " or " << STAT_TREES_type[STATE];
         error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE], (correction_message.str()).c_str());
      }
      else
         if ((nb_class < 2) || (nb_class >= _max_value.Int(variable)-_min_value.Int(variable)+1))
         {
            status= false;
            error.update(STAT_error[STATR_NB_CLASS]);
         }
   }

   if (status)
   {
      limit= new int[nb_class+1];
      limit[0]= _min_value.Int(variable);
      for(i= 1; i < nb_class; i++)
         limit[i]= ilimit[i-1];

      limit[nb_class]= _max_value.Int(variable) + 1;

      for(i= 0; i < nb_class; i++)
         if (limit[i] >= limit[i+1])
         {
            status= false;
            error.update(STAT_error[STATR_CLUSTER_LIMIT]);
         }

      if (status)
      {
         symbol= new int[_max_value.Int(variable)-_min_value.Int(variable)+1];

         i= 0;
         for(j= 0; j < nb_class; j++)
            for(k= limit[j]; k < limit[j+1]; k++)
               symbol[i++]= j;

         otrees= new observed_trees(_nb_trees, _type, trees);
         otrees->transcode(*this, variable, 0, nb_class-1, symbol);

         delete [] symbol;
      }
      delete [] limit;
   }
   return otrees;
}

/*****************************************************************
 *
 *  Cluster the values (e.g. symbols) of a given integral
 *  variable for Observed_trees, using a Format_error object,
 *  the variable index, the number of clusters and the cluster limits
 *
 **/
/*
template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::cluster(Format_error& error,
                                                  int variable,
                                                  int nb_class,
                                                  double_array ilimit) const
{
   bool status= true;
   register int i, j, k;
   int *symbol;
   double_array limit;
   observed_trees *otrees= NULL;

   error.init();

   if ((variable < 1) || (variable > _nb_float))
   {
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }
   else
   {
      variable--;

      if (_type[variable] != REAL_VALUE)
      {
         status= false;
         ostringstream correction_message;
         correction_message << STAT_TREES_type[REAL_VALUE];
         error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE], (correction_message.str()).c_str());
      }
      else
         if (nb_class < 2)
         {
            status= false;
            error.update(STAT_error[STATR_NB_CLASS]);
         }
   }

   if (status)
   {
      limit= new int[nb_class+1];
      limit[0]= _min_value.Double(variable);
      for(i= 1; i < nb_class; i++)
         limit[i]= ilimit[i-1];

      limit[nb_class]= _max_value.Double(variable) + 1;

      for(i= 0; i < nb_class; i++)
         if (limit[i] >= limit[i+1])
         {
            status= false;
            error.update(STAT_error[STATR_CLUSTER_LIMIT]);
         }

      if (status)
      {
         symbol= new int[_max_value.Int(variable)-_min_value.Int(variable)+1];

         i= 0;
         for(j= 0; j < nb_class; j++)
            for(k= limit[j]; k < limit[j+1]; k++)
               symbol[i++]= j;

         otrees= new observed_trees(_nb_trees, _type, trees);
         otrees->transcode(*this, variable, 0, nb_class-1, symbol);

         delete [] symbol;
      }
      delete [] limit;
   }
   return otrees;
}
*/

/*****************************************************************
 *
 *  Selection of Observed_trees on values taken by a given integral
 *  variable, using a Format_error object, the variable index,
 *  bounds on the values and a flag on keeping or rejecting
 *  the concerned trees. Trees having one node at least
 *  whose value is between the bounds are selected.
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::value_select(Format_error& error,
                                                       int variable,
                                                       int imin_value,
                                                       int imax_value,
                                                       bool keep) const
{
   bool status= true;
   register int t;
   int inb_trees, val;
   int_array index;
   observed_trees *otrees= NULL;
   vertex_iterator it, end;

   error.init();

   if ((variable < 1) || (variable > _nb_integral))
   {
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }
   if (_type[variable] == REAL_VALUE)
   {
      status= false;
      ostringstream correction_message;
      correction_message << STAT_TREES_type[REAL_VALUE];
      error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
   }
   else
   {
      variable--;

      if (_type[variable] != POSITION)
      {
         if ((imin_value > _max_value.Int(variable))
               || (imin_value > imax_value))
         {
            status= false;
            error.update(STAT_error[STATR_MIN_VALUE]);
         }
         if ((imax_value < _min_value.Int(variable))
               || (imax_value < imin_value)) // this second condition seems weird
         {
            status= false;
            error.update(STAT_error[STATR_MAX_VALUE]);
         }
      }
   }

   if (status)
   {
      // tree selection
      index= new int[_nb_trees];
      inb_trees= 0;

      for(t= 0; t < _nb_trees; t++)
      {
         tie(it, end)= trees[t]->vertices();
         while (it < end)
         {
            val= (trees[t]->get(*it)).Int(variable);
            if ((val >= imin_value) && (val <= imax_value))
            {
               if (keep)
                  index[inb_trees++]= t;
               break;
            }
            it++;
         }
         if ((!keep) && (*it == (int)trees[t]->get_size()))
            index[inb_trees++]= t;
      }

      if (inb_trees == 0)
      {
         status= false;
         error.update(STAT_error[STATR_EMPTY_SAMPLE]);
      }

      // copy of the trees
      if (status)
         otrees= new observed_trees(*this, inb_trees, index);

      delete [] index;
   }
   return otrees;
}

/*****************************************************************
 *
 *  Selection of Observed_trees on values taken by a given floating
 *  variable, using a Format_error object, the variable index,
 *  bounds on the values and a flag on keeping or rejecting
 *  the concerned trees. Trees having one node at least
 *  whose value is between the bounds are selected.
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::value_select(Format_error& error,
                                                       int variable,
                                                       double dmin_value,
                                                       double dmax_value,
                                                       bool keep) const
{
   bool status= true;
   register int t;
   int inb_trees;
   double val;
   int_array index;
   observed_trees *otrees= NULL;
   vertex_iterator it, end;

   error.init();

   if ((variable < 1) || (variable > _nb_float))
   {
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }
   if (_type[variable] != REAL_VALUE)
   {
      status= false;
      ostringstream correction_message;
      correction_message << STAT_TREES_type[INT_VALUE];
      error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE], (correction_message.str()).c_str());
   }
   else
   {
      variable--;

      if ((dmin_value > _max_value.Double(variable))
            || (dmin_value > dmax_value))
      {
         status= false;
         error.update(STAT_error[STATR_MIN_VALUE]);
      }
      if ((dmax_value < _min_value.Double(variable))
            || (dmax_value < dmin_value)) // this second condition seems weird
      {
         status= false;
         error.update(STAT_error[STATR_MAX_VALUE]);
      }
   }

   if (status)
   {
      // tree selection
      index= new int[_nb_trees];
      inb_trees= 0;

      for(t= 0; t < _nb_trees; t++)
      {
         tie(it, end)= trees[t]->vertices();
         while (it < end)
         {
            val= (trees[t]->get(*it)).Double(variable);
            if ((val >= this->imin_value) && (val <= this->imax_value))
            {
               if (keep)
                  index[inb_trees++]= t;
               break;
            }
            it++;
         }
         if ((!keep) && (*it == trees[t]->get_size()))
            index[inb_trees++]= t;
      }

      if (inb_trees == 0)
      {
         status= false;
         error.update(STAT_error[STATR_EMPTY_SAMPLE]);
      }

      // copy of the trees
      if (status)
         otrees= new observed_trees(*this, inb_trees, index);

      delete [] index;
   }
   return otrees;
}

/*****************************************************************
 *
 *  Observed_trees selection by identifier
 *  using a Format_error object, the number of trees, a tree identifier
 *  and a flag on keeping or rejecting the selected trees
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::select_individual(Format_error& error,
                                                            int inb_trees,
                                                            int_array iidentifier,
                                                            bool keep) const

{
   bool status= true , *selected_trees;
   register int t, j;
   int max_identifier;
   int_array index, identifier;
   observed_trees *otrees= NULL;

   error.init();

   if ((inb_trees < 1) || (inb_trees > (keep ? _nb_trees : _nb_trees-1)))
   {
      status= false;
      error.update(STAT_TREES_error[STATR_NB_TREES]);
   }
   else
   {
      max_identifier= 1;
      for(t= 0; t < inb_trees; t++)
         if (iidentifier[t] > max_identifier)
            max_identifier= iidentifier[t];

      selected_trees= new bool[max_identifier+1];
      identifier= new int[max_identifier+1]; // until identifier is properly defined
      for(t= 0; t <= max_identifier; t++)
      {
         selected_trees[t]= false;
         identifier[t]= t; // until identifier is properly defined
      }

      for(t= 0; t < inb_trees; t++)
      {
         for(j= 0; j < _nb_trees; j++)
             if (iidentifier[t] == identifier[j])
                break;

         if (j == _nb_trees)
         {
            status= false;
            ostringstream error_message;
            error_message << iidentifier[t] << ": " << STAT_TREES_error[STATR_TREE_IDENTIFIER];
            error.update((error_message.str()).c_str());
         }
         else
            if (selected_trees[iidentifier[t]])
            {
               status= false;
               ostringstream error_message;
               error_message << STAT_TREES_label[STATL_TREE] << " " << iidentifier[t] << " "
                             << STAT_error[STATR_ALREADY_SELECTED];
               error.update((error_message.str()).c_str());
            }
            else
               selected_trees[iidentifier[t]]= true;
      }
      delete [] selected_trees;
   }

   if (status)
   {
      index= identifier_select(_nb_trees, identifier, inb_trees, iidentifier, keep);
      otrees= new observed_trees(*this, (keep ? inb_trees : _nb_trees-inb_trees), index);
      delete [] index;
   }
   return otrees;
}

/*****************************************************************
 *
 *  Variables selection for Observed_trees
 *  using a Format_error object, the number of variables,
 *  the array of their indices and a flag on keeping or rejecting
 *  the selected variables
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::select_variable(Format_error& error,
                                                          int inb_variable,
                                                          int_array ivariable,
                                                          bool keep) const
{
   bool status= true, *selected_variable;
   register int var;
   int bnb_variable, nb_variable= _nb_integral+_nb_float,
       nb_integral= 0, nb_float= 0;
   int_array variable, itype;
   observed_trees *otrees= NULL;

   error.init();

   if ((inb_variable < 1) || (inb_variable > (keep ? nb_variable : nb_variable - 1)))
   {
      status= false;
      error.update(STAT_error[STATR_NB_SELECTED_VARIABLE]);
   }

   if (status)
   {
      selected_variable= new bool[nb_variable+1];
      for(var= 0; var < nb_variable; var++)
         selected_variable[var+1]= false;

      for(var= 0; var < inb_variable; var++)
      {
        // checking each variable designed by the array ivariable
        if ((ivariable[var] < 1) || (ivariable[var] > nb_variable))
        {
           status= false;
           ostringstream error_message;
           error_message << ivariable[var] << ": " << STAT_error[STATR_VARIABLE_INDEX];
           error.update((error_message.str()).c_str());
        }
        else
           if (selected_variable[ivariable[var]])
           {
              status= false;
              ostringstream error_message;
              error_message << STAT_label[STATL_VARIABLE] << " " << ivariable[var] << " "
                            << STAT_error[STATR_ALREADY_SELECTED];
              error.update((error_message.str()).c_str());
           }
           else
              selected_variable[ivariable[var]]= true;
      }
      delete [] selected_variable;
   }

   if (status)
   {
      variable= ::select_variable(nb_variable, inb_variable, ivariable, keep);

      bnb_variable= (keep ? inb_variable : nb_variable-inb_variable);

      itype= new int[bnb_variable];
      for(var= 0; var < bnb_variable; var++)
      {
         itype[var]= _type[variable[var]];
         if (itype[var] == REAL_VALUE)
            nb_float++;
         else
            nb_integral++;
      }

      if ((bnb_variable == 1)
            && ((itype[0] == TIME) || (itype[0] == POSITION)))
      {
         status= false;
         error.update(STAT_error[STATR_NB_SELECTED_VARIABLE]);
      }

      if (status)
      {
         // otrees= new observed_trees(nb_integral, nb_float, _type,
         //                            _nb_trees, false);
         otrees= new observed_trees(nb_integral, nb_float, itype,
                                    _nb_trees, false);
         // should include identifier field
         otrees->select_variable(*this, variable);
      }
      delete [] variable;
      delete [] itype;
   }
   return otrees;
}

/*****************************************************************
 *
 *  Variable concatenation for Observed_trees
 *  using a Format_error object, the number of Observed_trees objects,
 *  an array of pointers on these objects and a reference sample
 *  for the tree identifiers
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::merge_variable(Format_error& error,
                                                         int nb_sample,
                                                         const pt_observed_trees_array otrees,
                                                         int ref_sample) const
{
   bool status= true;
   ostringstream error_message, correction_message;
   register int s, var, t;
   int inb_integral, inb_float;
   // int_array iidentifier;
   observed_trees *res= NULL;
   const observed_trees **ptrees;
   value vsource, vdest;
   vertex_iterator it, end;
   children_iterator chs_it, chs_end, chd_it, chd_end;
   tree_type *cstree, *cdtree; // current source and destination trees
   Unlabelled_tree *utree= NULL;

   error.init();

   if ((_type[0] != INT_VALUE)
         && (_type[0] != REAL_VALUE) && (_type[0] != STATE))
   {
      status= false;
      error_message << STAT_label[STATL_SAMPLE] << " " << 1 << ": "
                    << STAT_TREES_error[STATR_VARIABLE_1_TYPE];
      correction_message << STAT_TREES_type[INT_VALUE] << " or " << STAT_TREES_type[REAL_VALUE]
                         << " or " << STAT_TREES_type[STATE];
      error.correction_update((error_message.str()).c_str(), (correction_message.str()).c_str());
   }

   for(s= 0; s < nb_sample; s++)
      if ((otrees[s]->_type[0] != INT_VALUE) && (otrees[s]->_type[0] != REAL_VALUE)
          && (otrees[s]->_type[0] != STATE))
      {
         status= false;
         error_message << STAT_label[STATL_SAMPLE] << " " << s+2 << ": "
                       << STAT_TREES_error[STATR_VARIABLE_1_TYPE];
         correction_message << STAT_TREES_type[INT_VALUE] << " or "
                            << STAT_TREES_type[REAL_VALUE] << " or " << STAT_TREES_type[STATE];
         error.correction_update((error_message.str()).c_str(), (correction_message.str()).c_str());
      }

   for(s= 0; s < nb_sample; s++)
   {
      if (otrees[s]->_nb_trees != _nb_trees)
      {
         status= false;
         error_message << STAT_label[STATL_SAMPLE] << " " << s+2 << ": "
                       << STAT_TREES_error[STATR_NB_TREES];
         error.update((error_message.str()).c_str());
      }
      if (status)
      {
         for(t= 0; t < _nb_trees; t++)
         {
            cstree= otrees[s]->trees[t];
            if ((status) && (cstree->get_size() != trees[t]->get_size()))
            {
               status= false;
               error_message << STAT_label[STATL_SAMPLE] << " " << s+2 << ": "
                             << STAT_TREES_label[STATL_TREE] << " " << t+1 << ": "
                             << STAT_TREES_error[STATR_TREE_SIZE];
               error.update((error_message.str()).c_str());
            }
            if (status)
            {
               tie(it, end)= cstree->vertices();
               while(it < end)
               {
                  if ((status) &&
                      (cstree->get_nb_children(*it) != trees[t]->get_nb_children(*it)))
                  {
                     status= false;
                     error_message << STAT_label[STATL_SAMPLE] << " " << s+2 << ": "
                                   << STAT_TREES_label[STATL_TREE] << " " << t+1 << ": "
                                   << STAT_TREES_label[STATL_VERTEX] << " " << *it << ": "
                                   << STAT_TREES_error[STATR_TREE_NB_CHILDREN];
                     error.update((error_message.str()).c_str());
                  }
                  if (status)
                  {
                     tie(chd_it, chd_end)= cstree->children(*it);
                     tie(chs_it, chs_end)= trees[t]->children(*it);
                     while (chs_it < chs_end)
                     {
                        if (*chs_it != *chd_it)
                        {
                           status= false;
                           error_message << STAT_label[STATL_SAMPLE] << " " << s+2 << ": "
                                         << STAT_TREES_label[STATL_TREE] << " " << t+1 << ": "
                                         << STAT_TREES_label[STATL_VERTEX] << " " << *it << ": "
                                         << STAT_label[STATL_CHILD] << " " << *chs_it << ": "
                                         << STAT_TREES_error[STATR_CHILD_ID];
                           error.update((error_message.str()).c_str());

                        }
                        chs_it++;
                        chd_it++;
                     }
                  }
                  it++;
               }
            }
         }
      }
   }

   if ((ref_sample != I_DEFAULT)
         && ((ref_sample < 1) || (ref_sample > nb_sample+1)))
   {
      status= false;
      error.update(STAT_error[STATR_SAMPLE_INDEX]);
   }

   if (status)
   {
      nb_sample++;
      ptrees= new const observed_trees*[nb_sample];
      // ptrees= new const (observed_trees*)[nb_sample];

      ptrees[0]= this;
      inb_integral= _nb_integral;
      inb_float= _nb_float; // nb_variable;
      for(s= 1; s < nb_sample; s++)
      {
         ptrees[s]= otrees[s-1];
         inb_integral+= otrees[s-1]->_nb_integral;
         inb_float+= otrees[s-1]->_nb_float;
      }

      // computation of the tree identifiers
/*
      if (ref_sample == I_DEFAULT)
      {
         for(t= 0; t < _nb_trees; t++)
         {
            for(s= 1; s < nb_sample; s++)
               if (ptrees[s]->identifier[t] != pseq[0]->identifier[t])
                 break;
            if (s < nb_sample)
              break;
         }

         if (t < _nb_trees)
            iidentifier= NULL;

         else
            iidentifier = pseq[0]->identifier;
      }
      else
         iidentifier= ptrees[--ref_sample]->identifier;
*/

      res= new observed_trees(inb_integral, inb_float, NULL, _nb_trees, true);

      if (hsize != NULL)
         res->hsize= new Histogram(*hsize);
      else
         res->hsize= NULL;
      if (hnb_children != NULL)
         res->hnb_children= new Histogram(*hnb_children);
      else
         res->hnb_children= NULL;

      vdest.reset(inb_integral, inb_float);

      // tree copy
      for(t= 0; t < _nb_trees; t++)
      {
         inb_integral= 0;
         inb_float= 0;
         utree= trees[t]->get_structure();
         res->trees[t]= new tree_type(*utree, vdest);
         delete utree;
         for(s= 0; s < nb_sample; s++)
         {
            cstree= ptrees[s]->trees[t];
            cdtree= res->trees[t];
#           ifdef DEBUG
            assert(cdtree != NULL);
            assert(cstree != NULL);
#           endif
            for(var= 0; var < ptrees[s]->_nb_integral; var++)
            {
               tie(it, end)= cstree->vertices();
               while (it < end)
               {
                  vsource= cstree->get(*it);
                  vdest= cdtree->get(*it);
                  vdest.Int(inb_integral)= vsource.Int(var);
                  cdtree->put(*it++, vdest);
               }
               inb_integral++;
            }

            for(var= 0; var < ptrees[s]->_nb_float; var++)
            {
               tie(it, end)= cstree->vertices();
               while (it < end)
               {
                  vsource= cstree->get(*it);
                  vdest= cdtree->get(*it);
                  vdest.Double(inb_float)= vsource.Double(var);
                  cdtree->put(*it++, vdest);
               }
               inb_float++;
            }
         }
      }

      inb_integral= 0;
      inb_float= 0;
      for(s= 0; s < nb_sample; s++)
      {
         for(var= 0; var < ptrees[s]->_nb_integral; var++)
         {
            res->_min_value.Int(inb_integral)= ptrees[s]->_min_value.Int(var);
            res->_max_value.Int(inb_integral)= ptrees[s]->_max_value.Int(var);
            if (ptrees[s]->characteristics[var] != NULL)
               res->characteristics[inb_integral]
                  = new Tree_characteristics(*(ptrees[s]->characteristics[var]));
            else
               res->build_characteristics(inb_integral);
            res->_type[inb_integral]= ptrees[s]->_type[var];
            inb_integral++;
         }

         for(var= 0; var < ptrees[s]->_nb_float; var++)
         {
            res->_min_value.Double(inb_float)= ptrees[s]->_min_value.Double(var);
            res->_max_value.Double(inb_float)= ptrees[s]->_max_value.Double(var);
            inb_float++;
         }
      }
      res->max_size_computation();
      res->max_depth_computation();

      delete [] ptrees;
   }

   return res;
}

/*****************************************************************
 *
 *  Observed_trees selection on a size criterion,
 *  using a Format_error object, size bounds and
 *  a flag on keeping or rejecting the selected trees
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::size_select(Format_error& error,
                                                      int imin_size,
                                                      int imax_size,
                                                      bool keep) const
{
   bool status= true;
   register int t;
   int inb_trees;
   int_array index;
   observed_trees *otrees= NULL;

   error.init();

   if ((imin_size < 1) || (imin_size >= hsize->nb_value)
         || (imin_size > imax_size))
   {
      status= false;
      error.update(STAT_TREES_error[STATR_MIN_TREE_SIZE]);
   }
   if ((imax_size < hsize->offset) || (imax_size < imin_size))
   {
      status= false;
      error.update(STAT_TREES_error[STATR_MAX_TREE_SIZE]);
   }

   if (status)
   {
     // tree selection
      index= new int[_nb_trees];
      inb_trees= 0;

      for(t= 0; t < _nb_trees; t++)
      {
         if (((int)(trees[t]->get_size()) >= imin_size)
               && ((int)(trees[t]->get_size()) <= imax_size))
         {
            if (keep)
               index[inb_trees++]= t;
         }
         else
            if (!keep)
               index[inb_trees++]= t;
      }

      if (inb_trees == 0)
      {
         status= false;
         error.update(STAT_error[STATR_EMPTY_SAMPLE]);
      }

      // copy of the selected trees
      if (status)
         otrees= new observed_trees(*this, inb_trees, index);

      delete [] index;
   }

   return otrees;
}

/*****************************************************************
 *
 *  Select subtrees of Observed_trees using a Format_error object,
 *  the subtree root index, the tree index (all trees if I_DEFAULT)
 *  and a flag on keeping or pruning the subtree.
 *  Pruning at root node results in a tree with one single node.
 *
 **/
 // should be improved by using real tree and vertex identifiers
 // or and array of pairs (tree, index)

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::select_subtrees(Format_error& error,
                                                          int iindex,
                                                          int itree,
                                                          bool keep) const
{
   bool status= true;
   register int t;
   observed_trees *otrees= NULL;
   tree_type **tree_set, *tmptree; // set of selected subtrees

   error.init();

   if ((itree != I_DEFAULT) && ((itree < 0) || (itree >= _nb_trees)))
   // one wants to select a particular tree with invalid index
   {
      status= false;
      ostringstream error_message;
      error_message << itree << ": " << STAT_TREES_error[STATR_TREE_IDENTIFIER];
      error.update((error_message.str()).c_str());
   }
   else
      for(t= 0; t < _nb_trees; t++)
         if ((itree == I_DEFAULT) || (itree == t))
            if ((iindex < 0) || (iindex > (int)(trees[t]->get_size())))
               {
                  status= false;
                  ostringstream error_message;
                  error_message << STAT_TREES_label[STATL_TREE] << " " << t+1 << ": "
                                << STAT_TREES_label[STATL_VERTEX] << " " << iindex << ": "
                                << STAT_TREES_error[STATR_VERTEX_ID];
                  error.update((error_message.str()).c_str());
               }

   if (status)
   {
      tree_set= new tree_type*[_nb_trees];
      for(t= 0; t < _nb_trees; t++)
      {
         if ((itree == t) || (itree == I_DEFAULT))
         {
            tmptree= select_subtree(*(trees[t]), iindex, keep);
            tree_set[t]= new tree_type(*tmptree);
            delete tmptree;
         }
         else
            tree_set[t]= new tree_type(*(trees[t]));
      }

      otrees= new observed_trees(_nb_trees, _type, tree_set);

      for(t= 0; t < _nb_trees; t++)
         delete tree_set[t];

      delete [] tree_set;
   }
   return otrees;

}

/*****************************************************************
 *
 *  Select homogeneous subtrees, using a Format_error object
 *  the variable for testing homogeneity,
 *  the number of values defining the homogeneous zones,
 *  the values that define the homogeneous zones
 *  and a flag on keeping or rejecting the subtrees
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::segmentation_extract(Format_error& error,
                                                               int variable,
                                                               int nb_value,
                                                               int *ivalue,
                                                               bool keep) const

{
   const int _nb_variables= _nb_integral + _nb_float;
   int i, j, k, var, cval, val, size, nb_present_value, nb_selected_value= 0;
       // pnb_values,
   int t, current_zone, nb_nodes;
   unsigned int ui;
   int *pfrequency= NULL, *itype= NULL;
   int *selected_value= NULL; // values to be kept in the segmentation
   bool status= true;
   bool *is_selected= NULL; // selection indicator
                            // for each value of the variable
   key v, vroot, cv;
   key *clist= NULL;
   value tvalue;
   value *ptvalue= new value(_nb_integral-1, _nb_float);
   Histogram *marginal= NULL;
   std::deque<std::deque<key> > zones;
   std::deque<key> node_list, zone_id;
   // set of the homogeneous zones
   std::deque<key> tmp_zone;
   children_iterator it, end;
   tree_type *ctree;
   observed_trees *otrees= NULL;
   pt_tree_type_array subtree_array= NULL;
   std::deque<tree_type*> subtree_vec;

   error.init();

   if ((variable < 1) || (variable > _nb_variables))
   {
      status= false;
      error.update(STAT_error[STATR_VARIABLE_INDEX]);
   }

   /*
   if ((_type[0] != INT_VALUE) && (_type[0] != STATE))
   {
     status= false;
     ostringstream correction_message;
     correction_message << STAT_TREES_word[INT_VALUE] << " or "
                        << STAT_TREES_word[STATE];
     error.correction_update(STAT_TREES_error[STATR_VARIABLE_1_TYPE] ,
                             (correction_message.str()).c_str());
   }
   */

   if (status)
   {
      variable--;

      if ((_type[variable] != INT_VALUE) && (_type[variable] != STATE))
      {
         status = false;
         ostringstream correction_message;
         correction_message << STAT_TREES_word[INT_VALUE] << " or " << STAT_TREES_word[STATE];
         error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE] , (correction_message.str()).c_str());
      }
      else
      {
         if (characteristics[variable] == NULL)
         {
            status= false;
            ostringstream error_message;
            error_message << STAT_label[STATL_VARIABLE] << " " << variable+1 << ": "
                          << STAT_error[STATR_MARGINAL_HISTOGRAM];
            error.update((error_message.str()).c_str());
         }
         else
         {
            marginal= characteristics[variable]->marginal;
            pfrequency= marginal->frequency + marginal->offset;
            nb_present_value= 0;
            // number of values that were observed at least once
            for(i= marginal->offset; i < marginal->nb_value; i++)
            {
               if (*pfrequency > 0)
                  nb_present_value++;
               pfrequency++;
            }
            // pnb_values= _max_value.Int(variable) - _min_value.Int(variable) + 1;
            // potential (maximal) number of values
            if ((nb_value < 1) || (nb_value > (keep ? nb_present_value : nb_present_value - 1)))
            {
               status= false;
               error.update(STAT_TREES_error[STATR_NB_SELECTED_VALUE]);
            }
            else
            {
               selected_value= new int[marginal->nb_value];
               // check that each value is not selected more than once
               for(i= marginal->offset; i < marginal->nb_value; i++)
                  selected_value[i]= false;

               for(i= 0; i < nb_value; i++)
               {
                  if ((ivalue[i] < marginal->offset) ||
                      (ivalue[i] >= marginal->nb_value) ||
                      (marginal->frequency[ivalue[i]] == 0))
                     // no tree containes the value ivalue[i]
                  {
                     status= false;
                     ostringstream error_message;
                     error_message << STAT_label[STATL_VARIABLE] << " " << variable + 1 << ": "
                                   << STAT_label[STATL_VALUE] << " " << ivalue[i] << " "
                                   << STAT_TREES_error[STATR_NOT_PRESENT];
                     error.update((error_message.str()).c_str());
                  }
                  else
                  {
                     if (selected_value[ivalue[i]])
                     {
                        status= false;
                        ostringstream error_message;
                        error_message << STAT_label[STATL_VALUE] << " " << ivalue[i] << " "
                                      << STAT_error[STATR_ALREADY_SELECTED];
                        error.update((error_message.str()).c_str());
                     }
                     else
                        selected_value[ivalue[i]]= true;
                  }
              }
              delete [] selected_value;
              selected_value= NULL;
            }
         }
      }
   }

   if (status)
   {
      itype= new int[_nb_variables - 1];
      for(i= 0; i < variable; i++)
         itype[i]= _type[i];
      for(i= variable + 1 ; i < _nb_variables; i++)
         itype[i-1]= _type[i];

      is_selected= new bool[marginal->nb_value];
      for(i=0; i < marginal->nb_value; i++)
         is_selected[i]= false;

      switch (keep)
      {
         case false :
         {
            nb_selected_value= nb_present_value - nb_value;
            selected_value= new int[nb_selected_value];

            pfrequency= marginal->frequency + marginal->offset;
            i= 0;

            for(j= marginal->offset; j < marginal->nb_value; j++)
            {
               if (*pfrequency > 0)
               {
                  for(k= 0; k < nb_value; k++)
                    if (ivalue[k] == j)
                       break;

                  if (k == nb_value)
                     selected_value[i++]= j;
               }
               pfrequency++;
            }
            break;
         }
         case true :
         {
            nb_selected_value= nb_value;
            selected_value= new int[nb_selected_value];
            for(i= 0; i < nb_selected_value; i++)
               selected_value[i]= ivalue[i];
            break;
         }
      }

      for(i= 0; i < nb_selected_value; i++)
         is_selected[selected_value[i]]= true;

      for(t= 0; t < _nb_trees; t++)
      {
#   ifdef DEBUG
         cerr << "Current tree: " << trees[t]->get_size() << " vertices." << endl;
         trees[t]->display(cout, trees[t]->root());
#   endif
         nb_nodes= 1;

         v= trees[t]->root();
         size= trees[t]->get_size();
         val= trees[t]->get(v).Int(variable);
         current_zone= 0;
         tmp_zone.push_back(v);
         zones.push_back(tmp_zone); // list of all homogeneous zones
         node_list.push_back(v);
         zone_id.push_back(current_zone);
         tmp_zone.pop_back();

         clist= new key[size];
         // correspondence between vids for entire tree and subtree
         for(i= 0; i < size; i++)
            clist[i]= size;

         while (!node_list.empty())
         {
            v= node_list.front();
            cval= (trees[t]->get(v)).Int(variable);
            // value within current zone
            current_zone= zone_id.front();
            tie(it, end)= trees[t]->children(v);
            while (it < end)
            // handle every child vertex
            {
               val= (trees[t]->get(*it)).Int(variable);
               if (val == cval)
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
            node_list.pop_front();
            zone_id.pop_front();
         }

         while (!zones.empty())
         {
            tmp_zone= zones.front(); // current zone
            vroot= tmp_zone.front(); // first vertex of zone
#   ifdef DEBUG
            cerr << "Handle zone rooted at: " << vroot << endl;
#   endif
            val= (trees[t]->get(vroot)).Int(variable);
            if (is_selected[val])
            {
               ctree= new tree_type(_nb_integral-1, _nb_float, 0, 0);
               while (!tmp_zone.empty())
               {
                  v= tmp_zone.front();
                  tvalue= trees[t]->get(v);
                  for(var= 0; var < variable; var++)
                     ptvalue->Int(var)= tvalue.Int(var);
                  for(var= variable+1; var < _nb_integral; var++)
                     ptvalue->Int(var-1)= tvalue.Int(var);
                  for(var= 0; var < _nb_float; var++)
                     ptvalue->Double(var-1)= tvalue.Double(var);
                  cv= ctree->add_vertex(*ptvalue);
                  clist[v]= cv;
                  if (v != vroot)
                     ctree->add_edge(clist[trees[t]->parent(v)], cv);

                  tmp_zone.pop_front();
               }
               subtree_vec.push_back(ctree);
            }
            else
            {
               while (!tmp_zone.empty())
                  tmp_zone.pop_front();
            }
            zones.pop_front();
#   ifdef DEBUG
            cerr << "Size of current tree: " << trees[t]->get_size() << endl;
#   endif
         }
         delete [] clist;
         clist= NULL;
         assert(zones.empty());
         assert(tmp_zone.empty());
         assert(node_list.empty());
         assert(zone_id.empty());
      } // each tree handled

      subtree_array= new tree_type*[subtree_vec.size()];
      for(ui= 0; ui < subtree_vec.size(); ui++)
         subtree_array[ui]= subtree_vec[ui];

      otrees= new Observed_trees(subtree_vec.size(), itype, subtree_array);

      for(ui= 0; ui < subtree_vec.size(); ui++)
      {
         delete subtree_array[ui];
         subtree_array[ui]= NULL;
      }
      delete [] subtree_array;
      subtree_array= NULL;
      delete [] is_selected;
      is_selected= NULL;
      delete [] selected_value;
      selected_value= NULL;
      delete ptvalue;
      ptvalue= NULL;
      delete [] itype;
      itype= NULL;

   } // end if status
   return otrees;
}

/*--------------------------------------------------------------*
 *
 *
 *--------------------------------------------------------------*/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>*
Observed_trees<Generic_Int_fl_container>::difference(Format_error &error,
                                                     int variable) const
{
   bool status= true;
   const int nb_variables= _nb_integral + _nb_float;
   int inb_variable= 0;
   register int t, i, j;
   key pvertex;
   int *ivariable= NULL;
   vertex_iterator it, end;
   value v, vp;
   Observed_trees *res= NULL, *tmp_res= NULL;

   error.init();

   if (variable != I_DEFAULT)
   {
      if ((variable < 1) || (variable > nb_variables))
      {
         status= false;
         error.update(STAT_error[STATR_VARIABLE_INDEX]);
      }
      else
      {
         variable--;

         if ((_type[variable] != INT_VALUE) && (_type[variable] != STATE))
         {
            status= false;
            ostringstream correction_message;
            correction_message << STAT_TREES_word[INT_VALUE] << " or "
                               << STAT_TREES_word[STATE];
            error.correction_update(STAT_TREES_error[STATR_VARIABLE_TYPE],
                                    (correction_message.str()).c_str());
         }
      }
   }

   if (status)
   {
      if (variable == I_DEFAULT)
         inb_variable= nb_variables;
      else
         inb_variable= 1;

      res= new Observed_trees(_nb_trees, _type, trees);

      for(t= 0; t < _nb_trees; t++)
      {
         tie(it, end)= (res->trees[t])->vertices();
         while (it < end)
         {
            if ((res->trees[t])->is_root(*it))
               it++;
            else
            {
               // differentiation
               v= trees[t]->get(*it);
               pvertex= trees[t]->parent(*it);
               vp= trees[t]->get(pvertex);
               for(i= 0; i < _nb_integral; i++)
                  if ((variable == I_DEFAULT) || (variable == i))
                     v.Int(i)= v.Int(i) - vp.Int(i);
               for(i= _nb_integral; i < nb_variables; i++)
                  if ((variable == I_DEFAULT) || (variable == i))
                     {
                        j= i - _nb_integral; // index of the floating variable
                        v.Double(j)= v.Double(j) - vp.Double(j);
                     }
               (res->trees[t])->put(*it++, v);
            }
         }
      }
      res->min_max_value_computation();
      for(i= 0; i < _nb_integral; i++)
         if ((variable == I_DEFAULT) || (variable == i))
            res->build_characteristics(i);

      if (variable != I_DEFAULT)
      {
         // keep only the variable of interest
         ivariable= new int[1];
         ivariable[0]= variable;
         tmp_res= select_variable(error, 1, ivariable);
         if (tmp_res != NULL)
         {
            delete res;
            res= NULL;
            res= tmp_res;
         }
         delete [] ivariable;
         ivariable= NULL;
      }
   }
   return res;
}

/*****************************************************************
 *
 *  Print an Observed_trees object on a single line
 *  using an output stream
 *
 **/

template<typename Generic_Int_fl_container>
std::ostream&
Observed_trees<Generic_Int_fl_container>::line_write(std::ostream& os) const
{
   const register int nb_variable= _nb_integral+_nb_float;
   const int cumul_size= cumul_size_computation();

   os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << "   "
      << _nb_trees << " " << STAT_TREES_label[_nb_trees == 1 ? STATL_TREE : STATL_TREES] << "   "
      << STAT_TREES_label[STATL_CUMULATIVE_SIZE] << ": " << cumul_size;

   return os;
}

/*****************************************************************
 *
 *  Writing of Observed_trees using an output stream,
 *  the format (line/column) and a flag on the level of detail.
 *
 **/

template<typename Generic_Int_fl_container>
ostream& Observed_trees<Generic_Int_fl_container>::ascii_data_write(ostream& os,
                                                                    char format,
                                                                    bool exhaustive) const
{
   ascii_write(os, exhaustive, false);
   ascii_print(os, format, false);

   return os;
}

/*****************************************************************
 *
 *  Writing of Observed_trees into a file
 *  using a Format_error object, the path,
 *  the format (line/column) and a flag on the level of detail.
 *
 **/

template<typename Generic_Int_fl_container>
bool Observed_trees<Generic_Int_fl_container>::ascii_data_write(Format_error& error,
                                                                const char * path,
                                                                char format,
                                                                bool exhaustive) const
{
   bool status= false;
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
      if (format != 'a')
         ascii_write(out_file, exhaustive, true);
      ascii_print(out_file, format, true);
   }
   return status;
}

/*****************************************************************
 *
 *  Gnuplot output of Trees using a Format_error
 *  a prefix for the files and the title of the figures
 *
 **/

template<typename Generic_Int_fl_container>
bool Observed_trees<Generic_Int_fl_container>::plot_data_write(Format_error& error,
                                                               const char * prefix,
                                                               const char * title) const
{
   bool status;
   register int i, j, k;
   int min_index, max_index, *pfrequency, *size_nb_trees;
   ostringstream *data_file_name;
   vertex_iterator it, end;
   // CPL for gcc 3.4
   register int t;
   vertex_iterator begin;

   error.init();

   if (_nb_trees > PLOT_NB_TREES)
   {
      status = false;
      error.update(STAT_error[STATR_NB_TREES]);
   }
   else
   {
      // writes the data files
      data_file_name= new ostringstream[hsize->nb_value];

      data_file_name[hsize->offset] << prefix << hsize->offset << ".dat";
      status= plot_print((data_file_name[hsize->offset].str()).c_str(), hsize->offset);

      if (!status)
         error.update(STAT_error[STATR_FILE_PREFIX]);
      else
      {
         pfrequency= hsize->frequency + hsize->offset + 1;
         for(i= hsize->offset+1; i < hsize->nb_value; i++)
            if (*pfrequency++ > 0)
            {
               data_file_name[i] << prefix << i << ".dat";
               plot_print((data_file_name[i].str()).c_str(), i);
            }

         size_nb_trees= new int[hsize->nb_value];

         if ((_type[0] == TIME) || (_type[0] == POSITION))
         {
            max_index= 0;
            for(i= 0; i < _nb_trees; i++)
            {
               tie(it, end)= trees[i]->vertices();
               if ((trees[i]->get(*(end-1))).Int(0) > max_index)
                  max_index= (trees[i]->get(*(end-1))).Int(0);
            }

            min_index= max_index;
            for(t= 0; t < this->_nb_trees; t++)
               if ((trees[i]->get(*begin)).Int(0) < min_index)
                  min_index= (trees[i]->get(*begin)).Int(0);
         }

         // writes the command and printing files

         for(i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
               case 0 :
                  file_name[0] << prefix << ".plot";
                  break;
               case 1 :
                  file_name[0] << prefix << ".print";
                  break;
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;
               file_name[1] << label(prefix) << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n";

            for(j= ((_type[0] == TIME) || (_type[0] == POSITION) ? 1 : 0); j < _nb_integral; j++)
            {
               for(k= 0; k < hsize->nb_value; k++)
                  size_nb_trees[k]= 0;

               out_file << "set title \"";
               if (title)
               {
                  out_file << title;
                  if (_nb_integral > ((_type[0] == TIME) || (_type[0] == POSITION) ? 2 : 1))
                     out_file << " - ";
               }

               if ((_type[0] == TIME) || (_type[0] == POSITION))
               {
                  if (_nb_integral > 2)
                     out_file << STAT_label[STATL_VARIABLE] << " " << j;

                  out_file << "\"\n\n";

                  if (max_index-min_index < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if (_max_value.Int(j)-_min_value(j) < TIC_THRESHOLD)
                     out_file << "set ytics " << MIN(_min_value.Int(j), 0) << ",1" << endl;

                  out_file << "plot [" << min_index << ":" << max_index << "] ["
                           << MIN(_min_value.Int(j), 0) << ":" << MAX(_max_value.Int(j), _min_value.Int(j)+1) << "] ";

                  for(k= 0; k < _nb_trees; k++)
                  {
                     out_file << "\"" << label((data_file_name[this->size[k]].str()).c_str()) << "\" using "
                              << size_nb_trees[trees[k]->get_size()]*_nb_integral+1 << " : "
                              << size_nb_trees[trees[k]->get_size()]*_nb_integral+j+1;
                     if (_nb_trees <= PLOT_TITLE_NB_TREES)
                        out_file << " title \"" << k // identifier[k]
                                 << "\" with linespoints";
                     else
                        out_file << " notitle with linespoints";

                     if (k < _nb_trees-1)
                        out_file << ",\\";

                     out_file << endl;
                     (size_nb_trees[trees[k]->get_size()])++;
                  }

                  if (max_index-min_index < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;

                  if (_max_value.Int(j)-_min_value.Int(j) < TIC_THRESHOLD)
                     out_file << "set ytics autofreq" << endl;
               }
               else // (_type[0] != TIME) && (_type[0] != POSITION)
               {
                  if (_nb_integral > 1)
                     out_file << STAT_label[STATL_VARIABLE] << " " << j + 1;

                  out_file << "\"\n\n";

                  if (this->max_size-1 < TIC_THRESHOLD)
                     out_file << "set xtics 0,1" << endl;

                  if (_max_value.Int(j)-_min_value.Int(j) < TIC_THRESHOLD)
                     out_file << "set ytics " << MIN(_min_value.Int(j), 0) << ",1" << endl;

                  out_file << "plot [0:" << this->max_size-1 << "] [" << MIN(_min_value.Int(j), 0)
                           << ":" << MAX(_max_value.Int(j), _min_value.Int(j)+1) << "] ";
                  for(k= 0; k < _nb_trees; k++)
                  {
                     out_file << "\"" << label((data_file_name[this->size[k]].str()).c_str()) << "\" using "
                              << size_nb_trees[trees[k]->get_size()]*_nb_trees+j+1;
                     if (_nb_trees <= PLOT_TITLE_NB_TREES)
                        out_file << " title \"" << k // identifier[k]
                                 << "\" with linespoints";
                     else
                        out_file << " notitle with linespoints";

                     if (k < _nb_trees-1)
                        out_file << ",\\";

                     out_file << endl;
                     (size_nb_trees[trees[k]->get_size()])++;
                  }

                  if (this->max_size-1 < TIC_THRESHOLD)
                     out_file << "set xtics autofreq" << endl;

                  if (_max_value.Int(j)-_min_value.Int(j)-1 < TIC_THRESHOLD)
                     out_file << "set ytics autofreq" << endl;
               }

               if ((i == 0) && (j < _nb_integral-1))
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;

               out_file << endl;
            }

            if (i == 1)
               out_file << "\nset terminal x11" << endl;

            out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
         }
         delete [] size_nb_trees;
      }
      delete [] data_file_name;
   }
   return status;
}


/*****************************************************************
 *
 *  Prints Observed_trees using an output stream and
 *  a flag on the detail level.
 *
 **/

template<typename Generic_Int_fl_container>
ostream& Observed_trees<Generic_Int_fl_container>::ascii_write(ostream& os,
                                                               bool exhaustive) const
{ return ascii_write(os, exhaustive, false); }

/*****************************************************************
 *
 *  Prints Observed_trees using an output stream,
 *  a flag on the detail level and one on the comment level.
 *
 **/

template<typename Generic_Int_fl_container>
ostream& Observed_trees<Generic_Int_fl_container>::ascii_write(ostream& os,
                                                               bool exhaustive,
                                                               bool comment_flag) const
{
   register int var, nb_variable= _nb_integral+_nb_float,
                cumul_size= cumul_size_computation(),
                cumul_children= cumul_nb_children_computation();
   // double mean, variance;

   // number of variables
   os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

   // integral variables
   for(var= 0; var < _nb_integral; var++)
   {
      os << "\n" << STAT_word[STATW_VARIABLE] << " " << var+1 << " : "
         << STAT_TREES_type[_type[var]] << "   ";
      if (comment_flag)
         os << "# ";

      if (characteristics[var] != NULL)
      {
         // cardinal number of the set of observed values
         if (characteristics[var]->marginal != NULL)
            os << characteristics[var]->marginal->nb_value << " ";

         switch (_type[var])
         {
            case STATE :
               os << STAT_label[characteristics[var]->marginal->nb_value == 1 ?
                                STATL_STATE : STATL_STATES] << endl;
            case INT_VALUE :
               os << STAT_label[characteristics[var]->marginal->nb_value == 1 ?
                                STATL_VALUE : STATL_VALUES] << endl;
         }

         os << "\n";
         if (comment_flag)
            os << "# ";

         os << STAT_label[_type[var] == STATE ? STATL_STATE : STATL_VALUE] << " "
            << STAT_label[STATL_HISTOGRAM] << " - ";

         (characteristics[var]->marginal)->ascii_characteristic_print(os, false, comment_flag);

         if ((characteristics[var]->marginal->nb_value <= ASCII_NB_VALUE) || (exhaustive))
         {
            os << "\n";
            if (comment_flag)
               os << "# ";

            os << "   | " << STAT_label[_type[var] == STATE ? STATL_STATE : STATL_VALUE] << " "
               << STAT_label[STATL_HISTOGRAM] << endl;
            characteristics[var]->marginal->ascii_print(os, comment_flag);
         }
         characteristics[var]->ascii_print(os, _type[var], *hsize,
                                           exhaustive, comment_flag);
      } // end if characteristics exist
   } // end "for each variable"

   // floating variables
   for(var= 0; var < _nb_float; var++)
   {
      os << "\n" << STAT_word[STATW_VARIABLE] << " " << _nb_integral+var+1 << " : "
         << STAT_TREES_type[REAL_VALUE];

      os << "   ";
      if (comment_flag)
         os << "# ";

      os << "(" << STAT_label[STATL_MIN_VALUE] << ": " << _min_value.Double(var) << ", "
         << STAT_label[STATL_MAX_VALUE] << ": " << _max_value.Double(var) << ")" << endl;

      os << "\n";
      if (comment_flag)
         os << "# ";
   }

   os << "\n";
   if (comment_flag)
      os << "# ";

   os << STAT_TREES_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
   hsize->ascii_characteristic_print(os, false, comment_flag);

   if (exhaustive)
   {
      os << "\n";
      if (comment_flag)
         os << "# ";

      os << "   | " << STAT_TREES_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      hsize->ascii_print(os, comment_flag);
   }

   os << "\n";
   if (comment_flag)
      os << "# ";

   os << STAT_TREES_label[STATL_CUMULATIVE_SIZE] << ": " << cumul_size << endl;

   os << "\n";
   if (comment_flag)
      os << "# ";

   os << STAT_TREES_label[STATL_TREE_CHILDREN] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
   hnb_children->ascii_characteristic_print(os, false, comment_flag);

   if (exhaustive)
   {
      os << "\n";
      if (comment_flag)
         os << "# ";

      os << "   | " << STAT_TREES_label[STATL_TREE_CHILDREN] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      hnb_children->ascii_print(os, comment_flag);
   }

   os << "\n";
   if (comment_flag)
      os << "# ";

   os << STAT_TREES_label[STATL_CUMULATIVE_CHILDREN] << ": " << cumul_children << endl;

   return os;
}

/*
template<typename Generic_Int_fl_container>
ostream& Observed_trees<Generic_Int_fl_container>::ascii_write(ostream& os,
                                                               bool exhaustive,
                                                               bool comment_flag) const
{
   register int var, nb_variable= _nb_integral+_nb_float,
                cumul_size= cumul_size_computation(),
                cumul_children= cumul_nb_children_computation();
   double mean, variance;

   os << nb_variable << " " << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

   // integral variables
   for(var= 0; var < _nb_integral; var++)
   {
      os << "\n" << STAT_word[STATW_VARIABLE] << " " << var+1 << " : "
         << STAT_TREES_type[_type[var]];

      if (_type[var] != POSITION)
      {
         os << "   ";
         if (comment_flag)
            os << "# ";

         os << "(" << STAT_label[STATL_MIN_VALUE] << ": " << _min_value.Int(var) << ", "
            << STAT_label[STATL_MAX_VALUE] << ": " << _max_value.Int(var) << ")" << endl;

         os << "\n";
         if (comment_flag)
            os << "# ";

         os << (_type[var] == TIME ? STAT_TREES_type[STATL_TIME] : STAT_label[STATL_VALUE]) << " "
            << STAT_label[STATL_HISTOGRAM] << " - ";

         if (characteristics[var]->marginal != NULL)
         {
            (characteristics[var]->marginal)->ascii_characteristic_print(os, exhaustive, comment_flag);

            if ((characteristics[var]->_max_value <= ASCII_NB_VALUE) || (exhaustive))
            {
               os << "\n";
               if (comment_flag)
                  os << "# ";

               os << "   | " << (_type[var] == TIME ? STAT_TREES_type[STATL_TIME] : STAT_label[STATL_VALUE]) << " "
                  << STAT_label[STATL_HISTOGRAM] << endl;
               characteristics[var]->marginal->ascii_print(os, comment_flag);
            }
         }
         else
            if (_type[var] == INT_VALUE)
            {
               os << STAT_label[STATL_SAMPLE_SIZE] << ": " << cumul_size << endl;

               mean= mean_computation(var);
               variance= variance_computation(var, mean);

               if (comment_flag)
                  os << "# ";

               os << STAT_label[STATL_MEAN] << ": " << mean << "   "
                  << STAT_label[STATL_VARIANCE] << ": " << variance << "   "
                  << STAT_label[STATL_STANDARD_DEVIATION] << ": " << sqrt(variance) << endl;

               if ((exhaustive) && (variance > 0.))
               {
                 if (comment_flag)
                    os << "# ";

                 os << STAT_label[STATL_SKEWNESS_COEFF] << ": " << skewness_computation(var , mean , variance) << "   "
                    << STAT_label[STATL_KURTOSIS_COEFF] << ": " << kurtosis_computation(var , mean , variance) << endl;
               }
            }
      }
   }

   // floating variables
   for(var= 0; var < _nb_float; var++)
   {
      os << "\n" << STAT_word[STATW_VARIABLE] << " " << _nb_integral+var+1 << " : "
         << STAT_TREES_type[REAL_VALUE];

      os << "   ";
      if (comment_flag)
         os << "# ";

      os << "(" << STAT_label[STATL_MIN_VALUE] << ": " << _min_value.Double(var) << ", "
         << STAT_label[STATL_MAX_VALUE] << ": " << _max_value.Double(var) << ")" << endl;

      os << "\n";
      if (comment_flag)
         os << "# ";
   }

   os << "\n";
   if (comment_flag)
      os << "# ";

   os << STAT_TREES_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
   hsize->ascii_characteristic_print(os, false, comment_flag);

   if (exhaustive)
   {
      os << "\n";
      if (comment_flag)
         os << "# ";

      os << "   | " << STAT_TREES_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      hsize->ascii_print(os, comment_flag);
   }

   os << "\n";
   if (comment_flag)
      os << "# ";

   os << STAT_TREES_label[STATL_CUMULATIVE_SIZE] << ": " << cumul_size << endl;

   os << STAT_TREES_label[STATL_TREE_CHILDREN] << " " << STAT_label[STATL_HISTOGRAM] << " - ";
   hnb_children->ascii_characteristic_print(os, false, comment_flag);

   if (exhaustive)
   {
      os << "\n";
      if (comment_flag)
         os << "# ";

      os << "   | " << STAT_TREES_label[STATL_TREE_CHILDREN] << " " << STAT_label[STATL_HISTOGRAM] << endl;
      hnb_children->ascii_print(os, comment_flag);
   }

   os << "\n";
   if (comment_flag)
      os << "# ";

   os << STAT_TREES_label[STATL_CUMULATIVE_CHILDREN] << ": " << cumul_children << endl;

   return os;
}*/

/*****************************************************************
 *
 *  Writes an Observed_trees object, including all the tree attributes,
 *  using and output stream, the format ('c'olumn, 'l'ine or 'a'array),
 *  a flag on the comments and the number of characters per line
 *
 **/

template<typename Generic_Int_fl_container>
ostream& Observed_trees<Generic_Int_fl_container>::ascii_print(ostream& os,
                                                               char format,
                                                               bool comment_flag,
                                                               int line_nb_character) const
{
#ifdef __GNUC__
#warning Observed_trees<Generic_Int_fl_container>::ascii_print  not implemented
#endif
   return os;
}

/*****************************************************************
 *
 *  Prints Observed_trees into a file using a Format_error object,
 *  the file path and a flag on the level of detail
 *
 **/

template<typename Generic_Int_fl_container>
bool Observed_trees<Generic_Int_fl_container>::ascii_write(Format_error& error,
                                                           const char * path,
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
      ascii_write(out_file, exhaustive, false);
   }

   return status;
}

/*****************************************************************
 *
 *  Prints Observed_trees in a spreadsheet fashion
 *  using a Format_error object and the file path
 *
 **/

template<typename Generic_Int_fl_container>
bool Observed_trees<Generic_Int_fl_container>::spreadsheet_write(Format_error& error,
                                                                 const char * path) const
{
   register int var, nb_variable= _nb_integral+_nb_float,
                cumul_size= cumul_size_computation(),
                cumul_children= cumul_nb_children_computation();
   double mean, variance;
   ofstream out_file(path);
   bool status;

   error.init();

   if (!out_file)
   {
      status= false;
      error.update(STAT_error[STATR_FILE_NAME]);
   }
   else
   {
      status = true;

      out_file << nb_variable << "\t" << STAT_word[nb_variable == 1 ? STATW_VARIABLE : STATW_VARIABLES] << endl;

      for(var= 0; var < _nb_integral; var++)
      {
         out_file << "\n" << STAT_word[STATW_VARIABLE] << "\t" << var+1 << "\t"
                  << STAT_TREES_type[_type[var]];

         if (_type[var] != POSITION)
         {
            out_file << "\t\t";
            out_file << STAT_label[STATL_MIN_VALUE] << "\t" << _min_value.Int(var) << "\t\t"
                     << STAT_label[STATL_MAX_VALUE] << "\t" << _max_value.Int(var) << endl;

            out_file << "\n";
            out_file << (_type[var] == TIME ? STAT_TREES_type[STATL_TIME] : STAT_label[STATL_VALUE]) << " "
                     << STAT_label[STATL_HISTOGRAM] << "\t";

            if (characteristics[var]->marginal != NULL)
            {
               (characteristics[var]->marginal)->spreadsheet_characteristic_print(out_file);

                out_file << "\n\t" << (_type[var] == TIME ? STAT_TREES_type[STATL_TIME] : STAT_label[STATL_VALUE]) << " "
                         << STAT_label[STATL_HISTOGRAM] << endl;
                (characteristics[var]->marginal)->spreadsheet_print(out_file);
            }
            else
               if (_type[var] == INT_VALUE)
               {
                  out_file << STAT_label[STATL_SAMPLE_SIZE] << "\t" << cumul_size << endl;

                  mean= mean_computation(var);
                  variance= variance_computation(var, mean);

                  out_file << STAT_label[STATL_MEAN] << "\t" << mean << "\t\t"
                           << STAT_label[STATL_VARIANCE] << "\t" << variance << "\t\t"
                           << STAT_label[STATL_STANDARD_DEVIATION] << "\t" << sqrt(variance) << endl;

                  if (variance > 0.)
                  {
                    out_file << STAT_label[STATL_SKEWNESS_COEFF] << "\t" << skewness_computation(var , mean , variance) << "\t\t"
                             << STAT_label[STATL_KURTOSIS_COEFF] << "\t" << kurtosis_computation(var , mean , variance) << endl;
                  }
               }
         }
      }
   }

   out_file << "\n";

   out_file << STAT_TREES_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
   hsize->spreadsheet_characteristic_print(out_file);

   out_file << "\n\t" << STAT_TREES_label[STATL_TREE_SIZE] << "\t" << STAT_label[STATL_HISTOGRAM] << endl;
   hsize->spreadsheet_print(out_file);

   out_file << "\n";

   out_file << STAT_TREES_label[STATL_CUMULATIVE_SIZE] << "\t" << cumul_size << endl;

   out_file << "\n";

   out_file << STAT_TREES_label[STATL_TREE_CHILDREN] << " " << STAT_label[STATL_HISTOGRAM] << "\t";
   hnb_children->spreadsheet_characteristic_print(out_file);

   out_file << "\n\t" << STAT_TREES_label[STATL_TREE_CHILDREN] << "\t" << STAT_label[STATL_HISTOGRAM] << endl;
   hnb_children->spreadsheet_print(out_file);

   out_file << "\n";

   out_file << STAT_TREES_label[STATL_CUMULATIVE_CHILDREN] << "\t" << cumul_children << endl;

   return status;
}

/*****************************************************************
 *
 *  Gnuplot output of Observed_trees
 *  using a Format_error object, a prefix for the files
 *  and the title of the output figures
 *
 **/

template<typename Generic_Int_fl_container>
bool Observed_trees<Generic_Int_fl_container>::plot_write(Format_error& error,
                                                          const char * prefix,
                                                          const char * title) const
{
   register int var;
   bool status= true, characteristics_computed;

   error.init();

   if (characteristics == NULL)
   {
      status= false;
      error.update(STAT_TREES_error[STATR_CHARACTERISTICS_NOT_COMPUTED]);
   }
   else
   {
      var= 0;
      while ((var < _nb_integral) && (status))
      {
         characteristics_computed=
            ((characteristics[var] != NULL) &&
               ((characteristics[var]->first_occurrence_root != NULL) ||
                (characteristics[var]->first_occurrence_leaves != NULL) ||
                (characteristics[var]->sojourn_size != NULL) ||
                (characteristics[var]->nb_zones != NULL) ||
                (characteristics[var]->nb_occurrences != NULL)));
         // check wether one characteristic at least is computed
         if (characteristics_computed)
            status= characteristics[var]->plot_print(prefix, title, var,
                                                     _nb_integral, _type[var] , *hsize);
         else
            status= plot_print(prefix, title, var, _nb_integral);
         var++;
      }
   }
   return status;
}

/*{
   bool status;
   register int var, i , j , k;
   int nb_histo;
   const Histogram **phisto; // = NULL;
   ostringstream data_file_name;

   error.init();

   // data file writing

   data_file_name << prefix << ".dat";

   phisto= new (const Histogram*)[_nb_integral+1];

   nb_histo= 0;
   for(var= 0; var < _nb_integral; var++)
      if ((_type[var] != POSITION) && (characteristics[var]->marginal != NULL))
         phisto[nb_histo++]= characteristics[var]->marginal;

   status= hsize->plot_print((data_file_name.str()).c_str(), nb_histo, phisto);
   // Is something similar required for hnb_children ?

   delete [] phisto;

   if (!status)
      error.update(STAT_error[STATR_FILE_PREFIX]);
   // writes the command and printing files
   else
   {
      for(i= 0; i < 2; i++)
      {
         j= 2;
         ostringstream file_name[2];

         switch (i)
         {
            case 0 :
              file_name[0] << prefix << ".plot";
              break;
            case 1 :
              file_name[0] << prefix << ".print";
              break;
         }

         ofstream out_file((file_name[0].str()).c_str());

         if (i == 1)
         {
            out_file << "set terminal postscript" << endl;
            file_name[1] << label(prefix) << ".ps";
            out_file << "set output \"" << file_name[1].str() << "\"\n\n";
         }

         out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                  << "set title";
         if (title)
            out_file << " \"" << title << "\"";

         out_file << "\n\n";

         for(k= 0; k < _nb_integral; k++)
         {
            if ((_type[k] != POSITION) && (characteristics[k]->marginal != NULL))
            {
               if ((characteristics[k]->marginal)->nb_value-1 < TIC_THRESHOLD)
                  out_file << "set xtics 0,1" << endl;

               if ((int)((characteristics[k]->marginal)->max * YSCALE) + 1 < TIC_THRESHOLD)
                  out_file << "set ytics 0,1" << endl;

               out_file << "plot [0:" << MAX((characteristics[k]->marginal)->nb_value - 1 , 1) << "] [0:"
                        << (int)((characteristics[k]->marginal)->max * YSCALE) + 1 << "] \""
                        << label((data_file_name.str()).c_str()) << "\" using " << j++ << " title \""
                        << STAT_label[STATL_VARIABLE] << " " << k + 1 << " - "
                        << (_type[k] == TIME ? STAT_TREES_type[STATL_TIME] : STAT_label[STATL_VALUE]) << " "
                        << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

               if ((characteristics[k]->marginal)->nb_value-1 < TIC_THRESHOLD)
                  out_file << "set xtics autofreq" << endl;

               if ((int)((characteristics[k]->marginal)->max * YSCALE)+1 < TIC_THRESHOLD)
                  out_file << "set ytics autofreq" << endl;


               if (i == 0)
                  out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;

               out_file << endl;
            }
        }

         if (hsize->nb_value-1 < TIC_THRESHOLD)
            out_file << "set xtics 0,1" << endl;

         if ((int)(hsize->max * YSCALE) + 1 < TIC_THRESHOLD)
            out_file << "set ytics 0,1" << endl;

         out_file << "plot [0:" << hsize->nb_value-1 << "] [0:"
                  << (int)(hsize->max * YSCALE)+1 << "] \""
                  << label((data_file_name.str()).c_str()) << "\" using 1 title \""
                  << STAT_TREES_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_HISTOGRAM]
                  << "\" with impulses" << endl;

         if (hsize->nb_value-1 < TIC_THRESHOLD)
            out_file << "set xtics autofreq" << endl;

         if ((int)(hsize->max * YSCALE)+1 < TIC_THRESHOLD)
            out_file << "set ytics autofreq" << endl;

         if (i == 1)
            out_file << "\nset terminal x11" << endl;

         out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
      }
  }
  return status;
}
*/

/*****************************************************************
 *
 *  Print the histograms corresponding to the tree size and
 *  children number for Observed_trees class
 *
 **/

template<typename Generic_Int_fl_container>
std::ostream& Observed_trees<Generic_Int_fl_container>::ascii_write_size_histogram(std::ostream &os,
                                                                                   bool exhaustive,
                                                                                   bool file_flag) const
{
   if (hsize != NULL)
      hsize->ascii_write(os, exhaustive, file_flag);
   return os;
}

template<typename Generic_Int_fl_container>
std::ostream& Observed_trees<Generic_Int_fl_container>::ascii_write_nb_children_histogram(std::ostream &os,
                                                                                          bool exhaustive,
                                                                                          bool file_flag) const
{
   if (hnb_children != NULL)
      hnb_children->ascii_write(os, exhaustive, file_flag);
   return os;
}

/*****************************************************************
 *
 *  Computation of the empirical mean of a given variable
 *  for Observed trees
 *
 **/

template<typename Generic_Int_fl_container>
double Observed_trees<Generic_Int_fl_container>::mean_computation(int variable) const
{
   typedef typename Int_fl_tree<Generic_Int_fl_container>::vertex_iterator vertex_iterator;
   typedef typename Int_fl_tree<Generic_Int_fl_container>::tree_type tree_type;

   register int t; // , j;
   double mean= D_INF;
   vertex_iterator it, end;
   tree_type *ctree; // current_tree;

   if (characteristics[variable]->marginal != NULL)
      mean= (characteristics[variable]->marginal)->mean;
   else
   {
      if ((_type[variable] == INT_VALUE) || (_type[variable] == NB_INTERNODE))
      {
         mean= 0.;
         for(t= 0; t < _nb_trees; t++)
         {
            ctree= trees[t];
            tie(it, end)= ctree->vertices();
            while (it < end)
               mean+= (ctree->get(*it++)).Int(variable);
         }
      }
      mean /= cumul_size_computation();
   }

   return mean;
}

/*****************************************************************
 *
 *  Computation of the (unbiased, or so-called estimated)
 *  empirical variance of a given variable for Observed trees
 *  using its empirical mean
 *
 **/

template<typename Generic_Int_fl_container>
double Observed_trees<Generic_Int_fl_container>::variance_computation(int variable,
                                                                      double mean) const
{
   typedef typename Int_fl_tree<Generic_Int_fl_container>::vertex_iterator vertex_iterator;
   typedef typename Int_fl_tree<Generic_Int_fl_container>::tree_type tree_type;

   register int t; // , j;
   int cumul_size= cumul_size_computation();
   double variance= D_DEFAULT, diff;
   vertex_iterator it, end;
   tree_type *ctree; // current_tree;

   if (characteristics[variable]->marginal != NULL)
      variance= (characteristics[variable]->marginal)->variance;
   else
   {
      if ((_type[variable] == INT_VALUE) || (_type[variable] == NB_INTERNODE))
      {
         variance= 0.;
         if (cumul_size > 1)
            for(t= 0; t < _nb_trees; t++)
            {
               ctree= trees[t];
               tie(it, end)= ctree->vertices();
               while (it < end)
               {
                  diff= (ctree->get(*it++)).Int(variable) - mean;
                  variance+= diff * diff;
               }
            }
      }
      variance /= (cumul_size - 1);
   }
   return variance;
}

/*****************************************************************
 *
 *  Computation of the empirical skewness of a given variable
 *  for Observed trees using its empirical mean and variance
 *
 **/

template<typename Generic_Int_fl_container>
double Observed_trees<Generic_Int_fl_container>::skewness_computation(int variable,
                                                                      double mean,
                                                                      double variance) const
{
   typedef typename Int_fl_tree<Generic_Int_fl_container>::vertex_iterator vertex_iterator;
   typedef typename Int_fl_tree<Generic_Int_fl_container>::tree_type tree_type;

   register int t; // , j;
   int cumul_size= cumul_size_computation();
   double skewness= D_DEFAULT, diff;
   vertex_iterator it, end;
   tree_type *ctree; // current_tree;

   if (characteristics[variable]->marginal != NULL)
      skewness= (characteristics[variable]->marginal)->skewness_computation();
   else
   {
      if ((_type[variable] == INT_VALUE) || (_type[variable] == NB_INTERNODE))
      {
         skewness= 0.;
         if ((cumul_size > 2) && (variance > 0.))
            for(t= 0; t < _nb_trees; t++)
            {
               ctree= trees[t];
               tie(it, end)= ctree->vertices();
               while (it < end)
               {
                  diff= (ctree->get(*it++)).Int(variable) - mean;
                  skewness+= diff * diff * diff;
               }
            }
      }
      skewness= skewness * cumul_size
                   / ((cumul_size - 1) * (double)(cumul_size - 2) * pow(variance, 1.5));
   }
   return skewness;
}

/*****************************************************************
 *
 *  Computation of the empirical kurtosis of a given variable
 *  for Observed trees using its empirical mean and variance
 *
 **/

template<typename Generic_Int_fl_container>
double Observed_trees<Generic_Int_fl_container>::kurtosis_computation(int variable,
                                                                      double mean, double variance) const
{
   typedef typename Int_fl_tree<Generic_Int_fl_container>::vertex_iterator vertex_iterator;
   typedef typename Int_fl_tree<Generic_Int_fl_container>::tree_type tree_type;

   register int t; // , j;
   int cumul_size= cumul_size_computation();
   double kurtosis= D_INF, diff;
   vertex_iterator it, end;
   tree_type *ctree; // current_tree;

   if (characteristics[variable]->marginal != NULL)
      kurtosis= (characteristics[variable]->marginal)->kurtosis_computation();
   else
   {
      if ((_type[variable] == INT_VALUE) || (_type[variable] == NB_INTERNODE))
      {
         kurtosis= 0.;
         if ((cumul_size < 2) || (variance == 0.))
            kurtosis= -2.;
         else
            for(t= 0; t < _nb_trees; t++)
            {
               ctree= trees[t];
               tie(it, end)= ctree->vertices();
               while (it < end)
               {
                  diff= (ctree->get(*it++)).Int(variable) - mean;
                  kurtosis+= diff * diff * diff * diff;
               }
            }
      }
      kurtosis= kurtosis / ((cumul_size - 1) * variance * variance) - 3;
   }
   return kurtosis;
}

template<typename Generic_Int_fl_container>
double Observed_trees<Generic_Int_fl_container>::iid_information_computation() const
{
   register int i;
   double information= 0.;


   for(i= (((_type[0] != STATE) || (_nb_integral == 1)) ? 0 : 1); i < _nb_integral; i++)
      if ((characteristics[i] != NULL) && (characteristics[i]->marginal != NULL))
         information+= characteristics[i]->marginal->information_computation();

   return information;
}


/*****************************************************************
 *
 *  Access to the number of integral variables
 *  and number of float variables
 *
 **/

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_nb_int() const
{ return _nb_integral; }

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_nb_float() const
{ return _nb_float; }

/*****************************************************************
 *
 *  Access to the types of variables, maximal and minimal values
 *
 **/

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_type(int variable) const
{ return _type[variable]; }

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_min_int_value(int variable) const
{ return _min_value.Int(variable); }

template<typename Generic_Int_fl_container>
double Observed_trees<Generic_Int_fl_container>::get_min_fl_value(int variable) const
{ return _min_value.Double(variable); }

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_max_int_value(int variable) const
{ return _max_value.Int(variable); }

template<typename Generic_Int_fl_container>
double Observed_trees<Generic_Int_fl_container>::get_max_fl_value(int variable) const
{ return _max_value.Double(variable); }

/*****************************************************************
 *
 *  Access to the maximal size and depth of the trees
 *  for Observed_trees
 *
 **/

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_max_size() const
{ return _max_size; }

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_max_depth() const
{ return _max_depth; }

/*****************************************************************
 *
 *  Return the size of a given tree of Observed_trees
 *
 **/

template<typename Generic_Int_fl_container>
unsigned int Observed_trees<Generic_Int_fl_container>::get_size(int index) const
{
   assert((index >= 0) and (index < _nb_trees));
   return trees[index]->get_size();
}

/*****************************************************************
 *
 *  Access to vertex identifiers of Observed_trees,
 *  presented as trees
 *
 **/

template<typename Generic_Int_fl_container>
pt_One_int_tree_array
Observed_trees<Generic_Int_fl_container>::get_identifier_trees()
{
   register int t;
   pt_One_int_tree_array res= new One_int_tree*[_nb_trees];

   for(t= 0; t < _nb_trees; t++)
      res[t]= get_identifier_tree(t);

   return res;
}

/*****************************************************************
 *
 *  Access to vertex identifiers of one tree of Observed_trees,
 *  presented as a tree
 *
 **/

template<typename Generic_Int_fl_container>
One_int_tree*
Observed_trees<Generic_Int_fl_container>::get_identifier_tree(int index)
{
   typedef One_int_tree::value value;
   typedef One_int_tree::vertex_iterator vertex_iterator;
   One_int_tree *res= NULL;
   Unlabelled_tree *utree= NULL;
   value v;
   vertex_iterator it, end;

   assert((index >= 0) && (index < _nb_trees));
   utree= trees[index]->get_structure();
   res= new One_int_tree(*utree, v);
   delete utree;
   tie(it, end)= res->vertices();
   while (it < end)
   {
      v.Int()= *it;
      res->put(*it++, v);
   }

   return res;
}

/*****************************************************************
 *
 *  Access to the histograms of topological indicators
 *  for Observed_trees
 *
 **/

template<typename Generic_Int_fl_container>
Distribution_data* Observed_trees<Generic_Int_fl_container>::extract_size() const
{  // access to the histogram of tree sizes
   Distribution_data *res= NULL;

   if (hsize != NULL)
      res= new Distribution_data(*hsize);

   return res;
}

template<typename Generic_Int_fl_container>
Distribution_data* Observed_trees<Generic_Int_fl_container>::extract_nb_children() const
{ // access to the histogram of tree sizes
   Distribution_data *res= NULL;

   if (hnb_children != NULL)
      res= new Distribution_data(*hnb_children);

   return res;
}

/*****************************************************************
 *
 * Return the maximal value for each variable of Observed_trees
 *
 **/

template<typename Generic_Int_fl_container> Generic_Int_fl_container
Observed_trees<Generic_Int_fl_container>::get_max_value() const
{
   Generic_Int_fl_container res= _max_value;
   return res;
}

/*****************************************************************
 *
 * Return the minimal value for each variable of Observed_trees
 *
 **/

template<typename Generic_Int_fl_container> Generic_Int_fl_container
Observed_trees<Generic_Int_fl_container>::get_min_value() const
{
   Generic_Int_fl_container res= _min_value;
   return res;
}

/*****************************************************************
 *
 *  Access to number of values of a finite variable
 *
 **/

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_nb_values(int variable) const
{
  assert(_type[variable] != REAL_VALUE);
  // the number of values for a float random variable does not make much sense
  return (_max_value.Int(variable)-_min_value.Int(variable)+1);
}

/*****************************************************************
 *
 *  Access to the histogram for the marginal distribution
 *
 **/

template<typename Generic_Int_fl_container>
Histogram* Observed_trees<Generic_Int_fl_container>::get_marginal(int variable) const
{
   Histogram *res;

   if (characteristics != NULL)
      if (characteristics[variable] != NULL)
         if (characteristics[variable]->get_marginal() != NULL)
            res= new Histogram(*(characteristics[variable]->get_marginal()));

   return res;
}

/*****************************************************************
 *
 *  Access to the number of observed trees and the trees
 *
 **/

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::get_nb_trees()
const { return _nb_trees; }

/*****************************************************************
 *
 *  Access to the characteristic quantity distributions
 *
 **/

template<typename Generic_Int_fl_container>
Tree_characteristics** Observed_trees<Generic_Int_fl_container>::get_characteristics() const
{
   Tree_characteristics **res;
   int var;

   if (characteristics != NULL)
   {
      res= new Tree_characteristics*[this->_nb_variable];

      for(var= 0; var < this->_nb_variable; var++)
      {
         if (characteristics[var] != NULL)
            res[var]= new Tree_characteristics(*(characteristics[var]));
         else
            res[var]= NULL;
      }
   }
   else
      res= NULL;

   return res;
}

template<typename Generic_Int_fl_container>
Tree_characteristics* Observed_trees<Generic_Int_fl_container>::get_characteristics(int variable) const
{
   Tree_characteristics *res;

   if (characteristics != NULL)
   {
      if (characteristics[variable] != NULL)
         res= new Tree_characteristics(*(characteristics[variable]));
      else
         res= NULL;
   }
   else
      res= NULL;

   return res; // characteristics[variable];
}

/*****************************************************************
 *
 *  Access to the individual histograms for the characteristic
 *  quantity distributions
 *
 **/

template<typename Generic_Int_fl_container>
Histogram* Observed_trees<Generic_Int_fl_container>::get_first_occurrence_root(int variable,
                                                                               int value) const
{
   Histogram *res;

   if (characteristics != NULL)
      if (characteristics[variable] != NULL)
         if (characteristics[variable]->get_first_occurrence_root() != NULL)
            res= new Histogram(*(characteristics[variable]->get_first_occurrence_root(value)));

   return res;
}

template<typename Generic_Int_fl_container>
Histogram* Observed_trees<Generic_Int_fl_container>::get_first_occurrence_leaves(int variable,
                                                                                 int value) const
{
   Histogram *res;

   if (characteristics != NULL)
      if (characteristics[variable] != NULL)
         if (characteristics[variable]->get_first_occurrence_leaves() != NULL)
            res= new Histogram(*(characteristics[variable]->get_first_occurrence_leaves(value)));

   return res;
}

template<typename Generic_Int_fl_container>
Histogram* Observed_trees<Generic_Int_fl_container>::get_sojourn_size(int variable,
                                                                      int value) const
{
   Histogram *res;

   if (characteristics != NULL)
      if (characteristics[variable] != NULL)
         if (characteristics[variable]->get_sojourn_size() != NULL)
            res= new Histogram(*(characteristics[variable]->get_sojourn_size(value)));

   return res;
}

template<typename Generic_Int_fl_container>
Histogram* Observed_trees<Generic_Int_fl_container>::get_nb_zones(int variable, int value) const
{
   Histogram *res;

   if (characteristics != NULL)
      if (characteristics[variable] != NULL)
         if (characteristics[variable]->get_sojourn_size() != NULL)
            res= new Histogram(*(characteristics[variable]->get_nb_zones(value)));

   return res;
}

template<typename Generic_Int_fl_container>
Histogram* Observed_trees<Generic_Int_fl_container>::get_nb_occurrence(int variable, int value) const
{
   Histogram *res;

   if (this->characteristics != NULL)
      if (this->characteristics[variable] != NULL)
         if (this->characteristics[variable]->get_nb_occurrence() != NULL)
            res= new Histogram(*(this->characteristics[variable]->get_nb_occurrence(value)));

   return res;
}

/*****************************************************************
 *
 *  Access to the individual observed_trees or
 *  to the set of observed_trees (or state trees)
 *
 **/

template<typename Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>** Observed_trees<Generic_Int_fl_container>::get_trees() const
{
   int t;
   Int_fl_tree<Generic_Int_fl_container> **res_trees;

   if (trees != NULL)
   {
      res_trees= new Int_fl_tree<Generic_Int_fl_container>*[_nb_trees];
      for(t= 0; t < _nb_trees; t++)
      {
         if (trees[t] != NULL)
            res_trees[t]= new Int_fl_tree<Generic_Int_fl_container>(*(trees[t]));
         else
            res_trees[t]= NULL;
      }
   }
   else
      res_trees= NULL;

   return res_trees;
}

template<typename Generic_Int_fl_container>
Int_fl_tree<Generic_Int_fl_container>* Observed_trees<Generic_Int_fl_container>::get_tree(int itree) const
{
   Int_fl_tree<Generic_Int_fl_container> *res_tree;

   assert(itree < _nb_trees);

   if (trees != NULL)
   {
      if (trees[itree] != NULL)
         res_tree= new Int_fl_tree<Generic_Int_fl_container>(*(trees[itree]));
      else
         res_tree= NULL;
   }
   else
      res_tree= NULL;

   return res_tree;
}

/*****************************************************************
 *
 *  Constructor of Observed_trees class using the number of
 *  integral and floating variables, the type of the variables,
 *  the number of trees and a flag on the field allocation
 *
 **/

template<typename Generic_Int_fl_container>
Observed_trees<Generic_Int_fl_container>::Observed_trees(int inb_integral,
                                                         int inb_float,
                                                         int_array itype,
                                                         int inb_trees,
                                                         bool init_flag)
 : _nb_integral(inb_integral)
 , _nb_float(inb_float)
 , _type(NULL)
 , _max_size(0)
 , _max_depth(0)
 , _nb_trees(inb_trees)
 , hsize(NULL)
 , hnb_children(NULL)
 , trees(NULL)
 , characteristics(NULL)
{ init(_nb_integral, _nb_float, itype, _nb_trees, init_flag); }

/*****************************************************************
 *
 *  Deallocation of the pointers for Observed_trees
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::remove()
{
   int var, t; //, j;

   if (trees != NULL)
   {
       for(t= 0; t < _nb_trees; t++)
       {
          if (trees[t]  != NULL)
          {
             delete trees[t];
             trees[t]= NULL;
          }
       }
       delete [] trees;
       trees= NULL;
   }

   if (characteristics != NULL)
   {
      for(var= 0; var < _nb_integral; var++)
      {
         if (characteristics[var] != NULL)
         {
            delete characteristics[var];
            characteristics[var]= NULL;
         }
      }
      delete [] characteristics;
      characteristics= NULL;
   }

   if (hnb_children != NULL)
   {
      delete hnb_children;
      hnb_children= NULL;
   }

   if (hsize != NULL)
   {
      delete hsize;
      hsize= NULL;
   }

   if (_type != NULL)
   {
      delete [] _type;
      _type= NULL;
   }

}

/*****************************************************************
 *
 *  Copy operator of Observed_trees class
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::copy(const Observed_trees& otrees)
{ // copy operator
   int var= 0, t; //nb_int= 0, j;

   _nb_integral= otrees._nb_integral;
   _nb_float= otrees._nb_float;
   _max_size= otrees._max_size;
   _max_depth= otrees._max_depth;
   _nb_trees= otrees._nb_trees;


   if (otrees.trees != NULL)
   {
      trees= new Int_fl_tree<Generic_Int_fl_container>*[_nb_trees];
      for(t= 0; t < _nb_trees; t++)
      {
          if (otrees.trees[t] != NULL)
             trees[t]= new Int_fl_tree<Generic_Int_fl_container>(*(otrees.trees[t]));
          else
             trees[t]= NULL;
      }
   }
   else
      trees= NULL;

   _type= new int[_nb_integral+_nb_float];

   for(var= 0; var < _nb_integral+_nb_float; var++)
      _type[var]= otrees._type[var];

   if (otrees.characteristics != NULL)
   {
      characteristics= new Tree_characteristics*[_nb_integral];
      for(var= 0; var < _nb_integral; var++)
      {
         if (otrees.characteristics[var] != NULL)
            characteristics[var]= new Tree_characteristics(*(otrees.characteristics[var]));
         else
            characteristics[var]= NULL;
      }
   }
   else
      characteristics= NULL;

   if (otrees.hsize != NULL)
      hsize= new Histogram(*(otrees.hsize));
   else
      hsize= NULL;

   if (otrees.hnb_children != NULL)
      hnb_children= new Histogram(*(otrees.hnb_children));
   else
      hnb_children= NULL;

   _min_value= otrees._min_value;
   _max_value= otrees._max_value;
}

/*****************************************************************
 *
 *  Field allocation for Observed_trees
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::init(int inb_integral,
                                                    int inb_float,
                                                    int_array itype,
                                                    int inb_trees,
                                                    bool init_flag)
{
   int nb_variables= inb_integral+inb_float;
   register int var, t;

   _nb_integral= inb_integral;
   _nb_float= inb_float;
   _nb_trees= inb_trees;

   _type= new int[nb_variables];
   if (itype != NULL)
      for(var= 0; var < nb_variables; var++)
         _type[var]= itype[var];
   else
   {
      for(var= _nb_integral; var < nb_variables; var++)
         _type[var]= REAL_VALUE;
   }

   _min_value.reset(_nb_integral, _nb_float);
   _max_value.reset(_nb_integral, _nb_float);


   if ((_nb_trees > 0) && (init_flag))
   {
      // trees= new (Int_fl_tree<Generic_Int_fl_container>*)[_nb_trees];
      trees= new tree_type*[_nb_trees];
      for(t= 0; t < _nb_trees; t++)
         trees[t]= NULL;

      if (_nb_integral > 0)
      {
         characteristics= new Tree_characteristics*[_nb_integral];
         for(var= 0; var < _nb_integral; var++)
            characteristics[var]= NULL;
      }
   }
}

/*****************************************************************
 *
 *  Cluster the values of a given integral variable
 *  for Observed_trees, using Observed_trees,
 *  the variable index and the clustering step
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::cluster(const Observed_trees& otrees,
                                                       int ivariable,
                                                       int step)
{
   register int t; //, k;
   vertex_iterator it, end;
   value v;

   for(t= 0; t < _nb_trees; t++)
   {
      tie(it, end)= trees[t]->vertices();
      // clustering of the values for variable ivariable
      while (it < end)
      {
         v= (otrees.trees[t])->get(*it);
         v.Int(ivariable)/= step;
         trees[t]->put(*it++, v);
      }
   }
   _min_value.Int(ivariable)= (otrees._min_value).Int(ivariable) / step;
   _max_value.Int(ivariable)= (otrees._max_value).Int(ivariable) / step;
   build_characteristics(ivariable);
}

/*****************************************************************
 *
 *  Transcoding of the values of a given integral variable
 *  for Observed_trees, using Observed_trees, the variable index,
 *  the smallest and greatest symbols and a transcoding table
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::transcode(const Observed_trees& otrees,
                                                         int ivariable,
                                                         int min_symbol,
                                                         int max_symbol,
                                                         int_array symbol)
{
   register int t;
   vertex_iterator it, end;
   value v;

   for(t= 0; t < _nb_trees; t++)
   {
      tie(it, end)= (otrees.trees[t])->vertices();
      while (it < end)
      {
         v= (otrees.trees[t])->get(*it);
         v.Int(ivariable)= symbol[v.Int(ivariable)-otrees._min_value.Int(ivariable)]
                             + min_symbol;
         trees[t]->put(*it++, v);
      }
   }

   _min_value.Int(ivariable)= min_symbol;
   _max_value.Int(ivariable)= max_symbol;
   build_characteristics(ivariable);
}

/*****************************************************************
 *
 *  Variable selection for Observed_trees using Observed_trees
 *  and the array of selected variables
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::select_variable(const Observed_trees& otrees,
                                                               int_array variable)
{
   register int t, var;
   value vsource, vdest;
   vertex_iterator it, end;
   tree_type *ctree; // current_tree
   Unlabelled_tree *utree= NULL;

   vdest.reset(_nb_integral, _nb_float);
   trees= new tree_type*[_nb_trees];
   for(t= 0; t < _nb_trees; t++)
   {
      ctree= otrees.trees[t];
      utree= ctree->get_structure();
      trees[t]= new tree_type(*utree, vdest);
      delete utree;
      utree= NULL;
      // should ensure that the number of variables is correct
      tie(it, end)= ctree->vertices();
      while(it < end)
      {
         vsource= ctree->get(*it);
         for(var= 0; var < _nb_integral; var++)
            vdest.Int(var)= vsource.Int(variable[var]);
         for(var= 0; var < _nb_float; var++)
            vdest.Double(var)= vsource.Double(variable[var]-_nb_float);
         trees[t]->put(*it++, vdest);
      }
   }

   if ((_nb_integral > 0) && (otrees.characteristics != NULL))
      characteristics= new Tree_characteristics*[_nb_integral];

   for(var= 0; var < _nb_integral; var++)
   {
      if (_type[var] != POSITION)
      {
        _min_value.Int(var)= otrees._min_value.Int(variable[var]);
        _max_value.Int(var)= otrees._max_value.Int(variable[var]);
        if (otrees.characteristics[variable[var]] != NULL)
           characteristics[var]= new Tree_characteristics(*(otrees.characteristics[variable[var]]));
           // new Tree_characterics(*(otrees.characteristics[variable[var]]));
      }
   }

   for(var= 0; var < _nb_float; var++)
   {
     _min_value.Double(var)= otrees._min_value.Double(variable[var]-_nb_float);
     _max_value.Double(var)= otrees._max_value.Double(variable[var]-_nb_float);
   }

   _max_size= otrees._max_size;
   _max_depth= otrees._max_depth;

   if (otrees.hsize != NULL)
      hsize= new Histogram(*otrees.hsize);

   if (otrees.hnb_children != NULL)
      hnb_children= new Histogram(*otrees.hnb_children);

}

/*****************************************************************
 *
 *  Gnuplot output of Observed_trees in the case of no other characteristic
 *  than the marginal histogram, using a prefix for the files,
 *  the title of output figures, the considered variable
 *  and the total number of variables
 *
 **/

template<typename Generic_Int_fl_container>
bool Observed_trees<Generic_Int_fl_container>::plot_print(const char * prefix,
                                                          const char * title,
                                                          int variable,
                                                          int nb_variables) const
{
   bool status= false;
   register int i;
   const Histogram *phisto[1];
   ostringstream data_file_name;

   // print the data file
   data_file_name << prefix << variable+1 << ".dat";

   phisto[0]= hsize;
   if (((characteristics != NULL) && (characteristics[variable] != NULL))
         && (characteristics[variable]->marginal != NULL))
            status= true;

   if (status)
   {
      status= characteristics[variable]->marginal->plot_print((data_file_name.str()).c_str(), 1, phisto);

      if (status)
      {
         for (i= 0; i < 2; i++)
         {
            ostringstream file_name[2];

            switch (i)
            {
               case 0 :
               {
                  if (nb_variables == 1)
                     file_name[0] << prefix << ".plot";
                  else
                     file_name[0] << prefix << variable+1 << ".plot";
                  break;
               }

               case 1 :
               {
                  if (nb_variables == 1)
                     file_name[0] << prefix << ".print";
                  else
                     file_name[0] << prefix << variable + 1 << ".print";
                  break;
               }
            }

            ofstream out_file((file_name[0].str()).c_str());

            if (i == 1)
            {
               out_file << "set terminal postscript" << endl;

               if (nb_variables == 1)
                  file_name[1] << label(prefix) << ".ps";
               else
                  file_name[1] << label(prefix) << variable+1 << ".ps";
               out_file << "set output \"" << file_name[1].str() << "\"\n\n";
            }

            out_file << "set border 15 lw 0\n" << "set tics out\n" << "set xtics nomirror\n"
                     << "set title";
            if (title != NULL)
               out_file << " \"" << title << "\"";

            out_file << "\n\n";

            if (characteristics[variable]->marginal->nb_value-1 < TIC_THRESHOLD)
               out_file << "set xtics 0,1" << endl;

            if ((int)(characteristics[variable]->marginal->max * YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics 0,1" << endl;

            out_file << "plot [0:" << MAX(characteristics[variable]->marginal->nb_value-1 , 1) << "] [0:"
                     << (int)(characteristics[variable]->marginal->max * YSCALE)+1 << "] \""
                     << label((data_file_name.str()).c_str()) << "\" using 1 title \""
                     << STAT_label[STATL_VARIABLE] << " " << variable+1 << " - "
                     << STAT_label[_type[variable] == STATE ? STATL_STATE : STATL_OUTPUT] << " "
                     << STAT_label[STATL_HISTOGRAM] << "\" with impulses" << endl;

            if (characteristics[variable]->marginal->nb_value-1 < TIC_THRESHOLD)
               out_file << "set xtics autofreq" << endl;

            if ((int)(characteristics[variable]->marginal->max*YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics autofreq" << endl;

            if (i == 0)
               out_file << "\npause -1 \"" << STAT_label[STATL_HIT_RETURN] << "\"" << endl;

            out_file << endl;

            if (hsize->nb_value-1 < TIC_THRESHOLD)
               out_file << "set xtics 0,1" << endl;

            if ((int)(hsize->max*YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics 0,1" << endl;

            out_file << "plot [0:" << hsize->nb_value-1 << "] [0:"
                     << (int)(hsize->max*YSCALE)+1 << "] \""
                     << label((data_file_name.str()).c_str()) << "\" using 2 title \""
                     << STAT_TREES_label[STATL_TREE_SIZE] << " " << STAT_label[STATL_HISTOGRAM]
                     << "\" with impulses" << endl;

            if (hsize->nb_value-1 < TIC_THRESHOLD)
               out_file << "set xtics autofreq" << endl;

            if ((int)(hsize->max*YSCALE)+1 < TIC_THRESHOLD)
               out_file << "set ytics autofreq" << endl;

            if (i == 1)
               out_file << "\nset terminal x11" << endl;

            out_file << "\npause 0 \"" << STAT_label[STATL_END] << "\"" << endl;
         }
      }
   }
   return status;
}

/*****************************************************************
 *
 *  Computation of the min and max value taken by an observed tree,
 *  for each variable (including the state variable)
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::min_max_value_computation()
{
    int t, variable;
    vertex_iterator it, end;
    value v;
#   ifdef DEBUG
    int nb_integral= get_nb_int();
    cerr << "Number of integer variables: " << nb_integral << endl;
#   endif

    for(variable= 0; variable < _nb_integral; variable++)
    {
        _min_value.Int(variable)= INT_MAX;
        _max_value.Int(variable)= INT_MIN;
    }

    for(variable= 0; variable < _nb_float; variable++)
    {
        _min_value.Double(variable)= FLT_MAX;
        _max_value.Double(variable)= -FLT_MAX;
    }

    for(t= 0; t < _nb_trees; t++)
    {
         tie(it, end)= trees[t]->vertices();
         while (it < end)
         {
                 v= trees[t]->get(*it);
                 for(variable= 0; variable < _nb_integral; variable++)
                 {
                     _min_value.Int(variable)= min(_min_value.Int(variable),
                                                   v.Int(variable));
                     _max_value.Int(variable)= max(_max_value.Int(variable),
                                                   v.Int(variable));
                 }

                 for(variable= 0; variable < _nb_float; variable++)
                 {
                      _min_value.Double(variable)= min(_min_value.Double(variable),
                                                       v.Double(variable));
                      _max_value.Double(variable)= max(_max_value.Double(variable),
                                                       v.Double(variable));
                 }
                 it++;
         }

    }
}

/*****************************************************************
 *
 *  Computation of the max size and depth of an observed tree
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::max_size_computation()
{
   _max_size= 0;
   int t;

   if (trees != NULL)
   {
      for(t= 0; t < _nb_trees; t++)
         if (trees[t] != NULL)
            _max_size= max(_max_size, trees[t]->get_size());
   }
}
template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::max_depth_computation()
{
   _max_depth= 0;
   int t;

   if (trees != NULL)
   {
      for(t= 0; t < _nb_trees; t++)
         if (trees[t] != NULL)
            _max_depth= max(_max_depth, trees[t]->get_depth());
   }
}

/*****************************************************************
 *
 *  Computation of the maximal number of children of an observed tree
 *
 **/

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::max_nb_children_computation()
{
   typedef typename Int_fl_tree<Generic_Int_fl_container>::vertex_iterator vertex_iterator;

   int t;
   unsigned int max_nb_children= 0;
   vertex_iterator it, end;

   if (trees != NULL)
   {
      for(t= 0; t < _nb_trees; t++)
         if (trees[t] != NULL)
         {
            tie(it, end)= trees[t]->vertices();
            while(it < end)
               max_nb_children= max(max_nb_children,
                                    trees[t]->get_nb_children(*it++));
         }
   }
   return max_nb_children;
}

/* template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::state_variable_init(int istate_variable)
{} */

/*****************************************************************
 *
 *  Computation of the histograms of the size and children number
 *  and their cumulated numbers for Observed_tree class
 *
 **/

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::cumul_size_computation() const
{  // computation of the cumulated tree sizes
   // i.e. the total number of vertices
   int res= 0;
   register int t;

   for(t= 0; t < _nb_trees; t++)
      res+= trees[t]->get_size();

   return res;
}

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::build_size_histogram()
{  // computation of the size histogram
   register int t;

   if (hsize != NULL)
      delete hsize;

   hsize= new Histogram(_max_size+1);
   hsize->nb_element= _nb_trees;
   for(t= 0; t < _nb_trees; t++)
      (hsize->frequency[trees[t]->get_size()])++;

   hsize->nb_value_computation();
   hsize->offset_computation();
   hsize->max_computation();
   hsize->mean_computation();
   hsize->variance_computation();
}

template<typename Generic_Int_fl_container>
int Observed_trees<Generic_Int_fl_container>::cumul_nb_children_computation() const
{  // computation of the cumulated number of children of each node
   int res= 0;
   register int t;
   vertex_iterator it, end;

   for(t= 0; t < _nb_trees; t++)
      if (trees[t] != NULL)
      {
         tie(it, end)= trees[t]->vertices();
         while (it < end)
            res+= trees[t]->get_nb_children(*it++);
      }
   return res;
}

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::build_nb_children_histogram()
{  // computation of histogram for the number of children
   register int t;
   int nb_vertices= cumul_size_computation();
   vertex_iterator it, end;

   if (hnb_children != NULL)
      delete hnb_children;

   hnb_children= new Histogram(max_nb_children_computation()+1);
   hnb_children->nb_element= nb_vertices;
   for(t= 0; t < _nb_trees; t++)
      if (trees[t] != NULL)
      {
         tie(it, end)= trees[t]->vertices();
         while (it < end)
            (hnb_children->frequency[trees[t]->get_nb_children(*it++)])++;
      }

   hnb_children->nb_value_computation();
   hnb_children->offset_computation();
   hnb_children->max_computation();
   hnb_children->mean_computation();
   hnb_children->variance_computation();
}

/*****************************************************************
 *
 *  Computation of the characteristic quantity histograms
 *  of Observed_tree class for one given discrete variable
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::build_characteristics(int variable)
{
   One_int_tree *otrees1= new One_int_tree[_nb_trees];
   int t; //i

   assert((variable < _nb_integral) && (characteristics != NULL));

   for(t= 0; t < _nb_trees; t++)
       otrees1[t]= *(trees[t]->select_int_variable(variable));

   if (characteristics[variable] != NULL)
   {
      delete characteristics[variable];
      characteristics[variable]= NULL;
   }

   characteristics[variable]= new Tree_characteristics(_min_value.Int(variable),
                                                       _max_value.Int(variable),
                                                       _max_size,
                                                       _max_depth,
                                                       _nb_trees,
                                                       otrees1,
                                                       variable);
   delete [] otrees1;
   otrees1= NULL;
}


/*****************************************************************
 *
 *  Computation of the characteristic quantity histograms
 *  of Observed_trees class for each discrete variable
 *
 **/

template<typename Generic_Int_fl_container>
void Observed_trees<Generic_Int_fl_container>::build_characteristics()
{
   int var; //i, t;

   _min_value.reset(_nb_integral, _nb_float);
   _max_value.reset(_nb_integral, _nb_float);
   min_max_value_computation();
   max_size_computation();
   max_depth_computation();

   // for(var= _nb_integral; var < _nb_integral; var++)
   //    build_value_histogram(var);

   if (_nb_integral > 0)
   {
       if (characteristics == NULL)
       {
          characteristics= new Tree_characteristics*[_nb_integral];
          for(var= 0; var < _nb_integral; var++)
              characteristics[var] = NULL;
       }

       for(var= 0; var < _nb_integral; var++)
          build_characteristics(var);
   }
   else
      characteristics= NULL;
}

/*****************************************************************
 *
 *  Left (bit) shift operator of Observed_trees
 *
 **/

template <typename Generic_Int_fl_container>
std::ostream& operator<<(std::ostream& os,
                         const Observed_trees<Generic_Int_fl_container>& otrees)
{ return otrees.ascii_write(os); }

#endif
