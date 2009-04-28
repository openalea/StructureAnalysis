/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Thomas cokelaer <Thomas.Cokelaer@cirad.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_compound.cpp 6168 2009-04-01 16:42:29Z cokelaer $
 *
 *-----------------------------------------------------------------------------*/
#include <stdio.h>
#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/markovian.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost::python;

class ChainWrap
{


};


void class_chain()
{
	class_< Chain >("_Chain", "Chain")
    .def(init<char, int , bool>())
    .def(init<char, int, int , bool>())
/*    .def_readonly("get_type", &Chain::get_type, "get type")
    .def_readonly("get_nb_state", &Chain::get_nb_state, "get nb state")
    .def_readonly("get_nb_row", &Chain::get_nb_row, "get nb_row")
    .def_readonly("get_nb_component", &Chain::get_nb_row, "get nb_component")
    .def("get_component_nb_state", &Chain::component_nb_state, "component nb state")
    .def("get_state_type", &Chain::get_state_type, "get state type")
    .def("get_initial", &Chain::get_initial, "get initial")
*/
//    .def()    ;
    ;

/*

bool get_accessibility(int state1 , int state2) const{ return accessibility[state1][state2]; }
int get_component(int icomponent , int index) const{ return component[icomponent][index]; }
double get_transition(int memory , int state) const{ return transition[memory][state]; } 
bool **accessibility;   // matrice d'accessibilite des etats
int *component_nb_state;  // nombre d'etats par classe
int **component;        // classes
char *state_type;       // types des etats ('r' : recurrent,
double *cumul_initial;  // fonction de repartition correspondant
double **transition;    // matrice des probabilites de transition
double **cumul_transition;  // fonctions de repartition correspondant aux lignes
void create_cumul();
void cumul_computation();
void remove_cumul();
void log_computation();
bool** logic_transition_computation() const;
bool connex_component_research(Format_error &error , bool **ilogic_transition = 0) const;
void graph_accessibility_computation(bool **ilogic_transition);
void probability_accessibility_computation();
void component_computation(bool **ilogic_transition = 0);
void thresholding(double min_probability);
int nb_parameter_computation(double min_probability = 0.) const;
double chi2_value_computation(const Chain_data &chain_data) const;
void init(bool left_right , double self_transition);
double likelihood_computation(const Chain_data &chain_data , bool initial_flag = true) const;
void chi2_fit(const Chain_data &chain_data , Test &test) const;

*/



}


