/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guedon <yann.guedon@cirad.fr>
 *                        Thomas cokelaer <Thomas.Cokelaer@cirad.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_chain.cpp 6168 2009-04-01 16:42:29Z cokelaer $
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
using namespace stat_tool;



class ChainWrap
{


};


void class_chain()
{
    class_< Chain >("_Chain", "Chain")
    .def(init<char, int , bool>())
    .def(init<char, int, int , bool>())
    .def(init< optional<char, int, bool> >())

    .def(self_ns::str(self)) //__str__

    .add_property("type", &Chain::type, "returns type")
    .add_property("nb_state", &Chain::nb_state, "returns nb_state")
    .add_property("nb_row", &Chain::nb_row, "returns nb_row")
    .add_property("nb_component", &Chain::nb_component, "returns nb_component")

/*
 *     .def("get_component_nb_state", &Chain::component_nb_state, "component nb state")
    .def("get_state_type", &Chain::get_state_type, "get state type")
    .def("get_initial", &Chain::get_initial, "get initial")
*/
//    .def()    ;
    ;

/*
    void init(bool left_right , double self_transition);

    double likelihood_computation(const Chain_data &chain_data , bool initial_flag = true) const;
    void chi2_fit(const Chain_data &chain_data , Test &test) const;

    bool **accessibility;   // matrice d'accessibilite des etats
    int *component_nb_state;  // nombre d'etats par classe
    int **component;        // classes
    char *state_type;       // types des etats ('r' : recurrent,
                            // 't' : transitoire, 'a' : absorbant)
    double *initial;        // probabilites initiales
    double *cumul_initial;  // fonction de repartition correspondant
                            // au probabilites initiales
    double **transition;    // matrice des probabilites de transition
    double **cumul_transition;  // fonctions de repartition correspondant aux lignes
                                // de la matrice des probabilites de transition

     std::ostream& ascii_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    void create_cumul();
    void cumul_computation();
    void remove_cumul();
    void log_computation();

    bool** logic_transition_computation() const;
    bool connex_component_research(StatError &error , bool **ilogic_transition = 0) const;
    void graph_accessibility_computation(bool **ilogic_transition);
    void probability_accessibility_computation();
    void component_computation(bool **ilogic_transition = 0);

    void thresholding(double min_probability, bool semi_markov);

    int nb_parameter_computation(double min_probability = 0.) const;
    double chi2_value_computation(const Chain_data &chain_data) const;

*/



}


