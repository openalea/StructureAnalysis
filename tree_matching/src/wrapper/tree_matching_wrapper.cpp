/* ---------------------------------------------------------------------------
 #
 #      tree_matching
 #
 #       Copyright 2003-2008 UMR LaBRI
 #
 #       File author(s): P. Ferraro (pascal.ferraro@labri.fr)
 #
 # ---------------------------------------------------------------------------
 #
 #                      GNU General Public Licence
 #
 #       This program is free software; you can redistribute it and/or
 #       modify it under the terms of the GNU General Public License as
 #       published by the Free Software Foundation; either version 2 of
 #       the License, or (at your option) any later version.
 #
 #       This program is distributed in the hope that it will be useful,
 #       but WITHOUT ANY WARRANTY; without even the implied warranty of
 #       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 #       GNU General Public License for more details.
 #
 #       You should have received a copy of the GNU General Public
 #       License along with this program; see the file COPYING. If not,
 #       write to the Free Software Foundation, Inc., 59
 #       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 #
 # ---------------------------------------------------------------------------
 */

#include <iostream>
#include <boost/python.hpp>
#include "export_tree_matching.h"


using namespace boost::python;

void py_factory_finalize();

BOOST_PYTHON_MODULE(__tree_matching__)
{
  export_TreeNode();
  export_TreeGraph();
  export_NodeCost();
  export_Matching();
  export_ExtMatching();
  export_MatchPath();
  export_ChoiceTable();
  export_GeneralMatchPath();

  Py_AtExit(&py_factory_finalize);
};
