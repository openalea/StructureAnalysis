/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 * 
 *       $Source$
 *       $Id: treenode.cpp 3258 2007-06-06 13:18:26Z dufourko $
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


#include "treenode.h"
using namespace std;


TreeNode::TreeNode(IdType id, IdType father):
  _id(id), _father(father), _values(), _depth(0)
{
}

TreeNode::~TreeNode()
{
}


void TreeNode::print() const
{
  cout<<"ID \t : "<<_id<<endl;
  cout<<"FATHER \t : "<<_father<<endl;
  cout<<"DEPTH \t : "<<_depth<<endl;
  cout<<"CHILD LIST \t : [ ";
  if(hasChild()){
	  for (size_t i=0; i < getChildNumber();i++){
		if (i>0) cout <<" , ";
		cout << getChild(i);
	  }
  }
  cout<<" ]"<<endl;

   for (int i=0;i<_values.size();i++)
     cout<<"ARG["<<i<<"] \t : "<<getValue(i)<<endl;
}



DistanceType TreeNode::getValue(int index) const
{
  if (!_values.empty())
  {
    return getTypedValue<DistanceType>(index);
    assert((index>=0)&&(index<_values.size()));
	return boost::any_cast<DistanceType>(_values[index]);
  }
  else  return(DIST_UNDEF);
}


TreeNode::IdType TreeNode::getChild(int id) const{
  assert ((id>=0) && (id<_childList.size()));
  return _childList[id];
}


/* ----------------------------------------------------------------------- */

static TreeNode::Factory * TREENODEFACTORY = NULL;

TreeNode::Factory::Factory(): __builder(NULL) { }


TreeNodePtr TreeNode::Factory::build(IdType id, IdType father_id)
{
	if(!__builder) return TreeNodePtr(new TreeNode(id, father_id));
	else return (*__builder)(id, father_id);
}

TreeNode::Factory& TreeNode::factory() {
	if(!TREENODEFACTORY) TREENODEFACTORY = new TreeNode::Factory();
	return *TREENODEFACTORY;
}
/* ----------------------------------------------------------------------- */
