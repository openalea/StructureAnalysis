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


#include "treenode.h"
using namespace std;


TreeNode::TreeNode(MTG& mtg ,int number, int depth, int father, int vertex, int complex)
{
  _number   = number;
  _depth    = depth;
  _father   = father;
  _vertex   = vertex;
  _complex  = complex;
  _mtg      = &mtg;
  _values   = (ValueVector*) NULL;
  _nb_value = 0;
}

TreeNode::~TreeNode()
{
  if (_values) delete (ValueVector*) _values;
}

void TreeNode::put(ostream& os)
{
        os<<"NUMBER   : "<<_number<<std::endl;
        os<<"NUMPOSTFIX: "<<_numPostfix<<std::endl;
        os<<"DEPTH    : "<<_depth<<std::endl;
        os<<"FATHER   : "<<_father<<std::endl;
        os<<"VERTEX   : "<<_vertex<<std::endl;
        os<<"COMPLEX  : "<<_complex;
        os<<"VALUE    : "<<_value;
        os<<"ORDER    : "<<_order;
        os<<"POSITION : "<<_position;
}

void TreeNode::print()
{
        cout<<"NUMBER    : "<<_number<<endl;
        cout<<"DEPTH     : "<<_depth<<endl;
        cout<<"FATHER    : "<<_father<<endl;
        cout<<"VERTEX    : "<<_vertex<<endl;
//         cout<<"COMPLEX   : "<<_complex<<endl;
//         cout<<"VALUE     : "<<_value<<endl;
//         cout<<"ORDER     : "<<_order<<endl;
//         cout<<"POSITION  : "<<_position<<endl;
//         cout<<" NB VALUE : "<<_nb_value<<endl;
//         if (_nb_value)
//         {
//           for (int i=0;i<_nb_value;i++)
//           {
//             cout<< " VALUE("<<i<<")="<<(*_values)[i]<<endl;
//           }
//         }
}

ostream& operator<<(ostream& os ,TreeNode node)
{
        node.put(os);
        return(os);
};

void TreeNode::resize(int new_size)
{
  if (_values)
  {
    _values->resize(new_size);
  }
   else
  {
    _values=new ValueVector(new_size,DIST_UNDEF);
  }
  _nb_value=new_size;
}

int TreeNode::getValueSize() const
{
  return(_nb_value);
}

DistanceType TreeNode::getValue(int index) const
{
  if (_values)
  {
    assert((index>=0)&&(index<_nb_value));
    return((*_values)[index]);
  }
  else
  {
    return(DIST_UNDEF);
  }
}

void TreeNode::putValue(int index, DistanceType new_value)
{
  assert((index>=0)&&(index<_nb_value));
  if (_values) { (*_values)[index]=new_value; }
}



