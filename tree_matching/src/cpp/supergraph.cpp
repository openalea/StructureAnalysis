/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.Ferraro (pascal.ferraro@cirad.fr)
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

#include "supergraph.h"

using namespace std;

// -------------
// Constructeur
// -------------
SuperGraph::SuperGraph(TreeGraph & input,TreeGraph & reference)
{
  T1=&input;
  T2=&reference;
  _size1 = T1->getNbVertex();
  _size2 = T2->getNbVertex();
}


// -------------
// Destructeur
// -------------
SuperGraph::~SuperGraph()
{
}

int SuperGraph::outdeg(node vertex)
{
  if (vertex <_size1)
    return (T1->getNbChild(vertex));
  else
    return (T2->getNbChild(vertex-_size1));
}

list<node> SuperGraph::getOutNodes(node vertex)
{
  list<node> children;
  NodeList* child_list;
  NodeList::const_iterator begin,end;
  if (vertex <_size1)
    {
      begin = (T1->sons(vertex))->begin();
      end   = (T1->sons(vertex))->end();
    }
  else
    {
      begin = (T2->sons(vertex-_size1))->begin();
      end   = (T2->sons(vertex-_size1))->end();
    }

  while (begin!=end)
    {
      if (vertex <_size1)
        children.push_back(*begin);
      else
        children.push_back(*begin+_size1);
      begin++;
    }
  return(children);
}


int SuperGraph::getNbVertex(int i)
{
  if (i==1)
    return(T1->getNbVertex());
  else
    return(T2->getNbVertex());
}
// <<<<<<< supergraph.cpp

VId SuperGraph::getVertexId(node vertex)
{
  if (vertex <_size1)
    return (T1->getNode(vertex))->getVertex();
  else
    return (T2->getNode(vertex-_size1))->getVertex();
}
// =======


// >>>>>>> 1.2
