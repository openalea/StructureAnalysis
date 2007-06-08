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


#include "amlnodefunction.h"
using namespace std;

NodeFunction::NodeFunction(string name,FNode* function,DistanceType default_value)
{
  _name=name;
  _function=function;
  _defaultValue=default_value;
}

DistanceType NodeFunction::operator() (VId vertex)
{
  assert(_function);

  DistanceType result=_defaultValue;

  AMObjVector argobjs(1); // Only one arg which is a vertex

  argobjs[0] = AMObj(AMObjType::VTX, vertex);

  assert(_function);

  AMObj result_obj = (*_function)(argobjs);

  switch(result_obj.tag())
    {
    case  AMObjType::INTEGER    : result = (DistanceType) result_obj.val.i;     break;
    case  AMObjType::REAL       : result = result_obj.val.r;                    break;
    default                     : result = _defaultValue;                       break;
    }
  return result;
}

DistanceType NodeFunction::operator() (VId vertex1, VId vertex2)
{
  assert(_function);

  DistanceType result=_defaultValue;

  AMObjVector argobjs(2); // two arg which is a vertex

  argobjs[0] = AMObj(AMObjType::VTX, vertex1);
  argobjs[1] = AMObj(AMObjType::VTX, vertex2);

  assert(_function);

  AMObj result_obj = (*_function)(argobjs);

  switch(result_obj.tag())
    {
    case  AMObjType::INTEGER    : result = (DistanceType) result_obj.val.i;     break;
    case  AMObjType::REAL       : result = result_obj.val.r;                    break;
    default                     : result = _defaultValue;                       break;
    }
  return result;
}

