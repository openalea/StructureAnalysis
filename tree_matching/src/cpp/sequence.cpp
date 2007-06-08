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


#include "sequence.h"

Sequence::Sequence()
{
  _size    = 0;
  _first   = (Relation*) NULL;
  _last    = (Relation*) NULL;
  _current = (Relation*) NULL;
  _used    = 0;
  _nbIns   = 0;
  _nbDel   = 0;
  _nbMat   = 0;
  _nbSub   = 0;
}

void Sequence::append(int imput_vertex,int reference_vertex,DistanceType cost)
{
  Relation* relation=new Relation(imput_vertex,reference_vertex,cost,0);
  if (_size)
    {
      _last->putNext(relation);
      _last=relation;
    }
  else
    {
      _first=relation;
      _current=relation;
      _last=relation;
    }
  _size++;
}


void Sequence::append(Relation* new_relation)
{
  new_relation->putNext(0);
  if (_size)
    {	
      _last->putNext(new_relation);
      _last=new_relation;
    }
  else
    {
      _first=new_relation;
      _current=new_relation;
      _last=new_relation;
    }
  _size++;
}	

void Sequence::clear()
{
  Relation* current;
  while ((_size)&&(_first))
    {
      current=_first->getNext();
      delete (Relation*) _first;
      _first=current;
      _size--;
    }
  _size=0;
  _first=0;
  _last=0;
  _current=0;
}

Sequence::~Sequence()
{
  if (_size)
  {
    do
    {
      Relation* current=_first->getNext();
      delete (Relation*) _first;
      _first=current;
      _size--;
    } while(_size);
  }
  _size=0;
  _first=0;
  _last=0;
  _current=0;
}

void Sequence::link(Sequence* sequence)
{
  if (sequence->getSize())
    {
      if (_size)
	{
	  _size=_size+sequence->getSize();
	  _last->putNext(sequence->getFirst());
	  _last=sequence->getLast();
	  _last->putNext(0);
	}
      else
	{
	  _size=sequence->getSize();
	  _first=sequence->getFirst();
	  _current=sequence->getCurrent();
	  _last=sequence->getLast();
	  _last->putNext(0);
	  _used=sequence->isUsed();
	}
    }
}

void Sequence::add(Sequence* sequence)
{
  if (sequence->getSize())
    {
      sequence->reset();
      do
	{
	  append(sequence->getCurrent()->getIV(),sequence->getCurrent()->getRV(),sequence->getCurrent()->getCost());
	}
      while(sequence->next());
    }
}

Relation* Sequence::getFirst() const
{
  return(_first);
}

Relation* Sequence::getLast() const
{
  return(_last);
}

int Sequence::getSize() const
{
  return(_size);
}

void Sequence::reset() 
{
  _current=_first;
}

int  Sequence::next()
{
  if (_current!=_last) 
    {
      _current=_current->getNext();
      return(1);
    }
  else
    {
      return(0);
    }
}

Relation* Sequence::getCurrent() const
{
  return(_current);
}

void Sequence::operator=(const Sequence sequence) 
{
  _size=sequence.getSize();
  _first=sequence.getFirst();
  _current=sequence.getCurrent();
  _last=sequence.getLast();
  _used=sequence.isUsed();
  /*Ajout des derniers champs de sequence */
  _nbIns=sequence.getNbIns();
  _nbDel=sequence.getNbDel();
  _nbMat=sequence.getNbMat();
  _nbSub=sequence.getNbSub();
}

void Sequence::print()
{
  reset();
  for (int i=0;i<_size;i++)
    {
      assert(_current);
      next();	
    }
}

	







