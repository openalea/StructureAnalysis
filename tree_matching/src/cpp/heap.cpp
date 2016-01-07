/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *
 *       $Source$
 *       $Id: heap.cpp 5532 2008-09-12 09:40:59Z boudon $
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


#include "heap.h"
using namespace std;

const int NO_FATHER=0;
const int NO_CHILD=0;

// Constructor
Heap::Heap()
{
  _order=10;
}

Heap::Heap(float heap_order)
{
  Item* FIRST=new Item(0,-1);

  _order=SINLT(heap_order);
  FIRST->putNumber(0);
  push_back(FIRST);
  //delete (Item*) FIRST;
  //at(0)->putNumber(0);
}

Heap::Heap(int heap_order,Item* first_elem)
{
  Item* FIRST=new Item(0,-1);

  _order=heap_order;
  FIRST->putNumber(0);
  push_back(FIRST);
  //at(0)->putNumber(0);
  first_elem->putNumber(1);
  push_back(first_elem);
  //delete (Item*) FIRST;
}

// Destructor
Heap::~Heap()
{
  //cout<<nbItem()<<endl;
  while (nbItem()>=0){
    delete (Item*) at(nbItem());
    //cout<<nbItem()<<endl;
    pop_back();
  }
  clear();
}

//-----------------------------------------------------------------------------//
// Useful functions

void Heap::replace(Item& old_item,const Item& new_item)
{
  old_item.putKey(new_item.getKey());
  old_item.putObject(new_item.getObject());
}


int Heap::father(int item_number) const
{
  int father_number=SINLT((float) (item_number-1)/((float) _order));
  assert(validIndex(item_number));
  assert(validIndex(father_number));
  if (father_number==0)
    {return(NO_FATHER);}
  else
    {
      return(father_number);
    }
}

int Heap::child(int item_number,int child_number) const
{
  int child_item_number=_order*(item_number-1)+child_number+1;

  assert(validIndex(item_number));
  if ((child_item_number>nbItem())||(child_number>_order))
    {
      return(NO_CHILD);
    }
  else
    {
      assert(child_item_number);
      return(child_item_number);
    }
}

int Heap::empty() const
{
  return(size()==1);
}

Item* Heap::at(size_type size) const
{
  vector<Item*>::const_iterator beg;
  beg = begin();
  for (size_type i=0;i<size;i++)
    beg++;
  return(*beg);
  //delete beg;
}

KeyType Heap::childKey(int item_number,int child_number) const
{
  return(at(child(item_number,child_number))->getKey());
}

KeyType Heap::fatherKey(int item_number) const
{
  return(at(father(item_number))->getKey());
}


int Heap::integrity() const
{
  int integrity=1;
  int i=1;

  while((!empty())&&(validIndex(i))&&(integrity))
    {
      if (at(i)->getNumber()>nbItem()){integrity=0;}
      if (!integrity) {cout<<"index : "<<i<< std::endl;at(i)->print();}
      i++;
    }
  return(integrity);
}

void Heap::print() const
{
  size_t index;
  for (index=1;index!=size();index++)
    {
      at(index)->print();
    }
}

int Heap::validIndex(int index) const
{
  return((index>=0)&&(index<(int)size()));
}

int Heap::position(ItemType object) const
{
  int pos=0;
  int index=1;
  while ((!pos)&&(validIndex(index)))
    {
      if (at(index)->getObject()==object) {pos=at(index)->getNumber();}
      index++;
    }
  return(pos);
}

int Heap::nbItem() const
{
  return(size()-1);
}

//-----------------------------------------------------------------------------//

// Special function
ItemType Heap::quickDeleteMin()
{
  ItemType min_object=at(1)->getObject();


  if(empty())
    {
      return(-1);
    }
  else
    {
      min_object=at(1)->getObject();
      if (nbItem()>1)
// <<<<<<< heap.cpp
	{               
	  siftDown(*at(nbItem()),1);
	}
      delete (Item*) at(nbItem());
// =======
//         {
//           siftDown(*at(nbItem()),1);
//         };
// >>>>>>> 1.3
      pop_back();
      return(min_object);
    }
}

bool Heap::order()
{
  int i;
  for (i=nbItem();i>0;i--)
    {
      siftDown(*at(i),i);
    }
  return(true);
}

bool Heap::quickInsert(KeyType key,ItemType object)
{
  Item* new_item=new Item(key,object);

  new_item->putNumber(nbItem()+1);
  push_back(new_item);
  //delete (Item*) new_item;
  return(false);
}


// Return the item of minimum key
Item* Heap::findmin() const
{
  if (empty())
    {return(at(0));}
  else
    {return(at(1));}
}

// Insert a new item
int Heap::insertItem(KeyType key,ItemType object)
{
  Item* new_item=new Item(key,object);
  int number_of_item=nbItem()+1;
  int position;

  new_item->putNumber(number_of_item);
  if (empty())
    {
      push_back(new_item);
      return(1);
    }
  else
    {
      push_back(new_item);
      position=siftUp(*new_item,number_of_item);
      return(position);
    }
  //delete (Item*) new_item;
}

// delete an item
int Heap::deleteItem(int item_to_delete_number)
{
  int last_item_number=nbItem();
  KeyType last_item_key;
  KeyType item_to_delete_key;
// <<<<<<< heap.cpp

// =======


// >>>>>>> 1.3
  assert(validIndex(item_to_delete_number));
  assert(validIndex(last_item_number));
  if (empty())
    {
      return(0);
    }
  else
    {
      last_item_key=at(last_item_number)->getKey();
      item_to_delete_key=at(item_to_delete_number)->getKey();
      if (item_to_delete_number!=last_item_number)
// <<<<<<< heap.cpp
	{
	  if (last_item_key <= item_to_delete_key)
	    {
	      siftUp(*at(last_item_number),item_to_delete_number);
	    }
	  else
	    {
	      siftDown(*at(last_item_number),item_to_delete_number);
	    }
	}
      delete (Item*) at(last_item_number);
// =======
//         {
//           if (last_item_key <= item_to_delete_key)
//             {
//               siftUp(*at(last_item_number),item_to_delete_number);
//             }
//           else
//             {
//               siftDown(*at(last_item_number),item_to_delete_number);
//             }
//         }
// >>>>>>> 1.3
      pop_back();
      return(1);
    }
}


// Delete the item of minum key
ItemType Heap::deleteMin()
{
  ItemType objet_min;

  if (empty())
    {
      return(-1);
    }
  else
    {
      objet_min=at(1)->getObject();
      deleteItem(1);
      return(objet_min);
    }
}

// Reorganize the heap during the insertion
int Heap::siftUp(Item& current_item,int current_number)
{
  KeyType father_Key=fatherKey(current_number);
  int father_number=father(current_number);
  KeyType item_Key=current_item.getKey();
  ItemType item_Object=current_item.getObject();

  assert(validIndex(current_number));
  while ((father(current_number))&&(father_Key > item_Key))
    {
      replace(*at(current_number),*at(father_number));
      current_number=father_number;
      assert(validIndex(current_number));
      father_number=father(current_number);
      assert(validIndex(father_number));
      father_Key=fatherKey(current_number);
    }
  at(current_number)->putKey(item_Key);
  at(current_number)->putObject(item_Object);

  return(current_number);
}



// Reorganize the heap during the deletion
int Heap::siftDown(Item& current_item,int current_number)
{
  int minchild_number=minChild(current_number);
  KeyType minchild_Key=at(minchild_number)->getKey();
  KeyType item_Key=(current_item).getKey();
  ItemType item_Object=current_item.getObject();

  assert(validIndex(current_number));
  while ((minchild_number!=0)&&(minchild_Key<item_Key))
    {
      replace(*at(current_number),*at(minchild_number));
      current_number=minchild_number;
      assert(validIndex(current_number));
      minchild_number=minChild(current_number);
      assert(validIndex(minchild_number));
      minchild_Key=at(minchild_number)->getKey();
    }
  at(current_number)->putKey(item_Key);
  at(current_number)->putObject(item_Object);

  return(current_number);
}


// Return the minimal key among the children's key
int Heap::minChild(int item_number)
{
  int child_index=1;
  KeyType child_Key=childKey(item_number,child_index);
  KeyType min_Key=child_Key;
  int min_child_number=1;

  assert(validIndex(item_number));
  while (child(item_number,child_index))
    {
      if (child_Key<min_Key) 
	{
	  min_Key=child_Key;
	  min_child_number=child_index;
	}
      child_index++;
      child_Key=childKey(item_number,child_index);
    }
  return(child(item_number,min_child_number));
}




