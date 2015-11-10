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


#ifndef SB_HEAP_HEADER
#define SB_HEAP_HEADER

#include "definitions.h"
#include "item.h"


/**
 *\class Heap
 *\brief Object used for optimizing the computation of the minimum 
 * cost maximum flow problem (Tarjan89)
 *\author Pascal ferraro
 *\date 1999
 */

class Heap  :public VectorOfItem
{

	public :

// Constructor
	Heap();
	Heap(float );
	Heap(int,Item*);

// Destructor
	~Heap();

// Renvoie l'item minimum
        Item* findmin() const;
// Insert a new item
        int insertItem(KeyType,ItemType);
// Delete an item
        int deleteItem(int);
// Delete the item of minimal key 
        ItemType deleteMin();
// Reorganize the heap during the insertion 
        int siftUp(Item&,int );
// Reorganize the heap during the deletion
        int siftDown(Item&,int );
// Return the number of object currently in the heap
	int nbItem() const;
// Print the whole heap
	void print() const;
// Return the position of an object in the heap
	int position(ItemType object) const;

  // Return ith Item
  Item* at(size_type size) const;

// Special functions for pathgraph use only 
        AmlBoolean quickInsert(KeyType,ItemType);
	ItemType quickDeleteMin();
	AmlBoolean order();	

	private :

// Order of the Heap (stands for d in Tarjan's book)
        int _order;
// Give the father ot anitem in the heap
	int father(int ) const ;
// Give the i-th child of an item in the heap
	int child(int ,int ) const ;
// Return the item of minimal key
        int minChild(int );
// Replace an item of the heap by an other item of the heap
	void replace(Item& ,const Item& );
// Tell if the heap is empty or not
	int empty() const;
// Verify that the heap is handled the right way
	int validIndex(int index) const;
	int integrity() const;
// Return the key of a child
 	KeyType childKey(int item_number,int child_number) const;
// Return the father's key
	KeyType fatherKey(int item_number) const;

};

#endif


	

                   	

