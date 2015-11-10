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


#include "item.h"
using namespace std;

Item::Item(KeyType key, ItemType object)
{
        _key=key;
        _object=object;
}

    Item::Item(const Item& item){ copy(item);}

void Item::copy(const Item& item)
{
        _key=item.getKey();
        _number=item.getNumber();
        _object=item.getObject();
}


const Item& Item::operator=(const Item& another_item)
{
        _key=another_item.getKey();
        _number=another_item.getNumber();
        _object=another_item.getObject();
        return *this;
}

Item::~Item()
{
};


 KeyType Item::getKey() const {
   return(_key);
 } 

int Item::getNumber() const {return(_number);}

ItemType Item::getObject() const {return(_object);}

void Item::putNumber(int new_number){ _number=new_number;};

void Item::putKey(const KeyType new_key)
{
        _key=new_key;
};

void Item::putObject(const ItemType new_object){ _object=new_object;};


int Item::operator==(const Item another_item) const
{
        return(_object==another_item.getObject());
}

int Item::operator!=(const Item another_item) const
{
        return(_object!=another_item.getObject());
}

int Item::operator<(const Item another_item) const
{
        return(_key>another_item.getKey());
}

int Item::operator>(const Item another_item) const
{
        return(_key<another_item.getKey());
}

void Item::print() const
{
        cout<<"DIJ :NUMBER : "<<getNumber()<<std::endl;
        cout<<"DIJ :KEY    : "<<getKey()<<std::endl;
        cout<<"DIJ :OBJECT : "<<getObject()<<std::endl;
}


