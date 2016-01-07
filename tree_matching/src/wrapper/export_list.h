/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       PlantGL: The Plant Graphic Library
 *
 *       Copyright 1995-2007 UMR CIRAD/INRIA/INRA DAP 
 *
 *       File author(s): F. Boudon et al.
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

#ifndef __export_list_h__
#define __export_list_h__

#include <boost/python.hpp>
#include <vector>


/* ----------------------------------------------------------------------- */
template <class T>
struct make_object {
	typedef T obj_type;

    const obj_type& __c_obj;

    make_object(const T& c_obj) : __c_obj(c_obj) {}
    boost::python::object operator()() const { return boost::python::object(__c_obj); }
};

template <class T, class Translator = make_object<typename T::value_type> >
struct list_converter {
	typedef T list_type;
    typedef typename T::const_iterator list_const_iterator;
    typedef typename T::value_type list_element_type;

    list_const_iterator __c_list_begin;
    list_const_iterator __c_list_end;

    list_converter(const T& c_list): __c_list_begin(c_list.begin()),__c_list_end(c_list.end()){}
    list_converter(const list_const_iterator& c_list_begin,
              const list_const_iterator& c_list_end): 
            __c_list_begin(c_list_begin),__c_list_end(c_list_end){}

    boost::python::object operator()() const {
        boost::python::list l;
        for (list_const_iterator it = __c_list_begin; it != __c_list_end; ++it)
            l.append(Translator(*it)());
        return l;
    }
};

template <class T, class Translator>
list_converter<T,Translator> make_list(const T& c_list) { return list_converter<T,Translator>(c_list); }

template <class T>
list_converter<T> make_list(const T& c_list) { return list_converter<T>(c_list); }

#endif
