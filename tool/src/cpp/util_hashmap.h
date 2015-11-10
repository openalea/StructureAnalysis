/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       PlantGL: Plant Graphic Library
 *
 *       Copyright 1995-2003 UMR Cirad/Inria/Inra Dap - Virtual Plant Team
 *
 *       File author(s): F. Boudon
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


#ifndef __util_tool_hashmap_h__
#define __util_tool_hashmap_h__

#include "config.h"

/*! \file util_hashmap.h
    \brief Utility for hashmap with std::string.
*/

#ifdef USING_UNORDERED_MAP

#ifdef GNU_TR1_STL_EXTENSION
	#include <tr1/unordered_map>
	#define tool_hash_map std::tr1::unordered_map

	#ifndef tool_hash
		#define tool_hash std::tr1::hash
	#endif

#else
	#include <unordered_map>
	#define tool_hash_map std::unordered_map

	#ifndef tool_hash
		#define tool_hash std::hash
	#endif
#endif



template <class T>
struct tool_hash_map_string : public tool_hash_map<std::string, T >{};

#else

#ifdef GNU_TR1_STL_EXTENSION
	#include <ext/hash_map>
#else
	#if defined(__GNUC__)
		#warning GNU STL Extension not activated ! Old GCC version used ?
	#endif
	#include <hash_map>
#endif

#ifndef tool_hash
#define tool_hash STDEXT::hash
#endif
#define tool_hash_map STDEXT::hash_map

#include "util_hash.h"

#ifndef WIN32_STL_EXTENSION

/**
   \class hash_map_string
   \brief Class used for parsing a stream.
*/

template <class T>
struct tool_hash_map_string : public tool_hash_map<std::string, T, tool_hashstr, tool_eqstr>
{};

#else

template <class T>
struct tool_hash_map_string : public tool_hash_map<std::string, T >
{};

#endif

#endif

#endif
