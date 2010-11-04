/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Ch. Pradal (christophe.pradal@cirad.fr)
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

/*! \file rw_cstring.h
    \brief File for Rogue Wave CString with the STL.
*/

#ifndef __rw_cstring_h__
#define __rw_cstring_h__

#include "config.h"

#ifdef RWOUT

/* ----------------------------------------------------------------------- */


#include <string>
#include <iostream>
#include "tools_namespace.h"
#include "rw_defs.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class RWString.

  \brief A string.

*/


class VPTOOLS_API RWString : public std::string
{
 public:
 enum stripType { leading, trailing, both };

 public:
  // Visual C++ bug
  typedef std::basic_string<char> stringM;

  RWString() : stringM() {}
  RWString( const RWString& s) : stringM(s) {}
  RWString( const std::string& s)  : stringM(s) {}
  RWString( char s)  : stringM(1,s) {}
  RWString( const char* s)  : stringM(s) {}
  RWString( const char* s, std::size_t n )  : stringM(s,n) {}
  virtual ~RWString() {}

  operator const char*() const { return c_str(); }

  std::istream& readLine(std::istream& s, bool _skipWhite= true);

  bool isNull() const { return this->empty(); }

  const char* data() const { return c_str(); }

  size_t first(char c) { return find(c); }
  size_t first(const char* c) { return find(c); }

  size_t last(char c) { return rfind(c); }
  size_t last(const char* c) { return rfind(c); }

  void remove( size_t pos ) { erase(pos); }
  void remove( size_t pos, size_t N ) { erase(pos,N); }

  size_t hash() const;

  char operator()(size_t i) const { return operator[](i); }
  char& operator()(size_t i) { return operator[](i); }

  RWString operator()( size_t start, size_t len )
  { return std::string(*this,start,len); }


  std::istream& readToDelim( std::istream& s, char delim='\n');

  RWString& prepend( const char* cs )
    {
    this->insert(0,cs);
    return *this;
    }

  RWString strip( stripType s= RWString::trailing, char c= ' ' );

  bool contains( const std::string& s ) const { return ( this->find(s) != npos ); }
}; // class RWString


/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

typedef VPTOOLS(RWString) RWCString;

/* ----------------------------------------------------------------------- */

#else

#include <rw/cstring.h>

#endif
// rw_cstring_h
#endif
