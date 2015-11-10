/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
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


#include "rw_cstring.h"
#include "util_hashmap.h"
using namespace std;

#ifdef RWOUT

/* ----------------------------------------------------------------------- */


VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */

istream& RWString::readLine(istream& s, bool skipWhite){
  if(skipWhite) ws(s);
  string st;
  getline( s, st );
  string::operator=(st);
  return s;
}

size_t RWString::hash() const
{
 tool_hash<std::string> my_hasher;
 return my_hasher(*this);
}

istream& RWString::readToDelim( istream& s, char delim)
{
  string st;
  getline( s, st, delim );
  string::operator=(st);
  return s;
}

RWString RWString::strip( stripType s, char c )
{
 size_t b= 0, e= size();
 if( s == leading || s == both )
   {
   b= this->find_first_not_of(c);
  if( b == npos ) b= 0;
   }
 else
  if( s == trailing || s == both )
    {
    e= this->find_last_not_of(c);
    if( e == npos ) e= size();
    }

 return substr(b,e-b);
}


/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#endif
