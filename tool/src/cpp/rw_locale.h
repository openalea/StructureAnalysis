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

/*! \file rw_locale.h
    \brief File for Rogue Wave locale with the STL.
*/

#ifndef __rw_locale_h__
#define __rw_locale_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */


#include <sstream>
#include <iostream>
//using namespace std;
#include <locale.h>

#include "rw_cstring.h"
#include "tools_namespace.h"
#include "util_string.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class Locale.

  \brief Handling cultural difference.

*/


class Locale
{
 public:

  enum dateOrder { DMY, MDY, YDM, YMD };
  dateOrder dateOrder_;

 public:

  Locale() : name() {}
  Locale( const std::string& localName ) : name(localName) {}
  virtual ~Locale() {}

  bool stringToNum( const RWString& s, double* fp ) const
    {
    if( fp == 0 ) return false;
        std::istringstream flow(s);
//    flow.imbue(std::locale(name));
    flow >> *fp;
    return true;
    }

  bool stringToNum( const RWString& s, long* ip ) const
    {
    if( ip == 0 ) return false;
        std::istringstream flow(s);
//    flow.imbue(std::locale(name));
    flow >> *ip;
    return true;
    }

  RWCString asString(long i) const { return RWCString(number(i));}
  RWCString asString(unsigned long ui) const { return RWCString(number(ui));}
  RWCString asString(double f ) const { return RWCString(number(f));}

  static const Locale* global( const Locale* loc ) { return loc; }

 private:
         std::string name;
}; // class Locale


/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

typedef VPTOOLS(Locale) RWLocale;
typedef VPTOOLS(Locale) RWLocaleSnapshot;

/* ----------------------------------------------------------------------- */

#else

#include <rw/locale.h>

#endif

// __rw_locale_h
#endif
