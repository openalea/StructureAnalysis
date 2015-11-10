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
/*! \file rw_tokenizer.h
    \brief File for Rogue Wave tokenizer with the STL.
*/

#ifndef __rw_tokenizer_h__
#define __rw_tokenizer_h__

#include "config.h"

#ifdef RWOUT

/* ----------------------------------------------------------------------- */


#include <string>
//using namespace std;

#include "tools_namespace.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class Tokenizer.

  \brief A tokenizer to lex a string.

*/


class Tokenizer
{

 public:

  /*!
    Construct a Tokeniser to lex the string s.
  */
	 Tokenizer( const std::string& s ) :
	  phrase(s), in(0), next(0)
	 {}

  virtual ~Tokenizer() {}

  /*!
    Advance to the next token and return it as a string. The tokens are
    considered to be delimated by any of the characters in s.
  */
  std::string operator() (const char* delimiters= " \t\n")
    {
    if(phrase.empty())
      return phrase;

    if( in == 0 )
      in= phrase.find_first_not_of(delimiters, next);

    if( in == std::string::npos )
      return std::string();

    next = phrase.find_first_of(delimiters, in);
    if( next == std::string::npos )
      next= phrase.length();

    std::string word(phrase, in, next - in);

    in= phrase.find_first_not_of(delimiters, next);
    return word;
    }


 private:
   std::string phrase;
   size_t in, next;

}; // class Tokenizer


/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

typedef VPTOOLS(Tokenizer) RWCTokenizer;

/* ----------------------------------------------------------------------- */

#else

#include <rw/ctoken.h>

#endif // RWOUT
// Tokenizer_h
#endif
