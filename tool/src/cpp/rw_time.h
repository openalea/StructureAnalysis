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
/*! \file rw_time.h
    \brief File for Rogue Wave time with the STL.
*/

#ifndef __rw_time_h__
#define __rw_time_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */


#include <ctime>
//using namespace std;

#include "tools_namespace.h"
#include "timer.h"


/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class RWTime

  \brief RW time wrapper.

*/


class CDate;

const unsigned long Epoch= 2177452800UL;

class VPTOOLS_API Time
{

// seconds 01/01/1901->01/01/1970

public:
  Time() { time(&s); }
  Time( const Time& t ) { s= t.s; }

  Time( unsigned hour, unsigned minute, unsigned second= 0 );

  Time( const CDate& d, unsigned hour= 0,
        unsigned minute= 0, unsigned second= 0 );

  Time( unsigned long time_s );

  virtual ~Time() {}

  void extract( struct tm* t ) const;

  unsigned long seconds() const
    { return (s+Epoch); }

  bool operator!=( const Time& t2 ) const { return s != t2.s; }
  bool operator==( const Time& t2 ) const { return s == t2.s; }
  bool operator<( const Time& t2 ) const { return s < t2.s; }

  bool isValid() const { return (s == (long)-1) ? false : true; }

public:
  time_t s;

};


/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#define RWTime VPTOOLS(Time)
#define RWTimer VPTOOLS(Timer)

/* ----------------------------------------------------------------------- */

#else

#include <rw/rwtime.h>

#endif

// __rw_time_h
#endif
