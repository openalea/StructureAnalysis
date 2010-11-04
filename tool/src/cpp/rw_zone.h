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
/*! \file rw_zone.h
    \brief File for Rogue Wave zone with the STL.
*/

#ifndef __rw_zone_h__
#define __rw_zone_h__

#include "config.h"

#ifdef RWOUT
/* ----------------------------------------------------------------------- */


#include <locale.h>
//using namespace std;

#include "tools_namespace.h"

/* ----------------------------------------------------------------------- */

VPTOOLS_BEGIN_NAMESPACE

/* ----------------------------------------------------------------------- */


/*!

  \class RWZone

  \brief RW zone wrapper.

*/


struct Zone
{
  static void locale( const Zone* z )
    { setlocale(LC_TIME,""); }

  enum DstRule{ NoDST, NoAm, WeEu };
  enum StdZone{ Europe, USEastern, USPacific };
};

/// Simple Zone
struct ZoneSimple : public Zone
{
  ZoneSimple( Zone::StdZone zone, Zone::DstRule ) {}
};

/* ----------------------------------------------------------------------- */

VPTOOLS_END_NAMESPACE

#define RWZone VPTOOLS(Zone)
#define RWZoneSimple VPTOOLS(ZoneSimple)

/* ----------------------------------------------------------------------- */

#else

#include <rw/zone.h>

#endif

// __rw_stack_h
#endif
