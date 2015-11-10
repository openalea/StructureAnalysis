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


#include "rw_time.h"
#include "rw_date.h"

#ifdef RWOUT

VPTOOLS_BEGIN_NAMESPACE

Time::Time( unsigned hour, unsigned minute, unsigned second )
{
    struct tm* timeinfo;
    time(&s);
    timeinfo= localtime(&s);
    timeinfo->tm_sec= second;
    timeinfo->tm_min= minute;
    timeinfo->tm_hour= hour;
    s= mktime(timeinfo);
}

Time::Time( const CDate& d, unsigned hour,
            unsigned minute, unsigned second )
{
  struct tm timeinfo;
  d.extract(&timeinfo);
  timeinfo.tm_sec= second;
  timeinfo.tm_min= minute;
  timeinfo.tm_hour= hour;
  s= mktime(&timeinfo);
}

Time::Time( unsigned long time_s )
{
    if( time_s < Epoch )
		s= (long)-1;
    else
		s= time_t(time_s-Epoch);
}

void Time::extract( struct tm* t ) const
{
    struct tm* timeinfo= localtime(&s);
	if(timeinfo!=NULL)*t = *timeinfo;
}

VPTOOLS_END_NAMESPACE

#endif
