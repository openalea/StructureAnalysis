/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): C. Nouguier & F. Boudon (frederic.boudon@cirad.fr) nouguier
 *
 *       $Source$
 *       $Id: all_tools.h 9869 2010-11-04 17:31:41Z dbarbeau $
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

/*  --------------------------------------------------------------------- */

/*! \file all_tools.h
    \brief File that include all the tools utility
*/

/// Configuration
#include "config.h"
#include "tools_namespace.h"

/// Utilities
#include "util_assert.h"
#include "util_string.h"

#include "util_types.h"

#include "util_math.h"

/// other tools
#include "dirnames.h"
#include "errormsg.h"
#include "rcobject.h"
#include "timer.h"

/// Parser files
/* 
Must be include with macro
#include "gparser.h"
#include "gscanner.h"
#include "gsmbtable.h"
*/
#include "readline.h"

/// STL
#include "std.h"

