/* -*-c++-*- 
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture 
 *
 *       Copyright 1995-2010 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Ch. Godin (christophe.godin@cirad.fr) 
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


/*! \file util_enviro.h
    \brief Function relative to the environment.
*/

#ifndef __util_vp_enviro_h__
#define __util_vp_enviro_h__

#include <string>
#include "tools_namespace.h"

VPTOOLS_BEGIN_NAMESPACE

/// Get the Home directory.
VPTOOLS_API std::string getHome();

/// Get the Current Working directory.
VPTOOLS_API std::string getCwd();

/// Get AMAPmod directory.
VPTOOLS_API std::string getAMAPmodDir();

/// Set AMAPmod directory.
VPTOOLS_API void setAMAPmodDir(const std::string&);

/// Get User name
VPTOOLS_API std::string getUserName();

/// Get Family Name of the Operating System [ex:Windows,Linux].
VPTOOLS_API std::string getOSFamily();

/// Get Name of the Operating System [ex:Windows NT].
VPTOOLS_API std::string getOSName();

/// Get Release of the Operating System [ex:5.1].
VPTOOLS_API std::string getOSRelease();

/// Get Version of the Operating System [ex:2840].
VPTOOLS_API std::string getOSVersion();

/// Get the machine name.
VPTOOLS_API std::string getMachineName();

/// Get the architecture of the machine.
VPTOOLS_API std::string getArchitecture();

/// Get the compiler name.
VPTOOLS_API std::string getCompilerName();

/// Get the compiler version.
VPTOOLS_API std::string getCompilerVersion();

/// Get the language.
VPTOOLS_API std::string getOSLanguage();
VPTOOLS_API std::string getLanguage();
VPTOOLS_API void setLanguage(const std::string&);

VPTOOLS_END_NAMESPACE

#endif
