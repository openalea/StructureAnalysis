/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture 
 *
 *       Copyright 1995-2010 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): C. Nouguier & F. Boudon (frederic.boudon@cirad.fr) nouguier 
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




#ifndef __util_math_h__
#define __util_math_h__

/*! \file util_math.h
    \brief File that contains some math utility.
*/

#include <cmath>     // Common math function.
#include <algorithm> // For min, max ...

/* ----------------------------------------------------------------------- */
#ifndef M_PI
/// value of pi 
#  define M_PI 3.14159265358979323846
#endif


/// Returns the cube value of \e v. 
inline float cb( const float& v ) {
  return v * v * v;
}

/// Returns the cube value of \e v. 
inline double cb( const double& v ) {
  return v * v * v;
}

/// Computes the square value of \e v.
inline float sq( const float& v ) {
  return v * v;
}

/// Computes the square value of \e v.
inline double sq( const double& v ) {
  return v * v;
}

inline bool isPowerOfTwo(int val){
   long i=0;
   while ((1<<i) < val) i++;
   if (val==(1<<i)) return true;
   return false;     
}

#ifdef _MSC_VER

/// On win32, redirect finite on _finite.
#define finite _finite

/// On win32, redefine cubic root.
inline double cbrt(double x){
	return pow(x,1/3);
}

/// On win32, redefine round double to int.
inline int rint(double x){
	int res = (int)x;
	if(x-res>0.5)res++;
	return res;
}

/// On win32, redefine round double to int.
#define round rint

/// On win32, redefine truncate double to int.
inline int trunc(double x){
	int res = (int)x;
	if(x-res>0.5)res++;
	return res;
}

#endif

/* ----------------------------------------------------------------------- */

//__util_math_h__
#endif 
