/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching
 *
 *       File author(s): F. Boudon (frederic.boudon@cirad.fr)
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
#ifndef __treematch_config_h__
#define __treematch_config_h__


/* ----------------------------------------------------------------------- */


/*! \def TREEMATCH_NODLL
    \brief Not creating dll

    Uncomment to use this functionnality
        Do nothing on other platform than windows
*/
/*! \def TREEMATCH_DLL
    \brief Using lib TreeMatching as a dll

    Uncomment to use this functionnality
        Do nothing on other platform than windows
*/
/*! \def TREEMATCH_MAKEDLL
    \brief Creating TreeMatching dll

    Uncomment to use this functionnality
        Do nothing on other platform than windows
*/


#if defined(_WIN32)
#if defined(TREEMATCH_NODLL)
#undef TREEMATCH_MAKEDLL
#undef TREEMATCH_DLL
#else
#ifndef TREEMATCH_DLL
#define TREEMATCH_DLL
#endif
#endif

#if defined(TREEMATCH_MAKEDLL)
#ifndef TREEMATCH_DLL
#define TREEMATCH_DLL
#endif
#endif

#ifdef TREEMATCH_DLL

#ifdef TREEMATCH_MAKEDLL             /* create a TREEMATCH DLL library */
#define TREEMATCH_API  __declspec(dllexport)
#else                                                   /* use a TREEMATCH DLL library */
#define TREEMATCH_API  __declspec(dllimport)
#endif

#endif

#else // OS != _WIN32

#undef TREEMATCH_MAKEDLL             /* ignore these for other platforms */
#undef TREEMATCH_DLL

#endif

#ifndef TREEMATCH_API
#define TREEMATCH_API
#endif


#if defined( _MSC_VER )
 // Make deprecated warning more important
 #   pragma warning ( 1 : 4996)
#endif

/// deprecated attribute definition
#ifdef __GNUC__
#define attribute_deprecated __attribute__((deprecated))
#elif defined( _MSC_VER )
#define attribute_deprecated __declspec(deprecated)
#else
#define attribute_deprecated
#endif

/* ----------------------------------------------------------------------- */

// __config_h__
#endif
