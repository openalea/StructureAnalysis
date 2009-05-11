/* -*-c++-*- 
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *
 *       $Source$
 *       $Id: definitions.h 3264 2007-06-06 14:22:22Z dufourko $
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


#ifndef SB_DEFINITIONS_HEADER
#define SB_DEFINITIONS_HEADER

#include "treematch_config.h"
#include <iostream>
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <fstream>
#include <list>
#include <vector>

typedef double 		DistanceType;
typedef int 		IntType;
typedef int 		ItemType;
typedef DistanceType 	KeyType;
typedef char 		Symbolic;
typedef double 		Numeric;

const DistanceType MAXDIST=1e+20;
const DistanceType MINDIST=0.0;

#define D_MAX(A,B) ((A)>(B)? (A):(B))
#define D_MIN(A,B) ((A)<(B)? (A):(B))
#define I_MAX(A,B) ((A)>(B)? (A):(B))
#define I_MIN(A,B) ((A)<(B)? (A):(B))


const DistanceType DIST_UNDEF=-1.0;
const int EMPTY_TREE=-1;
const int EMPTY_NODE=-1;


typedef enum { TM_RELATION=414,TM_SEQUENCE,TM_TREENODE,TM_TREEGRAPH,TM_TREEMATCH } TM_ClassIdent;


/**\par Tarjan's functions */

/** Return the smallest integer not less than val. */
TREEMATCH_API int SINLT(float );

/** Return the largest integer not greater than val */
TREEMATCH_API int LINGT(float);

/** Absolute value */
TREEMATCH_API DistanceType ABS(DistanceType );

TREEMATCH_API Symbolic readSymbol(const char* );

TREEMATCH_API Numeric readNumber(const char* );

#endif

