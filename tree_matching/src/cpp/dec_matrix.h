/* -*-c++-*- 
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture 
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr) 
 *
 *       $Source$
 *       $Id$
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

#ifndef __dec_matrix_h__
#define __dec_matrix_h__

#ifdef _WIN32

#define NEW_MAT(type,name,dim1,dim2) type ** name = new type *[dim2]; \
	int _i##name = 0; \
	for(_i##name = 0 ; _i##name < dim1 ; _i##name ++) name[_i##name] = new type[dim2]; \

#define DEL_MAT(name,dim1)  \
	int _j##name = 0; \
	for(_j##name = 0 ; _j##name < dim1 ; _j##name ++) delete name[_j##name]; \
	delete name; \

#else

#define NEW_MAT(type,name,dim1,dim2) type name[dim1][dim2];
#define DEL_MAT(name,dim1);
	
#endif

#endif

