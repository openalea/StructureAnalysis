 /*------------------------------------------------------------------------------
 *                                                                              
 *        VPlants.Stat_Tool : VPlants Statistics module
 *                                                                              
 *        Copyright 2006-2007 INRIA - CIRAD - INRA                      
 *                                                                              
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>         
 *                                                                              
 *        Distributed under the GPL 2.0 License.                               
 *        See accompanying file LICENSE.txt or copy at                          
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *                                                                              
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr                    
 *       
 *        $Id: export_base.h 8081 2010-02-18 09:32:17Z guedon $
 *                                                                       
 *-----------------------------------------------------------------------------*/

#ifndef __CLASS_STAT_TOOL_WBASE__
#define __CLASS_STAT_TOOL_WBASE__

// Boost.Python Wrapper export function
void class_constant();
void class_stat_error();


class StatInterfaceWrap;

void class_stat_interface();
void class_forward();


#endif
