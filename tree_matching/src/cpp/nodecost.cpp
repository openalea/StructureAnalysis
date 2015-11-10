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


#include"nodecost.h"

NodeCost::NodeCost( NodeCostType type)
{
 
  _type=type;
  _norm = L1;
  
  
}

NodeCost::NodeCost(char* file_name)
{
  _type=MATRIX;
  _norm = L1;
  char tmp [1000];
  char chaine [1000];
  char chaine1 [20];
  int cost;
  matrix_size = 0;
  FILE* input = fopen(file_name, "r");
  fgets(tmp, 1000, input);
  
  int k = 0;
  int i = 1;
  int j =0;
  while(tmp[i] != '\n'){
    if(tmp[i] != '\t'){
      chaine[j] = tmp[i];
      j++;
	}
    else{
      chaine[j] = '\n';
      sscanf(chaine,"%s\n%*s",&chaine1);
      //printf("%s\n",chaine1);
      symbols[chaine1] = k;
      k++;
      j =0;
    }
    i++;
  }

  int m = 0;
  int n =-1;
  while(fgets(tmp, 1000, input)){
    i = 0;
    j =0;
    n = -1;  
    while(tmp[i] != '\n' && tmp[i] != EOF){
      if(tmp[i] != '\t'){
	chaine[j] = tmp[i];
	j++;
      }
      else{
	chaine[j] = '\n';
	sscanf(chaine,"%s\n%*s",&chaine1);
	//printf("%s\n",chaine1);
	matrix[m][n] = atoi(chaine1);
	n++;
	j =0;
      }
      i++;
    }
    m++;
  }
  fclose(input);
}

NodeCost::NodeCost( NodeCostType type,Norm norm)
{
  _type=type;
  _norm = norm;

}

int NodeCost::getCost(const char* s1, const char* s2){
  //return matrix[symbols[s1]][symbols[s2]];
  if(s1=="-" || s2 == "-")
    return 0;
  else
    return 1;
}

DistanceType NodeCost::getInsertionCost(TreeNode* node)
{
  DistanceType cost;

  if (_type!= SCORE)
   {
     cost=node->getValue();
   }

  if (_type == SCORE)
   {
     cost=-node->getValue();
   }

  return(cost);
}

DistanceType NodeCost::getDeletionCost(TreeNode* node)
{
  DistanceType cost;
   if (_type!= SCORE)
   {
     cost=node->getValue();
   }
     if (_type== SCORE)
   {
     cost=-node->getValue();
   }
  return(cost);
}

DistanceType NodeCost::getChangingCost(TreeNode* i_node,TreeNode* r_node)
{
  DistanceType cost;
   if (_type!= SCORE)
   {
     cost=ABS(i_node->getValue()-r_node->getValue());
   }
   if (_type== SCORE)
   {
     cost=4-ABS(i_node->getValue()-r_node->getValue());
   }
  return(cost);
}


