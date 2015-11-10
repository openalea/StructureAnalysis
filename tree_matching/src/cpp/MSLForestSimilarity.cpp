#include "MSLSimilarity.h"
#include <iostream>
#include <strstream>
#include <cstdio>
#include <cmath>

//-----------------------------------------------------------------------------
// Classe permettant le calcul d'une valeur minimale parmi plusieurs 
//-----------------------------------------------------------------------------

class Max {
  float value;
public:
  Max(float v):value(v) {}
  Max &operator<<(float v) {
    if (v>value) value=v;
    return *this;
  }
  operator float() {
    return value;
  }
};


//-----------------------------------------------------------------------------
//Calcul de la distance entre deux forêts indexées par leurs racines en
//Remplissant récursivement le tableau fDist
//-----------------------------------------------------------------------------
void MSLSimilarity::computeForestsDistances(int input_vertex,int reference_vertex){

  int i,j,i1,j1;
  i= input_vertex;
  j= reference_vertex;   
  float max;
  //Calcul récursif de la distance entre les forêts T1[L1[i]..i] et T2[L2[j]..j]
  int fNum = fix(i,j);
  fDist[fNum][0][0] = 0;

  if(i>0){
    for (i1 = L1[i]; i1 <= i; i1++) {
      fDist[fNum][idx(L1[i],i1)][0] = 0;//fDist[fNum][idx(L1[i],i1-1)][0]    + delLocalNode(T1,i1);
    }
  }
  
  
  if(j>0){
    for (j1 = L2[j]; j1 <= j; j1++) {
      fDist[fNum][0][idx(L2[j],j1)] =  0;//fDist[fNum][0][idx(L2[j],j1-1)] + insLocalNode(T2,j1);
    }
 }

  if(i>0 || j>0){  
    for (i1 = L1[i]; i1 <= i; i1++) {
      for (j1 = L2[j]; j1 <= j; j1++) {
	if(L1[i1] == L1[i] && L2[j1] == L2[j]){
	  //foret[L1[i]...i1] et foret[L2[j]...j1] sont des arbres
	  
	  float dele = fDist[fNum][idx(L1[i],i1-1)][idx(L2[j],j1)] +  delLocalNode(T1,i1);
	  float inse = fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1-1)] +  insLocalNode(T2,j1);
	  float match = fDist[fNum][idx(L1[i],i1-1)][idx(L2[j],j1-1)] + matchLocalNode(i1,j1);
	  max = Max(dele) << inse << match << 0.;
	  
	  fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1)] = max;
	  _dst1[i1][j1] = max;
	  if(max>lsmax){
	    lsmax=max;
	    vmax=i1;
	    lvmax=L1[i];
	    wmax=j1;
	    lwmax=L2[j];
	  }
	}
	else{
	  
	  float dele = fDist[fNum][idx(L1[i],i1-1)][idx(L2[j],j1)] + delLocalNode(T1,i1);
	  float inse = fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1-1)] +  insLocalNode (T2,j1);
	  float match = fDist[fNum][idx(L1[i],L1[i1]-1)][idx(L2[j],L2[j1]-1)] + _dst1[i1][j1];
	  max = Max(dele) << inse << match << 0.;
	  
	  
	  fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1)] = max;
	 
	}
	/*if(max>lsmax){
	  lsmax=max;
	  vmax=i1;
	  lvmax=L1[i];
	  wmax=j1;
	  lwmax=L2[j];
	  }*/    
      }
    }
  }
}


//-----------------------------------------------------------------------------
//Calcul des indices de remplissage de la matrice fDist
//
//-----------------------------------------------------------------------------

int MSLSimilarity::fix(int a, int b){
  return a*sizeKeyroots2 +b;
}


int MSLSimilarity::idx(int a, int b){
  return (a>b)?0:(b-a+1);
}

