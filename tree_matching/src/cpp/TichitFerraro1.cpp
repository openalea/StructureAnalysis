//-*-c++-*-

#include "TichitFerraro1.h"
#include <iostream>
#include <strstream>
#include <cstdio>
#include <cmath>

using namespace std;

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

TichitFerraro1::TichitFerraro1(MTG& mtg,TreeGraph& input, TreeGraph & reference, int scale, int flag)
{
  _mtg = &mtg;
  superGraph = new SuperGraph(input,reference);
  sequence = new Sequence();
  _scale = scale; // représente l'echelle de comparaison : 0 = comparaison classique, 1 = cout indel utilise les complexes.
  _flag = flag; // pour l'affichage du temps de calcul
  this->init();
}

TichitFerraro1::~TichitFerraro1()
{}

const float TichitFerraro1::epsilon = 0.01;
const float TichitFerraro1::match_color = 0.0;
const float TichitFerraro1::diff_color = 1024.0;
const float TichitFerraro1::infiniteFloat = 1024.0;
const int TichitFerraro1::infiniteInt = 1024;

float TichitFerraro1::getEdgeValue(const edge& e)
{
  return 0.0;
}

float TichitFerraro1::getNodeValue(const node& n)
{
  return 0.0;
}

void TichitFerraro1::reset() {

}

static int totoro = 0;

/*
static float v[2][2] = {{0.5,-1},
                         {-1,0.5}};
static float h[2][2] = {{0.5,-1},
                         {-1,0.5}};
*/

float TichitFerraro1::v[2][2] = {{0.0,0.0},
                   {0.0,0.0}};


float TichitFerraro1::h[2][2] = {{0.0,0.0},
                   {0.0,0.0}};

float TichitFerraro1:: CostMatrix[NB_SYMBOLS][NB_SYMBOLS]
 = {
  { 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {-1,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  2},
};

Sequence* TichitFerraro1::getSequence() {
  return(sequence);
}


void TichitFerraro1::init() {
  bool cached, resultBool;
  string erreurMsg;
  sizeT1 = superGraph->getNbVertex(1);
  tree1.resize(sizeT1+1);
  tDist.resize(sizeT1+1);
  sizeT2 = superGraph->getNbVertex(2);
  tree2.resize(sizeT2+1);
  for(int i = 0; i<=sizeT1; i++) {
    tDist[i].resize(sizeT2+1);
  }

//   symbols["-"] = 0;
//   symbols["AU"] = 1;
//   symbols["UA"] = 2;
//   symbols["CG"] = 3;
//   symbols["GC"] = 4;
//   symbols["GU"] = 5;
//   symbols["UG"] = 6;
//   symbols["P"] = 7;
//   symbols["A"] = 8;
//   symbols["U"] = 9;
//   symbols["G"] = 10;
//   symbols["C"] = 11;
//   symbols["L"] = 12;
//   symbols["M"] = 13;
//   symbols["I"] = 14;
//   symbols["B5"] = 15;
//   symbols["B3"] = 16;
//   symbols["H"] = 17;
//   symbols["S"] = 18;
//   symbols["default"] = 19;
//   symbols["AG"] = 20;
//   symbols["GA"] = 21;
//   symbols["AC"] = 22;
//   symbols["CA"] = 23;
}


double TichitFerraro1::run()  {
  tree1[0] = 0;
  tree2[0] = 0;
  root1 = 0;
  root2 = root1 + sizeT1;
  storeNodes(root1,tree1,1);
  storeNodes(root2,tree2,1);
  preCompute();
  computeTreeDist();
  align(vmax,vmax,wmax,wmax);
  return tDist[sizeT1][sizeT2];
}

// int TichitFerraro1::storeNodes(node nd, vector<node>& tree, int i) {
//   Iterator<node> *itn = superGraph->getOutNodes(nd);
//   int j = i;
//   node n;
//   for(;itn->hasNext();) {
//     n = itn->next();
//     j = storeNodes(n, tree, j);
//   }
//   //cerr << "j = " << j << " nd = " << nd.id << " Outdeg = " << superGraph->outdeg(nd) << endl;
//   tree[j] = nd;
//   delete itn;
//   return j+1;
// }

int TichitFerraro1::storeNodes(node nd, vector<node>& tree, int i) {
  //list<node>::const_iterator *itn = superGraph->getOutNodes(nd);
  list<node> itn  = superGraph->getOutNodes(nd);
  list<node>::iterator begin,end;
  begin = itn.begin();
  end = itn.end();
  int j= i;
  node n;
  while (begin != end) {
    n = *(begin);
    j = storeNodes(n, tree, j);
    begin++;
  }
  //cerr << "j = " << j << " nd = " << nd<< " Outdeg = " << superGraph->outdeg(nd) << endl;
   tree[j] = nd;
   //   delete itn;
   return j+1;
 }

void TichitFerraro1::preCompute()  {
  L1.resize(sizeT1+1);
  keyroots1.resize(sizeT1);
  int i = 0;
  for(i = 1; i <= sizeT1; ++i) {
    if (superGraph->outdeg(tree1[i])==0) {
      //cerr << "L1[" << i << "] =" << i << endl;
      L1[i]=i;
    }
    else {
      list<node> itn  = superGraph->getOutNodes(tree1[i]);
      list<node>::iterator begin,end;
      begin = itn.begin();
      end = itn.end();
      node n = *(begin);
      //node n = (superGraph->getOutNodes(tree1[i]))->next();
      //cerr << n.id << endl;
      int j;
      for (j = 1; j<i;++j) {
        //cerr << "j = " << j << endl;
        //cerr << tree1[L1[j]].id << endl;
        if(tree1[j] == n) {
          //cerr << "L1[" << i << "] = " << "L1[" << j << "]" << endl;
          break;
        }
        else {
          //cerr << "crotte 1" << endl;
        }
      }
      L1[i]=L1[j];
      //L1[i]=L1[dfs1.comp_num(*(tree1[i].adj_nodes_begin()))];
    }
  }
  int idx = 0;
  for(i = 1; i <= sizeT1 - 1; ++i) {
    int minI = L1[i];
    int j = i+1;
    for(; j <= sizeT1; ++j) {
      if (L1[j] == minI)
        break;
    }
    if (j == sizeT1+1) {
      keyroots1[idx]=i;
      ++idx;
    }
  }
  keyroots1[idx] = sizeT1;
  sizeKeyroots1 = idx+1;
  keyroots1.resize(sizeKeyroots1);
  L2.resize(sizeT2+1);
  keyroots2.resize(sizeT2);
  for( i = 1; i <= sizeT2; ++i) {
    if (superGraph->outdeg(tree2[i])==0) {
      //cerr << "L2[" << i << "] =" << i << endl;
      L2[i]=i;
    }
    else {
      list<node> itn  = superGraph->getOutNodes(tree2[i]);
      list<node>::iterator begin,end;
      begin = itn.begin();
      end = itn.end();
      node n = *(begin);
      //      node n = (superGraph->getOutNodes(tree2[i]))->next();
      int j;
      for (j = 1; j<i;++j) {
        if(tree2[j] == n) {
          //cerr << "L2[" << i << "] = " << "L2[" << j << "]" << endl;
          break;
        }
        else {
          //cerr << "crotte 2" << endl;
        }
      }
      L2[i]=L2[j];
      //L2[i]=L2[dfs2.comp_num(*(tree2[i].adj_nodes_begin()))];
    }
  }
  idx = 0;
  for( i = 1; i <= sizeT2 - 1; ++i) {
    int minI = L2[i];
    int j = i+1;
    for(; j <= sizeT2; ++j) {
      if (L2[j] == minI)
        break;
    }
    if (j == sizeT2+1) {
      keyroots2[idx]=i;
      ++idx;
    }
  }
  keyroots2[idx] = sizeT2;
  sizeKeyroots2 = idx+1;
  keyroots2.resize(sizeKeyroots2);
  int sizeFD = sizeKeyroots1*sizeKeyroots2;

  fDist.resize(sizeFD);
  for( i = 0; i < sizeKeyroots1; ++i) {
    for (int j = 0; j < sizeKeyroots2; ++j) {
      int leni = keyroots1[i]-L1[keyroots1[i]]+1;
      fDist[i*sizeKeyroots2+j].resize(leni+1);
      for (int k = 0; k <= leni; ++k) {
        int lenj = keyroots2[j]-L2[keyroots2[j]]+1;
        fDist[i*sizeKeyroots2+j][k].resize(lenj+1);
      }
    }
  }

  rootF1.resize(sizeFD);
  for( i = 0; i < sizeKeyroots1; ++i) {
    for (int j = 0; j < sizeKeyroots2; ++j) {
      int leni = keyroots1[i]-L1[keyroots1[i]]+1;
      rootF1[i*sizeKeyroots2+j].resize(leni+1);
      for (int k = 0; k <= leni; ++k) {
        int lenj = keyroots2[j]-L2[keyroots2[j]]+1;
        rootF1[i*sizeKeyroots2+j][k].resize(lenj+1);
      }
    }
  }

  rootF2.resize(sizeFD);
  for( i = 0; i < sizeKeyroots1; ++i) {
    for (int j = 0; j < sizeKeyroots2; ++j) {
      int leni = keyroots1[i]-L1[keyroots1[i]]+1;
      rootF2[i*sizeKeyroots2+j].resize(leni+1);
      for (int k = 0; k <= leni; ++k) {
        int lenj = keyroots2[j]-L2[keyroots2[j]]+1;
        rootF2[i*sizeKeyroots2+j][k].resize(lenj+1);
      }
    }
  }
  rootT1.resize(sizeT1+1);
  for( i=0;i<=sizeT1;++i) {
    rootT1[i].resize(sizeT2+1);
  }
  rootT2.resize(sizeT1+1);
  for( i=0;i<=sizeT1;++i) {
    rootT2[i].resize(sizeT2+1);
  }
}

void TichitFerraro1::computeTreeDist() {
  lsmax = 0;
  for (int i=0; i<sizeKeyroots1; ++i) {
    if (_flag && (int(100.*i/sizeKeyroots1)%5 == 0))
      cerr << "\x0d" << "Already computed : "<<int(100.*i/sizeKeyroots1)<<" %   ...  " << flush;
    for (int j=0; j<sizeKeyroots2; ++j) {
      computeTreeDist(i,j);
    }
  }
  cerr << "\x0d" <<flush;
}


void TichitFerraro1::computeTreeDist(int ki, int kj)  {
  int i1;
  int j1;
  int i = keyroots1[ki];
  int j = keyroots2[kj];
  int fNum = fix(ki,kj);

  fDist[fNum][0][0] = 0;
  rootF1[fNum][0][0] = "";
  rootF2[fNum][0][0] = "";
  for (i1 = L1[i]; i1 <= i; i1++) {
    fDist[fNum][idx(L1[i], i1)][0] = 0;//fDist[fNum][idx(L1[i], i1-1)][0] +  del(tree1[i1]);
    rootF1[fNum][idx(L1[i], i1)][0] = setRootF1(fNum,i,i1,j,0);
    rootF2[fNum][idx(L1[i], i1)][0] = rootF2[fNum][idx(L1[i], i1-1)][0];
  }

  for (j1 = L2[j]; j1 <= j; j1++) {
    fDist[fNum][0][idx(L2[j], j1)] = 0;//fDist[fNum][0][idx(L2[j], j1-1)] +  ins(tree2[j1]);
    rootF1[fNum][0][idx(L2[j], j1)] = rootF1[fNum][0][idx(L2[j], j1-1)];
    rootF2[fNum][0][idx(L2[j], j1)] = setRootF2(fNum,i,0,j,j1);
  }

  for (i1=L1[i]; i1 <=i; i1++) {
    for (j1=L2[j]; j1<=j; j1++) {
      if (L1[i1] == L1[i] && L2[j1] == L2[j]) {
        //cerr << "delT :" << delTreeHoriz(fNum,i1,j1) + delTreeVert(fNum,i1,j1) <<endl;
        //cerr << "insT :" << insTreeHoriz(fNum,i1,j1) + insTreeVert(fNum,i1,j1)<<endl;
        //cerr << "matchT:" << matchTreeHoriz(fNum,i1,j1) +matchTreeVert(fNum,i1,j1)<<endl;
        float dele = fDist[fNum][idx(L1[i],i1-1)][idx(L2[j],j1)] +  del(tree1[i1]) + delTreeHoriz(fNum,i1,j1) + delTreeVert(fNum,i1,j1);
        float inse = fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1-1)] +  ins(tree2[j1]) + insTreeHoriz(fNum,i1,j1) + insTreeVert(fNum,i1,j1);;
        float match = fDist[fNum][idx(L1[i],i1-1)][idx(L2[j],j1-1)] + subst(tree1[i1],tree2[j1]) + matchTreeHoriz(fNum,i1,j1) +matchTreeVert(fNum,i1,j1);
        float maxi = Max(dele) << inse << match << 0.;

	if(maxi>lsmax){
	  lsmax = maxi;
	  vmax=i1;
	  wmax=j1;
	}
	
        fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1)] = maxi;
        if (maxi == dele) {
          rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = "0";
          rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = rootF2[fNum][idx(L1[i], i1-1)][idx(L2[j], j1)];
        }
        else if (maxi == inse) {
          rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = rootF1[fNum][idx(L1[i], i1)][idx(L2[j], j1-1)];
          rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = "0";
        }
        else  if (maxi == match) {
          rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = "1";
          rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = "1";
        }
	else{
	  rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = "0";
          rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = "0";
	}
        rootT1[i1][j1] = rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)];
        rootT2[i1][j1] = rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)];
        tDist[i1][j1] = fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1)];
      }
      else {
        //cerr << "del :" << delForestHoriz(fNum,i,i1,j,j1) + delForestVert(fNum,i,i1,j,j1) <<endl;
        //cerr << "ins :" << insForestHoriz(fNum,i,i1,j,j1) + insForestVert(fNum,i,i1,j,j1) <<endl;
        float dele = fDist[fNum][idx(L1[i],i1-1)][idx(L2[j],j1)] +  del(tree1[i1]) + delForestHoriz(fNum,i,i1,j,j1) + delForestVert(fNum,i,i1,j,j1);
        float inse = fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1-1)] +  ins(tree2[j1]) + insForestHoriz(fNum,i,i1,j,j1) + insForestVert(fNum,i,i1,j,j1);
        float match = fDist[fNum][idx(L1[i],L1[i1]-1)][idx(L2[j],L2[j1]-1)] + tDist[i1][j1];
        float maxi = Max(dele) << inse << match << 0.;
        fDist[fNum][idx(L1[i],i1)][idx(L2[j],j1)] = maxi;
        if (maxi == dele) {
          rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = setRootF1(fNum,i,i1,j,j1);
          rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = rootF2[fNum][idx(L1[i], i1-1)][idx(L2[j], j1)];
        }
        else if (maxi == inse) {
          rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = rootF1[fNum][idx(L1[i], i1)][idx(L2[j], j1-1)];
          rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = setRootF2(fNum,i,i1,j,j1);
        }
        else if (maxi == match) {
          rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = rootF1[fNum][idx(L1[i],L1[i1]-1)][idx(L2[j],L2[j1]-1)] + rootT1[i1][j1];
          rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = rootF2[fNum][idx(L1[i],L1[i1]-1)][idx(L2[j],L2[j1]-1)] + rootT2[i1][j1];
        }
	else{
	  rootF1[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = "0";
          rootF2[fNum][idx(L1[i],i1)][idx(L2[j], j1)] = "0";
	}
      }
    }
  }
}
/* exemple pour tree
 */
float TichitFerraro1::delTreeHoriz(int fNum, int i, int j) {
  float res = 0.0;

  string s = rootF1[fNum][idx(L1[i], 1-1)][idx(L2[j], j)];

  for (int a = 0; a <=1; a++) {
    for (int b = 0; b <=1; b++) {
      strstream tmp;
      string sab;
      tmp << a << b;
      tmp >> sab;
      res += h[a][b] * nbPatterns(sab,s);
    }
  }
  return res;
}

float TichitFerraro1::insTreeHoriz(int fNum, int i, int j) {
  float res = 0.0;

  string s = rootF2[fNum][idx(L1[i], 1)][idx(L2[j], j-1)];

  for (int a = 0; a <=1; a++) {
    for (int b = 0; b <=1; b++) {
      strstream tmp;
      string sab;
      tmp << a << b;
      tmp >> sab;
      res += h[a][b] * nbPatterns(sab,s);
    }
  }
  return res;
}
float TichitFerraro1::matchTreeHoriz(int fNum, int i, int j) {
  float res1 = 0.0;
  float res2 = 0.0;

  string s1 = rootF1[fNum][idx(L1[i], 1-1)][idx(L2[j], j-1)];
  string s2 = rootF2[fNum][idx(L1[i], 1-1)][idx(L2[j], j-1)];

  for (int a = 0; a <=1; a++) {
    for (int b = 0; b <=1; b++) {
      strstream tmp;
      string sab;
      tmp << a << b;
      tmp >> sab;

      res1 += h[a][b] * nbPatterns(sab,s1);
      res2 += h[a][b] * nbPatterns(sab,s2);
    }
  }
  return res1+res2;
}

float TichitFerraro1::matchTreeVert(int fNum, int i, int j) {
  float res1 = 0.0;
  float res2 = 0.0;

  string s1 = rootF1[fNum][idx(L1[i], 1-1)][idx(L2[j], j-1)];
  string s2 = rootF2[fNum][idx(L1[i], 1-1)][idx(L2[j], j-1)];

  for (int a = 0; a <=1; a++) {

      strstream tmp;
      string sa;
      tmp << a;
      tmp >> sa;

      res1 += v[a][1] * nbPatterns(sa,s1);
      res2 += v[a][1] * nbPatterns(sa,s2);
    }

  return res1+res2;
}

float TichitFerraro1::delTreeVert(int fNum, int i, int j) {
  float res = 0.0;

  string s = rootF1[fNum][idx(L1[i], 1-1)][idx(L2[j], j)];

  for (int a = 0; a <=1; a++) {

      strstream tmp;
      string sa;
      tmp << a;
      tmp >> sa;

      res += v[a][0] * nbPatterns(sa,s);
    }

  return res;
}

float TichitFerraro1::insTreeVert(int fNum, int i, int j) {
  float res = 0.0;

  string s = rootF2[fNum][idx(L1[i], 1)][idx(L2[j], j-1)];

  for (int a = 0; a <=1; a++) {
    strstream tmp;
      string sa;
      tmp << a;
      tmp >> sa;

      res += v[a][0] * nbPatterns(sa,s);

  }
  return res;
}


float TichitFerraro1::delForestHoriz(int fNum, int i, int i1, int j, int j1) {
  float res = 0.0;

  string prev = rootF1[fNum][idx(L1[i], i1-1)][idx(L2[j], j1)];
  string s = string(prev.end() - superGraph->outdeg(tree1[i1]), prev.end());

  for (int a = 0; a <=1; a++) {
    for (int b = 0; b <=1; b++) {
      strstream tmp;
      string sab;
      tmp << a << b;
      tmp >> sab;

      res += h[a][b] * nbPatterns(sab,s);
    }
  }
  return res;
}

float TichitFerraro1::insForestHoriz(int fNum, int i, int i1, int j, int j1) {
  float res = 0.0;

  string prev = rootF2[fNum][idx(L1[i], i1)][idx(L2[j], j1-1)];
  string s = string(prev.end() - superGraph->outdeg(tree1[i1]), prev.end());

  for (int a = 0; a <=1; a++) {
    for (int b = 0; b <=1; b++) {
      strstream tmp;
      string sab;
      tmp << a << b;
      tmp >> sab;

      res += h[a][b] * nbPatterns(sab,s);
    }
  }
  return res;
}

float TichitFerraro1::delForestVert(int fNum, int i, int i1, int j, int j1) {
  float res = 0.0;

  string prev = rootF1[fNum][idx(L1[i], i1-1)][idx(L2[j], j1)];
  string s = string(prev.end() - superGraph->outdeg(tree1[i1]), prev.end());

  for (int a = 0; a <=1; a++) {
    strstream tmp;
    string sa;
    tmp << a;
    tmp >> sa;

    res += v[a][0] * nbPatterns(sa,s);
  }
  return res;
}

float TichitFerraro1::insForestVert(int fNum, int i, int i1, int j, int j1) {
  float res = 0.0;

  string prev = rootF2[fNum][idx(L1[i], i1)][idx(L2[j], j1-1)];
  string s = string(prev.end() - superGraph->outdeg(tree1[i1]), prev.end());

  for (int a = 0; a <=1; a++) {
    strstream tmp;
    string sa;
    tmp << a;
    tmp >> sa;

    res += v[a][0] * nbPatterns(sa,s);
  }
  return res;
}

int TichitFerraro1::nbPatterns(string pattern, string str) {
  int count = 0;
  unsigned int i = 0;
  unsigned int maxBound = str.length() - pattern.length() + 1;
  while (i < str.length() - pattern.length() + 1) {
    unsigned int search = str.find(pattern,i);
    if(search < maxBound && search >= i) {
      count++;
      i = search+1;
    } else {
      break;
    }
  }
  return count;
}


float TichitFerraro1::del(const node& i) {
  //return CostMatrix[symbols[sProxy->getNodeValue(i)]][symbols["-"]];
  if (_scale)
    return del_complex(i);
  else
    return(-1.);
}

float TichitFerraro1::ins(const node& j) {
  //  return CostMatrix[symbols["-"]][symbols[sProxy->getNodeValue(j)]];
  if (_scale)
    return ins_complex(j);
  else
    return(-1.);
}

float TichitFerraro1::subst(const node& i, const node& j) {
  //  return CostMatrix[symbols[sProxy->getNodeValue(i)]][symbols[sProxy->getNodeValue(j)]];
  if (_scale)
    return subst_complex(i,j);
  else
    return(2.);
}

int TichitFerraro1::idx(int a, int b) const {
  return (a>b)?0:(b-a+1);
}


float TichitFerraro1::del_complex(const node& i) {
  //return CostMatrix[symbols[sProxy->getNodeValue(i)]][symbols["-"]];
  VId root1 = superGraph->getVertexId(i);
  TreeGraph *Tree1=new TreeGraph(*_mtg,root1,TOPO);
  return (-(Tree1->getNbVertex()));
}

float TichitFerraro1::ins_complex(const node& j) {
  //  return CostMatrix[symbols["-"]][symbols[sProxy->getNodeValue(j)]];
  VId root1 = superGraph->getVertexId(j);
  TreeGraph *Tree1=new TreeGraph(*_mtg,root1,TOPO);
  return(-(Tree1->getNbVertex()));
}

float TichitFerraro1::subst_complex(const node& i, const node& j) {
  //  return CostMatrix[symbols[sProxy->getNodeValue(i)]][symbols[sProxy->getNodeValue(j)]];
  VId root1 = superGraph->getVertexId(i);
  TreeGraph *Tree1=new TreeGraph(*_mtg,root1,TOPO);
  VId root2 = superGraph->getVertexId(j);
  TreeGraph *Tree2=new TreeGraph(*_mtg,root2,TOPO);
  TichitFerraro1* M=new TichitFerraro1(*_mtg,*Tree1,*Tree2,0,0);
  float D=M->run();
  delete (TichitFerraro1*) M;
  delete (TreeGraph*) Tree1;
  delete (TreeGraph*) Tree2;
  return(D);
}

/* Pour test
string TichitFerraro1::setRootF1(int fNum, int i, int i1, int j, int j1) {
  string prev = rootF1[fNum][idx(L1[i], i1-1)][idx(L2[j], j1)];
  if(L1[i1] <= L1[i1-1]) {
    return string(prev.begin(), prev.end() - tree1[i1].outdeg()) + "0";
  }
  return prev + "0";
}
*/

string TichitFerraro1::setRootF1(int fNum, int i, int i1, int j, int j1) {
  string prev = rootF1[fNum][idx(L1[i], i1-1)][idx(L2[j], j1)];
  return string(prev.begin(), prev.end() - superGraph->outdeg(tree1[i1])) + "0";
}

/*
string TichitFerraro1::setRootF2(int fNum, int i, int i1, int j, int j1) {
  string prev = rootF2[fNum][idx(L1[i], i1)][idx(L2[j], j1-1)];
  if(L2[j1] <= L2[j1-1]) {
    return string(prev.begin(), prev.end() - tree2[j1].outdeg()) + "0";
  }
  return prev + "0";
}
*/

string TichitFerraro1::setRootF2(int fNum, int i, int i1, int j, int j1) {
  string prev = rootF2[fNum][idx(L1[i], i1)][idx(L2[j], j1-1)];
  return string(prev.begin(), prev.end() - superGraph->outdeg(tree2[j1])) + "0";
}
/*
float TichitFerraro1::conex(int fNum, int i1, int j1) {

  string s1 = rootF1[fNum][idx(L1[i1], i1-1)][idx(L2[j1], j1-1)];
  string s2 = rootF2[fNum][idx(L1[i1], i1-1)][idx(L2[j1], j1-1)];

  int nb1 = 0;

  for(string::iterator it = s1.begin(); it != s1.end(); ++it) {
    if((*it) == '1') {
      nb1++;
    }
  }
  int nb2 = 0;
  for(string::iterator it = s2.begin(); it != s2.end(); ++it) {
    if((*it) == '1') {
      nb2++;
    }
  }
  return (2-nb1-nb2)*conex_cost;
}
*/

int TichitFerraro1::fix(int a, int b) const {
  return a*sizeKeyroots2 +b;
}


void TichitFerraro1::align(int t1, int i, int t2, int j) {
  //cerr << "ALIGN: t1: " << t1 << " i: " << i << " t2: " << t2 << " j: " << j << endl;
  if (i == 0 || j == 0) {
    return;
  }
  int lbi = L1[t1];
  int lbj = L2[t2];
  if (L1[i] < lbi || L2[j] < lbj) {
    return;
  }
  int fNum = getidx(t1,t2);
  float f =  fDist[fNum][idx(lbi,i)][idx(lbj,j)];
  float fi_1j_1 = fDist[fNum][idx(lbi,i-1)][idx(lbj,j-1)];
  float fi_1 = fDist[fNum][idx(lbi,i-1)][idx(lbj,j)];
  float fj_1 = fDist[fNum][idx(lbi,i)][idx(lbj,j-1)];
  /* si foret[lbi...i] et foret[lbj...j] sont des arbres */
  if(L1[i]==lbi && L2[j]==lbj) {
  if (fabs(f - fi_1 - del(tree1[i]) - delTreeHoriz(fNum,i,j) - delTreeVert(fNum,i,j)) < epsilon) {// Order is very important!!! test ins & del before match !!
      align(t1,i-1,t2,j);
      return;
    }
    if (fabs(f - fj_1 - ins(tree2[j]) - insTreeHoriz(fNum,i,j) - insTreeVert(fNum,i,j))< epsilon) {
      align(t1,i,t2,j-1);
      return;
    }
    if (fabs(f - fi_1j_1 - subst(tree1[i],tree2[j]) - matchTreeHoriz(fNum,i,j) - matchTreeVert(fNum,i,j)) < epsilon) {
      totoro++;
      //cerr << "prout prout" << endl;
      /*TreeNode* _node1 =  superGraph->T1->getNode(tree1[i]);
      TreeNode* _node2 =  superGraph->T2->getNode(tree2[j]-sizeT1);
      const Feature* f1 = superGraph->T1->getMTG()->si_feature(_node1->getVertex(),"base1id");
      const Feature* f2 = superGraph->T2->getMTG()->si_feature(_node2->getVertex(),"base1id");
      sequence->append( f1->i,f2->i, subst(i,j)); 
      const Feature* f3 = superGraph->T1->getMTG()->si_feature(_node1->getVertex(),"base2id");
      const Feature* f4 = superGraph->T2->getMTG()->si_feature(_node2->getVertex(),"base2id");
      sequence->append( f3->i,f4->i, subst(i,j));*/ 
      sequence->append(tree1[i],tree2[j]-sizeT1,subst(i,j));
      //metricProxy->setNodeValue(tree1[i],match_color);
      //metricProxy->setNodeValue(tree2[j],match_color);
      align(t1,i-1,t2,j-1);
      return;
    }
    //cerr << "ERROR5" << endl;
    return;
  }

  int newt1;
  int newt2;
  if (lbi == L1[i]) {
    newt1 = t1;
  }
  else {
    newt1 = getkey1(i);
  }
  if (lbj == L2[j]) {
    newt2 = t2;
  }
  else {
    newt2 = getkey2(j);
  }
  float t = tDist[i][j];
  float fprev = fDist[fNum][idx(lbi,L1[i]-1)][idx(lbj,L2[j]-1)];
//    cerr << "\tf: " << f << endl;
//    cerr << "\tfprev: " << fprev << endl;
//    cerr << "\tt: " << t << endl;
//    cerr << "\tfi_1: " << fi_1 << endl;
//    cerr << "\tfj_1: " << fj_1 << endl;
  if (fabs(f - fprev - t) < epsilon) {
    align(newt1,i,newt2,j);

    align(t1,L1[i]-1,t2,L2[j]-1);
    return;
  }
  if (fabs(f - fi_1 - del(tree1[i]) - delForestHoriz(fNum,t1,i,t2,j) - delForestVert(fNum,t1,i,t2,j)) < epsilon) {
    align(t1,i-1,t2,j);
    return;
  }
  if (fabs(f - fj_1 - ins(tree2[j]) - insForestHoriz(fNum,t1,i,t2,j) - insForestVert(fNum,t1,i,t2,j)) < epsilon) {
    align(t1,i,t2,j-1);
    return;
  }
  //cerr << "STRANGE" << endl;
  return;
}

int TichitFerraro1::getidx(int i, int j) {
  int i1,j1;
  for (i1 = sizeKeyroots1 - 1; keyroots1[i1] != i; --i1);
  for (j1 = sizeKeyroots2 - 1; keyroots2[j1] != j; --j1);
  return fix(i1,j1);
}

int TichitFerraro1::getkey1(int n) {
  int val = L1[n];
  int i;
  for (i = sizeT1; L1[i] != val; i--);
  return i;
}

int TichitFerraro1::getkey2(int n) {
  int val = L2[n];
  int i;
  for (i = sizeT2; L2[i] != val; i--);
  return i;
}



void TichitFerraro1::displayMatrix(vector<vector<float> > tDist) {
  int s1 = tDist.size();
  int s0 = tDist[0].size();
  cerr << "          " ;
  int j = 0;
  for ( j=1; j<s0; ++j) {
    cerr  << "  "<< tree2[j];
  }
  cerr << endl << endl << "  ";
  for ( j=0; j<s0; ++j) {
    fprintf(stderr," %2.1f",tDist[0][j]);
  }
  cerr << endl;
  for (int i=1; i<s1; ++i) {
    cerr << " " << tree1[i];
    int s2 = tDist[i].size();
    for ( j=0; j<s2; ++j) {
      fprintf(stderr," %2.1f",tDist[i][j]);
    }
    cerr << endl;
  }
  cerr << endl;
}


void TichitFerraro1::displayMatrix(vector<vector<string> > tDist) {
  int s1 = tDist.size();
  int s0 = tDist[0].size();
  cerr << "          " ;
  int j = 0;
  for ( j=1; j<s0; ++j) {
    cerr  << "  "<< tree2[j];
  }
  cerr << endl << endl << "  ";
  for ( j=0; j<s0; ++j) {
    cerr << " " << tDist[0][j];
  }
  cerr << endl;
  for (int i=1; i<s1; ++i) {
    cerr << " " << tree1[i];
    int s2 = tDist[i].size();
    for (int j=0; j<s2; ++j) {
      cerr << " " << tDist[i][j];
    }
    cerr << endl;
  }
  cerr << endl;
}

