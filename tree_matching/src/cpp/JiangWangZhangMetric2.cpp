/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): L. Tichit (P.Ferraro) (laurent.tichit@labri.fr)
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

#include "JiangWangZhangMetric2.h"

#include <iostream>
#include <strstream>
#include <cstdio>

using namespace std;

JiangWangZhangMetric2::JiangWangZhangMetric2(TreeGraph & input, TreeGraph & reference)
{
  superGraph = new SuperGraph(input,reference);
  sequence = new Sequence();
  this->init();
}

JiangWangZhangMetric2::~JiangWangZhangMetric2() {
}

const float JiangWangZhangMetric2::match_color = 0.0;
const float JiangWangZhangMetric2::diff_color = 1024.0;

float JiangWangZhangMetric2::getEdgeValue(const edge& e) {
  return 0.0;
}

float JiangWangZhangMetric2::getNodeValue(const node& n)
{
  return 0.0;
}

void JiangWangZhangMetric2::reset() {

}

Sequence* JiangWangZhangMetric2::getSequence() {
  return(sequence);
}

float JiangWangZhangMetric2:: CostMatrix[NB_SYMBOLS][NB_SYMBOLS]
 = {
  { 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
  {-1,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0, -2},
  {-1, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,  0},
};


void JiangWangZhangMetric2::init() {
  string erreurMsg;
  //Return the size of T1
  sizeT1 = superGraph->getNbVertex(1);
  tree1.resize(sizeT1+1);
  tDist.resize(sizeT1+1);
  //return the size of T2
  sizeT2 = superGraph->getNbVertex(2);
  tree2.resize(sizeT2+1);
  for(int i = 0; i<=sizeT1; i++) {
//     tree1[i]=i-1;
    tDist[i].resize(sizeT2+1);
  }
//   for(int j = 0; j<=sizeT2; j++) {
//     tree2[j]=j+sizeT1-1;
//   }

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
//   symbols["default"]=19;
//   symbols["AG"] = 20;
//   symbols["GA"] = 21;
//   symbols["AC"] = 22;
//   symbols["CA"] = 23;
 }


/* map a 2-coordinate set to a 1-coordinate set == hashfun
   size of the ordered parts of a set,
   s = number of sons of a tree => number of possible forests of its sons */
int JiangWangZhangMetric2::fSize(int s) {
  return (s*s +s)/2;
}

/* map a 2-coordinate set to a 1-coordinate set == hashfun */
int JiangWangZhangMetric2::id1(int i, int a, int b) const {
  if (a<=b) {
    return i*f1MaxSize + a - (b*b +b)/2;
  }
  return 0;
}

/* map a 2-coordinate set to a 1-coordinate set == hashfun */
int JiangWangZhangMetric2::id2(int j, int a, int b) const {
  if (a<=b) {
    return j*f2MaxSize + a - (b*b +b)/2;
  }
  return 0;
}

/* return the index of the forest of the sons of the node i  */
int JiangWangZhangMetric2::idx1(int i) const {
  if (i == 0) {
    return 0;
  }
  return id1(i, 1, superGraph->outdeg(tree1[i]));
}

int JiangWangZhangMetric2::idx2(int j) const {
  if (j == 0) {
    return 0;
  }
  return id2(j, 1, superGraph->outdeg(tree2[j]));
}

float JiangWangZhangMetric2::sumTree1(int i) const {
  float sum = 0.0;
  int max = superGraph->outdeg(tree1[i]);
  for (int ii = 1; ii <= max; ++ii) {
    sum += tDist[t1[i][ii]][0].dist;
  }
  return sum;
}

float JiangWangZhangMetric2::sumTree2(int j) const {
  float sum = 0.0;
  int max = superGraph->outdeg(tree2[j]);
  for (int jj = 1; jj <= max; ++jj) {
    sum += tDist[0][t2[j][jj]].dist;
  }
  return sum;
}

void JiangWangZhangMetric2::preCompute() {

  t1.resize(sizeT1+1);
  t2.resize(sizeT2+1);
  int i= 0;
  for ( i = 1; i <= sizeT1; ++i) {
    int size = superGraph->outdeg(tree1[i]);
    t1[i].resize(size + 1);
    t1[i][0] = 0;
    if (size != 0) {
      int j = 1;
      list<node> itn  = superGraph->getOutNodes(tree1[i]);
      list<node>::iterator begin,end;
      begin = itn.begin();
      end = itn.end();
      // DUR DUR
      for (; begin!=end; ++j) {
        node n = *(begin);
        begin++;
        //cerr << n.id << endl;
        int jj;
        for (jj = 1; jj<i;++jj) {
          if(tree1[jj] == n) {
            break;
          }
          else {
            //cerr << "crotte 1" << endl;
          }
        }
        t1[i][j] = jj;
        //cerr << "t1[" << i << "][" << j << "] = " << jj << endl;
      }
      //delete itn;
    }
  }
  int j= 1;
  for ( j = 1; j <= sizeT2; ++j) {
    int size = superGraph->outdeg(tree2[j]);
    t2[j].resize(size + 1);
    t2[j][0] = 0;
    if (size != 0) {
      int i = 1;
      list<node> itn  = superGraph->getOutNodes(tree2[j]);
      list<node>::iterator begin,end;
      begin = itn.begin();
      end = itn.end();
//       list<node>::const_iterator *itn;
      // DUR DUR
      for (; begin!=end; ++i) {
        node n = *begin;
        begin++;
        int ii;
        for (ii = 1; ii<j;++ii) {
          if(tree2[ii] == n) {
            break;
          }
          else {
            //cerr << "crotte 2" << endl;
          }
        }
        t2[j][i] = ii;
        //cerr << "t2[" << j << "][" << i << "] = " << ii << endl;
      }
 //      delete itn;
    }
  }

  /* look for the max degree of the trees*/
  int maxDeg1 = 0;
  for (i = 1; i <= sizeT1; ++i) {
    if (maxDeg1 < superGraph->outdeg(tree1[i])) {
      maxDeg1 = superGraph->outdeg(tree1[i]);
    }
  }
  int maxDeg2 = 0;
  for (j = 1; j <= sizeT2; ++j) {
    if (maxDeg2 < superGraph->outdeg(tree2[j])) {
      maxDeg2 = superGraph->outdeg(tree2[j]);
    }
  }
  f1MaxSize = fSize(maxDeg1);
  f2MaxSize = fSize(maxDeg2);
  int f1Size = f1MaxSize*sizeT1;
  int f2Size = f2MaxSize*sizeT2;
  fTmpDist.resize(f1Size+1);
  int k= 0;
  for (k = 0; k <= f1Size; k++) {
    fTmpDist[k].resize(f2Size+1);
  }
  fDist.resize(sizeT1+1);
  for (k = 0; k <= sizeT1; k++) {
    fDist[k].resize(sizeT2+1);
  }
}

void JiangWangZhangMetric2::displayMatrix(vector<vector<ftriplet> > tDist) {
  int s1 = tDist.size();
  int s0 = tDist[0].size();
  cerr << "     " ;
  int j= 1;
  for (j=1; j<s0; ++j) {
    cerr  << "  "<< tree2[j];
  }
  cerr << std::endl << std::endl << "   ";
  for (j=0; j<s0; ++j) {
    fprintf(stderr," %2.1f", tDist[0][j].dist);
  }
  cerr << endl;
  for (int i=1; i<s1; ++i) {
    cerr << "  " <<tree1[i];
    int s2 = tDist[i].size();
    for (j=0; j<s2; ++j) {
      fprintf(stderr," %2.1f", tDist[i][j].dist);
    }
    cerr << std::endl;
  }
  cerr << std::endl;
}

void JiangWangZhangMetric2::displayMatrix2(vector<vector<ftriplet> > tDist) {
  int s1 = tDist.size();
  int s0 = tDist[0].size();
  cerr << "     " ;
  int j= 1;
  for (j=1; j<s0; ++j) {
    cerr  << "  " << tree2[j];
  }
  cerr << endl << endl << "   ";
  for (j=0; j<s0; ++j) {
    cerr << tDist[0][j].howMuch;
  }
  cerr << endl;
  for (int i=1; i<s1; ++i) {
    cerr << "  " << tree1[i];
    int s2 = tDist[i].size();
    for (int j=0; j<s2; ++j) {
      cerr << tDist[i][j].howMuch;
    }
    cerr << endl;
  }
  cerr << endl;
}

void JiangWangZhangMetric2::displayMatrix3(vector<vector<ftriplet> > tDist) {
  int s1 = tDist.size();
  int s0 = tDist[0].size();
  cerr << "     " ;
  int j= 1;
  for (j=1; j<s0; ++j) {
    cerr  << "  "<< tree2[j];
  }
  cerr << endl << endl ;
  for (j=0; j<s0; ++j) {
    cerr << tDist[0][j].aligned;
  }
  cerr << endl;
  for (int i=1; i<s1; ++i) {
    cerr << "  " << tree1[i];
    int s2 = tDist[i].size();
    for (int j=0; j<s2; ++j) {
      cerr << tDist[i][j].aligned;
    }
    cerr << endl;
  }
  cerr << endl;
}

void JiangWangZhangMetric2::displayMatrix(vector<vector<float> > tDist) {
  int s1 = tDist.size();
  int s0 = tDist[0].size();
  cerr << "     " ;
  int j= 1;
  for (j=1; j<s0; ++j) {
//     cerr  << "  "<< sProxy->getNodeValue(tree2[j]);
  }
  cerr << endl << endl << "   ";
  for (j=0; j<s0; ++j) {
    fprintf(stderr," %2.1f", tDist[0][j]);
  }
  cerr << endl;
  for (int i=1; i<s1; ++i) {
    //    cerr << "  " << sProxy->getNodeValue(tree1[i]);
    int s2 = tDist[i].size();
    for (int j=0; j<s2; ++j) {
      fprintf(stderr," %2.1f", tDist[i][j]);
    }
    cerr << endl;
  }
  cerr << endl;
}

void JiangWangZhangMetric2::displayMatrix(vector<vector<string> > tDist) {
  int s1 = tDist.size();
  int s0 = tDist[0].size();
  cerr << "     " ;
  int j= 1;
  for (j=1; j<s0; ++j) {
    //cerr  << "   "<< sProxy->getNodeValue(tree2[j]);
  }
  cerr << endl << endl << "   ";
  for (j=0; j<s0; ++j) {
    cerr << " " << tDist[0][j];
  }
  cerr << endl;
  for (int i=1; i<s1; ++i) {
    //cerr << " " << sProxy->getNodeValue(tree1[i]);
    int s2 = tDist[i].size();
    for (int j=0; j<s2; ++j) {
      cerr << " " << tDist[i][j];
    }
    cerr << endl;
  }
  cerr << endl;
}


void JiangWangZhangMetric2::displayMatrix2(vector<vector<int> > tDist) {
  int s1 = tDist.size();
  for (int i=0; i<s1; ++i) {
    fprintf(stderr,"%2d:",i);
    int s2 = tDist[i].size();
    for (int j=0; j<s2; ++j) {
      fprintf(stderr," %2d", tDist[i][j]);
    }
    cerr << endl;
  }
  cerr << endl;
}

int JiangWangZhangMetric2::storeNodes(node nd, vector<node>& tree, int i) {
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

double JiangWangZhangMetric2::run()  {
  tree1[0] = 0;
  tree2[0] = 0;
  root1 = 0;
  root2 = root1 + sizeT1;
  storeNodes(root1,tree1,1);
  storeNodes(root2,tree2,1);
  preCompute();
  //cerr << "end precompute" << endl;

  /*cerr << "t1: " << endl;
  displayMatrix2(t1);
  cerr << "t2: " << endl;
  displayMatrix2(t2);
  */
  computeTreeDist();
  extractFDist();
  /*cerr << "fDist:" << endl;
  displayMatrix(fDist);
  cerr << "tDist:" << endl;
  displayMatrix(tDist);
  displayMatrix2(tDist);
  displayMatrix3(tDist);*/
  //cerr << "fTmpDist:" << endl;
  /*displayMatrix(fTmpDist);
  displayMatrix2(fTmpDist);
  displayMatrix3(fTmpDist);*/
  initNodesMetric();
  alignTree(new TreeI(sizeT1,sizeT2));
  return (tDist[sizeT1][sizeT2].dist);
}

void JiangWangZhangMetric2::initNodesMetric() const {
//   list<node>::const_iterator *itn;
//   for (itn = superGraph->getNodes();itn->hasNext();) {
//     metricProxy->setNodeValue(itn->next(),diff_color);
//   }
//   delete itn;
}

void JiangWangZhangMetric2::extractFDist() {
  for(int i = 0; i <= sizeT1; ++i) {
    for (int j = 0; j <= sizeT2; ++j) {
      fDist[i][j] = fTmpDist[idx1(i)][idx2(j)].dist;
    }
  }
}

void JiangWangZhangMetric2::alignTree(TreeI *tt){
  int i = tt->i;
  int j = tt->j;
  //cerr << "i= " << i << " j= " << j << " aligned= " << tDist[i][j].aligned << " howMuch=" << tDist[i][j].howMuch << endl;
  if(tDist[i][j].aligned == 1) {
    sequence->append(tree1[i],tree2[j]-sizeT1,subst(i,j));
//     metricProxy->setNodeValue(tree1[i],match_color);
//     metricProxy->setNodeValue(tree2[j],match_color);
  }
  if(tDist[i][j].howMuch == 1) {
    tDist[i][j].second->align(this);
  }
}

void JiangWangZhangMetric2::alignForest(ForestI *ft) {
  int i = ft->i;
  int j = ft->j;
  int p = ft->p;
  int q = ft->q;
  int s = ft->s;
  int t = ft->t;
  //cerr << "i= " << i << " j= " << j << " p= " << p << " q= " << q << " s= " << s << " t= " << t << " howMuch=" << fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch << endl;
  if(fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch == 1) {
    fTmpDist[id1(i,s,p)][id2(j,t,q)].second->align(this);
  }
  else if (fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch == 2) {
    fTmpDist[id1(i,s,p)][id2(j,t,q)].second->align(this);
    fTmpDist[id1(i,s,p)][id2(j,t,q)].third->align(this);
  }
}

void JiangWangZhangMetric2::computeTreeDist() {
  tDist[0][0].dist = 0;
  tDist[0][0].howMuch = 0;
  tDist[0][0].aligned = 0;
  fTmpDist[0][0].dist = 0;
  fTmpDist[0][0].howMuch = 0;
  fTmpDist[0][0].aligned = 0;
  int i= 1;
  for ( i=1; i<=sizeT1; ++i) {
    fTmpDist[idx1(i)][0].dist = sumTree1(i);
    fTmpDist[idx1(i)][0].howMuch = 0;
    fTmpDist[idx1(i)][0].aligned = 0;
    tDist[i][0].dist = fTmpDist[idx1(i)][0].dist + del(tree1[i]);
    tDist[i][0].howMuch = 0;
    tDist[i][0].aligned = 0;
  }
  int j= 1;
  for (j=1; j<=sizeT2; ++j) {
    fTmpDist[0][idx2(j)].dist = sumTree2(j);
    fTmpDist[0][idx2(j)].howMuch = 0;
    fTmpDist[0][idx2(j)].aligned = 0;
    tDist[0][j].dist = fTmpDist[0][idx2(j)].dist + ins(tree2[j]);
    tDist[0][j].howMuch = 0;
    tDist[0][j].aligned = 0;
  }
  for (i=1; i<=sizeT1; ++i) {
    int maxi = superGraph->outdeg(tree1[i]);
    for (int s = 1; s <= maxi; ++s) {
      for (int p = s; p <= maxi; ++p) {
        fTmpDist[id1(i,s,p)][0].dist = (fTmpDist[id1(i,s,p-1)][0]).dist+(tDist[t1[i][p]][0]).dist;
        fTmpDist[id1(i,s,p)][0].howMuch = 0;
        fTmpDist[id1(i,s,p)][0].aligned = 0;
      }
    }
  }
  for (j=1; j<=sizeT2; ++j) {
    int maxj = superGraph->outdeg(tree2[j]);
    for (int t = 1; t <= maxj; ++t) {
      for (int q = t; q <= maxj; ++q) {
        fTmpDist[0][id2(j,t,q)].dist = (fTmpDist[0][id2(j,t,q-1)]).dist+(tDist[0][t2[j][q]]).dist;
        (fTmpDist[0][id2(j,t,q)]).howMuch = 0;
        (fTmpDist[0][id2(j,t,q)]).aligned = 0;
      }
    }
  }

  for (i=1; i<=sizeT1; ++i) {
    int maxi = superGraph->outdeg(tree1[i]);
    for (int j=1; j<=sizeT2; ++j) {
      int maxj = superGraph->outdeg(tree2[j]);
      computeForestDist(i,1,j,1);
      for (int s = 2; s <= maxi; ++s) {
        computeForestDist(i,s,j,1);
      }
      int t= 2;
      for (t = 2; t <= maxj; ++t) {
        computeForestDist(i,1,j,t);
      }
      float minrj;
      Info *secondj = NULL;
      if (maxj == 0) {
        minrj = tDist[i][0].dist;
      }
      else {
        minrj = tDist[i][t2[j][1]].dist - tDist[0][t2[j][1]].dist;
        secondj = new TreeI(i,t2[j][1]);
      }
      int r= 2;
      for (r = 2; r <= maxj; ++r) {
        float tmp = tDist[i][t2[j][r]].dist - tDist[0][t2[j][r]].dist;
        Info *tmpSecondj = new TreeI(i,t2[j][r]);
        if (tmp > minrj) {
          minrj = tmp;
          Info * stmp = secondj;
          secondj = tmpSecondj;
          delete stmp;
        }
      }
      float minri;
      Info *secondi = NULL;
      if (maxi == 0) {
        minri = tDist[0][j].dist;
      }
      else {
        minri = tDist[t1[i][1]][j].dist - tDist[t1[i][1]][0].dist;
        secondi = new TreeI(t1[i][1],j);
      }
      for (r = 2; r <= maxi; ++r) {
        float tmp = tDist[t1[i][r]][j].dist - tDist[t1[i][r]][0].dist;
        Info *tmpSecondi = new TreeI(t1[i][r],j);
        if (tmp > minri) {
          minri = tmp;
          Info * stmp = secondi;
          secondi = tmpSecondi;
          delete stmp;
        }
      }
      float maxx = max(tDist[0][j].dist + minrj,
                           tDist[i][0].dist + minri,
                           (fTmpDist[idx1(i)][idx2(j)]).dist+ subst(tree1[i],tree2[j]));
      if (maxx == tDist[i][0].dist + minri) {
        tDist[i][j].second = secondi;
        if (tDist[i][j].second == NULL) {
          tDist[i][j].howMuch = 0;
        }
        else {
          tDist[i][j].howMuch = 1;
        }
        tDist[i][j].aligned = 0;
      }
      else if (maxx == tDist[0][j].dist + minrj) {
        tDist[i][j].second = secondj;
        if (tDist[i][j].second == NULL) {
          tDist[i][j].howMuch = 0;
        }
        else {
          tDist[i][j].howMuch = 1;
        }
        tDist[i][j].aligned = 0;
      }
      else {
        maxi = superGraph->outdeg(tree1[i]);
        maxj = superGraph->outdeg(tree2[j]);
        tDist[i][j].second = new ForestI(i,1,maxi,j,1,maxj);
        tDist[i][j].howMuch = 1;
        if (subst(tree1[i],tree2[j]) >= 0) {
          tDist[i][j].aligned = 1;
        }
        else {
          tDist[i][j].aligned = 0;
        }
      }
      tDist[i][j].dist = maxx;
      /*if (tDist[i][j].second != NULL) {
        tDist[i][j].second->display();
        }*/
    }
  }
}

void JiangWangZhangMetric2::computeForestDist(int i, int s, int j, int t) {
  int maxi = superGraph->outdeg(tree1[i]);
  int maxj = superGraph->outdeg(tree2[j]);

  for (int p = s; p <= maxi; ++p) {
    for (int q = t; q <= maxj; ++q) {
      //cerr << " i=" << i << " s=" << s << " p=" << p << " j=" << j << " t=" << t << " q=" << q << endl;
      float minks = fTmpDist[0][id2(j,t,q-1)].dist + fTmpDist[id1(i,s,p)][idx2(t2[j][q])].dist;
      //cerr << fTmpDist[0][id2(j,t,q-1)].dist << fTmpDist[id1(i,s,p)][idx2(t2[j][q])].dist << endl;
      Info *seconds = new ForestI(i,s,p,t2[j][q],1,superGraph->outdeg(tree2[t2[j][q]]));
      Info *thirds = NULL;

      int k= 2;
      for (k = 2; k < p; ++k) {
        float tmp = fTmpDist[id1(i,s,k-1)][id2(j,t,q-1)].dist + fTmpDist[id1(i,k,p)][idx2(t2[j][q])].dist;
        Info *tmpSeconds = new ForestI(i,k,p,t2[j][q],1,superGraph->outdeg(tree2[t2[j][q]]));
        Info *tmpThirds = new ForestI(i,s,k-1,j,t,q-1);
        if (tmp > minks) {
          minks = tmp;
          Info * stmp = seconds;
          Info * ttmp = thirds;
          seconds = tmpSeconds;
          thirds = tmpThirds;
          delete stmp;
          delete ttmp;
        }
      }
      float minkt = fTmpDist[id1(i,s,p-1)][0].dist + fTmpDist[idx1(t1[i][p])][id2(j,t,q)].dist;
      //cerr << fTmpDist[id1(i,s,p-1)][0].dist << fTmpDist[idx1(t1[i][p])][id2(j,t,q)].dist << endl;
      Info *secondt = new ForestI(t1[i][p],1,superGraph->outdeg(tree1[t1[i][p]]),j,t,q);
      Info *thirdt = NULL;

      for (k = t+1; k < q; ++k) {
        float tmp = fTmpDist[id1(i,s,p-1)][id2(j,t,k-1)].dist + fTmpDist[idx1(t1[i][p])][id2(j,k,q)].dist;
        Info * tmpSecondt = new ForestI(t1[i][p],1,superGraph->outdeg(tree1[t1[i][p]]),j,k,q);
        Info  * tmpThirdt = new ForestI(i,s,p-1,j,t,k-1);
        if (tmp > minks) {
          minks = tmp;
          Info * stmp = secondt;
          Info * ttmp = thirdt;
          secondt = tmpSecondt;
          thirdt = tmpThirdt;
          delete stmp;
          delete ttmp;
        }
      }
      float maxx = max(fTmpDist[id1(i,s,p-1)][id2(j,t,q)].dist + tDist[t1[i][p]][0].dist,
                        fTmpDist[id1(i,s,p)][id2(j,t,q-1)].dist + tDist[0][t2[j][q]].dist,
                        fTmpDist[id1(i,s,p-1)][id2(j,t,q-1)].dist + tDist[t1[i][p]][t2[j][q]].dist,
                        ins(tree2[t2[j][q]]) + minks,
                        del(tree1[t1[i][p]]) + minkt);
      if (maxx == fTmpDist[id1(i,s,p-1)][id2(j,t,q-1)].dist + tDist[t1[i][p]][t2[j][q]].dist) {
        fTmpDist[id1(i,s,p)][id2(j,t,q)].second = new ForestI(i,s,p-1,j,t,q-1);
        fTmpDist[id1(i,s,p)][id2(j,t,q)].third = new TreeI(t1[i][p],t2[j][q]);
        fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch = 2;
      }
      else if (maxx == ins(tree2[t2[j][q]]) + minks) {
        fTmpDist[id1(i,s,p)][id2(j,t,q)].second = seconds;
        fTmpDist[id1(i,s,p)][id2(j,t,q)].third = thirds;
        if (thirds == NULL) {
          fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch = 1;
          fTmpDist[id1(i,s,p)][id2(j,t,q)].third = NULL;
        }
        else {
          fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch = 2;
        }
      }
      else if (maxx == del(tree1[t1[i][p]]) + minkt) {
        fTmpDist[id1(i,s,p)][id2(j,t,q)].second = secondt;
        fTmpDist[id1(i,s,p)][id2(j,t,q)].third = thirdt;
        if (thirdt == NULL) {
          fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch = 1;
          fTmpDist[id1(i,s,p)][id2(j,t,q)].third = NULL;
        }
        else {
          fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch = 2;
        }
      }
      else if (maxx == fTmpDist[id1(i,s,p)][id2(j,t,q-1)].dist + tDist[0][t2[j][q]].dist) {
        fTmpDist[id1(i,s,p)][id2(j,t,q)].second = new ForestI(i,s,p,j,t,q-1);
        fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch = 1;
        fTmpDist[id1(i,s,p)][id2(j,t,q)].third = NULL;
      }
      else {
        fTmpDist[id1(i,s,p)][id2(j,t,q)].second = new ForestI(i,s,p-1,j,t,q);
        fTmpDist[id1(i,s,p)][id2(j,t,q)].howMuch = 1;
        fTmpDist[id1(i,s,p)][id2(j,t,q)].third = NULL;
      }
      fTmpDist[id1(i,s,p)][id2(j,t,q)].dist = maxx;
      /*fTmpDist[id1(i,s,p)][id2(j,t,q)].second->display();
      if (fTmpDist[id1(i,s,p)][id2(j,t,q)].third != NULL) {
        fTmpDist[id1(i,s,p)][id2(j,t,q)].third->display();
        } */
    }
  }
}

/* return the max of 5 values */
float JiangWangZhangMetric2::max(float f1, float f2, float f3, float f4, float f5) const {
  return max(max(f1,f2,f3),f4,f5);
}
/* return the max of 4 values */
float JiangWangZhangMetric2::max(float f1, float f2, float f3, float f4) const {
  return max(max(f1,f2),f3,f4);
}

/* return the max of 3 values */
float JiangWangZhangMetric2::max(float f1, float f2, float f3) const {
  return max(max(f1,f2),f3);
}

/* return the max of 2 values */
float JiangWangZhangMetric2::max(float f1, float f2) const {
  return (f1 > f2)?f1:f2;
}



float JiangWangZhangMetric2::del(const node& i) {
//   return CostMatrix[symbols[sProxy->getNodeValue(i)]][symbols["-"]];
  return(-1.);
}

float JiangWangZhangMetric2::ins(const node& j) {
//   return CostMatrix[symbols["-"]][symbols[sProxy->getNodeValue(j)]];
  return(-1.);
}


float JiangWangZhangMetric2::subst(const node& i, const node& j) {
//   return CostMatrix[symbols[sProxy->getNodeValue(i)]][symbols[sProxy->getNodeValue(j)]];
  return(0.0);
}

//--------------------------------
ForestI::ForestI(int ip,int sp,int pp,int jp,int tp,int qp): TreeI(ip,jp),p(pp),q(qp),s(sp),t(tp) {}

void ForestI::align(JiangWangZhangMetric2 *jwz) {
  //cerr << "ForestI::align" << endl;
  jwz->alignForest(this);
}

void ForestI::display() {
  cerr << "Forest: i=" << i << " s=" << s << " p=" << p << " j=" << j << " t=" << t << " q=" << q << endl;
}
//---------------------------------
TreeI::TreeI(int ip, int jp): i(ip), j(jp) {}

void TreeI::align(JiangWangZhangMetric2 *jwz) {
  //cerr << "TreeI::align" << endl;
  jwz->alignTree(this);
}

void TreeI::display() {
  cerr << "Tree: i=" << i << " j=" << j << endl;
}







