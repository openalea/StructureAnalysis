/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Tichit Laurent and Ferraro Pascal (laurent.tichit@labri.fr)
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




#ifndef _JIANGWANGZHANGMETRIC2_H
#define _JIANGWANGZHANGMETRIC2_H

#include "supergraph.h"
#include "treegraph.h"
#include "sequence.h"

///
#define NB_SYMBOLS 24

/**
 *\class JiangWangZhangMetric2
 *\brief Definition of teh distance metric between ordered tree graphs defined by Jiang et al.
 *\authors Laurent Tichit (and Pascal Ferraro)
 *\date 1st version 2000 (Tulip), updated 2003 (AMAPmod)
 */

class JiangWangZhangMetric2;

struct Info {
  virtual void display() = 0;
  virtual void align(JiangWangZhangMetric2 *) = 0;
};

class TreeI : public Info {
public:
  int i;
  int j;
  TreeI(int,int);
  void display();
  void align(JiangWangZhangMetric2 *);
};

class ForestI:public TreeI {
public:
  int p;
  int q;
  int s;
  int t;
  ForestI(int,int,int,int,int,int);
  void display();
  void align(JiangWangZhangMetric2 *);
};

struct ftriplet {
  float dist; // alignment distance
  int aligned;  // 1 == possibly aligned
  int howMuch;  // 0 1 2
  Info *second;
  Info *third;
};


class JiangWangZhangMetric2 {
 public:
  JiangWangZhangMetric2(TreeGraph &, TreeGraph &);
  ~JiangWangZhangMetric2();
  double  run();
  float getNodeValue(const node& n);
  float getEdgeValue(const edge& n);
  void  reset();
  bool  check(std::string &);
  void  alignTree(TreeI *);
  void  alignForest(ForestI *);
  Sequence* getSequence();

private:
  static float CostMatrix[NB_SYMBOLS][NB_SYMBOLS];
  //map<string, int> symbols;
  node root1;
  node root2;
  int sizeT1;
  int sizeT2;
  std::vector<node> tree1;
  std::vector<node> tree2;
  std::vector<std::vector<ftriplet> > tDist;

  int  storeNodes (node, std::vector<node>&, int);

  void initNodesMetric() const;

  std::vector<std::vector<float> > fDist;

  void preCompute();
  void computeTreeDist();
  int  fSize(int);
  void displayMatrix(std::vector<std::vector<ftriplet> >);
  void displayMatrix2(std::vector<std::vector<ftriplet> >);
  void displayMatrix3(std::vector<std::vector<ftriplet> >);
  void displayMatrix(std::vector<std::vector<float> >);
  void displayMatrix(std::vector<std::vector<std::string> >);

  int f1MaxSize;
  int f2MaxSize;
  std::vector<std::vector<ftriplet> > fTmpDist;
  std::vector<std::vector<int> > t1;
  std::vector<std::vector<int> > t2;
  void extractFDist();
  void computeForestDist(int, int, int, int);
  void displayMatrix2(std::vector<std::vector<int> >);
  float sumTree1(int) const;
  float sumTree2(int) const;
  int idx1(int) const;
  int idx2(int) const;
  int id1(int, int, int) const;
  int id2(int, int, int) const;
  static const float match_color;
  static const float diff_color;
  void  init();
  float subst(const node&, const node&);
  float ins(const node&);
  float del(const node&);
  float max(float, float, float, float, float) const;
  float max(float, float, float, float) const;
  float max(float, float, float) const;
  float max(float, float) const;
  SuperGraph   *superGraph;
  Sequence *sequence;
};

#endif







