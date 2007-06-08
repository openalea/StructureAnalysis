//-*-c++-*-
#ifndef _TICHITFERRARO1_H
#define _TICHITFERRARO1_H

#include "supergraph.h"
#include "treegraph.h"
#include "sequence.h"

#define NB_SYMBOLS 24

class TichitFerraro1 {
 public:
  TichitFerraro1(MTG&,TreeGraph &, TreeGraph &,int,int);
  ~TichitFerraro1();
  Sequence* getSequence();
  double run();
  float getNodeValue(const node& n);
  float getEdgeValue(const edge& n);
  void  reset();
  bool  check(std::string &);
  SuperGraph   *superGraph;
private:
  float lsmax;
  int vmax;
  int wmax;
  static float v[2][2];
  static float h[2][2];

  static float CostMatrix[NB_SYMBOLS][NB_SYMBOLS];
  //  map<string, int> symbols;
  node root1;
  node root2;
  int sizeT1;
  int sizeT2;
  std::vector<node> tree1;
  std::vector<node> tree2;
  std::vector<std::vector<float> > tDist;
  //  StringProxy *sProxy;

  std::vector<int> L1;
  std::vector<int> keyroots1;
  int sizeKeyroots1;
  std::vector<int> L2;
  std::vector<int> keyroots2;
  int sizeKeyroots2;

  int getkey1(int);
  int getkey2(int);
  void preCompute();
  void computeTreeDist(int, int);
  void computeTreeDist();
  int  idx(int, int) const;
  int fix(int, int) const;
  int           storeNodes (node, std::vector<node>&, int);
  void          initNodesMetric() const;
  void          init();
  void         displayMatrix(std::vector<std::vector<float> >);
  void         displayMatrix(std::vector<std::vector<std::string> >);


  std::vector<std::vector<std::vector<float> > > fDist;

  std::vector<std::vector<std::vector<std::string> > > rootF1;

  std::vector<std::vector<std::vector<std::string> > > rootF2;

  std::vector<std::vector<std::string> > rootT1;
  std::vector<std::vector<std::string> > rootT2;
  float conex(int, int, int);
  std::string setRootF1(int, int, int, int, int);
  std::string setRootF2(int, int, int, int, int);

  int getidx(int, int);
  void align(int, int, int, int);

  float insTreeVert(int, int, int) ;
  float delTreeVert(int, int, int) ;
  float matchTreeVert(int, int, int) ;

  float insTreeHoriz(int, int, int) ;
  float delTreeHoriz(int, int, int) ;
  float matchTreeHoriz(int, int, int) ;

  float insForestVert(int, int, int, int, int) ;
  float delForestVert(int, int, int, int, int) ;

  float insForestHoriz(int, int, int, int, int) ;
  float delForestHoriz(int, int, int, int, int) ;
  int nbPatterns(std::string, std::string) ;
  static const float epsilon;
  static const float match_color;
  static const float diff_color;
  static const float infiniteFloat;
  static const int infiniteInt;

  void setNbConexMatch(int, int);

  float       nbCompTree(int , int ) const;


  float       subst(const node&, const node&);
  float       ins(const node&);
  float       del(const node&);
  float       subst_complex(const node&, const node&);
  float       ins_complex(const node&);
  float       del_complex(const node&);

  Sequence *sequence;
  int _scale;
  MTG* _mtg;
  int _flag;
};

#endif



