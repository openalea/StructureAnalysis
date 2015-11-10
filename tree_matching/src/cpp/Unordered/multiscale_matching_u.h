
#ifndef SB_MULTISCALE_MATCHING_U_HEADER
#define SB_MULTISCALE_MATCHING_U_HEADER

#include "../definitions.h" 
#include "../matchpath.h" 
#include "../choicetable.h" 
#include "../mdtable.h"
#include "../inttable.h"
#include "../treegraph.h"
#include "../sequence.h" 
#include "../nodecost.h"
#include "../wnodecost.h"
#include "../mnodecost.h"
#include "matching_unordered.h"


/**
 *\class Multiscale Matching
 *\brief Algorithm for comparing two quotiented tree graph using the complex
 *\par Presentation
 * In order to compute a distance between two multiscale tree graph, we have extended 
 * the algorithm proposed by Zhang \cite{Zha93}. The distance is computed as the minimum 
 * cost of sequence of edit operations needed totransform one MTG into another one.
 * In fact we use only two scales of decomposition. We have added three new constraints in 
 * order to respect the MTG structures.
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance), this version take only into account
 * a topological cost.
 *\author Pascal ferraro
 *\date 10/2000
 */

class MultiscaleMatching_U : public Matching_U
{
	
  public :
    /** Default constructor.*/
    MultiscaleMatching_U() {}

    /** Constructs a MatchingWithComplex using two TreeGraphs defined at the same scale and a NodeCost. */
    MultiscaleMatching_U(TreeGraph& , TreeGraph& ,NodeCost& ) ;

    /** Constructs a MatchingWithComplex using two TreeGraphs defined at the same scale and a NodeCost. */
    void make(TreeGraph& , TreeGraph& ,NodeCost& ) ;

    /** Destructor. */
    ~MultiscaleMatching_U();
    
    /** Computes the distances between trees referenced by their roots.
     *  Five different types of distance are defined (see {fer00}) */
    DistanceType distanceBetweenTree(int ,int ) ; 
    DistanceType distanceBetweenTree_v_l(int ,int ) ; 
    DistanceType distanceBetweenTree_l_w(int ,int ) ; 
    DistanceType distanceBetweenTree_v_w(int ,int ) ; 
    DistanceType distanceBetweenTree_v_w_v_l(int ,int ) ; 
    DistanceType distanceBetweenTree_v_w_l_w(int ,int ) ; 
    DistanceType distanceBetweenTree_v_w_v_w(int ,int ) ; 
    DistanceType distanceBetweenTree_v_w_(int ,int ) ; 
    DistanceType distanceBetweenTree_l_l(int ,int ) ; 
    DistanceType distanceBetweenTree_l_l_v_l(int ,int ) ; 
    DistanceType distanceBetweenTree_l_l_l_w(int ,int ) ; 
    DistanceType distanceBetweenTree_l_l_v_w(int ,int ) ; 
    DistanceType distanceBetweenTree_l_l_l_l(int ,int ) ; 

    /** Computes the distances between forests referenced by their roots. */
    DistanceType distanceBetweenForest(int ,int ) ;
    DistanceType distanceBetweenForest_v_l(int ,int ) ; 
    DistanceType distanceBetweenForest_l_w(int ,int ) ; 
    DistanceType distanceBetweenForest_v_w(int ,int ) ; 
    DistanceType distanceBetweenForest_v_w_v_l(int ,int ) ; 
    DistanceType distanceBetweenForest_v_w_l_w(int ,int ) ; 
    DistanceType distanceBetweenForest_v_w_v_w(int ,int ) ; 
    DistanceType distanceBetweenForest_v_w_(int ,int ) ; 
    DistanceType distanceBetweenForest_l_l(int ,int ) ; 
    DistanceType distanceBetweenForest_l_l_v_l(int ,int ) ; 
    DistanceType distanceBetweenForest_l_l_l_w(int ,int ) ; 
    DistanceType distanceBetweenForest_l_l_v_w(int ,int ) ; 
    DistanceType distanceBetweenForest_l_l_l_l(int ,int ) ; 
    


    /**\par Operators */

    /** Returns the distance between two trees */
    DistanceType getDBT(int ,int ) const;

    /** Returns the distance between two forests */
    DistanceType getDBF(int ,int ) const;
 
   /** Computes recursively all the distance between subtrees and sub forests
    *  of the original tree graphs  */
    DistanceType  match();

   /** Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between trees which are referenced by their roots.  */
    void getList(int ,int ,Sequence*);

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined an image for the complex of \e v and no image for the complex of \e w.*/
    void ForestList(int ,int ,Sequence& );
    void ForestList_v_l(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined  the complex of \e v is the image  of the  complex of \e w.*/
    void ForestList_v_w_v_w(int ,int ,Sequence& );
    void ForestList_v_w_v_l(int ,int ,Sequence& );
    void ForestList_v_w_l_w(int ,int ,Sequence& );
    void ForestList_v_w_(int ,int ,Sequence& );
    void ForestList_v_w(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined an image for the complex of \e w and no image for the complex of \e v.*/
    void ForestList_l_w(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined no image for the complex of \e v and n the complex of \e w.*/
    void ForestList_l_l_l_l(int ,int ,Sequence& );
    void ForestList_l_l_v_w(int ,int ,Sequence& );
    void ForestList_l_l_v_l(int ,int ,Sequence& );
    void ForestList_l_l_l_w(int ,int ,Sequence& );
    //    void ForestList_l_l(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between trees  which are referenced by their roots.*/
    void TreeList(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between trees which are referenced by their roots which
    *  defined an image for the complex of \e v and no image for the complex of \e w.*/
    void TreeList_v_l(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between trees which are referenced by their roots in which
    *   the complex of \e v and the complex of \e w are image on each other .*/
    void TreeList_v_w_v_w(int ,int ,Sequence& );
    void TreeList_v_w_v_l(int ,int ,Sequence& );
    void TreeList_v_w_l_w(int ,int ,Sequence& );
    void TreeList_v_w_(int ,int ,Sequence& );
    void TreeList_v_w(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots which
    *  defined an image for the complex of \e w and no image for the complex of \e v.*/
   void TreeList_l_w(int ,int ,Sequence& );

    /* Returns the matching list (a pointer on Sequence) of the comparison 
    *  algorithm between forests which are referenced by their roots in which
    *  the complex and \e v and \e w are no image*/
    void TreeList_l_l_v_w(int ,int ,Sequence& );
    void TreeList_l_l_v_l(int ,int ,Sequence& );
    void TreeList_l_l_l_w(int ,int ,Sequence& );
    void TreeList_l_l_l_l(int ,int ,Sequence& );
    void TreeList_l_l(int ,int ,Sequence& );

    int operator()(int );

  protected :
    TreeGraph* T1;
    TreeGraph* T2;
    MatchingDistanceTable _d;  
    MatchingDistanceTable _d_v_w_v_w;
    MatchingDistanceTable _d_v_w_v_l;
    MatchingDistanceTable _d_v_w_l_w;
    MatchingDistanceTable _d_v_w_;
    MatchingDistanceTable _d_v_w;
    MatchingDistanceTable _d_v_l;
    MatchingDistanceTable _d_l_w;
    MatchingDistanceTable _d_l_l_v_w;
    MatchingDistanceTable _d_l_l_l_w;
    MatchingDistanceTable _d_l_l_v_l;
    MatchingDistanceTable _d_l_l_l_l;
    MatchingDistanceTable _d_l_l;
    MatchingDistanceTable _restrDistances_l_l_v_w;
    MatchingDistanceTable _restrDistances_v_w_v_w;
    MatchingDistanceTable _restrDistances_v_w_;
    MatchingDistanceTable _restrDistances;
    ChoiceTable _choices;
    ChoiceTable _choices_v_w_v_w;
    ChoiceTable _choices_v_w_v_l;
    ChoiceTable _choices_v_w_l_w;
    ChoiceTable _choices_v_w_;
    ChoiceTable _choices_v_w;
    ChoiceTable _choices_v_l;
    ChoiceTable _choices_l_w;
    ChoiceTable _choices_l_l_v_w;
    ChoiceTable _choices_l_l_v_l;
    ChoiceTable _choices_l_l_l_w;
    ChoiceTable _choices_l_l_l_l;
    ChoiceTable _choices_l_l;
    NodeCost* ND;
    MatchPath _restrMapp_v_w_v_w;
    MatchPath _restrMapp_v_w_;
    MatchPath _restrMapp;
    MatchPath _restrMapp_l_l_v_w;
    VertexVector _restrMappList_v_w_v_w;
    VertexVector _restrMappList;
    VertexVector _restrMappList_v_w_;
    VertexVector _restrMappList_l_l_v_w;
    int M(int,int);
    int M_v_l(int,int);
    int M_l_w(int,int);
    int M_v_w_v_w(int,int);
    int M_v_w_v_l(int,int);
    int M_v_w_l_w(int,int);
    int M_v_w_(int,int);
    int M_v_w(int,int);
    int M_l_l_v_w(int,int);
    int M_l_l_v_l(int,int);
    int M_l_l_l_w(int,int);
    int M_l_l_l_l(int,int);
    int M_l_l(int,int);
};


#endif
