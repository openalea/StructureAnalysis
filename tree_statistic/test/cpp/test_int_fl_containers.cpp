/****************************************************************
 *
 *  Test of the containers for multidimensional
 *  mixed data as defined in int_fl_containers.h
 */

#include "tree_statistic/base_int_fl_containers.h"
#include "tree_statistic/int_fl_containers.h"

#include "tool/rw_cstring.h"

using namespace std;
using namespace Stat_trees;

typedef Base_Int_fl_container<2,1> Base_ifc21;
typedef Base_Int_fl_container<0,1> Base_ifc01;
typedef Base_Int_fl_container<1,0> Base_ifc10;
typedef Base_Int_fl_container<3,0> Base_ifc30;
typedef Base_Int_fl_container<0,2> Base_ifc02;
typedef Base_Int_fl_container<0,0> Base_ifc00;

int main(void)
{
  DataSet d;
  IntDataSet id(1);
  DoubleDataSet dd(1);
  Base_ifc21 bifc21, bifc21b;
  Base_ifc01 bifc01, bifc01b;
  Base_ifc30 bifc30, bifc30b;
  Base_ifc10 bifc10, bifc10b;
  Base_ifc02 bifc02, bifc02b;
  Base_ifc00 bifc00, bifc00b;
  // i1f1container c;
  Int_fl_container ifc, ifcb;
  One_int_container o;

  RWCString buffer , token;


  cout<< "Default type for DataSet is " << d.type << endl;
  cout<< "Type for IntDataSet is " << id.type << endl;
  cout<< "Default value for IntDataSet elements is " << id.Int(0) << endl;
  cout<< "Type for DoubleDataSet is " << dd.type << endl;
  cout<< "Default value for DoubleDataSet elements is " << dd.Double(0) << endl;
  cout<< "Here is a container for " << bifc21.nb_int() <<
         " and " << bifc21.nb_float() << " double : " << bifc21 << endl;

  bifc21.Int(1-1)= 1;
  bifc21.Int(2-1)= 2;
  bifc21.Double(1-1)= 3.01E-06;

  cout<< "Each component of this container can be modified as follows : "
  << bifc21 << endl;

  bifc21b= bifc21;
  cout << "A copy of the 2-int-1-double container (assignement) : " << bifc21b << endl;

  cout << "Here is a container for " << bifc30.nb_int() << " int ( and "
       << bifc30.nb_float() << " double) : " << bifc30 << endl;

  // bifc30.Double(2-1)= 1;  should not work
  bifc30.Int(2-1)= 1; // should work


  cout << "The 2nd component of this container can be modified as follows : "
       << bifc30 << endl;

  bifc30b= bifc30;
  cout << "A copy of the 3-int container (assignement) : " << bifc30b << endl;

  cout << "Here is a container for " << bifc02.nb_float() << " double ( and "
       << bifc02.nb_int() << " int) : " << bifc02 << endl;

  bifc02.Double(2-1)= 1.145;

  cout << "The 2nd component of this container can be modified as follows : "
       << bifc02 << endl;

  cout << "Here is a container for " << bifc10.nb_int() << " int ( and "
       << bifc10.nb_float() << " double) : " << bifc10 << endl;

  // bifc30.Double(2-1)= 1;  should not work
  bifc10.Int(1-1)= 1; // should work


  cout<< "The component of this container can be modified as follows : "
  << bifc10 << endl;

  cout << "Here is a container for " <<  bifc01.nb_float() << " double ( and "
         << bifc01.nb_int() << " int) : " << bifc01 << endl;

  bifc01.Double(1-1)= 1.145;

  cout<< "The component of this container can be modified as follows : "
  << bifc01 << endl;

  cout << "An empty container (i.e. " << bifc00.nb_int() << " int and "
         << bifc00.nb_float() << " double) : " << bifc00 << endl;

  // Copy constructor
  Base_ifc10 bifc10a(bifc10);

  cout << "A copy of the one-int-container (constructor) : " << bifc10a << endl;

  bifc10a.Int(0)= 6;

  cout << "Container set to " << bifc10a << endl;
  // Assignement operator
  bifc10b= bifc10a;
  cout << "A copy of the one-int-container (assignement) : " << bifc10b << endl;

  // Copy constructor
  Base_ifc02 bifc02a(bifc02);

  cout << "A copy of the 2-float-container (constructor) : " << bifc02a << endl;

  bifc02a.Double(0)= 6;

  cout << "Container set to " << bifc02a << endl;
  // Assignement operator
  bifc02a= bifc02a;
  cout << "A second copy (assignement) : " << bifc02b << endl;

  // Copy constructor
  Base_ifc21 bifc21a(bifc21);

  cout << "A copy of the 2-int-1-float-container (constructor) : " << bifc21a << endl;

  bifc21a.Int(0)= 12;

  cout << "Container set to " << bifc21a << endl;
  // Assignement operator
  bifc21b= bifc21a;
  cout << "A second copy (assignement) : " << bifc21b << endl;

  // Copy constructor
  Base_ifc00 bifc00a(bifc00);

  cout << "A copy of the empty container (constructor) : " << bifc00a << endl;

  cout << "Container set to " << bifc00a << endl;
  // Assignement operator
  bifc00b= bifc00a;
  cout << "A second copy (assignement) : " << bifc00b << endl;

  // Copy constructor
  Base_ifc01 bifc01a(bifc01);

  cout << "A copy of the 1-float-container (constructor) : " << bifc01a << endl;

  bifc01a.Double(0)= 2.3456;

  cout << "Container set to " << bifc01a << endl;
  // Assignement operator
  bifc01b= bifc01a;
  cout << "A second copy (assignement) : " << bifc01b << endl;

  // Copy constructor
  Base_ifc30 bifc30a(bifc30);

  cout << "A copy of the 3-int-container (constructor) : " << bifc30a << endl;

  bifc30a.Int(2)= 8;

  cout << "Container set to " << bifc30a << endl;
  // Assignement operator
  bifc30b= bifc30a;
  cout << "A second copy (assignement) : " << bifc30b << endl;


  bifc10a= 3;

  // var_int= bifc10a;

  cout << "Container set to " << bifc10a << endl;

  bifc01b= bifc01;
  cout << "A copy of the one-double-container (assignement) : " << bifc01b << endl;

  bifc01= 3.1416;

  // var_int= bifc10a;

  cout << "Container set to " << bifc01 << endl;

  // test of the Int_fl_container class
  cout << endl;

  cout << "Test of the Int_fl_container class " << endl;

  cout<< "Here is a container for " << ifc.nb_int() <<
         " int and " << ifc.nb_float() << " double (default container) : " << ifc << endl;

  ifc.reset(1,1);
  ifc.Int(1-1)= 1;
  ifc.Double(1-1)= 3.01E-06;
  // ifc.Double(2-1)= 2; should not work

  cout<< "Each component of this container can be modified as follows : "
  << ifc << endl;

  ifcb= ifc;
  cout << "A copy of the 1-int-1-double container (assignement) : " << ifcb << endl;

  // Copy constructor
  Int_fl_container ifca(ifc);

  cout << "A copy of the 1-int-1-double container (constructor) : " << ifca << endl;

  ifca.Int(0)= 6;

  cout << "Container set to " << ifca << endl;
  // Assignement operator
  ifcb= ifca;
  cout << "A copy of the 1-int-1-double container (assignement) : " << ifcb << endl;



  ifc.reset(3,0);

  cout<< "Here is a container for " << ifc.nb_int() <<
  " int ( and " << ifc.nb_float() << " double) : " << ifc << endl;

  // bifc30.Double(2-1)= 1;  should not work
  ifc.Int(2-1)= 1; // should work


  cout<< "The 2nd component of this container can be modified as follows : "
  << ifc << endl;

  ifcb= ifc;
  cout << "A copy of the 3-int container (assignement) : " << ifcb << endl;

  ifc.reset(0,2);

  cout<< "Here is a container for " << ifc.nb_float() << " double ( and "
         << ifc.nb_int() << " int) : " << ifc << endl;

  ifc.Double(2-1)= 1.145;

  cout<< "The 2nd component of this container can be modified as follows : "
  << ifc << endl;

  ifc.reset(1,0);

  cout<< "Here is a container for " << ifc.nb_int() << " int ( and "
         << ifc.nb_float() << " double) : " << ifc << endl;

  // bifc30.Double(2-1)= 1;  should not work
  ifc.Int(1-1)= 1; // should work


  cout<< "The component of this container can be modified as follows : "
  << ifc << endl;

  ifc.reset(0,1);

  cout << "Here is a container for " << ifc.nb_float() << " double ( and "
         << ifc.nb_int() << " int) : " << ifc << endl;

  ifc.Double(1-1)= 1.145;

  cout<< "The component of this container can be modified as follows : "
  << ifc << endl;

  ifc.reset(0,0);

  cout << "An empty container (i.e. " << ifc.nb_int() << " int and "
         << ifc.nb_float() << " double) : " << ifc << endl;

  // test of the One_int_container class
  cout << endl;
  cout << "Test of the One_int_container class " << endl;


  cout << "A one-int-container : " << o << endl;
  o.Int()= 1;

  cout<< "Its component can be modified as follows : " << o << endl;

  // Copy constructor
  Int_fl_container oa(o), ob;

  ob= o;
  cout << "A copy of this container (assignement) : " << ob << endl;

  cout << "A copy of the 1-int-1-double container (constructor) : " << oa << endl;

  oa.Int(0)= 5;

  cout << "Container set to " << oa << endl;

  // oa.Int(1); should not work

  return 0;

}
