/****************************************************************
 *
 *  Test class Distribution
 */

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distribution.h"

using namespace stat_tool;

int main(void) {

  int iinf_bound = 1, isup_bound = I_DEFAULT, isimul = 100;
  double iparameter = 15, iprobability = 0.5;
  bool all;
  DiscreteParametric *d = NULL, *d_cp = NULL;
  DiscreteParametricModel *dm = NULL, *dm_cp = NULL;
  int var, i, j, nb_variable, nb_component;
  std::vector<int> simulated(isimul), simulated2(isimul);
  StatError error;


  d = new DiscreteParametric(NEGATIVE_BINOMIAL, iinf_bound, isup_bound, iparameter, iprobability);
  // d->init();
  // d->computation();
  d->ascii_print(cout);
  set_seed(0);
  cout << "Simulated: ";
  for (i = 0; i < isimul; i++) {
	  simulated[i] = d->simulation();
  	  if (i < 10)
  		  cout << simulated[i] << "; ";
  }
  all = true;
  set_seed(1);
  for (i = 0; i < isimul; i++) {
	  simulated2[i] = d->simulation();
	  if (simulated2[i] != simulated[i])
		  all = false;
  }
  if (all)
	  cout << "Both samples are equal" << endl;

  all = true;
  set_seed(0);
  for (i = 0; i < isimul; i++) {
	  simulated2[i] = d->simulation();
	  if (simulated2[i] != simulated[i])
		  all = false;
  }
  if (!all)
	  cout << "Both samples are different" << endl;

  delete d;
  d = NULL;

  return 0;
}
