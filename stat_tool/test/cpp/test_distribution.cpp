/****************************************************************
 *
 *  Test class Distribution
 */

#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distribution.h"

using namespace stat_tool;

int main(void) {

  int iinf_bound = 1, isup_bound = I_DEFAULT, isimul = 100;
  double iparameter = 15, iprobability = 0.5;
  DiscreteParametric *d = NULL, *d_cp = NULL;
  DiscreteParametricModel *dm = NULL, *dm_cp = NULL;
  int var, i, j, nb_variable, nb_component;
  std::vector<int> simulated(isimul);
  StatError error;
  const std::string nbinom_file = "./nbinom.dis", wrong_file = "./wrong_name";
  std::ofstream fout1(nbinom_file), fout2(wrong_file);
/*  const char * mixpath= "./tmp.mix", * spmixpath= "./tmp_sp.mix";
  const char * gnupath = "./tmp_mix", * gnu_datapath = "./tmp_mix_data";
  const char * margpath= "./marg_mix", * gnu_tmppath = "./tmp_mix_d";
  const char * np_modelpath= "./np_model.mix", * gnu_tmpnppath = "./tmp_mix_d";
  const char * entropy_data= "cluster_vectors.vec";
  const char * output_path = "tmp_entropy_output.vec";
  Distribution *marginal = NULL;
  DiscreteDistributionData *marginal_histo = NULL;*/
  Vectors *vec = NULL;


  d = new DiscreteParametric(NEGATIVE_BINOMIAL, iinf_bound, isup_bound, iparameter, iprobability);
  d->ascii_print(fout1); // will close automatically
  // d->init();
  // d->computation();
  d->ascii_print(cout);
  cout << "Simulated: ";
  for (i = 0; i < isimul; i++) {
	  simulated[i] = d->simulation();
  	  if (i < 10)
  		  cout << simulated[i] << "; ";
  }
  cout << endl;

  cout << "Copy distribution: " << endl;
  d_cp = new DiscreteParametric(*d);

  for (i = 0; i < d->nb_value; i++) {
	  assert(d->mass[i] == d_cp->mass[i]);
  }
  delete d;
  d = NULL;
  d_cp->ascii_print(cout);

  dm_cp = new DiscreteParametricModel();
  dm = dm_cp->ascii_read(error , nbinom_file);
  std::remove(nbinom_file.c_str());

  cout << "DiscreteParametricModel (from file): " << endl;
  dm->ascii_print(cout);

  delete dm;
  dm = NULL;

  dm = dm_cp->ascii_read(error , wrong_file);
  assert(dm == NULL);

  delete d_cp;
  d_cp = NULL;
  delete d;
  d = NULL;
  delete dm_cp;
  dm_cp = NULL;

  return 0;
}
