/****************************************************************
 *
 *  Test multivariate mixture models
 */

#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/multivariate_mixture.h"

using namespace stat_tool;

int main(void) {

  int var, i, j, nb_variable, nb_component;
  const int nb_vector = 500;
  double *pweight = NULL;
  bool *fparam = NULL;
  int* perm;
  StatError error;
  const char * mixpath= "./tmp.mix", * spmixpath= "./tmp_sp.mix";
  const char * gnupath = "./tmp_mix", * gnu_datapath = "./tmp_mix_data";
  const char * margpath= "./marg_mix", * gnu_tmppath = "./tmp_mix_d";
  const char * np_modelpath= "./np_model.mix", * gnu_tmpnppath = "./tmp_mix_d";
  const char * entropy_data= "cluster_vectors.vec";
  const char * output_path = "tmp_entropy_output.vec";
  Distribution *marginal = NULL;
  DiscreteDistributionData *marginal_histo = NULL;
  Vectors *vec = NULL;
  MultivariateMixture *mv1 = NULL, *mv_cp = NULL;
  MultivariateMixture *mv_np1 = NULL, *mv_np_estim = NULL;
  MultivariateMixture *mv_estim = NULL;
  MultivariateMixtureData *mv_data = NULL, *cluster = NULL;
  DiscreteParametric **dt1 = NULL, **dt2 = NULL;
  DiscreteParametricProcess **ppcomponent = NULL;


  // constructors of Mv_Mixture
  mv1 = new MultivariateMixture();

  // destructor of MultivariateMixture
  delete mv1;
  mv1= NULL;

  nb_variable = 2;
  nb_component = 3;

  dt1 = new DiscreteParametric*[nb_component];
  dt2 = new DiscreteParametric*[nb_component];
  pweight = new double[nb_component];
  ppcomponent = new DiscreteParametricProcess*[nb_variable];

  pweight[0] = 0.1;
  pweight[1] = 0.2;
  pweight[2] = 0.7;



  dt1[0] = new DiscreteParametric(0, BINOMIAL, 2, 12, D_DEFAULT, 0.1);
  dt1[1] = new DiscreteParametric(0, BINOMIAL, 0, 10, D_DEFAULT, 0.5);
  dt1[2] = new DiscreteParametric(0, BINOMIAL, 3, 10, D_DEFAULT, 0.8);

  dt2[0] = new DiscreteParametric(0, POISSON, 2, I_DEFAULT, 8.0, D_DEFAULT);
  dt2[1] = new DiscreteParametric(0, POISSON, 4, I_DEFAULT, 5.0, D_DEFAULT);
  dt2[2] = new DiscreteParametric(0, POISSON, 0, I_DEFAULT, 2.0, D_DEFAULT);

  cout << "Observation distributions for variable 1:" << endl;
  for (i = 0; i < nb_component; i++) {
    dt1[i]-> ascii_print(cout);
  }

  ppcomponent[0] = new DiscreteParametricProcess(nb_component, dt1);
  ppcomponent[1] = new DiscreteParametricProcess(nb_component, dt2);

  for (i = 0; i < nb_component; i++) {
    delete dt1[i];
    dt1[i] = NULL;
    delete dt2[i];
    dt2[i] = NULL;
  }

  cout << endl;

  mv1 = new MultivariateMixture(nb_component, pweight, nb_variable, ppcomponent, NULL);

  cout << "Mixture of " << nb_component << " components with " <<
    nb_variable << " variables:" << endl;

  mv1->ascii_write(cout, true);
  cout << endl;

  // copy
  mv_cp = new MultivariateMixture(*mv1);
  cout << "Copy constructor of MultivariateMixture: " << endl;
  mv_cp->ascii_write(cout);
  cout << endl;

  // destructor of MultivariateMixture
  delete mv_cp;
  mv_cp= NULL;

  delete mv1;
  mv1= NULL;

  cout << "MultivariateMixture_building (print into file " << mixpath << "): " << endl;

  mv1 = multivariate_mixture_building(error , nb_component ,
                                      nb_variable, pweight, ppcomponent, NULL);

  if (mv1 == NULL)
    cout << error;
  else {
    mv1->ascii_write(error, mixpath, true);
    delete mv1;
    mv1= NULL;
  }
  cout << endl;

  cout << "Read MultivariateMixture from file " << mixpath << ": " << endl;

  mv1 = multivariate_mixture_ascii_read(error , mixpath);

  if (mv1 == NULL) {
    cout << error;
    return 1;
  }
  else {
    mv1->ascii_write(cout, true);
  }
  cout << endl;

  mv_data = mv1->simulation(error, nb_vector);
  if (mv_data == NULL) {
    cout << error;
    return 1;
  }
  else {
    mv_data->ascii_write(cout, true);
  }
  cout << endl;

  cout << "Gnuplot output for MultivariateMixture ('" << gnupath << "' file)" << endl;
  mv1->plot_write(error, gnupath, "");
  cout << error << endl;

  cout << "Gnuplot output for MultivariateMixture_data ('" << gnu_datapath << "' file)" << endl;
  mv_data->plot_write(error, gnu_datapath, "");
  cout << error << endl;

  cout << "Extract marginal distribution for variable 1" << endl;
  marginal = mv1->extract_distribution(error, 1);
  marginal_histo = mv_data->extract_marginal(error, 1);

  if (marginal != NULL) {
    marginal->ascii_print(cout, false, false, false);
    delete marginal;
  }
  else
    cout << error;

  if (marginal_histo != NULL) {
    cout << "Gnuplot output for marginal ('" << margpath << "' file)" << endl;
    marginal_histo->plot_write(error, margpath);
    delete marginal_histo;
  }
  else
    cout << error;

  cout << "Estimate MultivariateMixture from initial model: " << endl;
  mv_estim = mv_data->mixture_estimation(error, cout, *mv1);

  if (mv_estim == NULL) {
    cout << error;
  }
  else {
    mv_estim->ascii_write(cout, true);
    mv_estim->plot_write(error, gnu_tmppath, "");
    mv_estim->spreadsheet_write(error, spmixpath);
    cout << error << endl;
    delete mv_estim;
    mv_estim = NULL;
  }
  cout << endl;

  fparam = new bool[2];
  fparam[0] = true;
  fparam[1] = false;

  cout << "Estimate MultivariateMixture from initial nb_component: " << endl;
  mv_estim = mv_data->mixture_estimation(error, cout, 3, I_DEFAULT, fparam);

  delete [] fparam;

  if (mv_estim == NULL) {
    cout << error;
  }
  else {
    mv_estim->ascii_write(cout, true);
    mv_estim->plot_write(error, gnu_tmppath, "");
    delete mv_estim;
    mv_estim = NULL;
  }
  cout << endl;

  delete mv1;
  mv1= NULL;
  delete mv_data;
  mv_data = NULL;

  for (var = 0; var < nb_variable; var++) {
    delete ppcomponent[var];
    ppcomponent[var] = NULL;
  }

  cout << "Read non parametric model: " << endl;

  mv_np1 = multivariate_mixture_ascii_read(error, np_modelpath);

  if (mv_np1 == NULL) {
    cout << error;
    return 1;
  }
  else {
    cout << "Value: " << endl;
    mv_np1->ascii_write(cout, false);
    mv_np1->spreadsheet_write(error, spmixpath);
    cout << error << endl;
  }


  // permutation of the states
  perm = new int[mv_np1->get_nb_component()];
  for(j= 0; j < mv_np1->get_nb_component(); j++)
     perm[mv_np1->get_nb_component()-j -1]=j;
  mv_np1->state_permutation(error, perm);
  cout << error;
  cout << "Permutation of the states: " << endl;
  mv_np1->ascii_write(cout, false);
  delete [] perm;
  perm = NULL;

  cout << endl;

  mv_data = mv_np1->simulation(error, nb_vector);

  if (mv_data == NULL) {
    cout << error;
    return 1;
  }

  cout << "Extract marginal distribution for variable 2" << endl;
  marginal = mv_np1->extract_distribution(error, 2);

  if (marginal != NULL) {
    marginal->ascii_print(cout, false, false, false);
    delete marginal;
  }
  else
    cout << error;

  cout << "Gnuplot output for MultivariateMixture_data ('" << gnu_tmpnppath << "' file)" << endl;
  mv_data->plot_write(error, gnu_tmpnppath, "");
  cout << error << endl;


  cout << "Estimate MultivariateMixture from initial nb_component: " << endl;
  mv_estim = mv_data->mixture_estimation(error, cout, 3);

  if (mv_estim == NULL) {
    cout << error;
    return 1;
  }
  else {
    mv_estim->ascii_write(cout, true);
    mv_estim->plot_write(error, gnu_tmpnppath, "");
  }

  // permutation of the states
  perm = new int[mv_estim->get_nb_component()];
  for(j= 0; j < mv_estim->get_nb_component(); j++)
     perm[mv_estim->get_nb_component()-j -1]=j;
  mv_estim->state_permutation(error, perm);
  cout << error;
  cout << "Permutation of the states: " << endl;
  mv_estim->ascii_write(cout, false);
  delete [] perm;
  perm = NULL;

  cout << endl;

  cout << "Cluster individuals: " << endl;
  cluster = mv_estim->cluster(error, *mv_data);

  if (cluster == NULL) {
    cout << error;
  }
  else {
    cout << "Estimated states: " << endl;
    for (i = 0; i < 10; i++)
       cout << cluster->get_int_vector(i, 0) << "\t";
    cout << "..." << endl;
    delete cluster;
    cluster = NULL;
  }

  cout << "State entropy: " << endl;
  cluster = mv_estim->cluster(error, *mv_data, VITERBI, true);

  if (cluster == NULL) {
    cout << error;
  }
  else {
    cout << "State entropies: " << endl;
    for (i = 0; i < 10; i++)
       cout << cluster->get_real_vector(i, 3) << "\t";
    cout << "..." << endl;
    delete mv_estim;
    mv_estim = NULL;
    delete cluster;
    cluster = NULL;
  }

  // compute state entropies on data
  cout << "Compute state entropies on data" << endl;

  vec = vectors_ascii_read(error, entropy_data);
  if (vec == NULL) {
    cout << error;
    return 1;
  }

  mv_estim = vec->mixture_estimation(error, cout, 3, 300);

  if (mv_estim == NULL) {
    cout << error;
    return 1;
  }
  else
    mv_estim->ascii_write(cout, true);

  cluster = mv_estim->cluster(error, *vec, VITERBI, true);

  if (cluster == NULL) {
    cout << error;
  }
  else {
    cluster->ascii_write(error, output_path);
  }

  delete mv_estim;
  delete cluster;

  delete vec;
  delete mv_data;
  delete mv_np1;

  delete [] dt1;
  delete [] dt2;
  delete [] pweight;
  delete [] ppcomponent;

  return 0;
}
