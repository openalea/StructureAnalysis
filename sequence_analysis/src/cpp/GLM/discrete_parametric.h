#ifndef DISCRETE_PARAMETRIC_H
#define DISCRETE_PARAMETRIC_H

#include <sstream>
#include <gsl/gsl_rng.h>

enum {
  NO_RANDOM,
  RANDOM
};

enum {
  IDENT,
  INV,
  LOG,
  LOGIT,
  PROBIT
};


class Parametric;

class Discrete_parametric : public Parametric {

  friend std::ostream& operator<<(std::ostream &os , const Discrete_parametric &disc);
  
 public:
  int nb_covariable; // nombre de covariables
  int nb_random; // nombre d'effets aléatoires
  int nb_cat; // nombre de categorie dans le cas multinomial
  int rando; // présence ou absences d'effets aléatoires
  int link; // fonction de lien
  double *random_variance; // variances des effets aleatoires de longuer nb_random
  double *regression; // parametres de regression de longueur nb_covariable (+nb_random?)
  double **regmult; //paramètres de regression de dimension nb_covariable * nb_cat
  double *paramult; // parametre de la loi multinomiale

  Discrete_parametric();

  Discrete_parametric(int inb_value, int inb_covariable, int inb_random, int inb_cat, double *irandom_variance, 
		      double *iregression, double **iregmult, int iident = POISSON, int ilink = LOG, int irando = NO_RANDOM); 

  Discrete_parametric(int inb_covariable, int inb_random, int inb_cat, double *irandom_variance, 
		      double *iregression, double **iregmult, int iident = POISSON, int ilink = LOG, int irando = NO_RANDOM);

  Discrete_parametric(double *covar, int inb_covariable, int inb_random, int inb_cat, double *irandom_variance, 
		      double *iregression, double **iregmult, int iident = POISSON, int ilink = LOG, int irando = NO_RANDOM,
		      int inb_value = NB_VALUE);

  Discrete_parametric(const Parametric &disc);
  Discrete_parametric(const Histogram &histo);
  Discrete_parametric(const Discrete_parametric &disc); 

  ~Discrete_parametric();

  void copy(const Discrete_parametric &disc);

  Discrete_parametric& operator=(const Discrete_parametric &disc);
 
  std::ostream& ascii_print(std::ostream &os) const;
  std::ostream& spreadsheet_print(std::ostream &os) const;

  int nb_parameter_computation();
  double parameter_computation(double *covar)const;
  double* parameter_computation_mult(double *covar)const;
  double* parameter_computation_mult_var(double *covar)const;

  double mass_computation(int observ, double *covar);

  double mean_computation(double *covar);
  double mean_computation_multi(int observ, double *covar);

  double variance_computation(double *covar);
  double variance_computation_multi(int observ, double *covar);


  int simulation(gsl_rng *r, double *covar); 



  // modifications des attributs


  void Set_regression(const double *iregression); 
  void Set_regmult(double **iregmult); 

  // private:
  // protected:
};

#endif
