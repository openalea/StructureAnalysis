#include <math.h>
#include <sstream>
#include "tool/rw_tokenizer.h"
#include "tool/rw_cstring.h"
#include "tool/rw_locale.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/stat_label.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "discrete_parametric.h"

using namespace std;


//-----------------------------------------------------
// Constructeur de la classe Discrete_parametric
//-----------------------------------------------------
Discrete_parametric::Discrete_parametric() : Parametric()
{
  rando = NO_RANDOM;
  link = LOG;
  nb_covariable = 0;
  nb_random = 0;
  nb_cat = 0;
  random_variance = 0;
  regression = 0;
  regmult = 0;
  paramult = 0;
}

//-----------------------------------------------------
// Constructeur de la classe Discrete_parametric
//-----------------------------------------------------
Discrete_parametric::Discrete_parametric(int inb_value, int inb_covariable, int inb_random, int inb_cat,
					 double *irandom_variance, double *iregression, double **iregmult,
					 int iident, int ilink, int irando)
 :Parametric(inb_value, iident, 0, 1)
{ 
  register int i,j;
  double *pregression, *prandom;
  
  ident = iident; 
  rando = irando;
  link = ilink;
  nb_covariable = inb_covariable;
  nb_random = inb_random;
  nb_cat = inb_cat;

    
  random_variance = new double[nb_random];
  prandom = random_variance;


  if(iregression){
    regmult = 0;
    paramult = 0;
    regression = new double[nb_covariable];
    pregression = regression;
    for (i = 0; i < nb_covariable ; i++) {
      *pregression++ = *iregression++;
    }
  }

  if(iregmult){
    regression = 0;
    regmult = new double*[nb_cat];
    paramult = new double[nb_cat];
    for (i = 0; i < nb_cat; i++){
      regmult[i] = new double[nb_covariable];
      paramult[i] = 0.;
      for (j = 0; j < nb_covariable; j++) {
	regmult[i][j] = iregmult[i][j];
      }
    }
  }

  for (i = 0; i < nb_random ; i++) {
    *prandom++ = *irandom_variance++;
  }
  
}

//-------------------------------------------------------
// Constructeur de la classe Discrete_parametric
//-------------------------------------------------------
Discrete_parametric::Discrete_parametric(int inb_covariable, int inb_random, int inb_cat, 
					 double *irandom_variance, double *iregression, double **iregmult,
					 int iident, int ilink, int irando)
  :Parametric(NB_VALUE, iident, 0,1)
{
  register int i, j;
  double *pregression, *prandom ;

  ident = iident; 
  rando = irando;
  link = ilink;
  nb_covariable = inb_covariable;
  nb_random = inb_random;
  nb_cat = inb_cat;
    
  random_variance = new double[nb_random];
  prandom = random_variance;

  if(iregression){
    regmult = 0;
    paramult = 0;
    regression = new double[nb_covariable];
    pregression = regression;
    for (i = 0; i < nb_covariable ; i++) {
      *pregression++ = *iregression++;
    }
  }
  if(iregmult){
    regression = 0;
    regmult = new double*[nb_cat];
    paramult = new double[nb_cat];
    for (i = 0; i < nb_cat; i++){
      regmult[i] = new double[nb_covariable];
      paramult[i] = 0.;
      for (j = 0; j < nb_covariable; j++) {
	regmult[i][j] = iregmult[i][j];
      }
    }
  }

  for (i = 0; i < nb_random ; i++) {
    *prandom++ = *irandom_variance++;
  }
    
}

//-------------------------------------------------------
// Constructeur de la classe Discrete_parametric
//-------------------------------------------------------
Discrete_parametric::Discrete_parametric(double *covar, int inb_covariable, int inb_random, int inb_cat, double *irandom_variance, 
					 double *iregression, double **iregmult, int iident, int ilink, int irando, int inb_value)
  : Parametric(inb_value, iident,0,1)// choix de inf_bound et de sup_bound????
{
  register int i,j ;
  double *pregression, *prandom;
  double tmp;

  ident = iident; 
  rando = irando;
  link = ilink;
  nb_covariable = inb_covariable;
  nb_random = inb_random;
  nb_cat = inb_cat;
    
  random_variance = new double[nb_random];
  prandom = random_variance;

  if(nb_cat==0){
    regmult = 0;
    paramult = 0;
    regression = new double[nb_covariable];
    pregression = regression;
    for (i = 0; i < nb_covariable ; i++) {
      *pregression++ = *iregression++;
    }
  }

  if(inb_cat > 0){
    regression = 0;
    regmult = new double*[nb_cat];
    paramult = new double[nb_cat];
    for (i = 0; i < nb_cat; i++){
      regmult[i] = new double[nb_covariable];
      paramult[i] = 0.;
      for (j = 0; j < nb_covariable; j++) {
	regmult[i][j] = iregmult[i][j];
      }
    }
  }
  
  for (i = 0; i < nb_random ; i++) {
    *prandom++ = *irandom_variance++;
  }
    
  parameter = 0.;
  probability = 0.;

  switch (ident) {
  case POISSON :
    for (i = 0;i < nb_covariable; i++) {
      parameter = parameter + covar[i] * regression[i];
    }
    switch (link) {
    case LOG:
      parameter = exp(parameter);
      break;
    case IDENT:
      parameter = parameter;
      break;
    }
    break;

  case BINOMIAL:
    for (i = 0;i < nb_covariable; i++) {
      probability = probability + covar[i] * regression[i];
    }
    switch (link){
    case LOGIT:
      probability = exp(probability)/ (1+ exp(probability));
      break;
    case IDENT:
      probability = probability;
      break;
    }
    break;

  case MULTINOMIAL:
    for (i = 0;i < nb_cat; i++) {
      for (j = 0; j < nb_covariable; j++){
	paramult[i] = paramult[i] + covar[j] * iregmult[i][j];
      }
    }
    switch (link){
    case LOGIT:
      tmp = 0.;
      for (i = 0 ; i < nb_cat; i++){
	paramult[i] = exp(paramult[i]);
      }
      for (i = 0; i < (nb_cat-1); i++){
	tmp = tmp + paramult[i];
      }
      for (i = 0; i < nb_cat; i++) {
	paramult[i] = paramult[i] / (1+tmp);
      }
      break;
    }
    break;
  }
  
  //  mean = parametric_mean_computation();
  // variance = parametric_variance_computation();

}

//-------------------------------------------------------
// Constructeur de la classe Discrete_parametric
//-------------------------------------------------------
Discrete_parametric::Discrete_parametric(const Parametric &disc)
  :Parametric(disc)
{
  rando = NO_RANDOM;
  link = LOG;
  nb_covariable = 0;
  nb_random = 0;
  nb_cat = 0;
  random_variance = 0;
  regression = 0;
  regmult = 0;
  paramult = 0;
}

//-------------------------------------------------------
// Constructeur de la classe Discrete_parametric
//-------------------------------------------------------  
Discrete_parametric::Discrete_parametric(const Histogram &histo)
  :Parametric(histo)
{
  rando = NO_RANDOM;
  link = LOG;
  nb_covariable = 0;
  nb_random = 0;
  nb_cat = 0;
  random_variance = 0;
  regression = 0;
  regmult = 0;
  paramult = 0;
}

//-------------------------------------------------------
// Constructeur de la classe Discrete_parametric
//-------------------------------------------------------  
Discrete_parametric::Discrete_parametric(const Discrete_parametric &disc)
  : Parametric(disc)
{
  copy(disc);
}


//-----------------------------------------------------
// Destructeur de la classe Discrete_parametric
//-----------------------------------------------------
Discrete_parametric::~Discrete_parametric()
{
  register int i;

  if (random_variance) {
    delete [] random_variance;
  }

  if (regression) {
    delete [] regression;
  }

  if (paramult) {
    delete [] paramult;
  }

  if(regmult){
    for ( i = 0; i < nb_cat; i++){
      delete regmult[i];
    }
    delete [] regmult;
  }
}

//---------------------------------------------------------------
// Copie d'un objet Discrete_parametric
//---------------------------------------------------------------
void Discrete_parametric::copy(const Discrete_parametric &disc)
{
  register int i,j;
  double *pregression, *prandom, *cregression, *crandom, *pparamult, *cparamult;

  ident = disc.ident;
  rando = disc.rando;
  link = disc.link;
  nb_covariable = disc.nb_covariable;
  nb_random = disc.nb_random;
  nb_cat = disc.nb_cat;

  random_variance = new double[nb_random];
  prandom = random_variance;
  crandom = disc.random_variance;
  for (i = 0; i < nb_random ; i++) {
    *prandom++ = *crandom++;
  }


  if(nb_cat > 0){
    paramult = new double[nb_cat];
    pparamult = paramult;
    cparamult = disc.paramult;
    regmult = new double*[nb_cat];
    for (i = 0; i < nb_cat ; i++) {
      *pparamult++ = *cparamult++;
      regmult[i] = new double[nb_covariable];
      for (j = 0; j < nb_covariable; j++){
	regmult[i][j] = disc.regmult[i][j];
      }
    }
  }
  else {
    regression = new double[nb_covariable];
    pregression = regression;
    cregression = disc.regression;
    for (i = 0; i < nb_covariable ; i++) {
      *pregression++ = *cregression++;
    }
  }

}

//--------------------------------------------------------------
// Operateur d'assignement de la classe Discrete_parametric.
//--------------------------------------------------------------
Discrete_parametric& Discrete_parametric::operator=(const Discrete_parametric &disc)
{
  if (&disc != this) {
    delete [] mass;
    delete [] cumul;

    Parametric::copy(disc);
    copy(disc);
  }

  return *this;
}

//--------------------------------------------------------------
// Ecriture des parametres d'une loi discrète
//--------------------------------------------------------------
ostream& Discrete_parametric::ascii_print(ostream &os) const
{
  register int i,j ;

  os << "IDENT: ";
  switch (ident){
  case BINOMIAL:
    os << " BINOMIAL " << "   " << endl;
    break;  
  case POISSON:
    os << " POISSON " << "   " << endl;
    break;
  case MULTINOMIAL:
    os << " MULTINOMIAL " << "   " << endl;
    break;
  }

  os << "LINK: ";
  switch (link){
  case LOG:
    os << " LOG " << "   " << endl;
    break;  
  case LOGIT:
    os << " LOGIT " << "   " << endl;
    break;
  case PROBIT:
    os << " PROBIT " << "   " << endl;
    break;
  }


  if (!rando)
    os << "NO_RANDOM " << "   " << endl;
  else
    os << "RANDOM " << "   " << endl;
  
  if (nb_covariable != 0) {
    os << "Number of covariables: " << nb_covariable << "   ";
  }
  if (nb_random != 0) {
    os << "Number of random effects: " << nb_random << "   ";
  }
  if (nb_cat != 0) {
    os << "Number of category: " << nb_cat << "   ";
  }
  os << endl;
  
  if (nb_cat == 0) {
    for (i = 0; i < nb_covariable; i++) {
      os << "Covariable " << i << " : " << regression[i] << " ; "<<endl;
    }
    os << endl;
  }  
  
  if (nb_cat > 0) {
    for (i = 0; i < nb_cat; i++){
      for (j = 0; j < nb_covariable; j++) {
	os << "Category " << i <<", Covariable " << j << " : " << regmult[i][j] << " ; "<<endl;
      }
      os << "Category " << i <<", probability : " << paramult[i] << " ; "<<endl;
      os << endl;
    }
  }

  if (random_variance) {
    for (i = 0; i < nb_random; i++) {
      os << "Random effect variance " << i << " : " << random_variance[i] << " ; " <<endl;
    }
    os << endl;
  }   


  return os;
}

//--------------------------------------------------------------
// Ecriture des parametres d'une loi continue au format tableur.
//--------------------------------------------------------------
ostream& Discrete_parametric::spreadsheet_print(ostream &os) const
{
  register int i, j;

  os << "IDENT: ";
  switch (ident){
  case BINOMIAL:
    os << " BINOMIAL " << "   " << endl;
    break;  
  case POISSON:
    os << " POISSON " << "   " << endl;
    break;
  case MULTINOMIAL:
    os << " MULTINOMIAL " << "   " << endl;
    break;
  }

  os << "LINK: ";
  switch (link){
  case LOG:
    os << " LOG " << "   " << endl;
    break;  
  case LOGIT:
    os << " LOGIT " << "   " << endl;
    break;
  case PROBIT:
    os << " PROBIT " << "   " << endl;
    break;
  }

  if (!rando)
    os << "NO_RANDOM " << "   " << endl;
  else
    os << "RANDOM " << "   " << endl;

  if (nb_covariable != 0) {
    os << "Number of covariables " << "\t" << nb_covariable << "\t";
  }
  if (nb_random != 0) {
    os << "Number of random effects " <<"\t" << nb_random << "\t";
  }
  if (nb_cat != 0) {
    os << "Number of category: " << nb_cat << "   ";
  }
  os << endl;

  if (regression) {
    for (i = 0; i < nb_covariable; i++) {
      os << "Covariable " << i << "\t" << regression[i] << "\t" << endl;
    }
    os << endl;
  }  

  if (regmult) {
    for (i = 0; i < nb_cat; i++){
      for (j = 0; j < nb_covariable; j++) {
	os << "Category " << i <<", Covariable " << j << " : " << regmult[i][j] << " ; "<<endl;
      }
      os << "Category " << i <<", probability : " << paramult[i] << " ; "<<endl;
      os << endl;
    }
  }

  if (random_variance) {
    for (i = 0; i < nb_random; i++) {
      os << "Random effect variance " << i << "\t" << random_variance[i] << "\t" << endl;
    }
    os << endl;
  }   


  return os;
}

//--------------------------------------------------------------
// Visualisation d'une loi discrete parametrique.
//--------------------------------------------------------------
ostream& operator<<(ostream &os , const Discrete_parametric &disc)
{
  os.precision(5);

  os << endl;
  disc.ascii_print(os);
  if(disc.nb_cat == 0){
    disc.Parametric::ascii_print(os);
  }

  os.precision(6);

  return os;
}

//--------------------------------------------------------------
// Calcul du nombre de parametres d'une loi.
//--------------------------------------------------------------
int Discrete_parametric::nb_parameter_computation()
{
  int nb_parameter;

  switch (rando) {
  case NO_RANDOM :
    nb_parameter = nb_covariable + 1;
    switch (ident) {
    case MULTINOMIAL:
      nb_parameter = nb_covariable * (nb_cat-1) + 1; 
      break;
    }
    break;
  case RANDOM :
    nb_parameter = nb_covariable + nb_random + 1;
    break; 
  }

  return nb_parameter;
}

//--------------------------------------------------------------
// Calcul de la valeur des paramètres pour chaque loi
//--------------------------------------------------------------
double Discrete_parametric::parameter_computation(double *covar)const
{
  register int i;

  double param = 0.;

  switch (ident) {
  case POISSON :
    for (i = 0;i < nb_covariable; i++) {
      param = param + covar[i] * regression[i];
    }
    switch (link) {
    case LOG:
      param = exp(param);
      break;
    case IDENT:
      param = param;
      break;
    }
    break;
  case BINOMIAL:
    for (i = 0;i < nb_covariable; i++) {
      param = param + covar[i] * regression[i];
    }
    switch (link){
    case LOGIT:
      param = exp(param)/ (1.+ exp(param));
      break;
    case IDENT:
      param = param;
      break;
    }
    break;
  }

  return param;
}

//--------------------------------------------------------------
// Calcul de la valeur des paramètres pour chaque loi
//--------------------------------------------------------------
double* Discrete_parametric::parameter_computation_mult(double *covar)const
{
  register int i,j;
  double tmp;
  double *param;

  param = new double[nb_cat];

  switch (ident) {
  case MULTINOMIAL:
    for (i = 0;i < nb_cat; i++) {
      param[i] = 0.;
      for (j = 0; j < nb_covariable; j++){
	param[i] = param[i] + covar[j] * regmult[i][j];
      }
    }
    switch (link){
    case LOGIT:
      tmp = 0.;
      for (i = 0 ; i < nb_cat; i++){
	param[i] = exp(param[i]);
      }
      for (i = 0; i < (nb_cat-1); i++){
	tmp = tmp + param[i];
      }
      for (i = 0; i < nb_cat; i++) {
	param[i] = param[i] / (1+tmp);
      }
      break;
    }
    break;
  }

  return param;
}

//--------------------------------------------------------------
// Calcul de la valeur des paramètres pour chaque loi
//--------------------------------------------------------------
double* Discrete_parametric::parameter_computation_mult_var(double *covar)const
{
  register int i,j;
  double tmp;
  double *param;

  param = new double[nb_cat];

  switch (ident) {
  case MULTINOMIAL:
    for (i = 0;i < nb_cat; i++) {
      param[i] = 0.;
      for (j = 0; j < nb_covariable; j++){
	param[i] = param[i] + covar[j] * regmult[i][j];
      }
    }
    switch (link){
    case LOGIT:
      tmp = 0.;
      for (i = 0 ; i < nb_cat; i++){
	param[i] = exp(param[i]);
      }
      for (i = 0; i < (nb_cat-1); i++){
	tmp = tmp + param[i];
      }
      for (i = 0; i < nb_cat; i++) {
	param[i] = param[i]/(1.+tmp); /* / ((1.+tmp)*(1.+tmp));*/
      }
      break;
    }
    break;
  }

  return param;
}



//-----------------------------------------------------------
// Calcul de la probabilité d'une observation 
//-----------------------------------------------------------
double Discrete_parametric::mass_computation(int observ, double *covar)
{
  register int i;

  double pmass;
  double param;
  double fact1 = 1.;
  double fact2 = 1.;

  double *parama; 

 

  switch (ident) {
  case POISSON :
    param = (*this).parameter_computation(covar);
    for (i=1; i<=observ; i++){
      fact1 = fact1 * i;
      fact2 = fact2 * param;
    }
    pmass = (fact2/fact1) * exp(-param); 
    break;

  case BINOMIAL :
    param = (*this).parameter_computation(covar);
    if(observ==1){
      pmass = param;
    }
    else {
      pmass = 1.-param;
    }
    break;

  case MULTINOMIAL:
    parama = new double[nb_cat];
    parama = (*this).parameter_computation_mult(covar);
    for (i = 0; i < nb_cat; i++){
      if (observ == i){
	pmass = parama[i];
      }
    }
    break;
  }
  
  return pmass;
}



//-----------------------------------------------------------
// Tirage aléatoire selon une distribution discrete paramètrique
//-----------------------------------------------------------
int Discrete_parametric::simulation(gsl_rng *r, double *covar)
{
  int temp;
  double param;


  param = (*this).parameter_computation(covar);
    
  switch (ident) {
  case POISSON :
    temp = gsl_ran_poisson(r, param);
    break;
  case BINOMIAL :
    temp = gsl_ran_bernoulli(r, param);
    break;
  default : 
    break;
  }
  
  return temp;
}


// //------------------------------------------------------------
// // Mise à jour de la variance des effets aléatoires
// //------------------------------------------------------------
// void Continuous_parametric::Set_random_variance(const double *irandom_variance)
// { 
//   register int i;

//   double *prandom_variance;
//   random_variance = new double[nb_random];

//   prandom_variance = random_variance;

//   for (i = 0; i < nb_covariable ; i++) {
//     *prandom_variance++ = *irandom_variance++;
//   }

//   this->variance = parametric_variance_computation();
//   this->min_value_computation();
//   this->max_value_computation();
//   this->nb_value_computation();
//   //  this->density_cont_computation();
//   // this->cumul_cont_computation();
// }

//-------------------------------------------------------------
// Calcul de l'élément mu_ajt
//-------------------------------------------------------------
double Discrete_parametric::mean_computation(double *covar)
{
  double moy;

  switch (ident){
  case POISSON:
    moy = parameter_computation(covar);
    break;
  case BINOMIAL:
    moy = parameter_computation(covar);
    break;
  case MULTINOMIAL:
    break;
  }
  
  return moy;
}

//-------------------------------------------------------------
// Calcul de l'élément mu_ajt dans le cas multinomial
//-------------------------------------------------------------
double Discrete_parametric::mean_computation_multi(int observ, double *covar)
{
  register int i;
  double moy;
  double *moyi;
  moyi = new double[nb_cat];

  switch (ident){
  case MULTINOMIAL:
    moyi = parameter_computation_mult(covar);
    break;
  default:
    break;
  }
  
  for (i = 0; i < nb_cat; i++){
    if (observ == i){
      moy = moyi[i];
    }
  }

  return moy;
}

//--------------------------------------------------------------
// Calcul de l'élément W_ajt
//--------------------------------------------------------------
double Discrete_parametric::variance_computation(double *covar)
{
  double var;

  switch (ident){
  case POISSON:
    var = parameter_computation(covar);
    break;
  case BINOMIAL:
    var = parameter_computation(covar);
    var = var*(1.- var );
    break; 
  case MULTINOMIAL:
    break;
  }
  
  return var;
}

//--------------------------------------------------------------
// Calcul de l'élément W_acatjt dans le cas multinomial
//--------------------------------------------------------------
double Discrete_parametric::variance_computation_multi(int observ, double *covar)
{
  register int i;
  double var;
  double *vari;
  vari = new double[nb_cat];

  switch (ident){
  case MULTINOMIAL:
    vari = parameter_computation_mult_var(covar);
    for (i = 0; i < nb_cat; i++){
      vari[i] = vari[i]*(1.- vari[i] );
    }  
    
    break; 
  default:
    break;
  }
  
  for (i = 0; i < nb_cat; i++){
    if (observ == i){
      var = vari[i];
    }
  }

  return var;
}


//-------------------------------------------------------------
// Mise à jour des paramètres de régression fixes
//-------------------------------------------------------------
void Discrete_parametric::Set_regression(const double *iregression) 
{ 
  register int i;
  double *pregression;
  regression = new double[nb_covariable];

  pregression = regression;

  for (i = 0; i < nb_covariable ; i++) {
    *pregression++ = *iregression++;
  }

}

//----------------------------------------------------------------
// Mise à jour des paramètres de régression fixes, cas multinomial
//----------------------------------------------------------------
void Discrete_parametric::Set_regmult(double **iregmult) 
{ 
  register int i,j;

  regmult = new double*[nb_cat];
  for (i = 0; i < nb_cat; i++) {
    regmult[i] = new double[nb_covariable];
  }

  for (i = 0; i < nb_cat ; i++) {
    for (j = 0; j < nb_covariable; j++){
    regmult[i][j] = iregmult[i][j];
    }
  }

}
