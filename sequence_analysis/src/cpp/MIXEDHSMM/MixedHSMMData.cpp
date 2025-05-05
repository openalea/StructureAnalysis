#include "MixedHSMMData.h"



// Constructeur nul par défaul (0 séquence produite)
MixedHSMMData::MixedHSMMData ()
{
    SemiMarkovData::SemiMarkovData();

    covariate_are_dynamic = false;
    covariate = NULL;
}


// Constructeur à partir d'un objet SemiMarkovData pour les données observées, et d'un objet Sequences pour les covariables
MixedHSMMData::MixedHSMMData (bool _covariate_are_dynamic, SemiMarkovData* observed_data, Sequences* covariate_data)
{
    SemiMarkovData::SemiMarkovData();

    covariate_are_dynamic = _covariate_are_dynamic;
    covariate = covariate_data;
}


































/*--------------------------------------------------------------*/
/**
 *  \brief Construction of a SemiMarkovData object from a MarkovianSequences object.
 *
 *  \param[in] seq              reference on a MarkovianSequences object,
 *  \param[in] transform        type of transform (SEQUENCE_COPY/ADD_STATE_VARIABLE),
 *  \param[in] initial_run_flag addition/removing of the initial run length frequency distributions.
 */
/*--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData(const MarkovianSequences &seq , sequence_transformation transform ,
    bool initial_run_flag)
:MarkovianSequences(seq , transform , (initial_run_flag ? ADD_INITIAL_RUN : REMOVE_INITIAL_RUN))

{
semi_markov = NULL;
chain_data = NULL;

likelihood = D_INF;
restoration_likelihood = D_INF;
sample_entropy = D_DEFAULT;

posterior_probability = NULL;
posterior_state_probability = NULL;
entropy = NULL;
nb_state_sequence = NULL;
}




































////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the SemiMarkovData class.
 *
 *  \param[in] ilength_distribution sequence length frequency distribution,
 *  \param[in] inb_variable         number of variables,
 *  \param[in] itype                variable types,
 *  \param[in] init_flag            flag initialization.
 */
/*--------------------------------------------------------------*/

SemiMarkovData::SemiMarkovData(const FrequencyDistribution &ilength_distribution , int inb_variable ,
    variable_nature *itype , bool init_flag)
:MarkovianSequences(ilength_distribution , inb_variable , itype , init_flag)

{
semi_markov = NULL;
chain_data = NULL;

likelihood = D_INF;
restoration_likelihood = D_INF;
sample_entropy = D_DEFAULT;

posterior_probability = NULL;
posterior_state_probability = NULL;
entropy = NULL;
nb_state_sequence = NULL;
}


/*--------------------------------------------------------------*/
/**
*  \brief Construction of a MarkovianSequences exactly as a Sequences object
*
*  \param[in] ilength_distribution sequence length frequency distribution,
*  \param[in] inb_variable         number of variables,
*  \param[in] itype                variable types,
*  \param[in] init_flag            flag initialization.
*/
/*--------------------------------------------------------------*/

MarkovianSequences::MarkovianSequences(const stat_tool::FrequencyDistribution &ilength_distribution , int inb_variable ,
                         stat_tool::variable_nature *itype , bool init_flag)
: min_interval(NULL),
self_transition(NULL),
observation_distribution(NULL),
observation_histogram(NULL),
characteristics(NULL),
Sequences(ilength_distribution , inb_variable , itype , init_flag)
{
init();
}



/*--------------------------------------------------------------*/
/**
 *  \brief Constructor of the Sequences class.
 *
 *  \param[in] ilength_distribution sequence length frequency distribution,
 *  \param[in] inb_variable         number of variables,
 *  \param[in] itype                variable types,
 *  \param[in] init_flag            flag initialization.
 */
/*--------------------------------------------------------------*/

Sequences::Sequences(const FrequencyDistribution &ilength_distribution ,
    int inb_variable , variable_nature *itype , bool init_flag)

{
int i , j , k;
int *plength;


nb_sequence = ilength_distribution.nb_element;

identifier = new int[nb_sequence];
for (i = 0;i < nb_sequence;i++) {
identifier[i] = i + 1;
}

length = new int[nb_sequence];
plength = length;
for (i = ilength_distribution.offset;i < ilength_distribution.nb_value;i++) {
for (j = 0;j < ilength_distribution.frequency[i];j++) {
*plength++ = i;
}
}

max_length = ilength_distribution.nb_value - 1;
cumul_length_computation();
length_distribution = new FrequencyDistribution(ilength_distribution);

vertex_identifier = NULL;

index_param_type = IMPLICIT_TYPE;
index_parameter_distribution = NULL;
index_interval = NULL;
index_parameter = NULL;

nb_variable = inb_variable;

type = new variable_nature[nb_variable];
min_value = new double[nb_variable];
max_value = new double[nb_variable];
marginal_distribution = new FrequencyDistribution*[nb_variable];
marginal_histogram = new Histogram*[nb_variable];

if (itype) {
for (i = 0;i < nb_variable;i++) {
type[i] = itype[i];
}
}

else {
type[0] = STATE;
for (i = 1;i < nb_variable;i++) {
type[i] = INT_VALUE;
}
}

for (i = 0;i < nb_variable;i++) {
min_value[i] = 0.;
max_value[i] = 0.;
marginal_distribution[i] = NULL;
marginal_histogram[i] = NULL;
}

int_sequence = new int**[nb_sequence];
real_sequence = new double**[nb_sequence];
for (i = 0;i < nb_sequence;i++) {
int_sequence[i] = new int*[nb_variable];
real_sequence[i] = new double*[nb_variable];
for (j = 0;j < nb_variable;j++) {
if (type[j] != REAL_VALUE) {
int_sequence[i][j] = new int[length[i]];
real_sequence[i][j] = NULL;

if (init_flag) {
for (k = 0;k < length[i];k++) {
int_sequence[i][j][k] = 0;
}
}
}

else {
int_sequence[i][j] = NULL;
real_sequence[i][j] = new double[length[i]];

if (init_flag) {
for (k = 0;k < length[i];k++) {
real_sequence[i][j][k] = 0.;
}
}
}
}
}
}


