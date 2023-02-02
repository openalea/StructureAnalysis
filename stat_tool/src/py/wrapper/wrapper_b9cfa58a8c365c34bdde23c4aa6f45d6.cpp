#include "_stat_tool.h"


void wrapper_b9cfa58a8c365c34bdde23c4aa6f45d6()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::stat_tool_label > enum_b9cfa58a8c365c34bdde23c4aa6f45d6("stat_tool_label");

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_LINE", ::stat_tool::STATL_LINE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_WORD", ::stat_tool::STATL_WORD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_HIT_RETURN", ::stat_tool::STATL_HIT_RETURN);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_END", ::stat_tool::STATL_END);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MODEL", ::stat_tool::STATL_MODEL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MEAN", ::stat_tool::STATL_MEAN);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MEDIAN", ::stat_tool::STATL_MEDIAN);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MODE", ::stat_tool::STATL_MODE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VARIANCE", ::stat_tool::STATL_VARIANCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_STANDARD_DEVIATION", ::stat_tool::STATL_STANDARD_DEVIATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_LOWER_QUARTILE", ::stat_tool::STATL_LOWER_QUARTILE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_UPPER_QUARTILE", ::stat_tool::STATL_UPPER_QUARTILE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MEAN_ABSOLUTE_DEVIATION", ::stat_tool::STATL_MEAN_ABSOLUTE_DEVIATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VARIATION_COEFF", ::stat_tool::STATL_VARIATION_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VARIANCE_MEAN_RATIO", ::stat_tool::STATL_VARIANCE_MEAN_RATIO);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SKEWNESS_COEFF", ::stat_tool::STATL_SKEWNESS_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_KURTOSIS_COEFF", ::stat_tool::STATL_KURTOSIS_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_INFORMATION", ::stat_tool::STATL_INFORMATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CONCENTRATION_COEFF", ::stat_tool::STATL_CONCENTRATION_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SMOOTHNESS", ::stat_tool::STATL_SMOOTHNESS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MEAN_DIRECTION", ::stat_tool::STATL_MEAN_DIRECTION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MEAN_RESULTANT_LENGTH", ::stat_tool::STATL_MEAN_RESULTANT_LENGTH);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CIRCULAR_STANDARD_DEVIATION", ::stat_tool::STATL_CIRCULAR_STANDARD_DEVIATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MEAN_CONFIDENCE_INTERVAL", ::stat_tool::STATL_MEAN_CONFIDENCE_INTERVAL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_ONE_SIDED", ::stat_tool::STATL_ONE_SIDED);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_TWO_SIDED", ::stat_tool::STATL_TWO_SIDED);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CHI2_TEST", ::stat_tool::STATL_CHI2_TEST);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_F_TEST", ::stat_tool::STATL_F_TEST);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_T_TEST", ::stat_tool::STATL_T_TEST);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_WILCOXON_MANN_WHITNEY_TEST", ::stat_tool::STATL_WILCOXON_MANN_WHITNEY_TEST);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_KRUSKAL_WALLIS_TEST", ::stat_tool::STATL_KRUSKAL_WALLIS_TEST);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FREEDOM_DEGREE", ::stat_tool::STATL_FREEDOM_DEGREE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FREEDOM_DEGREES", ::stat_tool::STATL_FREEDOM_DEGREES);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_STANDARD_NORMAL_VALUE", ::stat_tool::STATL_STANDARD_NORMAL_VALUE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CHI2_VALUE", ::stat_tool::STATL_CHI2_VALUE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_F_VALUE", ::stat_tool::STATL_F_VALUE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_T_VALUE", ::stat_tool::STATL_T_VALUE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CRITICAL_PROBABILITY", ::stat_tool::STATL_CRITICAL_PROBABILITY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_REFERENCE", ::stat_tool::STATL_REFERENCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_TEST", ::stat_tool::STATL_TEST);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MANN_WHITNEY_INFERIOR_PROBABILITY", ::stat_tool::STATL_MANN_WHITNEY_INFERIOR_PROBABILITY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MANN_WHITNEY_EQUAL_PROBABILITY", ::stat_tool::STATL_MANN_WHITNEY_EQUAL_PROBABILITY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MANN_WHITNEY_SUPERIOR_PROBABILITY", ::stat_tool::STATL_MANN_WHITNEY_SUPERIOR_PROBABILITY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FIT", ::stat_tool::STATL_FIT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_LIKELIHOOD", ::stat_tool::STATL_LIKELIHOOD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_NORMALIZED", ::stat_tool::STATL_NORMALIZED);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MAX_LIKELIHOOD", ::stat_tool::STATL_MAX_LIKELIHOOD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DEVIANCE", ::stat_tool::STATL_DEVIANCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CLASSIFICATION_LIKELIHOOD", ::stat_tool::STATL_CLASSIFICATION_LIKELIHOOD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MAX_CLASSIFICATION_LIKELIHOOD", ::stat_tool::STATL_MAX_CLASSIFICATION_LIKELIHOOD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CLASSIFICATION_ENTROPY", ::stat_tool::STATL_CLASSIFICATION_ENTROPY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FREE_PARAMETER", ::stat_tool::STATL_FREE_PARAMETER);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FREE_PARAMETERS", ::stat_tool::STATL_FREE_PARAMETERS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_PENALIZED_LIKELIHOOD", ::stat_tool::STATL_PENALIZED_LIKELIHOOD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_ITERATION", ::stat_tool::STATL_ITERATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_ITERATIONS", ::stat_tool::STATL_ITERATIONS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DISTRIBUTION", ::stat_tool::STATL_DISTRIBUTION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DISTRIBUTIONS", ::stat_tool::STATL_DISTRIBUTIONS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_UNPROPER", ::stat_tool::STATL_UNPROPER);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_COMPLEMENTARY_PROBABILITY", ::stat_tool::STATL_COMPLEMENTARY_PROBABILITY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CUMULATIVE", ::stat_tool::STATL_CUMULATIVE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SURVIVOR", ::stat_tool::STATL_SURVIVOR);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FUNCTION", ::stat_tool::STATL_FUNCTION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MATCHING", ::stat_tool::STATL_MATCHING);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CONCENTRATION", ::stat_tool::STATL_CONCENTRATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CURVE", ::stat_tool::STATL_CURVE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MIXTURE", ::stat_tool::STATL_MIXTURE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_WEIGHT", ::stat_tool::STATL_WEIGHT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_COMPONENT", ::stat_tool::STATL_COMPONENT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CATEGORICAL", ::stat_tool::STATL_CATEGORICAL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DISCRETE_PARAMETRIC", ::stat_tool::STATL_DISCRETE_PARAMETRIC);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_COMPONENT_DISTANCE", ::stat_tool::STATL_COMPONENT_DISTANCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_POSTERIOR_ASSIGNMENT_PROBABILITY", ::stat_tool::STATL_POSTERIOR_ASSIGNMENT_PROBABILITY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_ASSIGNMENT_ENTROPY", ::stat_tool::STATL_ASSIGNMENT_ENTROPY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CONVOLUTION", ::stat_tool::STATL_CONVOLUTION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_COMPOUND", ::stat_tool::STATL_COMPOUND);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SUM", ::stat_tool::STATL_SUM);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_ELEMENTARY", ::stat_tool::STATL_ELEMENTARY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DEATH_PROBABILITY", ::stat_tool::STATL_DEATH_PROBABILITY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SURVIVAL_PROBABILITY", ::stat_tool::STATL_SURVIVAL_PROBABILITY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FREQUENCY", ::stat_tool::STATL_FREQUENCY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_BACKWARD", ::stat_tool::STATL_BACKWARD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FORWARD", ::stat_tool::STATL_FORWARD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_OBSERVATION_INTER_EVENT", ::stat_tool::STATL_OBSERVATION_INTER_EVENT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SOJOURN_TIME", ::stat_tool::STATL_SOJOURN_TIME);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FREQUENCY_DISTRIBUTION", ::stat_tool::STATL_FREQUENCY_DISTRIBUTION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_FREQUENCY_DISTRIBUTIONS", ::stat_tool::STATL_FREQUENCY_DISTRIBUTIONS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_HISTOGRAM", ::stat_tool::STATL_HISTOGRAM);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_BIN_WIDTH", ::stat_tool::STATL_BIN_WIDTH);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SAMPLE", ::stat_tool::STATL_SAMPLE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SAMPLE_SIZE", ::stat_tool::STATL_SAMPLE_SIZE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DISSIMILARITIES", ::stat_tool::STATL_DISSIMILARITIES);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_INFORMATION_RATIO", ::stat_tool::STATL_INFORMATION_RATIO);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CLUSTERING_STEP", ::stat_tool::STATL_CLUSTERING_STEP);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CATEGORY", ::stat_tool::STATL_CATEGORY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VECTOR", ::stat_tool::STATL_VECTOR);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VECTORS", ::stat_tool::STATL_VECTORS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_NB_VECTOR", ::stat_tool::STATL_NB_VECTOR);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VARIABLE", ::stat_tool::STATL_VARIABLE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_IDENTIFIER", ::stat_tool::STATL_IDENTIFIER);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MIN_VALUE", ::stat_tool::STATL_MIN_VALUE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MAX_VALUE", ::stat_tool::STATL_MAX_VALUE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MARGINAL", ::stat_tool::STATL_MARGINAL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VARIANCE_COVARIANCE_MATRIX", ::stat_tool::STATL_VARIANCE_COVARIANCE_MATRIX);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CORRELATION_MATRIX", ::stat_tool::STATL_CORRELATION_MATRIX);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SPEARMAN_RANK_CORRELATION_MATRIX", ::stat_tool::STATL_SPEARMAN_RANK_CORRELATION_MATRIX);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_KENDALL_RANK_CORRELATION_MATRIX", ::stat_tool::STATL_KENDALL_RANK_CORRELATION_MATRIX);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_LIMIT_CORRELATION_COEFF", ::stat_tool::STATL_LIMIT_CORRELATION_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SPEARMAN_LIMIT_RANK_CORRELATION_COEFF", ::stat_tool::STATL_SPEARMAN_LIMIT_RANK_CORRELATION_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_KENDALL_LIMIT_RANK_CORRELATION_COEFF", ::stat_tool::STATL_KENDALL_LIMIT_RANK_CORRELATION_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SUP_NORM_DISTANCE", ::stat_tool::STATL_SUP_NORM_DISTANCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_OVERLAP", ::stat_tool::STATL_OVERLAP);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CONTINGENCY_TABLE", ::stat_tool::STATL_CONTINGENCY_TABLE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DEVIATION_TABLE", ::stat_tool::STATL_DEVIATION_TABLE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CHI2_CONTRBUTION_TABLE", ::stat_tool::STATL_CHI2_CONTRBUTION_TABLE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VARIANCE_ANALYSIS", ::stat_tool::STATL_VARIANCE_ANALYSIS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VARIATION_SOURCE", ::stat_tool::STATL_VARIATION_SOURCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_BETWEEN_SAMPLES", ::stat_tool::STATL_BETWEEN_SAMPLES);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_WITHIN_SAMPLES", ::stat_tool::STATL_WITHIN_SAMPLES);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_TOTAL", ::stat_tool::STATL_TOTAL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SQUARE_SUM", ::stat_tool::STATL_SQUARE_SUM);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MEAN_SQUARE", ::stat_tool::STATL_MEAN_SQUARE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_INTERCEPT", ::stat_tool::STATL_INTERCEPT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SLOPE", ::stat_tool::STATL_SLOPE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_EXPLANATORY_VARIABLE", ::stat_tool::STATL_EXPLANATORY_VARIABLE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_RESPONSE_VARIABLE", ::stat_tool::STATL_RESPONSE_VARIABLE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CORRELATION_COEFF", ::stat_tool::STATL_CORRELATION_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_R_SQUARED", ::stat_tool::STATL_R_SQUARED);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DETERMINATION_COEFF", ::stat_tool::STATL_DETERMINATION_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_NULL_AUTOREGRESSIVE_COEFF_95_CONFIDENCE_LIMIT", ::stat_tool::STATL_NULL_AUTOREGRESSIVE_COEFF_95_CONFIDENCE_LIMIT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_RESIDUAL", ::stat_tool::STATL_RESIDUAL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_STANDARDIZED_RESIDUAL", ::stat_tool::STATL_STANDARDIZED_RESIDUAL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_LINEAR", ::stat_tool::STATL_LINEAR);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_LINEAR_MODEL", ::stat_tool::STATL_LINEAR_MODEL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_NONPARAMETRIC", ::stat_tool::STATL_NONPARAMETRIC);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_REGRESSION", ::stat_tool::STATL_REGRESSION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_ESTIMATION", ::stat_tool::STATL_ESTIMATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DISTANCE", ::stat_tool::STATL_DISTANCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CUMUL_DISTANCE", ::stat_tool::STATL_CUMUL_DISTANCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_LENGTH", ::stat_tool::STATL_LENGTH);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DELETION", ::stat_tool::STATL_DELETION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_INSERTION", ::stat_tool::STATL_INSERTION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_INDEL", ::stat_tool::STATL_INDEL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MATCH", ::stat_tool::STATL_MATCH);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SUBSTITUTION", ::stat_tool::STATL_SUBSTITUTION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_TRANSPOSITION", ::stat_tool::STATL_TRANSPOSITION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_RATE", ::stat_tool::STATL_RATE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MATRIX", ::stat_tool::STATL_MATRIX);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_NB_ROW", ::stat_tool::STATL_NB_ROW);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_NB_COLUMN", ::stat_tool::STATL_NB_COLUMN);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SYMMETRY_RATE", ::stat_tool::STATL_SYMMETRY_RATE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_TRIANGLE_INEQUALITY_RATE", ::stat_tool::STATL_TRIANGLE_INEQUALITY_RATE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CLUSTER", ::stat_tool::STATL_CLUSTER);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CLUSTERS", ::stat_tool::STATL_CLUSTERS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_NEIGHBOR", ::stat_tool::STATL_NEIGHBOR);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MOST_DISTANT", ::stat_tool::STATL_MOST_DISTANT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DIAMETER", ::stat_tool::STATL_DIAMETER);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_SEPARATION", ::stat_tool::STATL_SEPARATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_PATTERN_LEVEL", ::stat_tool::STATL_PATTERN_LEVEL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_ISOLATED", ::stat_tool::STATL_ISOLATED);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_NON_ISOLATED", ::stat_tool::STATL_NON_ISOLATED);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_WITHIN", ::stat_tool::STATL_WITHIN);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_BETWEEN", ::stat_tool::STATL_BETWEEN);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_RATIO", ::stat_tool::STATL_RATIO);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_PROTOTYPE", ::stat_tool::STATL_PROTOTYPE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_STEP", ::stat_tool::STATL_STEP);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CHILD", ::stat_tool::STATL_CHILD);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_COMPOSITION", ::stat_tool::STATL_COMPOSITION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CHILD_CLUSTER_DISTANCE_COEFF", ::stat_tool::STATL_CHILD_CLUSTER_DISTANCE_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DIAMETER_COEFF", ::stat_tool::STATL_DIAMETER_COEFF);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_DENDROGRAM_SCALE", ::stat_tool::STATL_DENDROGRAM_SCALE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CLASS", ::stat_tool::STATL_CLASS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_STATE", ::stat_tool::STATL_STATE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_STATES", ::stat_tool::STATL_STATES);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_TRANSIENT", ::stat_tool::STATL_TRANSIENT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_RECURRENT", ::stat_tool::STATL_RECURRENT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_ABSORBING", ::stat_tool::STATL_ABSORBING);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_MEMORY", ::stat_tool::STATL_MEMORY);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_STATIONARY_PROBABILITIES", ::stat_tool::STATL_STATIONARY_PROBABILITIES);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VALUE", ::stat_tool::STATL_VALUE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_VALUES", ::stat_tool::STATL_VALUES);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_OUTPUT", ::stat_tool::STATL_OUTPUT);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_OUTPUT_PROCESS", ::stat_tool::STATL_OUTPUT_PROCESS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_OBSERVATION", ::stat_tool::STATL_OBSERVATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_OBSERVATION_PROBABILITIY_MATRIX", ::stat_tool::STATL_OBSERVATION_PROBABILITIY_MATRIX);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_THEORETICAL", ::stat_tool::STATL_THEORETICAL);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_RESTORATION", ::stat_tool::STATL_RESTORATION);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_WEIGHTS", ::stat_tool::STATL_WEIGHTS);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_OBSERVATION_DISTRIBUTION_DISTANCE", ::stat_tool::STATL_OBSERVATION_DISTRIBUTION_DISTANCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE", ::stat_tool::STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE);

    enum_b9cfa58a8c365c34bdde23c4aa6f45d6.value("STATL_Q_Q_PLOT", ::stat_tool::STATL_Q_Q_PLOT);

}