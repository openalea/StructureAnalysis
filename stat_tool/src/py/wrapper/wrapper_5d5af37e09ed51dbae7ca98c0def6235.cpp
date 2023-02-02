#include "_stat_tool.h"


void wrapper_5d5af37e09ed51dbae7ca98c0def6235()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::stat_tool_error > enum_5d5af37e09ed51dbae7ca98c0def6235("stat_tool_error");

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_FILE_NAME", ::stat_tool::STATR_FILE_NAME);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_FILE_PREFIX", ::stat_tool::STATR_FILE_PREFIX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_EMPTY_SAMPLE", ::stat_tool::STATR_EMPTY_SAMPLE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SAMPLE_SIZE", ::stat_tool::STATR_SAMPLE_SIZE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SAMPLE_SIZES", ::stat_tool::STATR_SAMPLE_SIZES);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NON_EXISTING_DISTRIBUTION", ::stat_tool::STATR_NON_EXISTING_DISTRIBUTION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION", ::stat_tool::STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NO_DATA", ::stat_tool::STATR_NO_DATA);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_PLOT_NB_DISTRIBUTION", ::stat_tool::STATR_PLOT_NB_DISTRIBUTION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_PLOT_NB_HISTOGRAM", ::stat_tool::STATR_PLOT_NB_HISTOGRAM);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_PLOT_NULL_VARIANCE", ::stat_tool::STATR_PLOT_NULL_VARIANCE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_VALUE_RANGE", ::stat_tool::STATR_VALUE_RANGE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MIN_INF_BOUND", ::stat_tool::STATR_MIN_INF_BOUND);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_ESTIMATION_FAILURE", ::stat_tool::STATR_ESTIMATION_FAILURE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_COMPLETE_INTERVAL_TOO_SMALL", ::stat_tool::STATR_NB_COMPLETE_INTERVAL_TOO_SMALL);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_COMPLETE_MIN_VALUE", ::stat_tool::STATR_COMPLETE_MIN_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_FORWARD_MIN_VALUE", ::stat_tool::STATR_FORWARD_MIN_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NO_EVENT_MIN_VALUE", ::stat_tool::STATR_NO_EVENT_MIN_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MEAN_ESTIMATION", ::stat_tool::STATR_MEAN_ESTIMATION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_INTER_EVENT_SUPPORT", ::stat_tool::STATR_INTER_EVENT_SUPPORT);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_DISTRIBUTION_INDEX", ::stat_tool::STATR_DISTRIBUTION_INDEX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_FREQUENCY_DISTRIBUTION_INDEX", ::stat_tool::STATR_FREQUENCY_DISTRIBUTION_INDEX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_DISTRIBUTION", ::stat_tool::STATR_NB_DISTRIBUTION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_WEIGHT_STEP", ::stat_tool::STATR_WEIGHT_STEP);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_DISTRIBUTION_FLAG", ::stat_tool::STATR_DISTRIBUTION_FLAG);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NOT_PRESENT", ::stat_tool::STATR_NOT_PRESENT);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_OUTPUT_PROCESS_INDEX", ::stat_tool::STATR_OUTPUT_PROCESS_INDEX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_OUTPUT_PROCESS_TYPE", ::stat_tool::STATR_OUTPUT_PROCESS_TYPE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_OPTIMAL_ASSIGNMENT", ::stat_tool::STATR_OPTIMAL_ASSIGNMENT);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_VARIABLE_NB_VALUE", ::stat_tool::STATR_VARIABLE_NB_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_OUTPUT_PROCESS", ::stat_tool::STATR_NB_OUTPUT_PROCESS);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_OUTPUT", ::stat_tool::STATR_NB_OUTPUT);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_POSITIVE_MIN_VALUE", ::stat_tool::STATR_POSITIVE_MIN_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_KNOWN_DISTRIBUTION_MIN_VALUE", ::stat_tool::STATR_KNOWN_DISTRIBUTION_MIN_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_KNOWN_DISTRIBUTION_MEAN", ::stat_tool::STATR_KNOWN_DISTRIBUTION_MEAN);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_ITERATION", ::stat_tool::STATR_NB_ITERATION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MIN_NB_ASSIGNMENT", ::stat_tool::STATR_MIN_NB_ASSIGNMENT);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_PENALTY_WEIGHT", ::stat_tool::STATR_PENALTY_WEIGHT);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_CLUSTERING_STEP", ::stat_tool::STATR_CLUSTERING_STEP);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_CLASS", ::stat_tool::STATR_NB_CLASS);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SHIFT_VALUE", ::stat_tool::STATR_SHIFT_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_THRESHOLD_VALUE", ::stat_tool::STATR_THRESHOLD_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_ROUNDED_VALUE", ::stat_tool::STATR_ROUNDED_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SMALLER_THAN", ::stat_tool::STATR_SMALLER_THAN);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_GREATER_THAN", ::stat_tool::STATR_GREATER_THAN);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NOT_ALLOWED", ::stat_tool::STATR_NOT_ALLOWED);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_CATEGORY", ::stat_tool::STATR_NB_CATEGORY);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NON_CONSECUTIVE_CATEGORIES", ::stat_tool::STATR_NON_CONSECUTIVE_CATEGORIES);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MISSING_CATEGORY", ::stat_tool::STATR_MISSING_CATEGORY);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_CLUSTER_LIMIT", ::stat_tool::STATR_CLUSTER_LIMIT);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_INFORMATION_RATIO", ::stat_tool::STATR_INFORMATION_RATIO);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NULL_INFORMATION", ::stat_tool::STATR_NULL_INFORMATION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SCALING_COEFF", ::stat_tool::STATR_SCALING_COEFF);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_VALUE", ::stat_tool::STATR_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MIN_VALUE", ::stat_tool::STATR_MIN_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MAX_VALUE", ::stat_tool::STATR_MAX_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MARGINAL_FREQUENCY_DISTRIBUTION", ::stat_tool::STATR_MARGINAL_FREQUENCY_DISTRIBUTION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MARGINAL_HISTOGRAM", ::stat_tool::STATR_MARGINAL_HISTOGRAM);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_HISTOGRAM_BIN_WIDTH", ::stat_tool::STATR_HISTOGRAM_BIN_WIDTH);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_HISTOGRAM_MIN_VALUE", ::stat_tool::STATR_HISTOGRAM_MIN_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_BAD", ::stat_tool::STATR_BAD);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_VECTOR", ::stat_tool::STATR_NB_VECTOR);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_VARIABLE", ::stat_tool::STATR_NB_VARIABLE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_SELECTED_VARIABLE", ::stat_tool::STATR_NB_SELECTED_VARIABLE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_VARIABLE_TYPE", ::stat_tool::STATR_VARIABLE_TYPE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_VARIABLE_INDEX", ::stat_tool::STATR_VARIABLE_INDEX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_VARIABLE_INDICES", ::stat_tool::STATR_VARIABLE_INDICES);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_ALREADY_USED", ::stat_tool::STATR_ALREADY_USED);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_ALREADY_SELECTED", ::stat_tool::STATR_ALREADY_SELECTED);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SAMPLE_INDEX", ::stat_tool::STATR_SAMPLE_INDEX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_RANK_CORRELATION_COMPUTATION", ::stat_tool::STATR_RANK_CORRELATION_COMPUTATION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SHIFTED_SCALED", ::stat_tool::STATR_SHIFTED_SCALED);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_VALUE_PERIOD", ::stat_tool::STATR_NB_VALUE_PERIOD);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_VALUE", ::stat_tool::STATR_NB_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MISSING_VALUE", ::stat_tool::STATR_MISSING_VALUE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_LEAST_SQUARE_ALGORITHM", ::stat_tool::STATR_LEAST_SQUARE_ALGORITHM);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_REGRESSION_FAILURE", ::stat_tool::STATR_REGRESSION_FAILURE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MATRIX_DIMENSIONS", ::stat_tool::STATR_MATRIX_DIMENSIONS);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_INFINITE_DISTANCES", ::stat_tool::STATR_INFINITE_DISTANCES);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SYMMETRICAL_MATRIX", ::stat_tool::STATR_SYMMETRICAL_MATRIX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_UNSYMMETRICAL_MATRIX", ::stat_tool::STATR_UNSYMMETRICAL_MATRIX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_UNNORMALIZED_DISSIMILARITY_MEASURES", ::stat_tool::STATR_UNNORMALIZED_DISSIMILARITY_MEASURES);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SYMMETRY", ::stat_tool::STATR_SYMMETRY);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_TRIANGLE_INEQUALITY", ::stat_tool::STATR_TRIANGLE_INEQUALITY);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_MATRIX_STRUCTURE", ::stat_tool::STATR_MATRIX_STRUCTURE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SQUARE_MATRIX", ::stat_tool::STATR_SQUARE_MATRIX);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_SINGLE_ELEMENT_CLUSTERS", ::stat_tool::STATR_SINGLE_ELEMENT_CLUSTERS);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NB_CLUSTER", ::stat_tool::STATR_NB_CLUSTER);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_PATTERN_TYPE", ::stat_tool::STATR_PATTERN_TYPE);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_PROTOTYPE_IDENTIFIER", ::stat_tool::STATR_PROTOTYPE_IDENTIFIER);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NUMBER", ::stat_tool::STATR_NUMBER);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NULL_INITIAL_PROBABILITY", ::stat_tool::STATR_NULL_INITIAL_PROBABILITY);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_ODD", ::stat_tool::STATR_ODD);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NON_SYMMETRICAL_DISTRIBUTION", ::stat_tool::STATR_NON_SYMMETRICAL_DISTRIBUTION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_UNPROPER_DISTRIBUTION", ::stat_tool::STATR_UNPROPER_DISTRIBUTION);

    enum_5d5af37e09ed51dbae7ca98c0def6235.value("STATR_NO_PERMUTATION", ::stat_tool::STATR_NO_PERMUTATION);

}