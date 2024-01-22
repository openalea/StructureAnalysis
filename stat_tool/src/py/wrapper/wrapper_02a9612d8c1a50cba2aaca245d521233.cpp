#include "_stat_tool.h"


void wrapper_02a9612d8c1a50cba2aaca245d521233()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::stat_tool_keyword > enum_02a9612d8c1a50cba2aaca245d521233("stat_tool_keyword");
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_INF_BOUND", ::stat_tool::STATW_INF_BOUND);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_SUP_BOUND", ::stat_tool::STATW_SUP_BOUND);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_PARAMETER", ::stat_tool::STATW_PARAMETER);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_PROBABILITY", ::stat_tool::STATW_PROBABILITY);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_NO_SEGMENT", ::stat_tool::STATW_NO_SEGMENT);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_SEQUENCE_LENGTH", ::stat_tool::STATW_SEQUENCE_LENGTH);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_SHAPE", ::stat_tool::STATW_SHAPE);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_SCALE", ::stat_tool::STATW_SCALE);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_ZERO_PROBABILITY", ::stat_tool::STATW_ZERO_PROBABILITY);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_MEAN", ::stat_tool::STATW_MEAN);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_STANDARD_DEVIATION", ::stat_tool::STATW_STANDARD_DEVIATION);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_MEAN_DIRECTION", ::stat_tool::STATW_MEAN_DIRECTION);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_CONCENTRATION", ::stat_tool::STATW_CONCENTRATION);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_INTERCEPT", ::stat_tool::STATW_INTERCEPT);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_SLOPE", ::stat_tool::STATW_SLOPE);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_AUTOREGRESSIVE_COEFF", ::stat_tool::STATW_AUTOREGRESSIVE_COEFF);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_MIXTURE", ::stat_tool::STATW_MIXTURE);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_CONVOLUTION", ::stat_tool::STATW_CONVOLUTION);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_DISTRIBUTIONS", ::stat_tool::STATW_DISTRIBUTIONS);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_DISTRIBUTION", ::stat_tool::STATW_DISTRIBUTION);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_WEIGHT", ::stat_tool::STATW_WEIGHT);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_WEIGHTS", ::stat_tool::STATW_WEIGHTS);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_COMPONENTS", ::stat_tool::STATW_COMPONENTS);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_COMPONENT", ::stat_tool::STATW_COMPONENT);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_COMPOUND", ::stat_tool::STATW_COMPOUND);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_SUM", ::stat_tool::STATW_SUM);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_ELEMENTARY", ::stat_tool::STATW_ELEMENTARY);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_VARIABLE", ::stat_tool::STATW_VARIABLE);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_VARIABLES", ::stat_tool::STATW_VARIABLES);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_DISTANCE", ::stat_tool::STATW_DISTANCE);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_CATEGORIES", ::stat_tool::STATW_CATEGORIES);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_PERIOD", ::stat_tool::STATW_PERIOD);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_STATES", ::stat_tool::STATW_STATES);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_INITIAL_PROBABILITIES", ::stat_tool::STATW_INITIAL_PROBABILITIES);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_TRANSITION_PROBABILITIES", ::stat_tool::STATW_TRANSITION_PROBABILITIES);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_STATE", ::stat_tool::STATW_STATE);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_FUNCTION", ::stat_tool::STATW_FUNCTION);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_OUTPUT_PROCESS", ::stat_tool::STATW_OUTPUT_PROCESS);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_OUTPUT_PROCESSES", ::stat_tool::STATW_OUTPUT_PROCESSES);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_NONPARAMETRIC", ::stat_tool::STATW_NONPARAMETRIC);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_CATEGORICAL", ::stat_tool::STATW_CATEGORICAL);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_PARAMETRIC", ::stat_tool::STATW_PARAMETRIC);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_DISCRETE_PARAMETRIC", ::stat_tool::STATW_DISCRETE_PARAMETRIC);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_CONTINUOUS_PARAMETRIC", ::stat_tool::STATW_CONTINUOUS_PARAMETRIC);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_OBSERVATION_DISTRIBUTION", ::stat_tool::STATW_OBSERVATION_DISTRIBUTION);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_OBSERVATION_MODEL", ::stat_tool::STATW_OBSERVATION_MODEL);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_OBSERVATION_PROBABILITIES", ::stat_tool::STATW_OBSERVATION_PROBABILITIES);
    enum_02a9612d8c1a50cba2aaca245d521233.value("STATW_OUTPUT", ::stat_tool::STATW_OUTPUT);

}