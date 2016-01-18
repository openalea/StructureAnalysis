#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_model_selection_criterion()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::model_selection_criterion >("model_selection_criterion")
            .value("AIC", ::stat_tool::model_selection_criterion::AIC)
            .value("AI_CC", ::stat_tool::model_selection_criterion::AICc)
            .value("BIC", ::stat_tool::model_selection_criterion::BIC)
            .value("BI_CC", ::stat_tool::model_selection_criterion::BICc)
            .value("ICL", ::stat_tool::model_selection_criterion::ICL)
            .value("IC_LC", ::stat_tool::model_selection_criterion::ICLc)
            .value("M_BIC", ::stat_tool::model_selection_criterion::mBIC)
            .value("LIKELIHOOD_SLOPE", ::stat_tool::model_selection_criterion::LIKELIHOOD_SLOPE)
            .value("SEGMENTATION_LIKELIHOOD_SLOPE", ::stat_tool::model_selection_criterion::SEGMENTATION_LIKELIHOOD_SLOPE);
}