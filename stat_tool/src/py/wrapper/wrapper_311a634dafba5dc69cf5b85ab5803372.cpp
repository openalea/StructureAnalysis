#include "_stat_tool.h"


void wrapper_311a634dafba5dc69cf5b85ab5803372()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::model_selection_criterion > enum_311a634dafba5dc69cf5b85ab5803372("model_selection_criterion");
    enum_311a634dafba5dc69cf5b85ab5803372.value("AIC", ::stat_tool::AIC);
    enum_311a634dafba5dc69cf5b85ab5803372.value("AI_CC", ::stat_tool::AICc);
    enum_311a634dafba5dc69cf5b85ab5803372.value("BIC", ::stat_tool::BIC);
    enum_311a634dafba5dc69cf5b85ab5803372.value("BI_CC", ::stat_tool::BICc);
    enum_311a634dafba5dc69cf5b85ab5803372.value("ICL", ::stat_tool::ICL);
    enum_311a634dafba5dc69cf5b85ab5803372.value("IC_LC", ::stat_tool::ICLc);
    enum_311a634dafba5dc69cf5b85ab5803372.value("M_BIC", ::stat_tool::mBIC);
    enum_311a634dafba5dc69cf5b85ab5803372.value("LIKELIHOOD_SLOPE", ::stat_tool::LIKELIHOOD_SLOPE);
    enum_311a634dafba5dc69cf5b85ab5803372.value("DIMENSION_JUMP", ::stat_tool::DIMENSION_JUMP);
    enum_311a634dafba5dc69cf5b85ab5803372.value("SEGMENTATION_LIKELIHOOD_SLOPE", ::stat_tool::SEGMENTATION_LIKELIHOOD_SLOPE);
    enum_311a634dafba5dc69cf5b85ab5803372.value("SEGMENTATION_DIMENSION_JUMP", ::stat_tool::SEGMENTATION_DIMENSION_JUMP);

}