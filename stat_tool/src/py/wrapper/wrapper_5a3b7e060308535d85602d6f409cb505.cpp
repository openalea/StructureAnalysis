#include "_stat_tool.h"


void wrapper_5a3b7e060308535d85602d6f409cb505()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::matrix_transform > enum_5a3b7e060308535d85602d6f409cb505("matrix_transform");
    enum_5a3b7e060308535d85602d6f409cb505.value("COPY", ::stat_tool::COPY);
    enum_5a3b7e060308535d85602d6f409cb505.value("SYMMETRIZATION", ::stat_tool::SYMMETRIZATION);
    enum_5a3b7e060308535d85602d6f409cb505.value("UNNORMALIZATION", ::stat_tool::UNNORMALIZATION);

}