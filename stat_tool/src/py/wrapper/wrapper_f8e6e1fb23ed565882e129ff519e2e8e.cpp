#include "_stat_tool.h"


void wrapper_f8e6e1fb23ed565882e129ff519e2e8e()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::discrete_parametric > enum_f8e6e1fb23ed565882e129ff519e2e8e("discrete_parametric");
    enum_f8e6e1fb23ed565882e129ff519e2e8e.value("CATEGORICAL", ::stat_tool::CATEGORICAL);
    enum_f8e6e1fb23ed565882e129ff519e2e8e.value("BINOMIAL", ::stat_tool::BINOMIAL);
    enum_f8e6e1fb23ed565882e129ff519e2e8e.value("POISSON", ::stat_tool::POISSON);
    enum_f8e6e1fb23ed565882e129ff519e2e8e.value("NEGATIVE_BINOMIAL", ::stat_tool::NEGATIVE_BINOMIAL);
    enum_f8e6e1fb23ed565882e129ff519e2e8e.value("POISSON_GEOMETRIC", ::stat_tool::POISSON_GEOMETRIC);
    enum_f8e6e1fb23ed565882e129ff519e2e8e.value("UNIFORM", ::stat_tool::UNIFORM);
    enum_f8e6e1fb23ed565882e129ff519e2e8e.value("PRIOR_SEGMENT_LENGTH", ::stat_tool::PRIOR_SEGMENT_LENGTH);
    enum_f8e6e1fb23ed565882e129ff519e2e8e.value("MULTINOMIAL", ::stat_tool::MULTINOMIAL);

}