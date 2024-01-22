#include "_stat_tool.h"


void wrapper_ed1db0c35af1528792e2875aed1f9caf()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::enum_< enum ::stat_tool::output_format > enum_ed1db0c35af1528792e2875aed1f9caf("output_format");
    enum_ed1db0c35af1528792e2875aed1f9caf.value("ASCII", ::stat_tool::ASCII);
    enum_ed1db0c35af1528792e2875aed1f9caf.value("SPREADSHEET", ::stat_tool::SPREADSHEET);
    enum_ed1db0c35af1528792e2875aed1f9caf.value("GNUPLOT", ::stat_tool::GNUPLOT);
    enum_ed1db0c35af1528792e2875aed1f9caf.value("PLOT", ::stat_tool::PLOT);

}