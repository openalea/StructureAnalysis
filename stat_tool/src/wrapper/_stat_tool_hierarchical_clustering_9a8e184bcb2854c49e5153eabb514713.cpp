#include <boost/python.hpp>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_hierarchical_clustering_9a8e184bcb2854c49e5153eabb514713()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_49ccb2f8dd8750fc9322c21bd8c7e470)(class ::stat_tool::DistanceMatrix const &, class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, enum ::stat_tool::hierarchical_strategy, enum ::stat_tool::linkage, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, enum ::stat_tool::output_format) = ::stat_tool::hierarchical_clustering_9a8e184bcb2854c49e5153eabb514713;
        boost::python::def("hierarchical_clustering_9a8e184bcb2854c49e5153eabb514713", function_pointer_49ccb2f8dd8750fc9322c21bd8c7e470);
}