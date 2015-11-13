#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_rank_correlation_computation_00d124e24d345c9f927ea5d72a2a405e()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_67c0d3c47c89561cb1d38aab0e04cfc3)(class ::stat_tool::Vectors const &, class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, enum ::stat_tool::correlation_type, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = ::stat_tool::rank_correlation_computation_00d124e24d345c9f927ea5d72a2a405e;
        boost::python::def("rank_correlation_computation_00d124e24d345c9f927ea5d72a2a405e", function_pointer_67c0d3c47c89561cb1d38aab0e04cfc3);
}