#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_ascii_read_543f048128b55b8ca39ac22bf2f802ce()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::DiscreteMixture * (*function_pointer_69948f1ac8225c6d95fae6dbcb3396a9)(class ::stat_tool::StatError &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, double) = ::stat_tool::ascii_read_543f048128b55b8ca39ac22bf2f802ce;
        boost::python::def("ascii_read_543f048128b55b8ca39ac22bf2f802ce", function_pointer_69948f1ac8225c6d95fae6dbcb3396a9, boost::python::return_value_policy< boost::python::reference_existing_object >());
}