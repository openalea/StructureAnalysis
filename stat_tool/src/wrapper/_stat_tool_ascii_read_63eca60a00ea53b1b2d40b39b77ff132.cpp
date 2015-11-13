#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_ascii_read_63eca60a00ea53b1b2d40b39b77ff132()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::Convolution * (*function_pointer_145f6199c2e455fc8de8a2084f41ed85)(class ::stat_tool::StatError &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, double) = ::stat_tool::ascii_read_63eca60a00ea53b1b2d40b39b77ff132;
        boost::python::def("ascii_read_63eca60a00ea53b1b2d40b39b77ff132", function_pointer_145f6199c2e455fc8de8a2084f41ed85, boost::python::return_value_policy< boost::python::reference_existing_object >());
}