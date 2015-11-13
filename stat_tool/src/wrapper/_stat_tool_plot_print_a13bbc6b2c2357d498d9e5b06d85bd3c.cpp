#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/char_ptr_to_string.h>
#include <stat_tool/stat_tools.h>

void _stat_tool_plot_print_a13bbc6b2c2357d498d9e5b06d85bd3c()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        bool (*function_pointer_cd90cb498b375f8ba96ab81519a86146)(class ::stat_tool::ContinuousParametric &, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &, class ::stat_tool::Histogram const *, class ::stat_tool::FrequencyDistribution const *) = ::stat_tool::plot_print_a13bbc6b2c2357d498d9e5b06d85bd3c;
        boost::python::def("plot_print_a13bbc6b2c2357d498d9e5b06d85bd3c", function_pointer_cd90cb498b375f8ba96ab81519a86146);
}