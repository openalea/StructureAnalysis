#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_histogram()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::Histogram::*method_pointer_c75330da123c5028b9eab03621ce51a6)(class ::stat_tool::Histogram const &) = &::stat_tool::Histogram::copy;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Histogram::*method_pointer_299e7eadecda52b981d1c38b42254355)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::Histogram::ascii_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Histogram::*method_pointer_5861a387e27753ca9c5b1bd46211b9d4)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::Histogram::spreadsheet_print;
        void (::stat_tool::Histogram::*method_pointer_9fdb0ce985a755cfbd7d4382c1b686ab)(class ::stat_tool::SinglePlot &) const = &::stat_tool::Histogram::plotable_write;
        void (::stat_tool::Histogram::*method_pointer_44692df90a8a558b9f0b9e588c6a1cc9)() = &::stat_tool::Histogram::max_computation;
        boost::python::class_< class ::stat_tool::Histogram, std::shared_ptr< class ::stat_tool::Histogram > >("Histogram", boost::python::no_init)
            .def(boost::python::init< int, bool >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >())
            .def(boost::python::init< class ::stat_tool::Histogram const & >())
            .def("copy", method_pointer_c75330da123c5028b9eab03621ce51a6)
            .def("ascii_print", method_pointer_299e7eadecda52b981d1c38b42254355, boost::python::return_internal_reference<>())
            .def("spreadsheet_print", method_pointer_5861a387e27753ca9c5b1bd46211b9d4, boost::python::return_internal_reference<>())
            .def("plotable_write", method_pointer_9fdb0ce985a755cfbd7d4382c1b686ab)
            .def("max_computation", method_pointer_44692df90a8a558b9f0b9e588c6a1cc9)
            .def_readwrite("nb_element", &::stat_tool::Histogram::nb_element)
            .def_readwrite("nb_category", &::stat_tool::Histogram::nb_category)
            .def_readwrite("step", &::stat_tool::Histogram::step)
            .def_readwrite("max", &::stat_tool::Histogram::max)
            .def_readwrite("type", &::stat_tool::Histogram::type)
            .def_readwrite("min_value", &::stat_tool::Histogram::min_value)
            .def_readwrite("max_value", &::stat_tool::Histogram::max_value);
}