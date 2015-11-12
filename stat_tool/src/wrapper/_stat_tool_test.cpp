#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_test()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::Test::*method_pointer_5d1e82f49af15250bdfc9ed3c007426f)(class ::stat_tool::Test const &) = &::stat_tool::Test::copy;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Test::*method_pointer_4da67ea9bb6f59dc8ce44b976fb440ff)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool, bool) const = &::stat_tool::Test::ascii_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Test::*method_pointer_2e7f62768c54586baba69036ba2b6ff4)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::Test::spreadsheet_print;
        void (::stat_tool::Test::*method_pointer_fe1275e4546a5da09641c670d5b6177d)() = &::stat_tool::Test::standard_normal_critical_probability_computation;
        void (::stat_tool::Test::*method_pointer_1828d9c2d2415ae68b778f6d26c05036)() = &::stat_tool::Test::standard_normal_value_computation;
        void (::stat_tool::Test::*method_pointer_a2f30cdc9f8656f6a5b9b12990f0187c)() = &::stat_tool::Test::chi2_critical_probability_computation;
        void (::stat_tool::Test::*method_pointer_e2b31f57c62459d2b55949f85f609ab5)() = &::stat_tool::Test::chi2_value_computation;
        void (::stat_tool::Test::*method_pointer_0496cfc256345921a65582f307882888)() = &::stat_tool::Test::F_critical_probability_computation;
        void (::stat_tool::Test::*method_pointer_395c5b3a9b045414a12c16407e5beb09)() = &::stat_tool::Test::F_value_computation;
        void (::stat_tool::Test::*method_pointer_56006b70c92b5107b03c64bf8d3a993a)() = &::stat_tool::Test::t_critical_probability_computation;
        void (::stat_tool::Test::*method_pointer_6dc9bf3a118e5863904b7e8efab67296)() = &::stat_tool::Test::t_value_computation;
        boost::python::class_< class ::stat_tool::Test, std::shared_ptr< class ::stat_tool::Test > >("Test", boost::python::no_init)
            .def(boost::python::init< enum ::stat_tool::test_distribution, bool >())
            .def(boost::python::init< enum ::stat_tool::test_distribution, bool, int, int, double >())
            .def(boost::python::init< enum ::stat_tool::test_distribution, bool, int, int, double, double >())
            .def(boost::python::init< class ::stat_tool::Test const & >())
            .def("copy", method_pointer_5d1e82f49af15250bdfc9ed3c007426f)
            .def("ascii_print", method_pointer_4da67ea9bb6f59dc8ce44b976fb440ff, boost::python::return_internal_reference<>())
            .def("spreadsheet_print", method_pointer_2e7f62768c54586baba69036ba2b6ff4, boost::python::return_internal_reference<>())
            .def("standard_normal_critical_probability_computation", method_pointer_fe1275e4546a5da09641c670d5b6177d)
            .def("standard_normal_value_computation", method_pointer_1828d9c2d2415ae68b778f6d26c05036)
            .def("chi2_critical_probability_computation", method_pointer_a2f30cdc9f8656f6a5b9b12990f0187c)
            .def("chi2_value_computation", method_pointer_e2b31f57c62459d2b55949f85f609ab5)
            .def("f__critical_probability_computation", method_pointer_0496cfc256345921a65582f307882888)
            .def("f__value_computation", method_pointer_395c5b3a9b045414a12c16407e5beb09)
            .def("t_critical_probability_computation", method_pointer_56006b70c92b5107b03c64bf8d3a993a)
            .def("t_value_computation", method_pointer_6dc9bf3a118e5863904b7e8efab67296)
            .def_readwrite("ident", &::stat_tool::Test::ident)
            .def_readwrite("one_side", &::stat_tool::Test::one_side)
            .def_readwrite("df1", &::stat_tool::Test::df1)
            .def_readwrite("df2", &::stat_tool::Test::df2)
            .def_readwrite("value", &::stat_tool::Test::value)
            .def_readwrite("critical_probability", &::stat_tool::Test::critical_probability);
}