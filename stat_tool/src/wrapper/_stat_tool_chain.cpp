#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>
#include <stat_tool/markovian.h>

void _stat_tool_chain()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::Chain::*method_pointer_a255fb643996510ca503e0599fd445c5)(class ::stat_tool::Chain const &) = &::stat_tool::Chain::parameter_copy;
        void (::stat_tool::Chain::*method_pointer_3e1e9bbe45495cac8fe305ca8064a73d)(class ::stat_tool::Chain const &) = &::stat_tool::Chain::copy;
        void (::stat_tool::Chain::*method_pointer_06eebefd01815b2cbf83a87c5b7f4349)() = &::stat_tool::Chain::remove;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Chain::*method_pointer_f7574889cdbb5fd381197d07f04e7af0)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::Chain::ascii_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Chain::*method_pointer_cac61d75eba05a709779f52328e898f9)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::Chain::spreadsheet_print;
        void (::stat_tool::Chain::*method_pointer_a7a1ab89213756b282b6ca868005bb8b)() = &::stat_tool::Chain::create_cumul;
        void (::stat_tool::Chain::*method_pointer_4d2ce28bf1025608a64163e69ce7beb3)() = &::stat_tool::Chain::cumul_computation;
        void (::stat_tool::Chain::*method_pointer_e4b964206af254b7aa6dc1c1df876b0e)() = &::stat_tool::Chain::remove_cumul;
        void (::stat_tool::Chain::*method_pointer_479008557ee351fa94dd837219f04e24)() = &::stat_tool::Chain::log_computation;
        void (::stat_tool::Chain::*method_pointer_a6a8144bc965575cb12cb1bcd1b31328)() = &::stat_tool::Chain::probability_accessibility_computation;
        void (::stat_tool::Chain::*method_pointer_a1c3a26e57da572fbbe0a5fdf925e646)(double, bool) = &::stat_tool::Chain::thresholding;
        int (::stat_tool::Chain::*method_pointer_169d063bc4985ff49208c4e648dc7b11)(double) const = &::stat_tool::Chain::nb_parameter_computation;
        double (::stat_tool::Chain::*method_pointer_31c042fd128258f39ed5c70ff948be31)(class ::stat_tool::ChainData const &) const = &::stat_tool::Chain::chi2_value_computation;
        void (::stat_tool::Chain::*method_pointer_ad88b036fac55065aae22812dc3a4b4e)(bool, double) = &::stat_tool::Chain::init;
        double (::stat_tool::Chain::*method_pointer_e303c8fae3345709a2ee6b3a84cd16d5)(class ::stat_tool::ChainData const &, bool) const = &::stat_tool::Chain::likelihood_computation;
        void (::stat_tool::Chain::*method_pointer_895ffb83a01d57c99eeb11a71fef736d)(class ::stat_tool::ChainData const &, class ::stat_tool::Test &) const = &::stat_tool::Chain::chi2_fit;
        boost::python::class_< class ::stat_tool::Chain, std::shared_ptr< class ::stat_tool::Chain > >("Chain", boost::python::no_init)
            .def(boost::python::init< enum ::stat_tool::process_type, int, bool >())
            .def(boost::python::init< enum ::stat_tool::process_type, int, int, bool >())
            .def(boost::python::init< class ::stat_tool::Chain const & >())
            .def("parameter_copy", method_pointer_a255fb643996510ca503e0599fd445c5)
            .def("copy", method_pointer_3e1e9bbe45495cac8fe305ca8064a73d)
            .def("remove", method_pointer_06eebefd01815b2cbf83a87c5b7f4349)
            .def("ascii_print", method_pointer_f7574889cdbb5fd381197d07f04e7af0, boost::python::return_internal_reference<>())
            .def("spreadsheet_print", method_pointer_cac61d75eba05a709779f52328e898f9, boost::python::return_internal_reference<>())
            .def("create_cumul", method_pointer_a7a1ab89213756b282b6ca868005bb8b)
            .def("cumul_computation", method_pointer_4d2ce28bf1025608a64163e69ce7beb3)
            .def("remove_cumul", method_pointer_e4b964206af254b7aa6dc1c1df876b0e)
            .def("log_computation", method_pointer_479008557ee351fa94dd837219f04e24)
            .def("probability_accessibility_computation", method_pointer_a6a8144bc965575cb12cb1bcd1b31328)
            .def("thresholding", method_pointer_a1c3a26e57da572fbbe0a5fdf925e646)
            .def("nb_parameter_computation", method_pointer_169d063bc4985ff49208c4e648dc7b11)
            .def("chi2_value_computation", method_pointer_31c042fd128258f39ed5c70ff948be31)
            .def("init", method_pointer_ad88b036fac55065aae22812dc3a4b4e)
            .def("likelihood_computation", method_pointer_e303c8fae3345709a2ee6b3a84cd16d5)
            .def("chi2_fit", method_pointer_895ffb83a01d57c99eeb11a71fef736d)
            .def_readwrite("type", &::stat_tool::Chain::type)
            .def_readwrite("nb_state", &::stat_tool::Chain::nb_state)
            .def_readwrite("nb_row", &::stat_tool::Chain::nb_row)
            .def_readwrite("nb_component", &::stat_tool::Chain::nb_component);
}