#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/discrete_mixture.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>

void _stat_tool_continuous_parametric()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::continuous_parametric >("continuous_parametric")
            .value("GAMMA", ::stat_tool::continuous_parametric::GAMMA)
            .value("INVERSE_GAUSSIAN", ::stat_tool::continuous_parametric::INVERSE_GAUSSIAN)
            .value("GAUSSIAN", ::stat_tool::continuous_parametric::GAUSSIAN)
            .value("VON_MISES", ::stat_tool::continuous_parametric::VON_MISES)
            .value("ZERO_INFLATED_GAMMA", ::stat_tool::continuous_parametric::ZERO_INFLATED_GAMMA)
            .value("LINEAR_MODEL", ::stat_tool::continuous_parametric::LINEAR_MODEL);
        void (::stat_tool::ContinuousParametric::*method_pointer_51d3f6aaff1c56f9958e931bf18034cc)(class ::stat_tool::ContinuousParametric const &) = &::stat_tool::ContinuousParametric::copy;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ContinuousParametric::*method_pointer_f33e653ce3c45484aaa449282e51e100)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::ContinuousParametric::ascii_parameter_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ContinuousParametric::*method_pointer_094becabde1558ffa4de95bf41b96b5a)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::ContinuousParametric::ascii_characteristic_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ContinuousParametric::*method_pointer_5d741de66ca05becb1d7e36a0888e789)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool, bool, class ::stat_tool::Histogram const *, class ::stat_tool::FrequencyDistribution const *) = &::stat_tool::ContinuousParametric::ascii_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ContinuousParametric::*method_pointer_13b5f743aa025178a29482819de99c1b)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::ContinuousParametric::spreadsheet_parameter_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ContinuousParametric::*method_pointer_6479e0bbce00577caf356e5efe797bcb)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool, class ::stat_tool::Histogram const *, class ::stat_tool::FrequencyDistribution const *) = &::stat_tool::ContinuousParametric::spreadsheet_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ContinuousParametric::*method_pointer_38e87bcc468751c9b329303b0fd5b1ce)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::ContinuousParametric::spreadsheet_characteristic_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::ContinuousParametric::*method_pointer_838f5b7384dd529d9a3cac3e2e017e66)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::ContinuousParametric::plot_title_print;
        void (::stat_tool::ContinuousParametric::*method_pointer_e07a6588990458909e46971a75350dae)(class ::stat_tool::SinglePlot &, class ::stat_tool::Histogram const *, class ::stat_tool::FrequencyDistribution const *) = &::stat_tool::ContinuousParametric::plotable_write;
        int (::stat_tool::ContinuousParametric::*method_pointer_eb05dd7eac87515c9bc1748a0004daa7)() const = &::stat_tool::ContinuousParametric::nb_parameter_computation;
        double (::stat_tool::ContinuousParametric::*method_pointer_2912c7228ded5fe0ad4edc7b170b41fd)(double, double) const = &::stat_tool::ContinuousParametric::von_mises_mass_computation;
        double (::stat_tool::ContinuousParametric::*method_pointer_34ac7fb3a6e95c50b999ebc59d6f7bbe)(double, double) const = &::stat_tool::ContinuousParametric::mass_computation;
        void (::stat_tool::ContinuousParametric::*method_pointer_2753680aa4bf56e888b2b587bf22bfd9)() = &::stat_tool::ContinuousParametric::von_mises_cumul_computation;
        double (::stat_tool::ContinuousParametric::*method_pointer_744741da83c85b21b16946091ba3b49c)(class ::stat_tool::ContinuousParametric &) = &::stat_tool::ContinuousParametric::sup_norm_distance_computation;
        double (::stat_tool::ContinuousParametric::*method_pointer_2a908d5857785e3ba04ff72f253dd84b)(class ::stat_tool::FrequencyDistribution const &, int) const = &::stat_tool::ContinuousParametric::likelihood_computation;
        double (::stat_tool::ContinuousParametric::*method_pointer_b440105ffde95546868447c3762fa258)() = &::stat_tool::ContinuousParametric::simulation;
        boost::python::class_< class ::stat_tool::ContinuousParametric, std::shared_ptr< class ::stat_tool::ContinuousParametric > >("ContinuousParametric", boost::python::no_init)
            .def(boost::python::init< enum ::stat_tool::continuous_parametric, double, double, double, enum ::stat_tool::angle_unit >())
            .def(boost::python::init< class ::stat_tool::ContinuousParametric const & >())
            .def("copy", method_pointer_51d3f6aaff1c56f9958e931bf18034cc)
            .def("ascii_parameter_print", method_pointer_f33e653ce3c45484aaa449282e51e100, boost::python::return_internal_reference<>())
            .def("ascii_characteristic_print", method_pointer_094becabde1558ffa4de95bf41b96b5a, boost::python::return_internal_reference<>())
            .def("ascii_print", method_pointer_5d741de66ca05becb1d7e36a0888e789, boost::python::return_internal_reference<>())
            .def("spreadsheet_parameter_print", method_pointer_13b5f743aa025178a29482819de99c1b, boost::python::return_internal_reference<>())
            .def("spreadsheet_print", method_pointer_6479e0bbce00577caf356e5efe797bcb, boost::python::return_internal_reference<>())
            .def("spreadsheet_characteristic_print", method_pointer_38e87bcc468751c9b329303b0fd5b1ce, boost::python::return_internal_reference<>())
            .def("plot_title_print", method_pointer_838f5b7384dd529d9a3cac3e2e017e66, boost::python::return_internal_reference<>())
            .def("plotable_write", method_pointer_e07a6588990458909e46971a75350dae)
            .def("nb_parameter_computation", method_pointer_eb05dd7eac87515c9bc1748a0004daa7)
            .def("von_mises_mass_computation", method_pointer_2912c7228ded5fe0ad4edc7b170b41fd)
            .def("mass_computation", method_pointer_34ac7fb3a6e95c50b999ebc59d6f7bbe)
            .def("von_mises_cumul_computation", method_pointer_2753680aa4bf56e888b2b587bf22bfd9)
            .def("sup_norm_distance_computation", method_pointer_744741da83c85b21b16946091ba3b49c)
            .def("likelihood_computation", method_pointer_2a908d5857785e3ba04ff72f253dd84b)
            .def("simulation", method_pointer_b440105ffde95546868447c3762fa258)
            .def_readwrite("ident", &::stat_tool::ContinuousParametric::ident)
            .def_readwrite("min_value", &::stat_tool::ContinuousParametric::min_value)
            .def_readwrite("max_value", &::stat_tool::ContinuousParametric::max_value)
            .def_readwrite("unit", &::stat_tool::ContinuousParametric::unit)
            .def_readwrite("slope_standard_deviation", &::stat_tool::ContinuousParametric::slope_standard_deviation)
            .def_readwrite("sample_size", &::stat_tool::ContinuousParametric::sample_size)
            .def_readwrite("correlation", &::stat_tool::ContinuousParametric::correlation);
}