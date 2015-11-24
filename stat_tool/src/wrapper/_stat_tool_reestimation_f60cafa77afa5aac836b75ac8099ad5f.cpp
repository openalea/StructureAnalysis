#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_reestimation_f60cafa77afa5aac836b75ac8099ad5f()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::Reestimation<double>::*method_pointer_f18ecadec1f75cb8a1ee6731ad08f8de)(int) = &::stat_tool::Reestimation<double>::init;
        void (::stat_tool::Reestimation<double>::*method_pointer_6cd3e0e967ad5d2586d5a0bdd730bfe9)(class ::stat_tool::Reestimation<double> const &) = &::stat_tool::Reestimation<double>::copy;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Reestimation<double>::*method_pointer_3083678d01fa53ce9491d5a6ef584d22)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool, bool) const = &::stat_tool::Reestimation<double>::ascii_characteristic_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Reestimation<double>::*method_pointer_aa58fe2872665642b1e0159a046f1c7d)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::Reestimation<double>::ascii_circular_characteristic_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::Reestimation<double>::*method_pointer_243ea5296a695aa1aad9278141aba833)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::Reestimation<double>::print;
        void (::stat_tool::Reestimation<double>::*method_pointer_8f96e69a86e651149c66ff25d886ecec)() = &::stat_tool::Reestimation<double>::nb_value_computation;
        void (::stat_tool::Reestimation<double>::*method_pointer_9015730ccf2f5c73b7fbaef9b33fcbaf)() = &::stat_tool::Reestimation<double>::offset_computation;
        void (::stat_tool::Reestimation<double>::*method_pointer_eec7acd160fd5cd5932ab4b434022872)() = &::stat_tool::Reestimation<double>::nb_element_computation;
        void (::stat_tool::Reestimation<double>::*method_pointer_bc7651ab78ac577baf8f7dc8a126da2a)() = &::stat_tool::Reestimation<double>::max_computation;
        void (::stat_tool::Reestimation<double>::*method_pointer_9745dfce6f47520791befcf24b636841)() = &::stat_tool::Reestimation<double>::mean_computation;
        void (::stat_tool::Reestimation<double>::*method_pointer_48a7c9b195de56e2b01a549b0acaad39)(bool) = &::stat_tool::Reestimation<double>::variance_computation;
        double (::stat_tool::Reestimation<double>::*method_pointer_cf4f2348814e5db8bc50c283b7049ea6)() const = &::stat_tool::Reestimation<double>::mean_absolute_deviation_computation;
        double (::stat_tool::Reestimation<double>::*method_pointer_984bf8b493be55a9b6acaf787ec4c375)() const = &::stat_tool::Reestimation<double>::log_geometric_mean_computation;
        double (::stat_tool::Reestimation<double>::*method_pointer_2d5cbd0adfb954bdab4734f8b3377682)() const = &::stat_tool::Reestimation<double>::skewness_computation;
        double (::stat_tool::Reestimation<double>::*method_pointer_a5820387ce135ff9943e4e7d5439be11)() const = &::stat_tool::Reestimation<double>::kurtosis_computation;
        double (::stat_tool::Reestimation<double>::*method_pointer_0207716e29015917be3b157e781cec98)() const = &::stat_tool::Reestimation<double>::information_computation;
        double (::stat_tool::Reestimation<double>::*method_pointer_1ec5601b58d55865a9a00608bd9c8f0b)(class ::stat_tool::Distribution const &) const = &::stat_tool::Reestimation<double>::likelihood_computation;
        void (::stat_tool::Reestimation<double>::*method_pointer_e104d96f06ae53ea9ccd57ca42d0d60a)(class ::stat_tool::Distribution *) const = &::stat_tool::Reestimation<double>::distribution_estimation;
        double (::stat_tool::Reestimation<double>::*method_pointer_c8cbb95a801a5f8dafdab61cc3c93620)(class ::stat_tool::DiscreteParametric *, int, bool) const = &::stat_tool::Reestimation<double>::binomial_estimation;
        double (::stat_tool::Reestimation<double>::*method_pointer_6c44cee83702546aa24e1ae2aef0f2b4)(class ::stat_tool::DiscreteParametric *, int, bool, double) const = &::stat_tool::Reestimation<double>::poisson_estimation;
        double (::stat_tool::Reestimation<double>::*method_pointer_156c045ee4065701bdcb7ba92fc21e93)(class ::stat_tool::DiscreteParametric *, int, bool, double) const = &::stat_tool::Reestimation<double>::negative_binomial_estimation;
        double (::stat_tool::Reestimation<double>::*method_pointer_8ed1fe0f85f55060a4e1da40752dd13e)(class ::stat_tool::DiscreteParametric *, int, bool, double) const = &::stat_tool::Reestimation<double>::parametric_estimation;
        double (::stat_tool::Reestimation<double>::*method_pointer_90725ca276e45ffe8d63772e875b8f10)(class ::stat_tool::DiscreteParametric *, int, bool, double) const = &::stat_tool::Reestimation<double>::type_parametric_estimation;
        class ::stat_tool::DiscreteParametric * (::stat_tool::Reestimation<double>::*method_pointer_ec141b3c67a85c9cb9e0d3c14437642f)(int, bool, double) const = &::stat_tool::Reestimation<double>::type_parametric_estimation;
        void (::stat_tool::Reestimation<double>::*method_pointer_781241ffbbeb578aa34cdf9df9db7e16)(class ::stat_tool::Reestimation<double> const *, double) = &::stat_tool::Reestimation<double>::equilibrium_process_combination;
        void (::stat_tool::Reestimation<double>::*method_pointer_5012426ad7c55128a0e06e1f283ca721)(class ::stat_tool::Reestimation<double> const *, class ::stat_tool::Distribution *, double) const = &::stat_tool::Reestimation<double>::equilibrium_process_estimation;
        void (::stat_tool::Reestimation<double>::*method_pointer_902f74d547d5527b88d49fcb64e42dcb)(class ::stat_tool::ContinuousParametric *, int) const = &::stat_tool::Reestimation<double>::gamma_estimation;
        void (::stat_tool::Reestimation<double>::*method_pointer_1f3650165f6c546ea589f46946e5f312)(class ::stat_tool::ContinuousParametric *, int) const = &::stat_tool::Reestimation<double>::zero_inflated_gamma_estimation;
        boost::python::class_< class ::stat_tool::Reestimation<double>, std::shared_ptr< class ::stat_tool::Reestimation<double> > >("_Reestimation_f60cafa77afa5aac836b75ac8099ad5f", boost::python::no_init)
            .def(boost::python::init< int >())
            .def(boost::python::init< class ::stat_tool::Reestimation<double> const & >())
            .def(boost::python::init< class ::std::vector<stat_tool::Reestimation<double> *, std::allocator<stat_tool::Reestimation<double> *> > const & >())
            .def("init", method_pointer_f18ecadec1f75cb8a1ee6731ad08f8de)
            .def("copy", method_pointer_6cd3e0e967ad5d2586d5a0bdd730bfe9)
            .def("ascii_characteristic_print", method_pointer_3083678d01fa53ce9491d5a6ef584d22, boost::python::return_internal_reference<>())
            .def("ascii_circular_characteristic_print", method_pointer_aa58fe2872665642b1e0159a046f1c7d, boost::python::return_internal_reference<>())
            .def("print", method_pointer_243ea5296a695aa1aad9278141aba833, boost::python::return_internal_reference<>())
            .def("nb_value_computation", method_pointer_8f96e69a86e651149c66ff25d886ecec)
            .def("offset_computation", method_pointer_9015730ccf2f5c73b7fbaef9b33fcbaf)
            .def("nb_element_computation", method_pointer_eec7acd160fd5cd5932ab4b434022872)
            .def("max_computation", method_pointer_bc7651ab78ac577baf8f7dc8a126da2a)
            .def("mean_computation", method_pointer_9745dfce6f47520791befcf24b636841)
            .def("variance_computation", method_pointer_48a7c9b195de56e2b01a549b0acaad39)
            .def("mean_absolute_deviation_computation", method_pointer_cf4f2348814e5db8bc50c283b7049ea6)
            .def("log_geometric_mean_computation", method_pointer_984bf8b493be55a9b6acaf787ec4c375)
            .def("skewness_computation", method_pointer_2d5cbd0adfb954bdab4734f8b3377682)
            .def("kurtosis_computation", method_pointer_a5820387ce135ff9943e4e7d5439be11)
            .def("information_computation", method_pointer_0207716e29015917be3b157e781cec98)
            .def("likelihood_computation", method_pointer_1ec5601b58d55865a9a00608bd9c8f0b)
            .def("distribution_estimation", method_pointer_e104d96f06ae53ea9ccd57ca42d0d60a)
            .def("binomial_estimation", method_pointer_c8cbb95a801a5f8dafdab61cc3c93620)
            .def("poisson_estimation", method_pointer_6c44cee83702546aa24e1ae2aef0f2b4)
            .def("negative_binomial_estimation", method_pointer_156c045ee4065701bdcb7ba92fc21e93)
            .def("parametric_estimation", method_pointer_8ed1fe0f85f55060a4e1da40752dd13e)
            .def("type_parametric_estimation", method_pointer_90725ca276e45ffe8d63772e875b8f10)
            .def("type_parametric_estimation", method_pointer_ec141b3c67a85c9cb9e0d3c14437642f, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("equilibrium_process_combination", method_pointer_781241ffbbeb578aa34cdf9df9db7e16)
            .def("equilibrium_process_estimation", method_pointer_5012426ad7c55128a0e06e1f283ca721)
            .def("gamma_estimation", method_pointer_902f74d547d5527b88d49fcb64e42dcb)
            .def("zero_inflated_gamma_estimation", method_pointer_1f3650165f6c546ea589f46946e5f312)
            .def_readwrite("nb_value", &::stat_tool::Reestimation<double>::nb_value)
            .def_readwrite("alloc_nb_value", &::stat_tool::Reestimation<double>::alloc_nb_value)
            .def_readwrite("offset", &::stat_tool::Reestimation<double>::offset)
            .def_readwrite("nb_element", &::stat_tool::Reestimation<double>::nb_element)
            .def_readwrite("max", &::stat_tool::Reestimation<double>::max)
            .def_readwrite("mean", &::stat_tool::Reestimation<double>::mean)
            .def_readwrite("variance", &::stat_tool::Reestimation<double>::variance);
}