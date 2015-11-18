#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/convolution.h>
#include <stat_tool/discrete_mixture.h>

void _stat_tool_discrete_parametric()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::enum_< enum ::stat_tool::discrete_parametric >("discrete_parametric")
            .value("CATEGORICAL", ::stat_tool::discrete_parametric::CATEGORICAL)
            .value("BINOMIAL", ::stat_tool::discrete_parametric::BINOMIAL)
            .value("POISSON", ::stat_tool::discrete_parametric::POISSON)
            .value("NEGATIVE_BINOMIAL", ::stat_tool::discrete_parametric::NEGATIVE_BINOMIAL)
            .value("UNIFORM", ::stat_tool::discrete_parametric::UNIFORM)
            .value("MULTINOMIAL", ::stat_tool::discrete_parametric::MULTINOMIAL);
        void (::stat_tool::DiscreteParametric::*method_pointer_abfba932b25d52f7bfd2f2ea5b43da88)(int, int, double, double) = &::stat_tool::DiscreteParametric::init;
        void (::stat_tool::DiscreteParametric::*method_pointer_d60dea7c925a54c0bd2a855bbe1a5f58)(enum ::stat_tool::discrete_parametric, int, int, double, double) = &::stat_tool::DiscreteParametric::init;
        void (::stat_tool::DiscreteParametric::*method_pointer_ae8dca5e086c5c93a61df11a32ad57c0)(class ::stat_tool::DiscreteParametric const &) = &::stat_tool::DiscreteParametric::copy;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteParametric::*method_pointer_42a9ef98ce235c25ba78c8418e2cff9c)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::DiscreteParametric::ascii_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteParametric::*method_pointer_1e2470ecea845c5ebc32d2bcfbc98ee4)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool, bool) const = &::stat_tool::DiscreteParametric::ascii_parametric_characteristic_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteParametric::*method_pointer_bea7fd4010685dd38f61049f78e7197f)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::DiscreteParametric::spreadsheet_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteParametric::*method_pointer_26ddfbf705055ef4800204a666cf0e0b)(class ::std::basic_ostream<char, std::char_traits<char> > &, bool) const = &::stat_tool::DiscreteParametric::spreadsheet_parametric_characteristic_print;
        class ::std::basic_ostream<char, std::char_traits<char> > & (::stat_tool::DiscreteParametric::*method_pointer_44c56fc88c6856aba7893009d7e78e2c)(class ::std::basic_ostream<char, std::char_traits<char> > &) const = &::stat_tool::DiscreteParametric::plot_title_print;
        int (*method_pointer_60c76b854c1d563496f805ccd43c71d0)(enum ::stat_tool::discrete_parametric, int, int, double, double, double) = ::stat_tool::DiscreteParametric::nb_value_computation;
        int (::stat_tool::DiscreteParametric::*method_pointer_a24cc9da2b445b8fa9fbd5e1fa58f5f1)() = &::stat_tool::DiscreteParametric::nb_parameter_computation;
        void (::stat_tool::DiscreteParametric::*method_pointer_4e8d96996526545bb344532cc0ebec2d)() = &::stat_tool::DiscreteParametric::nb_parameter_update;
        double (::stat_tool::DiscreteParametric::*method_pointer_3ba4fb2776bf5bbda2feb6b1d71eb37e)() const = &::stat_tool::DiscreteParametric::parametric_mean_computation;
        double (::stat_tool::DiscreteParametric::*method_pointer_82ad6ecc57b75a5382f3c3a1215c6a9f)() const = &::stat_tool::DiscreteParametric::parametric_variance_computation;
        double (::stat_tool::DiscreteParametric::*method_pointer_2007520c52a35b8dae3f03347d06460b)() const = &::stat_tool::DiscreteParametric::parametric_skewness_computation;
        double (::stat_tool::DiscreteParametric::*method_pointer_35685c0d655e595c80326e55d42db7aa)() const = &::stat_tool::DiscreteParametric::parametric_kurtosis_computation;
        double (::stat_tool::DiscreteParametric::*method_pointer_e93c764042a8541886aef2c7b94ef8de)(class ::stat_tool::DiscreteParametric const &) const = &::stat_tool::DiscreteParametric::sup_norm_distance_computation;
        void (::stat_tool::DiscreteParametric::*method_pointer_f08fd4cf5057564bb24d4fb8eaae5a5e)(int, enum ::stat_tool::distribution_computation) = &::stat_tool::DiscreteParametric::binomial_computation;
        void (::stat_tool::DiscreteParametric::*method_pointer_e2982f3ab18058e589f27cc1c9eb2e31)(int, double, enum ::stat_tool::distribution_computation) = &::stat_tool::DiscreteParametric::poisson_computation;
        void (::stat_tool::DiscreteParametric::*method_pointer_88f3a00219075c28988cca75b0b86470)(int, double, enum ::stat_tool::distribution_computation) = &::stat_tool::DiscreteParametric::negative_binomial_computation;
        void (::stat_tool::DiscreteParametric::*method_pointer_51178caa41f9529b9e99d27f1cc17698)() = &::stat_tool::DiscreteParametric::uniform_computation;
        void (::stat_tool::DiscreteParametric::*method_pointer_57ac64387476563898463f45cfd5f06d)(int, double) = &::stat_tool::DiscreteParametric::computation;
        int (::stat_tool::DiscreteParametric::*method_pointer_f92a8082d2cf58e8acbe4f9d734c54eb)() const = &::stat_tool::DiscreteParametric::simulation;
        double (::stat_tool::DiscreteParametric::*method_pointer_f4e78197979c52a0a7542d0ff9dfaf99)(class ::stat_tool::Forward const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const *) const = &::stat_tool::DiscreteParametric::renewal_likelihood_computation;
        void (::stat_tool::DiscreteParametric::*method_pointer_fb81b674097f501c9e2fe702b1ce1772)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const *, class ::stat_tool::Reestimation<double> *, class ::stat_tool::Reestimation<double> *, int) const = &::stat_tool::DiscreteParametric::expectation_step;
        void (::stat_tool::DiscreteParametric::*method_pointer_d6c5d94497175d61bdacfc03605de75a)(class ::stat_tool::Reestimation<double> const *, int) = &::stat_tool::DiscreteParametric::reestimation;
        double (::stat_tool::DiscreteParametric::*method_pointer_ea4f3357a8be53688998cb7f27f3b865)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::DiscreteParametric::state_occupancy_likelihood_computation;
        double (::stat_tool::DiscreteParametric::*method_pointer_15a05927895c516cad12850712c13c6b)(class ::stat_tool::Forward const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::DiscreteParametric::state_occupancy_likelihood_computation;
        void (::stat_tool::DiscreteParametric::*method_pointer_2695736208d05343809912a569eace13)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Reestimation<double> *, int) const = &::stat_tool::DiscreteParametric::expectation_step;
        void (::stat_tool::DiscreteParametric::*method_pointer_def459f4cd2a5beab74b4f631c0bd319)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Reestimation<double> *, class ::stat_tool::Reestimation<double> *, int, bool, enum ::stat_tool::duration_distribution_mean_estimator) const = &::stat_tool::DiscreteParametric::expectation_step;
        boost::python::class_< class ::stat_tool::DiscreteParametric, std::shared_ptr< class ::stat_tool::DiscreteParametric >, boost::python::bases< class ::stat_tool::Distribution > >("DiscreteParametric", boost::python::no_init)
            .def(boost::python::init< int, enum ::stat_tool::discrete_parametric, int, int, double, double >())
            .def(boost::python::init< enum ::stat_tool::discrete_parametric, int, int, double, double, double >())
            .def(boost::python::init< class ::stat_tool::Distribution const &, int >())
            .def(boost::python::init< class ::stat_tool::Distribution const &, double >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametric const &, double >())
            .def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >())
            .def(boost::python::init< class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::distribution_transformation, int >())
            .def("init", method_pointer_abfba932b25d52f7bfd2f2ea5b43da88)
            .def("init", method_pointer_d60dea7c925a54c0bd2a855bbe1a5f58)
            .def("copy", method_pointer_ae8dca5e086c5c93a61df11a32ad57c0)
            .def("ascii_print", method_pointer_42a9ef98ce235c25ba78c8418e2cff9c, boost::python::return_internal_reference<>())
            .def("ascii_parametric_characteristic_print", method_pointer_1e2470ecea845c5ebc32d2bcfbc98ee4, boost::python::return_internal_reference<>())
            .def("spreadsheet_print", method_pointer_bea7fd4010685dd38f61049f78e7197f, boost::python::return_internal_reference<>())
            .def("spreadsheet_parametric_characteristic_print", method_pointer_26ddfbf705055ef4800204a666cf0e0b, boost::python::return_internal_reference<>())
            .def("plot_title_print", method_pointer_44c56fc88c6856aba7893009d7e78e2c, boost::python::return_internal_reference<>())
            .def("nb_value_computation", method_pointer_60c76b854c1d563496f805ccd43c71d0)
            .def("nb_parameter_computation", method_pointer_a24cc9da2b445b8fa9fbd5e1fa58f5f1)
            .def("nb_parameter_update", method_pointer_4e8d96996526545bb344532cc0ebec2d)
            .def("parametric_mean_computation", method_pointer_3ba4fb2776bf5bbda2feb6b1d71eb37e)
            .def("parametric_variance_computation", method_pointer_82ad6ecc57b75a5382f3c3a1215c6a9f)
            .def("parametric_skewness_computation", method_pointer_2007520c52a35b8dae3f03347d06460b)
            .def("parametric_kurtosis_computation", method_pointer_35685c0d655e595c80326e55d42db7aa)
            .def("sup_norm_distance_computation", method_pointer_e93c764042a8541886aef2c7b94ef8de)
            .def("binomial_computation", method_pointer_f08fd4cf5057564bb24d4fb8eaae5a5e)
            .def("poisson_computation", method_pointer_e2982f3ab18058e589f27cc1c9eb2e31)
            .def("negative_binomial_computation", method_pointer_88f3a00219075c28988cca75b0b86470)
            .def("uniform_computation", method_pointer_51178caa41f9529b9e99d27f1cc17698)
            .def("computation", method_pointer_57ac64387476563898463f45cfd5f06d)
            .def("simulation", method_pointer_f92a8082d2cf58e8acbe4f9d734c54eb)
            .def("renewal_likelihood_computation", method_pointer_f4e78197979c52a0a7542d0ff9dfaf99)
            .def("expectation_step", method_pointer_fb81b674097f501c9e2fe702b1ce1772)
            .def("reestimation", method_pointer_d6c5d94497175d61bdacfc03605de75a)
            .def("state_occupancy_likelihood_computation", method_pointer_ea4f3357a8be53688998cb7f27f3b865)
            .def("state_occupancy_likelihood_computation", method_pointer_15a05927895c516cad12850712c13c6b)
            .def("expectation_step", method_pointer_2695736208d05343809912a569eace13)
            .def("expectation_step", method_pointer_def459f4cd2a5beab74b4f631c0bd319)
            .staticmethod("nb_value_computation")
            .def_readwrite("ident", &::stat_tool::DiscreteParametric::ident)
            .def_readwrite("inf_bound", &::stat_tool::DiscreteParametric::inf_bound)
            .def_readwrite("sup_bound", &::stat_tool::DiscreteParametric::sup_bound)
            .def_readwrite("parameter", &::stat_tool::DiscreteParametric::parameter)
            .def_readwrite("probability", &::stat_tool::DiscreteParametric::probability);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::DiscreteParametric >, std::shared_ptr< class ::stat_tool::Distribution > >();
}