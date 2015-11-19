#include <boost/python.hpp>
#include <stat_tool/compound.h>
#include <stat_tool/multivariate_mixture.h>
#include <stat_tool/regression.h>
#include <stat_tool/convolution.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/mixture.h>

void _stat_tool_vectors()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::Vectors::*method_pointer_7325afaa1c79508aa61a5ea7a3eb1802)(int) = &::stat_tool::Vectors::min_interval_computation;
        class ::stat_tool::DiscreteDistributionData * (::stat_tool::Vectors::*method_pointer_0a64824f80f95d7c95e677c5900996a9)(class ::stat_tool::StatError &, int) const = &::stat_tool::Vectors::extract;
        bool (::stat_tool::Vectors::*method_pointer_d4a0c706785b560180824a6675b5a254)(class ::stat_tool::StatError &) = &::stat_tool::Vectors::check;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_7e36c2e63fa851f3873c2c91109940a5)(class ::stat_tool::StatError &, int, int) const = &::stat_tool::Vectors::shift;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_b01568ce752456ab8cd8af6fece78b51)(class ::stat_tool::StatError &, int, double) const = &::stat_tool::Vectors::shift;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_95853d148357577baa33cff6d5a50044)(class ::stat_tool::StatError &, int, int, enum ::stat_tool::threshold_direction) const = &::stat_tool::Vectors::thresholding;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_87a4345382b053478c72fb72d1c44ec1)(class ::stat_tool::StatError &, int, double, enum ::stat_tool::threshold_direction) const = &::stat_tool::Vectors::thresholding;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_cefde5dc98065d0586c4cac2dc7cccd0)(class ::stat_tool::StatError &, int, int, enum ::stat_tool::rounding) const = &::stat_tool::Vectors::cluster;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_687cdb26391b5a42b246993209abc6f6)(class ::stat_tool::StatError &, int, int) const = &::stat_tool::Vectors::scaling;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_38fd0e7a6bfb5cdcb6ecc3570ac1fba6)(class ::stat_tool::StatError &, int, double) const = &::stat_tool::Vectors::scaling;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_a1c95f46a3ce50d9a07a6f3c086fe120)(class ::stat_tool::StatError &, int, enum ::stat_tool::rounding) const = &::stat_tool::Vectors::round;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_3851b5934dd1541ea4e5a99df6faa717)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, int, int, int, bool) const = &::stat_tool::Vectors::value_select;
        class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_bd307c6b891f584f9831b5fa2a5df28f)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, int, double, double, bool) const = &::stat_tool::Vectors::value_select;
        bool (::stat_tool::Vectors::*method_pointer_9eaaa5f3017e581299fb0f28342bc777)(class ::stat_tool::StatError &, int, double, double) = &::stat_tool::Vectors::select_step;
        double (::stat_tool::Vectors::*method_pointer_cf538072d5af5e2dbe926ceab38d9e00)(int) const = &::stat_tool::Vectors::mean_absolute_deviation_computation;
        double (::stat_tool::Vectors::*method_pointer_d906dd83782c5681a7f3165911f1ea0f)(int) const = &::stat_tool::Vectors::mean_absolute_difference_computation;
        double (::stat_tool::Vectors::*method_pointer_1944b7da5bdc599aac220549d5f4e647)(int) const = &::stat_tool::Vectors::skewness_computation;
        double (::stat_tool::Vectors::*method_pointer_22883ffe79065e8c9ce97373fad77b65)(int) const = &::stat_tool::Vectors::kurtosis_computation;
        double (::stat_tool::Vectors::*method_pointer_e1746706f96a55f2a52200107c7278a8)() const = &::stat_tool::Vectors::spearman_rank_single_correlation_computation;
        double (::stat_tool::Vectors::*method_pointer_55efe34adebb5ed6bcb06f3a1dacb480)() const = &::stat_tool::Vectors::kendall_rank_single_correlation_computation;
        class ::stat_tool::DistanceMatrix * (::stat_tool::Vectors::*method_pointer_dfdcf634f5695d78afda61471d4aae6b)(class ::stat_tool::StatError &, class ::stat_tool::VectorDistance const &, bool) const = &::stat_tool::Vectors::comparison;
        class ::stat_tool::Regression * (::stat_tool::Vectors::*method_pointer_0d0afb76c3bc5f4d80335e57a185ead2)(class ::stat_tool::StatError &, int, int) const = &::stat_tool::Vectors::linear_regression;
        class ::stat_tool::Regression * (::stat_tool::Vectors::*method_pointer_86d0ca671ac2548595161412d53ed757)(class ::stat_tool::StatError &, int, int, class ::stat_tool::Distribution const &, enum ::stat_tool::moving_average_method) const = &::stat_tool::Vectors::moving_average;
        class ::stat_tool::Regression * (::stat_tool::Vectors::*method_pointer_46b74a48e9645b688643826906daf90c)(class ::stat_tool::StatError &, int, int, double, bool) const = &::stat_tool::Vectors::nearest_neighbor_smoother;
        class ::stat_tool::Mixture * (::stat_tool::Vectors::*method_pointer_52f6ab03d2235abdac1eb4a9149dbf5e)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::Mixture const &, bool, bool, enum ::stat_tool::tying_rule, bool, int) const = &::stat_tool::Vectors::mixture_estimation;
        class ::stat_tool::Mixture * (::stat_tool::Vectors::*method_pointer_f2e551884a7b5eb1a5fc786b6751a45b)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, int, int, double, double, bool, enum ::stat_tool::tying_rule, bool, int) const = &::stat_tool::Vectors::mixture_estimation;
        class ::stat_tool::Mixture * (::stat_tool::Vectors::*method_pointer_48c22638372c50c5a09a496acb892022)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, class ::stat_tool::Mixture const &, bool, bool, enum ::stat_tool::tying_rule, int, int, double, bool, int) const = &::stat_tool::Vectors::mixture_stochastic_estimation;
        class ::stat_tool::Mixture * (::stat_tool::Vectors::*method_pointer_6893d541ef945da9be04cbab7fac6dfb)(class ::stat_tool::StatError &, class ::std::basic_ostream<char, std::char_traits<char> > &, int, int, double, double, bool, enum ::stat_tool::tying_rule, int, int, double, bool, int) const = &::stat_tool::Vectors::mixture_stochastic_estimation;
        int (::stat_tool::Vectors::*method_pointer_99bc1d3ca6a7569cb81ab8af2e6327da)() const = &::stat_tool::Vectors::get_nb_vector;
        int (::stat_tool::Vectors::*method_pointer_4bf15c5573385ba787ef95b022e982a0)(int) const = &::stat_tool::Vectors::get_identifier;
        int (::stat_tool::Vectors::*method_pointer_e773df74b82353cdac053cb8a4c9b1b4)() const = &::stat_tool::Vectors::get_nb_variable;
        enum ::stat_tool::variable_nature (::stat_tool::Vectors::*method_pointer_c8d2df673930525a970860ee1ec19b0a)(int) const = &::stat_tool::Vectors::get_type;
        double (::stat_tool::Vectors::*method_pointer_fdd06df80540534e83deee2a175b9242)(int) const = &::stat_tool::Vectors::get_min_value;
        double (::stat_tool::Vectors::*method_pointer_962d04415e265ca9855000b21c89ee14)(int) const = &::stat_tool::Vectors::get_max_value;
        class ::stat_tool::FrequencyDistribution * (::stat_tool::Vectors::*method_pointer_1021e618f6165332b14acb62f85c5d25)(int) const = &::stat_tool::Vectors::get_marginal_distribution;
        class ::stat_tool::Histogram * (::stat_tool::Vectors::*method_pointer_7f7c9bf629565f34931a0f605a4acb01)(int) const = &::stat_tool::Vectors::get_marginal_histogram;
        double (::stat_tool::Vectors::*method_pointer_ee71d48f3565591b931ca8b751a71252)(int) const = &::stat_tool::Vectors::get_mean;
        double (::stat_tool::Vectors::*method_pointer_8701e1e47b125c18921a4992c06dc1f4)(int, int) const = &::stat_tool::Vectors::get_covariance;
        int (::stat_tool::Vectors::*method_pointer_662805d7237158d187f7c504808c38a5)(int, int) const = &::stat_tool::Vectors::get_int_vector;
        double (::stat_tool::Vectors::*method_pointer_fe52e7dfa3fa5c2894a06029ba3bf75e)(int, int) const = &::stat_tool::Vectors::get_real_vector;
        boost::python::class_< class ::stat_tool::Vectors, std::shared_ptr< class ::stat_tool::Vectors >, boost::python::bases< class ::stat_tool::StatInterface > >("Vectors", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::Vectors const &, int, enum ::stat_tool::variable_nature >())
            .def(boost::python::init< class ::stat_tool::Vectors const &, enum ::stat_tool::vector_transformation >())
            .def("min_interval_computation", method_pointer_7325afaa1c79508aa61a5ea7a3eb1802)
            .def("extract", method_pointer_0a64824f80f95d7c95e677c5900996a9, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("check", method_pointer_d4a0c706785b560180824a6675b5a254)
            .def("shift", method_pointer_7e36c2e63fa851f3873c2c91109940a5, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("shift", method_pointer_b01568ce752456ab8cd8af6fece78b51, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("thresholding", method_pointer_95853d148357577baa33cff6d5a50044, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("thresholding", method_pointer_87a4345382b053478c72fb72d1c44ec1, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("cluster", method_pointer_cefde5dc98065d0586c4cac2dc7cccd0, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("scaling", method_pointer_687cdb26391b5a42b246993209abc6f6, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("scaling", method_pointer_38fd0e7a6bfb5cdcb6ecc3570ac1fba6, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("round", method_pointer_a1c95f46a3ce50d9a07a6f3c086fe120, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("value_select", method_pointer_3851b5934dd1541ea4e5a99df6faa717, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("value_select", method_pointer_bd307c6b891f584f9831b5fa2a5df28f, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("select_step", method_pointer_9eaaa5f3017e581299fb0f28342bc777)
            .def("mean_absolute_deviation_computation", method_pointer_cf538072d5af5e2dbe926ceab38d9e00)
            .def("mean_absolute_difference_computation", method_pointer_d906dd83782c5681a7f3165911f1ea0f)
            .def("skewness_computation", method_pointer_1944b7da5bdc599aac220549d5f4e647)
            .def("kurtosis_computation", method_pointer_22883ffe79065e8c9ce97373fad77b65)
            .def("spearman_rank_single_correlation_computation", method_pointer_e1746706f96a55f2a52200107c7278a8)
            .def("kendall_rank_single_correlation_computation", method_pointer_55efe34adebb5ed6bcb06f3a1dacb480)
            .def("comparison", method_pointer_dfdcf634f5695d78afda61471d4aae6b, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("linear_regression", method_pointer_0d0afb76c3bc5f4d80335e57a185ead2, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("moving_average", method_pointer_86d0ca671ac2548595161412d53ed757, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("nearest_neighbor_smoother", method_pointer_46b74a48e9645b688643826906daf90c, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("mixture_estimation", method_pointer_52f6ab03d2235abdac1eb4a9149dbf5e, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("mixture_estimation", method_pointer_f2e551884a7b5eb1a5fc786b6751a45b, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("mixture_stochastic_estimation", method_pointer_48c22638372c50c5a09a496acb892022, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("mixture_stochastic_estimation", method_pointer_6893d541ef945da9be04cbab7fac6dfb, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_nb_vector", method_pointer_99bc1d3ca6a7569cb81ab8af2e6327da)
            .def("get_identifier", method_pointer_4bf15c5573385ba787ef95b022e982a0)
            .def("get_nb_variable", method_pointer_e773df74b82353cdac053cb8a4c9b1b4)
            .def("get_type", method_pointer_c8d2df673930525a970860ee1ec19b0a)
            .def("get_min_value", method_pointer_fdd06df80540534e83deee2a175b9242)
            .def("get_max_value", method_pointer_962d04415e265ca9855000b21c89ee14)
            .def("get_marginal_distribution", method_pointer_1021e618f6165332b14acb62f85c5d25, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_marginal_histogram", method_pointer_7f7c9bf629565f34931a0f605a4acb01, boost::python::return_value_policy< boost::python::reference_existing_object >())
            .def("get_mean", method_pointer_ee71d48f3565591b931ca8b751a71252)
            .def("get_covariance", method_pointer_8701e1e47b125c18921a4992c06dc1f4)
            .def("get_int_vector", method_pointer_662805d7237158d187f7c504808c38a5)
            .def("get_real_vector", method_pointer_fe52e7dfa3fa5c2894a06029ba3bf75e);
        boost::python::implicitly_convertible< std::shared_ptr< class ::stat_tool::Vectors >, std::shared_ptr< class ::stat_tool::StatInterface > >();
}