#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Vectors const volatile * get_pointer<class ::stat_tool::Vectors const volatile >(class ::stat_tool::Vectors const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_e00cc711eccf595f84c6952a456d7688()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Vectors::*method_pointer_7325afaa1c79508aa61a5ea7a3eb1802)(int ) = &::stat_tool::Vectors::min_interval_computation;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::Vectors::*method_pointer_0a64824f80f95d7c95e677c5900996a9)(class ::stat_tool::StatError &, int ) const = &::stat_tool::Vectors::extract;
    bool  (::stat_tool::Vectors::*method_pointer_d4a0c706785b560180824a6675b5a254)(class ::stat_tool::StatError &) = &::stat_tool::Vectors::check;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_7e36c2e63fa851f3873c2c91109940a5)(class ::stat_tool::StatError &, int , int ) const = &::stat_tool::Vectors::shift;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_b01568ce752456ab8cd8af6fece78b51)(class ::stat_tool::StatError &, int , double ) const = &::stat_tool::Vectors::shift;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_95853d148357577baa33cff6d5a50044)(class ::stat_tool::StatError &, int , int , enum ::stat_tool::threshold_direction ) const = &::stat_tool::Vectors::thresholding;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_87a4345382b053478c72fb72d1c44ec1)(class ::stat_tool::StatError &, int , double , enum ::stat_tool::threshold_direction ) const = &::stat_tool::Vectors::thresholding;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_efc5029fa8fb5fdd951f2b02cd00e5e5)(class ::stat_tool::StatError &, int , class ::std::vector< int, class ::std::allocator< int > > ) const = &::stat_tool::Vectors::transcode;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_cefde5dc98065d0586c4cac2dc7cccd0)(class ::stat_tool::StatError &, int , int , enum ::stat_tool::rounding ) const = &::stat_tool::Vectors::cluster;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_d643a00055dc585ba022b70b5659e416)(class ::stat_tool::StatError &, int , int , class ::std::vector< int, class ::std::allocator< int > > ) const = &::stat_tool::Vectors::cluster;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_1c8bd1bb6c65517887b5a28a67f5edb7)(class ::stat_tool::StatError &, int , int , class ::std::vector< double, class ::std::allocator< double > > ) const = &::stat_tool::Vectors::cluster;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_687cdb26391b5a42b246993209abc6f6)(class ::stat_tool::StatError &, int , int ) const = &::stat_tool::Vectors::scaling;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_38fd0e7a6bfb5cdcb6ecc3570ac1fba6)(class ::stat_tool::StatError &, int , double ) const = &::stat_tool::Vectors::scaling;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_a1c95f46a3ce50d9a07a6f3c086fe120)(class ::stat_tool::StatError &, int , enum ::stat_tool::rounding ) const = &::stat_tool::Vectors::round;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_7d1da63ddd5c57f8a5d531d6de79078d)(class ::stat_tool::StatError &, bool , int , int , int , bool ) const = &::stat_tool::Vectors::value_select;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_d797bc5036e2503cad783720b1088de5)(class ::stat_tool::StatError &, bool , int , double , double , bool ) const = &::stat_tool::Vectors::value_select;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_de250270454c51d790e48877c46f86d9)(class ::stat_tool::StatError &, int , class ::std::vector< int, class ::std::allocator< int > > , bool ) const = &::stat_tool::Vectors::select_individual;
    class ::stat_tool::Vectors * (::stat_tool::Vectors::*method_pointer_d761498f88915e81986ca3c4386af241)(class ::stat_tool::StatError &, int , class ::std::vector< int, class ::std::allocator< int > > , bool ) const = &::stat_tool::Vectors::select_variable;
    class ::stat_tool::Vectors * (*method_pointer_ff2504e6b0f153a39cb10ab6a6311f23)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) = ::stat_tool::Vectors::ascii_read;
    bool  (::stat_tool::Vectors::*method_pointer_e65a013b8a23512da704fb416d54dc32)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, bool ) const = &::stat_tool::Vectors::ascii_data_write;
    bool  (::stat_tool::Vectors::*method_pointer_ecc04a018f6a5b208ee872074d8f8b75)(class ::stat_tool::StatError &, int , double , double ) = &::stat_tool::Vectors::select_bin_width;
    double  (::stat_tool::Vectors::*method_pointer_26c98b2e5f5c5bb086bea4965a5b78a2)(int , double ) const = &::stat_tool::Vectors::mean_absolute_deviation_computation;
    double  (::stat_tool::Vectors::*method_pointer_d906dd83782c5681a7f3165911f1ea0f)(int ) const = &::stat_tool::Vectors::mean_absolute_difference_computation;
    double  (::stat_tool::Vectors::*method_pointer_1944b7da5bdc599aac220549d5f4e647)(int ) const = &::stat_tool::Vectors::skewness_computation;
    double  (::stat_tool::Vectors::*method_pointer_22883ffe79065e8c9ce97373fad77b65)(int ) const = &::stat_tool::Vectors::kurtosis_computation;
    double  (::stat_tool::Vectors::*method_pointer_e1746706f96a55f2a52200107c7278a8)() const = &::stat_tool::Vectors::spearman_rank_single_correlation_computation;
    double  (::stat_tool::Vectors::*method_pointer_55efe34adebb5ed6bcb06f3a1dacb480)() const = &::stat_tool::Vectors::kendall_rank_single_correlation_computation;
    bool  (::stat_tool::Vectors::*method_pointer_94a5e269c48e54fdbdb94c6846d20e64)(class ::stat_tool::StatError &, bool , enum ::stat_tool::correlation_type , class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) const = &::stat_tool::Vectors::rank_correlation_computation;
    class ::stat_tool::DistanceMatrix * (::stat_tool::Vectors::*method_pointer_dfdcf634f5695d78afda61471d4aae6b)(class ::stat_tool::StatError &, class ::stat_tool::VectorDistance const &, bool ) const = &::stat_tool::Vectors::comparison;
    bool  (::stat_tool::Vectors::*method_pointer_a84fe68ae7ac5af2adc050f704321bb9)(class ::stat_tool::StatError &, bool , int , int , class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, enum ::stat_tool::output_format ) const = &::stat_tool::Vectors::contingency_table;
    bool  (::stat_tool::Vectors::*method_pointer_b5852633431e5a0d987ebaa9360a389b)(class ::stat_tool::StatError &, bool , int , int , int , class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, enum ::stat_tool::output_format ) const = &::stat_tool::Vectors::variance_analysis;
    bool  (::stat_tool::Vectors::*method_pointer_07b80bea06045d9d9ea624eb00268514)(class ::stat_tool::StatError &, bool , class ::stat_tool::Vectors const &) const = &::stat_tool::Vectors::sup_norm_distance;
    class ::stat_tool::Regression * (::stat_tool::Vectors::*method_pointer_0d0afb76c3bc5f4d80335e57a185ead2)(class ::stat_tool::StatError &, int , int ) const = &::stat_tool::Vectors::linear_regression;
    class ::stat_tool::Regression * (::stat_tool::Vectors::*method_pointer_86d0ca671ac2548595161412d53ed757)(class ::stat_tool::StatError &, int , int , class ::stat_tool::Distribution const &, enum ::stat_tool::moving_average_method ) const = &::stat_tool::Vectors::moving_average;
    class ::stat_tool::Regression * (::stat_tool::Vectors::*method_pointer_46b74a48e9645b688643826906daf90c)(class ::stat_tool::StatError &, int , int , double , bool ) const = &::stat_tool::Vectors::nearest_neighbor_smoother;
    class ::stat_tool::Mixture * (::stat_tool::Vectors::*method_pointer_fed11c312f67586fb2c5044a9d8a1c22)(class ::stat_tool::StatError &, bool , class ::stat_tool::Mixture const &, bool , bool , enum ::stat_tool::tying_rule , bool , int ) const = &::stat_tool::Vectors::mixture_estimation;
    class ::stat_tool::Mixture * (::stat_tool::Vectors::*method_pointer_d0677aba3f785bcbbe86a1ad850b225e)(class ::stat_tool::StatError &, bool , int , int , double , double , bool , enum ::stat_tool::tying_rule , bool , int ) const = &::stat_tool::Vectors::mixture_estimation;
    class ::stat_tool::Mixture * (::stat_tool::Vectors::*method_pointer_3e8aef16c1c052eaa8c877b0562b3459)(class ::stat_tool::StatError &, bool , class ::stat_tool::Mixture const &, bool , bool , enum ::stat_tool::tying_rule , int , int , double , bool , int ) const = &::stat_tool::Vectors::mixture_stochastic_estimation;
    class ::stat_tool::Mixture * (::stat_tool::Vectors::*method_pointer_f6218f54c6745899ab677bfe4385f6fb)(class ::stat_tool::StatError &, bool , int , int , double , double , bool , enum ::stat_tool::tying_rule , int , int , double , bool , int ) const = &::stat_tool::Vectors::mixture_stochastic_estimation;
    int  (::stat_tool::Vectors::*method_pointer_99bc1d3ca6a7569cb81ab8af2e6327da)() const = &::stat_tool::Vectors::get_nb_vector;
    int  (::stat_tool::Vectors::*method_pointer_4bf15c5573385ba787ef95b022e982a0)(int ) const = &::stat_tool::Vectors::get_identifier;
    int  (::stat_tool::Vectors::*method_pointer_e773df74b82353cdac053cb8a4c9b1b4)() const = &::stat_tool::Vectors::get_nb_variable;
    enum ::stat_tool::variable_nature  (::stat_tool::Vectors::*method_pointer_c8d2df673930525a970860ee1ec19b0a)(int ) const = &::stat_tool::Vectors::get_type;
    double  (::stat_tool::Vectors::*method_pointer_fdd06df80540534e83deee2a175b9242)(int ) const = &::stat_tool::Vectors::get_min_value;
    double  (::stat_tool::Vectors::*method_pointer_962d04415e265ca9855000b21c89ee14)(int ) const = &::stat_tool::Vectors::get_max_value;
    class ::stat_tool::FrequencyDistribution * (::stat_tool::Vectors::*method_pointer_1021e618f6165332b14acb62f85c5d25)(int ) const = &::stat_tool::Vectors::get_marginal_distribution;
    class ::stat_tool::Histogram * (::stat_tool::Vectors::*method_pointer_7f7c9bf629565f34931a0f605a4acb01)(int ) const = &::stat_tool::Vectors::get_marginal_histogram;
    double  (::stat_tool::Vectors::*method_pointer_ee71d48f3565591b931ca8b751a71252)(int ) const = &::stat_tool::Vectors::get_mean;
    double  (::stat_tool::Vectors::*method_pointer_8701e1e47b125c18921a4992c06dc1f4)(int , int ) const = &::stat_tool::Vectors::get_covariance;
    int  (::stat_tool::Vectors::*method_pointer_662805d7237158d187f7c504808c38a5)(int , int ) const = &::stat_tool::Vectors::get_int_vector;
    double  (::stat_tool::Vectors::*method_pointer_fe52e7dfa3fa5c2894a06029ba3bf75e)(int , int ) const = &::stat_tool::Vectors::get_real_vector;
    boost::python::class_< class ::stat_tool::Vectors, autowig::Held< class ::stat_tool::Vectors >::Type, boost::python::bases< class ::stat_tool::StatInterface > > class_e00cc711eccf595f84c6952a456d7688("Vectors", "Vectors with integer- and real-valued variables\n\n", boost::python::no_init);
    class_e00cc711eccf595f84c6952a456d7688.def(boost::python::init<  >(""));
    class_e00cc711eccf595f84c6952a456d7688.def(boost::python::init< class ::stat_tool::Vectors const &, int , enum ::stat_tool::variable_nature  >(""));
    class_e00cc711eccf595f84c6952a456d7688.def(boost::python::init< class ::stat_tool::Vectors const &, enum ::stat_tool::vector_transformation  >(""));
    class_e00cc711eccf595f84c6952a456d7688.def("min_interval_computation", method_pointer_7325afaa1c79508aa61a5ea7a3eb1802, "");
    class_e00cc711eccf595f84c6952a456d7688.def("extract", method_pointer_0a64824f80f95d7c95e677c5900996a9, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("check", method_pointer_d4a0c706785b560180824a6675b5a254, "");
    class_e00cc711eccf595f84c6952a456d7688.def("shift", method_pointer_7e36c2e63fa851f3873c2c91109940a5, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("shift", method_pointer_b01568ce752456ab8cd8af6fece78b51, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("thresholding", method_pointer_95853d148357577baa33cff6d5a50044, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("thresholding", method_pointer_87a4345382b053478c72fb72d1c44ec1, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("transcode", method_pointer_efc5029fa8fb5fdd951f2b02cd00e5e5, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("cluster", method_pointer_cefde5dc98065d0586c4cac2dc7cccd0, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("cluster", method_pointer_d643a00055dc585ba022b70b5659e416, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("cluster", method_pointer_1c8bd1bb6c65517887b5a28a67f5edb7, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("scaling", method_pointer_687cdb26391b5a42b246993209abc6f6, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("scaling", method_pointer_38fd0e7a6bfb5cdcb6ecc3570ac1fba6, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("round", method_pointer_a1c95f46a3ce50d9a07a6f3c086fe120, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("value_select", method_pointer_7d1da63ddd5c57f8a5d531d6de79078d, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("value_select", method_pointer_d797bc5036e2503cad783720b1088de5, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("select_individual", method_pointer_de250270454c51d790e48877c46f86d9, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("select_variable", method_pointer_d761498f88915e81986ca3c4386af241, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("ascii_read", method_pointer_ff2504e6b0f153a39cb10ab6a6311f23, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("ascii_data_write", method_pointer_e65a013b8a23512da704fb416d54dc32, "");
    class_e00cc711eccf595f84c6952a456d7688.def("select_bin_width", method_pointer_ecc04a018f6a5b208ee872074d8f8b75, "");
    class_e00cc711eccf595f84c6952a456d7688.def("mean_absolute_deviation_computation", method_pointer_26c98b2e5f5c5bb086bea4965a5b78a2, "");
    class_e00cc711eccf595f84c6952a456d7688.def("mean_absolute_difference_computation", method_pointer_d906dd83782c5681a7f3165911f1ea0f, "");
    class_e00cc711eccf595f84c6952a456d7688.def("skewness_computation", method_pointer_1944b7da5bdc599aac220549d5f4e647, "");
    class_e00cc711eccf595f84c6952a456d7688.def("kurtosis_computation", method_pointer_22883ffe79065e8c9ce97373fad77b65, "");
    class_e00cc711eccf595f84c6952a456d7688.def("spearman_rank_single_correlation_computation", method_pointer_e1746706f96a55f2a52200107c7278a8, "");
    class_e00cc711eccf595f84c6952a456d7688.def("kendall_rank_single_correlation_computation", method_pointer_55efe34adebb5ed6bcb06f3a1dacb480, "");
    class_e00cc711eccf595f84c6952a456d7688.def("rank_correlation_computation", method_pointer_94a5e269c48e54fdbdb94c6846d20e64, "");
    class_e00cc711eccf595f84c6952a456d7688.def("comparison", method_pointer_dfdcf634f5695d78afda61471d4aae6b, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("contingency_table", method_pointer_a84fe68ae7ac5af2adc050f704321bb9, "");
    class_e00cc711eccf595f84c6952a456d7688.def("variance_analysis", method_pointer_b5852633431e5a0d987ebaa9360a389b, "");
    class_e00cc711eccf595f84c6952a456d7688.def("sup_norm_distance", method_pointer_07b80bea06045d9d9ea624eb00268514, "");
    class_e00cc711eccf595f84c6952a456d7688.def("linear_regression", method_pointer_0d0afb76c3bc5f4d80335e57a185ead2, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("moving_average", method_pointer_86d0ca671ac2548595161412d53ed757, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("nearest_neighbor_smoother", method_pointer_46b74a48e9645b688643826906daf90c, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("mixture_estimation", method_pointer_fed11c312f67586fb2c5044a9d8a1c22, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("mixture_estimation", method_pointer_d0677aba3f785bcbbe86a1ad850b225e, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("mixture_stochastic_estimation", method_pointer_3e8aef16c1c052eaa8c877b0562b3459, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("mixture_stochastic_estimation", method_pointer_f6218f54c6745899ab677bfe4385f6fb, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_nb_vector", method_pointer_99bc1d3ca6a7569cb81ab8af2e6327da, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_identifier", method_pointer_4bf15c5573385ba787ef95b022e982a0, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_nb_variable", method_pointer_e773df74b82353cdac053cb8a4c9b1b4, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_type", method_pointer_c8d2df673930525a970860ee1ec19b0a, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_min_value", method_pointer_fdd06df80540534e83deee2a175b9242, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_max_value", method_pointer_962d04415e265ca9855000b21c89ee14, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_marginal_distribution", method_pointer_1021e618f6165332b14acb62f85c5d25, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_marginal_histogram", method_pointer_7f7c9bf629565f34931a0f605a4acb01, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_mean", method_pointer_ee71d48f3565591b931ca8b751a71252, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_covariance", method_pointer_8701e1e47b125c18921a4992c06dc1f4, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_int_vector", method_pointer_662805d7237158d187f7c504808c38a5, "");
    class_e00cc711eccf595f84c6952a456d7688.def("get_real_vector", method_pointer_fe52e7dfa3fa5c2894a06029ba3bf75e, "");
    class_e00cc711eccf595f84c6952a456d7688.staticmethod("ascii_read");

    if(autowig::Held< class ::stat_tool::Vectors >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::Vectors >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        //boost::python::objects::class_value_wrapper< autowig::Held< class ::stat_tool::Vectors >::Type, boost::python::objects::make_ptr_instance< class ::stat_tool::Vectors, boost::python::objects::pointer_holder< autowig::Held< class ::stat_tool::Vectors >::Type, class ::stat_tool::Vectors > > >();
    }

}