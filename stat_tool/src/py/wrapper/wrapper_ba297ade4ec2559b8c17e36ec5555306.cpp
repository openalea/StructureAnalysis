#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Distribution const volatile * get_pointer<class ::stat_tool::Distribution const volatile >(class ::stat_tool::Distribution const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_ba297ade4ec2559b8c17e36ec5555306()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Distribution::*method_pointer_20e207dfe4fa5b2babd36ed8759d2f50)(class ::stat_tool::Distribution const &, int ) = &::stat_tool::Distribution::mass_copy;
    void  (::stat_tool::Distribution::*method_pointer_8241fcfc781c5c28b73bfb49118270eb)(class ::stat_tool::Distribution const &) = &::stat_tool::Distribution::equal_size_copy;
    void  (::stat_tool::Distribution::*method_pointer_7861e1796ec65c5ea4e6e3d3662e93a0)(int ) = &::stat_tool::Distribution::init;
    void  (::stat_tool::Distribution::*method_pointer_1a02716b6a1e5933a2d1c01e5bdc9787)(class ::stat_tool::Distribution const &, int ) = &::stat_tool::Distribution::copy;
    void  (::stat_tool::Distribution::*method_pointer_e1a38282387e5b098062256c24c48250)(class ::stat_tool::Distribution const &) = &::stat_tool::Distribution::normalization_copy;
    bool  (::stat_tool::Distribution::*method_pointer_5de0c8c69ef25a18bfb2dca9f5476dbe)(class ::stat_tool::Distribution const &) const = &::stat_tool::Distribution::operator==;
    bool  (::stat_tool::Distribution::*method_pointer_6f4fd7c01d5952e6ae30a1c2c26a9986)(class ::stat_tool::Distribution const &) const = &::stat_tool::Distribution::operator!=;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Distribution::*method_pointer_f51b0b6dc37f55b3a3b05517aa8cf59f)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , bool ) const = &::stat_tool::Distribution::ascii_characteristic_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Distribution::*method_pointer_277767f0416e53a0b5a065106bdb18f6)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , bool , bool , class ::stat_tool::FrequencyDistribution const *) const = &::stat_tool::Distribution::ascii_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Distribution::*method_pointer_d634a4260bda5d7e970eb0982b72a64a)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::Distribution::print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Distribution::*method_pointer_7a6a3e16802d50c58fdec748505c3fd3)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::Distribution::spreadsheet_characteristic_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Distribution::*method_pointer_362ca0374f375980ad00b2f4884b669b)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , bool , bool , class ::stat_tool::FrequencyDistribution const *) const = &::stat_tool::Distribution::spreadsheet_print;
    int  (::stat_tool::Distribution::*method_pointer_bf5ecf127084544b8020df64a53b8af6)(class ::stat_tool::FrequencyDistribution const *) const = &::stat_tool::Distribution::plot_nb_value_computation;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Distribution::*method_pointer_c148df8ed35e591d89cc61426bde8f80)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::Distribution::plot_title_print;
    void  (::stat_tool::Distribution::*method_pointer_1611ec32975b592ab99017cd3e9155dc)(class ::stat_tool::SinglePlot &, double ) const = &::stat_tool::Distribution::plotable_mass_write;
    void  (::stat_tool::Distribution::*method_pointer_cd7d6da0b27c5c28a19645b64dcf3d4b)(class ::stat_tool::SinglePlot &) const = &::stat_tool::Distribution::plotable_cumul_write;
    void  (::stat_tool::Distribution::*method_pointer_024ea3f353af53b49d19e8e16f5c81d6)(class ::stat_tool::SinglePlot &, class ::stat_tool::Distribution const &) const = &::stat_tool::Distribution::plotable_cumul_matching_write;
    void  (::stat_tool::Distribution::*method_pointer_7fb5ba57b0bd557294ba67c9df8632a1)(class ::stat_tool::SinglePlot &) const = &::stat_tool::Distribution::plotable_concentration_write;
    void  (::stat_tool::Distribution::*method_pointer_639913b5ea6450be9ffe34d7cd00519d)(class ::stat_tool::SinglePlot &) const = &::stat_tool::Distribution::plotable_survivor_write;
    ::stat_tool::MultiPlotSet * (::stat_tool::Distribution::*method_pointer_3152bf3e508658348d0279533cb2e149)() const = &::stat_tool::Distribution::get_plotable;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Distribution::*method_pointer_97bb1de239bf5a11bff913b13f26f1ab)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::Distribution::survival_ascii_write;
    class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > >  (::stat_tool::Distribution::*method_pointer_02ad13e5a4df5dd4979523e393819cf5)() const = &::stat_tool::Distribution::survival_ascii_write;
    bool  (::stat_tool::Distribution::*method_pointer_29d745aa02ad5c6c9dc87e670cd2bc49)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) const = &::stat_tool::Distribution::survival_ascii_write;
    bool  (::stat_tool::Distribution::*method_pointer_f0655161c09e5108bac7ffac96d0d04e)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) const = &::stat_tool::Distribution::survival_spreadsheet_write;
    ::stat_tool::MultiPlotSet * (::stat_tool::Distribution::*method_pointer_d4f02c25f9d75770a8a65d9a6542b45a)(class ::stat_tool::StatError &) const = &::stat_tool::Distribution::survival_get_plotable;
    void  (::stat_tool::Distribution::*method_pointer_953ee629873855a290beb8855d09104c)() = &::stat_tool::Distribution::max_computation;
    double  (::stat_tool::Distribution::*method_pointer_51ad43af0274525bbe63fee07acc29c5)() const = &::stat_tool::Distribution::mode_computation;
    void  (::stat_tool::Distribution::*method_pointer_a537bbe5a6ae55b6a8016e8b9067acb3)() = &::stat_tool::Distribution::mean_computation;
    double  (::stat_tool::Distribution::*method_pointer_22ebb6c8e915553c8a18a6a52eec23f5)(double ) const = &::stat_tool::Distribution::quantile_computation;
    void  (::stat_tool::Distribution::*method_pointer_57ba5f4689b45904b44c587e93dc15c0)() = &::stat_tool::Distribution::variance_computation;
    void  (::stat_tool::Distribution::*method_pointer_2001ad55df8d52adb76fc96030206f1e)() = &::stat_tool::Distribution::nb_value_computation;
    void  (::stat_tool::Distribution::*method_pointer_ad1911b42c4d59a9bcd8b807891de938)() = &::stat_tool::Distribution::offset_computation;
    double  (::stat_tool::Distribution::*method_pointer_6dbaddf1f1155bb1b8bc3fed76147286)() const = &::stat_tool::Distribution::concentration_computation;
    double  (::stat_tool::Distribution::*method_pointer_abccde0ab0b65b6f9141021ff8800c37)(double ) const = &::stat_tool::Distribution::mean_absolute_deviation_computation;
    double  (::stat_tool::Distribution::*method_pointer_4deb6d8102835bf2a36830a3a0f5def6)() const = &::stat_tool::Distribution::skewness_computation;
    double  (::stat_tool::Distribution::*method_pointer_8e61b763068c50deb5d19ba5f2436015)() const = &::stat_tool::Distribution::kurtosis_computation;
    double  (::stat_tool::Distribution::*method_pointer_54dd5bfbf7475a8c9fc4a6cec745a4b9)() const = &::stat_tool::Distribution::information_computation;
    double  (::stat_tool::Distribution::*method_pointer_9b858f2cc01d5741a6a964db52a07be8)() const = &::stat_tool::Distribution::first_difference_norm_computation;
    double  (::stat_tool::Distribution::*method_pointer_1d3db396096a50529a2f6089d240114a)() const = &::stat_tool::Distribution::second_difference_norm_computation;
    void  (::stat_tool::Distribution::*method_pointer_942b4b2dbe0150daa8b9f041f9977451)() = &::stat_tool::Distribution::cumul_computation;
    double  (::stat_tool::Distribution::*method_pointer_6e70e9f9f8145bb887d3170890f185df)(class ::stat_tool::Distribution const &) const = &::stat_tool::Distribution::overlap_distance_computation;
    void  (::stat_tool::Distribution::*method_pointer_5550dd3da9575a6a9d3ea3a25be23a33)() = &::stat_tool::Distribution::log_computation;
    double  (::stat_tool::Distribution::*method_pointer_09855248fce95cbab1ab68194f1aa161)(class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::Distribution::survivor_likelihood_computation;
    double  (::stat_tool::Distribution::*method_pointer_943efd6d629457b0b3517b8554783f5f)(class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::Distribution::chi2_value_computation;
    void  (::stat_tool::Distribution::*method_pointer_51a1252546c85dbfb77553dad869271e)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Test &) const = &::stat_tool::Distribution::chi2_degree_of_freedom;
    double  (::stat_tool::Distribution::*method_pointer_ef0baa45f93f5c24874234cea8b1123a)(class ::stat_tool::Reestimation< int > const &) const = &::stat_tool::Distribution::likelihood_computation;
    double  (::stat_tool::Distribution::*method_pointer_e57d2d73cda251c0b53bdab57f19a016)(class ::stat_tool::Reestimation< double > const &) const = &::stat_tool::Distribution::likelihood_computation;
    void  (::stat_tool::Distribution::*method_pointer_3b525dc4ae96537eb00d73abeb8ec6b4)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Test &) const = &::stat_tool::Distribution::chi2_fit;
    void  (::stat_tool::Distribution::*method_pointer_99e6e2ea1dba5aa1bd8782a37e564723)(class ::stat_tool::Distribution &, class ::stat_tool::Distribution &, int ) = &::stat_tool::Distribution::convolution;
    int  (::stat_tool::Distribution::*method_pointer_8baa16edfb3f5100ad7a4104e8ec0140)() const = &::stat_tool::Distribution::simulation;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::Distribution::*method_pointer_9c6f4f3ac7e255d8981b418b2aae0b49)(class ::stat_tool::StatError &, int ) const = &::stat_tool::Distribution::truncate;
    boost::python::class_< class ::stat_tool::Distribution, autowig::Held< class ::stat_tool::Distribution >::Type > class_ba297ade4ec2559b8c17e36ec5555306("Distribution", "Discrete distribution\n\n", boost::python::no_init);
    class_ba297ade4ec2559b8c17e36ec5555306.def(boost::python::init< int  >(""));
    class_ba297ade4ec2559b8c17e36ec5555306.def(boost::python::init< class ::stat_tool::Distribution const &, double  >(""));
    class_ba297ade4ec2559b8c17e36ec5555306.def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >(""));
    class_ba297ade4ec2559b8c17e36ec5555306.def(boost::python::init< class ::stat_tool::Distribution const &, enum ::stat_tool::distribution_transformation , int  >(""));
    class_ba297ade4ec2559b8c17e36ec5555306.def("mass_copy", method_pointer_20e207dfe4fa5b2babd36ed8759d2f50, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("equal_size_copy", method_pointer_8241fcfc781c5c28b73bfb49118270eb, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("init", method_pointer_7861e1796ec65c5ea4e6e3d3662e93a0, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("copy", method_pointer_1a02716b6a1e5933a2d1c01e5bdc9787, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("normalization_copy", method_pointer_e1a38282387e5b098062256c24c48250, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("__eq__", method_pointer_5de0c8c69ef25a18bfb2dca9f5476dbe, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("__neq__", method_pointer_6f4fd7c01d5952e6ae30a1c2c26a9986, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("ascii_characteristic_print", method_pointer_f51b0b6dc37f55b3a3b05517aa8cf59f, boost::python::return_internal_reference<>(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("ascii_print", method_pointer_277767f0416e53a0b5a065106bdb18f6, boost::python::return_internal_reference<>(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("print", method_pointer_d634a4260bda5d7e970eb0982b72a64a, boost::python::return_internal_reference<>(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("spreadsheet_characteristic_print", method_pointer_7a6a3e16802d50c58fdec748505c3fd3, boost::python::return_internal_reference<>(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("spreadsheet_print", method_pointer_362ca0374f375980ad00b2f4884b669b, boost::python::return_internal_reference<>(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("plot_nb_value_computation", method_pointer_bf5ecf127084544b8020df64a53b8af6, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("plot_title_print", method_pointer_c148df8ed35e591d89cc61426bde8f80, boost::python::return_internal_reference<>(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("plotable_mass_write", method_pointer_1611ec32975b592ab99017cd3e9155dc, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("plotable_cumul_write", method_pointer_cd7d6da0b27c5c28a19645b64dcf3d4b, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("plotable_cumul_matching_write", method_pointer_024ea3f353af53b49d19e8e16f5c81d6, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("plotable_concentration_write", method_pointer_7fb5ba57b0bd557294ba67c9df8632a1, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("plotable_survivor_write", method_pointer_639913b5ea6450be9ffe34d7cd00519d, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("get_plotable", method_pointer_3152bf3e508658348d0279533cb2e149, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("survival_ascii_write", method_pointer_97bb1de239bf5a11bff913b13f26f1ab, boost::python::return_internal_reference<>(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("survival_ascii_write", method_pointer_02ad13e5a4df5dd4979523e393819cf5, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("survival_ascii_write", method_pointer_29d745aa02ad5c6c9dc87e670cd2bc49, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("survival_spreadsheet_write", method_pointer_f0655161c09e5108bac7ffac96d0d04e, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("survival_get_plotable", method_pointer_d4f02c25f9d75770a8a65d9a6542b45a, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("max_computation", method_pointer_953ee629873855a290beb8855d09104c, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("mode_computation", method_pointer_51ad43af0274525bbe63fee07acc29c5, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("mean_computation", method_pointer_a537bbe5a6ae55b6a8016e8b9067acb3, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("quantile_computation", method_pointer_22ebb6c8e915553c8a18a6a52eec23f5, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("variance_computation", method_pointer_57ba5f4689b45904b44c587e93dc15c0, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("nb_value_computation", method_pointer_2001ad55df8d52adb76fc96030206f1e, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("offset_computation", method_pointer_ad1911b42c4d59a9bcd8b807891de938, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("concentration_computation", method_pointer_6dbaddf1f1155bb1b8bc3fed76147286, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("mean_absolute_deviation_computation", method_pointer_abccde0ab0b65b6f9141021ff8800c37, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("skewness_computation", method_pointer_4deb6d8102835bf2a36830a3a0f5def6, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("kurtosis_computation", method_pointer_8e61b763068c50deb5d19ba5f2436015, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("information_computation", method_pointer_54dd5bfbf7475a8c9fc4a6cec745a4b9, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("first_difference_norm_computation", method_pointer_9b858f2cc01d5741a6a964db52a07be8, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("second_difference_norm_computation", method_pointer_1d3db396096a50529a2f6089d240114a, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("cumul_computation", method_pointer_942b4b2dbe0150daa8b9f041f9977451, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("overlap_distance_computation", method_pointer_6e70e9f9f8145bb887d3170890f185df, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("log_computation", method_pointer_5550dd3da9575a6a9d3ea3a25be23a33, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("survivor_likelihood_computation", method_pointer_09855248fce95cbab1ab68194f1aa161, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("chi_2__value_computation", method_pointer_943efd6d629457b0b3517b8554783f5f, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("chi_2__degree_of_freedom", method_pointer_51a1252546c85dbfb77553dad869271e, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("likelihood_computation", method_pointer_ef0baa45f93f5c24874234cea8b1123a, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("likelihood_computation", method_pointer_e57d2d73cda251c0b53bdab57f19a016, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("chi_2__fit", method_pointer_3b525dc4ae96537eb00d73abeb8ec6b4, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("convolution", method_pointer_99e6e2ea1dba5aa1bd8782a37e564723, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("simulation", method_pointer_8baa16edfb3f5100ad7a4104e8ec0140, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def("truncate", method_pointer_9c6f4f3ac7e255d8981b418b2aae0b49, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_ba297ade4ec2559b8c17e36ec5555306.def_readwrite("nb_value", &::stat_tool::Distribution::nb_value, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def_readwrite("alloc_nb_value", &::stat_tool::Distribution::alloc_nb_value, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def_readwrite("offset", &::stat_tool::Distribution::offset, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def_readwrite("max", &::stat_tool::Distribution::max, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def_readwrite("complement", &::stat_tool::Distribution::complement, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def_readwrite("mean", &::stat_tool::Distribution::mean, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def_readwrite("variance", &::stat_tool::Distribution::variance, "");
    class_ba297ade4ec2559b8c17e36ec5555306.def_readwrite("nb_parameter", &::stat_tool::Distribution::nb_parameter, "");

}