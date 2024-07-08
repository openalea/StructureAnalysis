#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::DiscreteParametric const volatile * get_pointer<class ::stat_tool::DiscreteParametric const volatile >(class ::stat_tool::DiscreteParametric const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_82d14e3aa0215f91b29ab0b4c2b24290()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::DiscreteParametric::*method_pointer_abfba932b25d52f7bfd2f2ea5b43da88)(int , int , double , double ) = &::stat_tool::DiscreteParametric::init;
    void  (::stat_tool::DiscreteParametric::*method_pointer_d60dea7c925a54c0bd2a855bbe1a5f58)(enum ::stat_tool::discrete_parametric , int , int , double , double ) = &::stat_tool::DiscreteParametric::init;
    void  (::stat_tool::DiscreteParametric::*method_pointer_ae8dca5e086c5c93a61df11a32ad57c0)(class ::stat_tool::DiscreteParametric const &) = &::stat_tool::DiscreteParametric::copy;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::DiscreteParametric::*method_pointer_a7856a6f123a5467b6f843bedf375a0a)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::DiscreteParametric::ascii_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::DiscreteParametric::*method_pointer_eb1ae7cc3992584a8ca9ce58502b9b3e)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , bool ) const = &::stat_tool::DiscreteParametric::ascii_parametric_characteristic_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::DiscreteParametric::*method_pointer_f758963d54b55ec4aed4030060db4181)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::DiscreteParametric::spreadsheet_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::DiscreteParametric::*method_pointer_5f9ad6e464ef56df911b721da26566b8)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::DiscreteParametric::spreadsheet_parametric_characteristic_print;
    int  (*method_pointer_60c76b854c1d563496f805ccd43c71d0)(enum ::stat_tool::discrete_parametric , int , int , double , double , double ) = ::stat_tool::DiscreteParametric::nb_value_computation;
    int  (::stat_tool::DiscreteParametric::*method_pointer_a24cc9da2b445b8fa9fbd5e1fa58f5f1)() = &::stat_tool::DiscreteParametric::nb_parameter_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_4e8d96996526545bb344532cc0ebec2d)() = &::stat_tool::DiscreteParametric::nb_parameter_update;
    double  (::stat_tool::DiscreteParametric::*method_pointer_3ba4fb2776bf5bbda2feb6b1d71eb37e)() const = &::stat_tool::DiscreteParametric::parametric_mean_computation;
    double  (::stat_tool::DiscreteParametric::*method_pointer_82ad6ecc57b75a5382f3c3a1215c6a9f)() const = &::stat_tool::DiscreteParametric::parametric_variance_computation;
    double  (::stat_tool::DiscreteParametric::*method_pointer_2007520c52a35b8dae3f03347d06460b)() const = &::stat_tool::DiscreteParametric::parametric_skewness_computation;
    double  (::stat_tool::DiscreteParametric::*method_pointer_35685c0d655e595c80326e55d42db7aa)() const = &::stat_tool::DiscreteParametric::parametric_kurtosis_computation;
    double  (::stat_tool::DiscreteParametric::*method_pointer_e93c764042a8541886aef2c7b94ef8de)(class ::stat_tool::DiscreteParametric const &) const = &::stat_tool::DiscreteParametric::sup_norm_distance_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_f08fd4cf5057564bb24d4fb8eaae5a5e)(int , enum ::stat_tool::distribution_computation ) = &::stat_tool::DiscreteParametric::binomial_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_e2982f3ab18058e589f27cc1c9eb2e31)(int , double , enum ::stat_tool::distribution_computation ) = &::stat_tool::DiscreteParametric::poisson_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_88f3a00219075c28988cca75b0b86470)(int , double , enum ::stat_tool::distribution_computation ) = &::stat_tool::DiscreteParametric::negative_binomial_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_cb32ab393d84536f8de01edeb15c25f8)(int , double ) = &::stat_tool::DiscreteParametric::poisson_geometric_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_51178caa41f9529b9e99d27f1cc17698)() = &::stat_tool::DiscreteParametric::uniform_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_45e2e6f032ec5d659f20afe34bb3c959)() = &::stat_tool::DiscreteParametric::prior_segment_length_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_57ac64387476563898463f45cfd5f06d)(int , double ) = &::stat_tool::DiscreteParametric::computation;
    int  (::stat_tool::DiscreteParametric::*method_pointer_f92a8082d2cf58e8acbe4f9d734c54eb)() const = &::stat_tool::DiscreteParametric::simulation;
    double  (::stat_tool::DiscreteParametric::*method_pointer_f4e78197979c52a0a7542d0ff9dfaf99)(class ::stat_tool::Forward const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const *) const = &::stat_tool::DiscreteParametric::renewal_likelihood_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_fb81b674097f501c9e2fe702b1ce1772)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const *, class ::stat_tool::Reestimation< double > *, class ::stat_tool::Reestimation< double > *, int ) const = &::stat_tool::DiscreteParametric::expectation_step;
    double  (::stat_tool::DiscreteParametric::*method_pointer_ea4f3357a8be53688998cb7f27f3b865)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::DiscreteParametric::state_occupancy_likelihood_computation;
    double  (::stat_tool::DiscreteParametric::*method_pointer_15a05927895c516cad12850712c13c6b)(class ::stat_tool::Forward const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::DiscreteParametric::state_occupancy_likelihood_computation;
    void  (::stat_tool::DiscreteParametric::*method_pointer_2695736208d05343809912a569eace13)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Reestimation< double > *, int ) const = &::stat_tool::DiscreteParametric::expectation_step;
    void  (::stat_tool::DiscreteParametric::*method_pointer_def459f4cd2a5beab74b4f631c0bd319)(class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Reestimation< double > *, class ::stat_tool::Reestimation< double > *, int , bool , enum ::stat_tool::duration_distribution_mean_estimator ) const = &::stat_tool::DiscreteParametric::expectation_step;
    boost::python::class_< class ::stat_tool::DiscreteParametric, autowig::Held< class ::stat_tool::DiscreteParametric >::Type, boost::python::bases< class ::stat_tool::Distribution > > class_82d14e3aa0215f91b29ab0b4c2b24290("DiscreteParametric", "Discrete parametric distribution\n\n", boost::python::no_init);
    class_82d14e3aa0215f91b29ab0b4c2b24290.def(boost::python::init< int , enum ::stat_tool::discrete_parametric , int , int , double , double  >(""));
    class_82d14e3aa0215f91b29ab0b4c2b24290.def(boost::python::init< enum ::stat_tool::discrete_parametric , int , int , double , double , double  >(""));
    class_82d14e3aa0215f91b29ab0b4c2b24290.def(boost::python::init< class ::stat_tool::Distribution const &, int  >(""));
    class_82d14e3aa0215f91b29ab0b4c2b24290.def(boost::python::init< class ::stat_tool::Distribution const &, double  >(""));
    class_82d14e3aa0215f91b29ab0b4c2b24290.def(boost::python::init< class ::stat_tool::DiscreteParametric const &, double  >(""));
    class_82d14e3aa0215f91b29ab0b4c2b24290.def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >(""));
    class_82d14e3aa0215f91b29ab0b4c2b24290.def(boost::python::init< class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::distribution_transformation , int  >(""));
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("init", method_pointer_abfba932b25d52f7bfd2f2ea5b43da88, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("init", method_pointer_d60dea7c925a54c0bd2a855bbe1a5f58, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("copy", method_pointer_ae8dca5e086c5c93a61df11a32ad57c0, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("ascii_print", method_pointer_a7856a6f123a5467b6f843bedf375a0a, boost::python::return_internal_reference<>(), "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("ascii_parametric_characteristic_print", method_pointer_eb1ae7cc3992584a8ca9ce58502b9b3e, boost::python::return_internal_reference<>(), "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("spreadsheet_print", method_pointer_f758963d54b55ec4aed4030060db4181, boost::python::return_internal_reference<>(), "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("spreadsheet_parametric_characteristic_print", method_pointer_5f9ad6e464ef56df911b721da26566b8, boost::python::return_internal_reference<>(), "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("nb_value_computation", method_pointer_60c76b854c1d563496f805ccd43c71d0, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("nb_parameter_computation", method_pointer_a24cc9da2b445b8fa9fbd5e1fa58f5f1, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("nb_parameter_update", method_pointer_4e8d96996526545bb344532cc0ebec2d, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("parametric_mean_computation", method_pointer_3ba4fb2776bf5bbda2feb6b1d71eb37e, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("parametric_variance_computation", method_pointer_82ad6ecc57b75a5382f3c3a1215c6a9f, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("parametric_skewness_computation", method_pointer_2007520c52a35b8dae3f03347d06460b, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("parametric_kurtosis_computation", method_pointer_35685c0d655e595c80326e55d42db7aa, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("sup_norm_distance_computation", method_pointer_e93c764042a8541886aef2c7b94ef8de, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("binomial_computation", method_pointer_f08fd4cf5057564bb24d4fb8eaae5a5e, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("poisson_computation", method_pointer_e2982f3ab18058e589f27cc1c9eb2e31, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("negative_binomial_computation", method_pointer_88f3a00219075c28988cca75b0b86470, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("poisson_geometric_computation", method_pointer_cb32ab393d84536f8de01edeb15c25f8, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("uniform_computation", method_pointer_51178caa41f9529b9e99d27f1cc17698, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("prior_segment_length_computation", method_pointer_45e2e6f032ec5d659f20afe34bb3c959, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("computation", method_pointer_57ac64387476563898463f45cfd5f06d, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("simulation", method_pointer_f92a8082d2cf58e8acbe4f9d734c54eb, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("renewal_likelihood_computation", method_pointer_f4e78197979c52a0a7542d0ff9dfaf99, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("expectation_step", method_pointer_fb81b674097f501c9e2fe702b1ce1772, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("state_occupancy_likelihood_computation", method_pointer_ea4f3357a8be53688998cb7f27f3b865, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("state_occupancy_likelihood_computation", method_pointer_15a05927895c516cad12850712c13c6b, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("expectation_step", method_pointer_2695736208d05343809912a569eace13, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def("expectation_step", method_pointer_def459f4cd2a5beab74b4f631c0bd319, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.staticmethod("nb_value_computation");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def_readwrite("ident", &::stat_tool::DiscreteParametric::ident, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def_readwrite("parameter", &::stat_tool::DiscreteParametric::parameter, "");
    class_82d14e3aa0215f91b29ab0b4c2b24290.def_readwrite("probability", &::stat_tool::DiscreteParametric::probability, "");

    if(autowig::Held< class ::stat_tool::DiscreteParametric >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteParametric >::Type, autowig::Held< class ::stat_tool::Distribution >::Type >();
    }

}