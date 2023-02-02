#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::FrequencyDistribution const volatile * get_pointer<class ::stat_tool::FrequencyDistribution const volatile >(class ::stat_tool::FrequencyDistribution const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_f7aa7e5718795bbe854edb6369819dd9()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    bool  (::stat_tool::FrequencyDistribution::*method_pointer_1e33eb1d4cac5053a5297cdd62f11f86)(class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::operator==;
    bool  (::stat_tool::FrequencyDistribution::*method_pointer_27baba88c70154d1a16dd9895e36ee28)(class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::operator!=;
    void  (::stat_tool::FrequencyDistribution::*method_pointer_9c4ca4152b24578185b165f5eada378f)(class ::stat_tool::SinglePlot &) const = &::stat_tool::FrequencyDistribution::plotable_frequency_write;
    void  (::stat_tool::FrequencyDistribution::*method_pointer_cdf19052baa0516db9e6c59e9ad197cd)(class ::stat_tool::SinglePlot &) const = &::stat_tool::FrequencyDistribution::plotable_mass_write;
    void  (::stat_tool::FrequencyDistribution::*method_pointer_8e060e16db8b56009d6d4f42547088bd)(class ::stat_tool::SinglePlot &) const = &::stat_tool::FrequencyDistribution::plotable_survivor_write;
    double  (::stat_tool::FrequencyDistribution::*method_pointer_0b8608a35c3550a99ac8a0958d802fd4)() const = &::stat_tool::FrequencyDistribution::concentration_computation;
    void  (::stat_tool::FrequencyDistribution::*method_pointer_9fe9a16abbbb56babbba6a548924f162)(class ::stat_tool::Reestimation< double > const *, int ) = &::stat_tool::FrequencyDistribution::update;
    class ::stat_tool::FrequencyDistribution * (::stat_tool::FrequencyDistribution::*method_pointer_5882a600801d5cc5ac78f111a990f6b4)(int ) const = &::stat_tool::FrequencyDistribution::frequency_scale;
    int  (::stat_tool::FrequencyDistribution::*method_pointer_2b8a299165215e44ad3ec5cb8e6c85b0)() const = &::stat_tool::FrequencyDistribution::min_interval_computation;
    class ::stat_tool::DiscreteParametric * (::stat_tool::FrequencyDistribution::*method_pointer_5378bc5deac05aedb57ce3e0fba4344b)(enum ::stat_tool::discrete_parametric , int , bool , double ) const = &::stat_tool::FrequencyDistribution::parametric_estimation;
    double  (::stat_tool::FrequencyDistribution::*method_pointer_b0a0fd68ce7952998ac6c6e0c971980f)(class ::stat_tool::ContinuousParametric const &, int ) const = &::stat_tool::FrequencyDistribution::likelihood_computation;
    void  (::stat_tool::FrequencyDistribution::*method_pointer_01abb08d00865de783d5f214e2fb750c)(class ::stat_tool::FrequencyDistribution const &, int ) = &::stat_tool::FrequencyDistribution::shift;
    void  (::stat_tool::FrequencyDistribution::*method_pointer_2010daaf8405536eb6244ab6f6ae7e30)(class ::stat_tool::FrequencyDistribution const &, int , enum ::stat_tool::rounding ) = &::stat_tool::FrequencyDistribution::cluster;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_258945c32f1c5ddfafe44231fd2ce510)(class ::stat_tool::StatError &, int ) const = &::stat_tool::FrequencyDistribution::shift;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_8d4e878e460e5fe7ad1180daa4d9591e)(class ::stat_tool::StatError &, int , enum ::stat_tool::rounding ) const = &::stat_tool::FrequencyDistribution::cluster;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_0306613ca8c659f78e816cd4421b3eef)(class ::stat_tool::StatError &, double , bool ) const = &::stat_tool::FrequencyDistribution::cluster;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_f435e82c2ba25da89607e27fb6ff69ca)(class ::stat_tool::StatError &, int , class ::std::vector< int, class ::std::allocator< int > > ) const = &::stat_tool::FrequencyDistribution::cluster;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_edef132f745e5822956701457993c716)(class ::stat_tool::StatError &, class ::std::vector< int, class ::std::allocator< int > > ) const = &::stat_tool::FrequencyDistribution::transcode;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::FrequencyDistribution::*method_pointer_e924ee1e64a45a2ebf02d8b85778ec21)(class ::stat_tool::StatError &, int , int , bool ) const = &::stat_tool::FrequencyDistribution::value_select;
    bool  (::stat_tool::FrequencyDistribution::*method_pointer_1cd26df01a8e55b2bdd42cfcbf3172c4)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) const = &::stat_tool::FrequencyDistribution::ascii_write;
    ::stat_tool::MultiPlotSet * (::stat_tool::FrequencyDistribution::*method_pointer_22f6f58186ea54b6adfad7525b30c70c)() const = &::stat_tool::FrequencyDistribution::get_plotable;
    bool  (::stat_tool::FrequencyDistribution::*method_pointer_a08d3756f3945bef95dbe223fd480b5f)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) const = &::stat_tool::FrequencyDistribution::survival_ascii_write;
    bool  (::stat_tool::FrequencyDistribution::*method_pointer_b00a5fd2f81b59f6a9a98778e5dd1ab5)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) const = &::stat_tool::FrequencyDistribution::survival_spreadsheet_write;
    ::stat_tool::MultiPlotSet * (::stat_tool::FrequencyDistribution::*method_pointer_bed5046815d5536db150188d88618460)(class ::stat_tool::StatError &) const = &::stat_tool::FrequencyDistribution::survival_get_plotable;
    void  (::stat_tool::FrequencyDistribution::*method_pointer_2c92ad0bc8e3590bad08a18cc201e9f1)(bool , class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::F_comparison;
    void  (::stat_tool::FrequencyDistribution::*method_pointer_bf6cbe53aef2532db17fd29f79846fc5)(bool , class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::t_comparison;
    bool  (::stat_tool::FrequencyDistribution::*method_pointer_8feb4c63d9275090a65f51f743982f4c)(class ::stat_tool::StatError &, bool , class ::stat_tool::FrequencyDistribution const &) const = &::stat_tool::FrequencyDistribution::wilcoxon_mann_whitney_comparison;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_0e6e66d328825e58ac42478d950343a9)(class ::stat_tool::StatError &, class ::stat_tool::DiscreteParametric const &) const = &::stat_tool::FrequencyDistribution::fit;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_4443f5df74fd5265b08286162c0e8df5)(class ::stat_tool::StatError &, enum ::stat_tool::discrete_parametric , int , bool , double ) const = &::stat_tool::FrequencyDistribution::parametric_estimation;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_93be0d5f583551e48f38113575b8cdcb)(class ::stat_tool::StatError &, int , bool , double ) const = &::stat_tool::FrequencyDistribution::type_parametric_estimation;
    class ::stat_tool::DiscreteMixture * (::stat_tool::FrequencyDistribution::*method_pointer_3c7877e75df35d599ee7ea6ab878440e)(class ::stat_tool::StatError &, class ::stat_tool::DiscreteMixture const &, int , bool , bool , double ) const = &::stat_tool::FrequencyDistribution::discrete_mixture_estimation;
    class ::stat_tool::DiscreteMixture * (::stat_tool::FrequencyDistribution::*method_pointer_e9ca12dfcfc85ef889bdd92cd00ab762)(class ::stat_tool::StatError &, int , enum ::stat_tool::discrete_parametric *, int , bool , bool , double ) const = &::stat_tool::FrequencyDistribution::discrete_mixture_estimation;
    class ::stat_tool::DiscreteMixture * (::stat_tool::FrequencyDistribution::*method_pointer_7db0ac831cc25d709beee6e0ce93cbac)(class ::stat_tool::StatError &, bool , int , int , enum ::stat_tool::discrete_parametric *, int , bool , bool , enum ::stat_tool::model_selection_criterion , double ) const = &::stat_tool::FrequencyDistribution::discrete_mixture_estimation;
    class ::stat_tool::Convolution * (::stat_tool::FrequencyDistribution::*method_pointer_68d32044c568566696531ad7fadf075f)(class ::stat_tool::StatError &, bool , class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::estimation_criterion , int , double , enum ::stat_tool::penalty_type , enum ::stat_tool::side_effect ) const = &::stat_tool::FrequencyDistribution::convolution_estimation;
    class ::stat_tool::Convolution * (::stat_tool::FrequencyDistribution::*method_pointer_3afd01d5515b5255aefd56e874a626d3)(class ::stat_tool::StatError &, bool , class ::stat_tool::DiscreteParametric const &, int , enum ::stat_tool::estimation_criterion , int , double , enum ::stat_tool::penalty_type , enum ::stat_tool::side_effect ) const = &::stat_tool::FrequencyDistribution::convolution_estimation;
    class ::stat_tool::Compound * (::stat_tool::FrequencyDistribution::*method_pointer_8102f5a2d92d5f7eb734002ea9df6099)(class ::stat_tool::StatError &, bool , class ::stat_tool::DiscreteParametric const &, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::compound_distribution , enum ::stat_tool::estimation_criterion , int , double , enum ::stat_tool::penalty_type , enum ::stat_tool::side_effect ) const = &::stat_tool::FrequencyDistribution::compound_estimation;
    class ::stat_tool::Compound * (::stat_tool::FrequencyDistribution::*method_pointer_be92445f015a554aaa2e1109dc236eda)(class ::stat_tool::StatError &, bool , class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::compound_distribution , int , enum ::stat_tool::estimation_criterion , int , double , enum ::stat_tool::penalty_type , enum ::stat_tool::side_effect ) const = &::stat_tool::FrequencyDistribution::compound_estimation;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_212bcd08409c5e84b1e7fa08776369c4)(class ::stat_tool::StatError &, bool , class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const *, class ::stat_tool::DiscreteParametric const &, enum ::stat_tool::estimation_criterion , int , enum ::stat_tool::duration_distribution_mean_estimator , double , enum ::stat_tool::penalty_type , enum ::stat_tool::side_effect , double ) const = &::stat_tool::FrequencyDistribution::estimation;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::FrequencyDistribution::*method_pointer_5af91847cb1c58d39e20f54bcd09e8d5)(class ::stat_tool::StatError &, bool , class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::FrequencyDistribution const *, enum ::stat_tool::estimation_criterion , int , enum ::stat_tool::duration_distribution_mean_estimator , double , enum ::stat_tool::penalty_type , enum ::stat_tool::side_effect ) const = &::stat_tool::FrequencyDistribution::estimation;
    boost::python::class_< class ::stat_tool::FrequencyDistribution, autowig::Held< class ::stat_tool::FrequencyDistribution >::Type, boost::python::bases< class ::stat_tool::Reestimation< int > > > class_f7aa7e5718795bbe854edb6369819dd9("FrequencyDistribution", "Frequency distribution\n\n", boost::python::no_init);
    class_f7aa7e5718795bbe854edb6369819dd9.def(boost::python::init< int  >(""));
    class_f7aa7e5718795bbe854edb6369819dd9.def(boost::python::init< class ::stat_tool::Distribution const & >(""));
    class_f7aa7e5718795bbe854edb6369819dd9.def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, enum ::stat_tool::frequency_distribution_transformation , int , enum ::stat_tool::rounding  >(""));
    class_f7aa7e5718795bbe854edb6369819dd9.def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >(""));
    class_f7aa7e5718795bbe854edb6369819dd9.def("__eq__", method_pointer_1e33eb1d4cac5053a5297cdd62f11f86, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("__neq__", method_pointer_27baba88c70154d1a16dd9895e36ee28, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("plotable_frequency_write", method_pointer_9c4ca4152b24578185b165f5eada378f, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("plotable_mass_write", method_pointer_cdf19052baa0516db9e6c59e9ad197cd, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("plotable_survivor_write", method_pointer_8e060e16db8b56009d6d4f42547088bd, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("concentration_computation", method_pointer_0b8608a35c3550a99ac8a0958d802fd4, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("update", method_pointer_9fe9a16abbbb56babbba6a548924f162, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("frequency_scale", method_pointer_5882a600801d5cc5ac78f111a990f6b4, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("min_interval_computation", method_pointer_2b8a299165215e44ad3ec5cb8e6c85b0, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("parametric_estimation", method_pointer_5378bc5deac05aedb57ce3e0fba4344b, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("likelihood_computation", method_pointer_b0a0fd68ce7952998ac6c6e0c971980f, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("shift", method_pointer_01abb08d00865de783d5f214e2fb750c, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("cluster", method_pointer_2010daaf8405536eb6244ab6f6ae7e30, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("shift", method_pointer_258945c32f1c5ddfafe44231fd2ce510, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("cluster", method_pointer_8d4e878e460e5fe7ad1180daa4d9591e, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("cluster", method_pointer_0306613ca8c659f78e816cd4421b3eef, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("cluster", method_pointer_f435e82c2ba25da89607e27fb6ff69ca, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("transcode", method_pointer_edef132f745e5822956701457993c716, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("value_select", method_pointer_e924ee1e64a45a2ebf02d8b85778ec21, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("ascii_write", method_pointer_1cd26df01a8e55b2bdd42cfcbf3172c4, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("get_plotable", method_pointer_22f6f58186ea54b6adfad7525b30c70c, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("survival_ascii_write", method_pointer_a08d3756f3945bef95dbe223fd480b5f, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("survival_spreadsheet_write", method_pointer_b00a5fd2f81b59f6a9a98778e5dd1ab5, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("survival_get_plotable", method_pointer_bed5046815d5536db150188d88618460, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("f__comparison", method_pointer_2c92ad0bc8e3590bad08a18cc201e9f1, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("t_comparison", method_pointer_bf6cbe53aef2532db17fd29f79846fc5, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("wilcoxon_mann_whitney_comparison", method_pointer_8feb4c63d9275090a65f51f743982f4c, "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("fit", method_pointer_0e6e66d328825e58ac42478d950343a9, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("parametric_estimation", method_pointer_4443f5df74fd5265b08286162c0e8df5, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("type_parametric_estimation", method_pointer_93be0d5f583551e48f38113575b8cdcb, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("discrete_mixture_estimation", method_pointer_3c7877e75df35d599ee7ea6ab878440e, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("discrete_mixture_estimation", method_pointer_e9ca12dfcfc85ef889bdd92cd00ab762, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("discrete_mixture_estimation", method_pointer_7db0ac831cc25d709beee6e0ce93cbac, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("convolution_estimation", method_pointer_68d32044c568566696531ad7fadf075f, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("convolution_estimation", method_pointer_3afd01d5515b5255aefd56e874a626d3, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("compound_estimation", method_pointer_8102f5a2d92d5f7eb734002ea9df6099, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("compound_estimation", method_pointer_be92445f015a554aaa2e1109dc236eda, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("estimation", method_pointer_212bcd08409c5e84b1e7fa08776369c4, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f7aa7e5718795bbe854edb6369819dd9.def("estimation", method_pointer_5af91847cb1c58d39e20f54bcd09e8d5, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");

    if(autowig::Held< class ::stat_tool::FrequencyDistribution >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::FrequencyDistribution >::Type, autowig::Held< class ::stat_tool::Reestimation< int > >::Type >();
        //boost::python::objects::class_value_wrapper< autowig::Held< class ::stat_tool::FrequencyDistribution >::Type, boost::python::objects::make_ptr_instance< class ::stat_tool::FrequencyDistribution, boost::python::objects::pointer_holder< autowig::Held< class ::stat_tool::FrequencyDistribution >::Type, class ::stat_tool::FrequencyDistribution > > >();
    }

}