#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::MixtureData const volatile * get_pointer<class ::stat_tool::MixtureData const volatile >(class ::stat_tool::MixtureData const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_0bd75707ac98530d926ac17eb4fe47d3()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::MixtureData::*method_pointer_22b0a6e9b8985815a70d8c14b0968354)(class ::stat_tool::StatError &, int , int ) const = &::stat_tool::MixtureData::extract;
    class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > >  (::stat_tool::MixtureData::*method_pointer_9b1cb4d742805bfcb50ab9953cec9468)(bool ) const = &::stat_tool::MixtureData::ascii_data_write;
    void  (::stat_tool::MixtureData::*method_pointer_65b10f92a726577cbf6564411c83c993)(enum ::stat_tool::variable_nature ) = &::stat_tool::MixtureData::state_variable_init;
    double  (::stat_tool::MixtureData::*method_pointer_33ca36d613b75e3ebc8bbd312abb4ab5)() const = &::stat_tool::MixtureData::classification_information_computation;
    double  (::stat_tool::MixtureData::*method_pointer_500bb41fb5f55d4e8a20e25372760385)() const = &::stat_tool::MixtureData::information_computation;
    void  (::stat_tool::MixtureData::*method_pointer_e3d6ff48072e511caff022cc29ffe844)(int ) = &::stat_tool::MixtureData::build_observation_frequency_distribution;
    void  (::stat_tool::MixtureData::*method_pointer_acc63112985659a790213c952e1d2bdf)(int , int , double ) = &::stat_tool::MixtureData::build_observation_histogram;
    void  (::stat_tool::MixtureData::*method_pointer_d2c29ecee265561ca35772e3fb810044)(int ) = &::stat_tool::MixtureData::build_observation_histogram;
    bool  (::stat_tool::MixtureData::*method_pointer_cefb5df5edf75356bb3720d082eabc58)(class ::stat_tool::StatError &, int , double , double ) = &::stat_tool::MixtureData::select_bin_width;
    class ::stat_tool::Mixture * (::stat_tool::MixtureData::*method_pointer_7f8d731d35af5419bd666194a383fb26)() const = &::stat_tool::MixtureData::get_mixture;
    class ::stat_tool::FrequencyDistribution * (::stat_tool::MixtureData::*method_pointer_a03eb68b67645232b08139d24f6872bc)(int , int ) const = &::stat_tool::MixtureData::get_observation_distribution;
    class ::stat_tool::Histogram * (::stat_tool::MixtureData::*method_pointer_7a2d10c5ce5952b68902366170c32ad0)(int , int ) const = &::stat_tool::MixtureData::get_observation_histogram;
    double  (::stat_tool::MixtureData::*method_pointer_fb0825b198655748bbb63293623f2952)() const = &::stat_tool::MixtureData::get_likelihood;
    double  (::stat_tool::MixtureData::*method_pointer_297cd32e02845d91946c609b7dabe7a0)() const = &::stat_tool::MixtureData::get_restoration_likelihood;
    double  (::stat_tool::MixtureData::*method_pointer_d994a13faeee5eecbd9fd6c88cba2457)() const = &::stat_tool::MixtureData::get_sample_entropy;
    double  (::stat_tool::MixtureData::*method_pointer_ba471cb1da16567d8a412133323a293d)(int ) const = &::stat_tool::MixtureData::get_posterior_probability;
    double  (::stat_tool::MixtureData::*method_pointer_c2cf941a4bdb513dac8d83e9192ce107)(int ) const = &::stat_tool::MixtureData::get_entropy;
    boost::python::class_< class ::stat_tool::MixtureData, autowig::Held< class ::stat_tool::MixtureData >::Type, boost::python::bases< class ::stat_tool::Vectors > > class_0bd75707ac98530d926ac17eb4fe47d3("MixtureData", "Data structure corresponding to a multivariate mixture of distributions\n\n", boost::python::no_init);
    class_0bd75707ac98530d926ac17eb4fe47d3.def(boost::python::init<  >(""));
    class_0bd75707ac98530d926ac17eb4fe47d3.def(boost::python::init< int , int , enum ::stat_tool::variable_nature *, bool  >(""));
    class_0bd75707ac98530d926ac17eb4fe47d3.def(boost::python::init< class ::stat_tool::Vectors const &, enum ::stat_tool::vector_transformation  >(""));
    class_0bd75707ac98530d926ac17eb4fe47d3.def(boost::python::init< class ::stat_tool::MixtureData const &, bool , enum ::stat_tool::vector_transformation  >(""));
    class_0bd75707ac98530d926ac17eb4fe47d3.def("extract", method_pointer_22b0a6e9b8985815a70d8c14b0968354, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("ascii_data_write", method_pointer_9b1cb4d742805bfcb50ab9953cec9468, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("state_variable_init", method_pointer_65b10f92a726577cbf6564411c83c993, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("classification_information_computation", method_pointer_33ca36d613b75e3ebc8bbd312abb4ab5, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("information_computation", method_pointer_500bb41fb5f55d4e8a20e25372760385, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("build_observation_frequency_distribution", method_pointer_e3d6ff48072e511caff022cc29ffe844, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("build_observation_histogram", method_pointer_acc63112985659a790213c952e1d2bdf, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("build_observation_histogram", method_pointer_d2c29ecee265561ca35772e3fb810044, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("select_bin_width", method_pointer_cefb5df5edf75356bb3720d082eabc58, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("get_mixture", method_pointer_7f8d731d35af5419bd666194a383fb26, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("get_observation_distribution", method_pointer_a03eb68b67645232b08139d24f6872bc, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("get_observation_histogram", method_pointer_7a2d10c5ce5952b68902366170c32ad0, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("get_likelihood", method_pointer_fb0825b198655748bbb63293623f2952, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("get_restoration_likelihood", method_pointer_297cd32e02845d91946c609b7dabe7a0, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("get_sample_entropy", method_pointer_d994a13faeee5eecbd9fd6c88cba2457, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("get_posterior_probability", method_pointer_ba471cb1da16567d8a412133323a293d, "");
    class_0bd75707ac98530d926ac17eb4fe47d3.def("get_entropy", method_pointer_c2cf941a4bdb513dac8d83e9192ce107, "");

    if(autowig::Held< class ::stat_tool::MixtureData >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::MixtureData >::Type, autowig::Held< class ::stat_tool::Vectors >::Type >();
    }

}