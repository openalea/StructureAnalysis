#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::ConvolutionData const volatile * get_pointer<class ::stat_tool::ConvolutionData const volatile >(class ::stat_tool::ConvolutionData const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_acd6178a21305dc0a048b18ff459a4b3()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::ConvolutionData::*method_pointer_9095edf2ab535d02b6e5dedd623f48d2)(class ::stat_tool::StatError &, int ) const = &::stat_tool::ConvolutionData::extract;
    class ::stat_tool::Convolution * (::stat_tool::ConvolutionData::*method_pointer_3d736d0b87aa50ea9a7b6f951a72aee0)() const = &::stat_tool::ConvolutionData::get_convolution;
    int  (::stat_tool::ConvolutionData::*method_pointer_88f77df9758959e5af42130b5e24340a)() const = &::stat_tool::ConvolutionData::get_nb_distribution;
    class ::stat_tool::FrequencyDistribution * (::stat_tool::ConvolutionData::*method_pointer_407c3928cc115eb7b55d1ae35f3c8cc9)(int ) const = &::stat_tool::ConvolutionData::get_frequency_distribution;
    boost::python::class_< class ::stat_tool::ConvolutionData, autowig::Held< class ::stat_tool::ConvolutionData >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::FrequencyDistribution > > class_acd6178a21305dc0a048b18ff459a4b3("ConvolutionData", "Data structure corresponding to a convolution of discrete distributions\n\n", boost::python::no_init);
    class_acd6178a21305dc0a048b18ff459a4b3.def(boost::python::init<  >(""));
    class_acd6178a21305dc0a048b18ff459a4b3.def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, int  >(""));
    class_acd6178a21305dc0a048b18ff459a4b3.def(boost::python::init< class ::stat_tool::Convolution const & >(""));
    class_acd6178a21305dc0a048b18ff459a4b3.def(boost::python::init< class ::stat_tool::ConvolutionData const &, bool  >(""));
    class_acd6178a21305dc0a048b18ff459a4b3.def("extract", method_pointer_9095edf2ab535d02b6e5dedd623f48d2, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_acd6178a21305dc0a048b18ff459a4b3.def("get_convolution", method_pointer_3d736d0b87aa50ea9a7b6f951a72aee0, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_acd6178a21305dc0a048b18ff459a4b3.def("get_nb_distribution", method_pointer_88f77df9758959e5af42130b5e24340a, "");
    class_acd6178a21305dc0a048b18ff459a4b3.def("get_frequency_distribution", method_pointer_407c3928cc115eb7b55d1ae35f3c8cc9, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");

    if(autowig::Held< class ::stat_tool::ConvolutionData >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::ConvolutionData >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::ConvolutionData >::Type, autowig::Held< class ::stat_tool::FrequencyDistribution >::Type >();
    }

}