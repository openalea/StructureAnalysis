#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::DiscreteDistributionData const volatile * get_pointer<class ::stat_tool::DiscreteDistributionData const volatile >(class ::stat_tool::DiscreteDistributionData const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_699cdb7e19d058c5bb8d30066fdc55ac()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteParametricModel * (::stat_tool::DiscreteDistributionData::*method_pointer_09f4b4034c225520bfc427bc9432f05a)(class ::stat_tool::StatError &) const = &::stat_tool::DiscreteDistributionData::extract_model;
    class ::stat_tool::DiscreteDistributionData * (*method_pointer_59d66a506b4f5e198e7598f67a1365fd)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) = ::stat_tool::DiscreteDistributionData::ascii_read;
    class ::stat_tool::DiscreteParametric * (::stat_tool::DiscreteDistributionData::*method_pointer_3dc769f178ff5cacaf81003605ec49dc)() const = &::stat_tool::DiscreteDistributionData::get_distribution;
    boost::python::class_< class ::stat_tool::DiscreteDistributionData, autowig::Held< class ::stat_tool::DiscreteDistributionData >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::FrequencyDistribution > > class_699cdb7e19d058c5bb8d30066fdc55ac("DiscreteDistributionData", "", boost::python::no_init);
    class_699cdb7e19d058c5bb8d30066fdc55ac.def(boost::python::init< int  >(""));
    class_699cdb7e19d058c5bb8d30066fdc55ac.def(boost::python::init< class ::stat_tool::Distribution const & >(""));
    class_699cdb7e19d058c5bb8d30066fdc55ac.def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >(""));
    class_699cdb7e19d058c5bb8d30066fdc55ac.def(boost::python::init< class ::std::vector< int, class ::std::allocator< int > > const >(""));
    class_699cdb7e19d058c5bb8d30066fdc55ac.def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, enum ::stat_tool::frequency_distribution_transformation , int , enum ::stat_tool::rounding  >(""));
    class_699cdb7e19d058c5bb8d30066fdc55ac.def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Distribution const * >(""));
    class_699cdb7e19d058c5bb8d30066fdc55ac.def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::DiscreteParametric const * >(""));
    class_699cdb7e19d058c5bb8d30066fdc55ac.def(boost::python::init< class ::stat_tool::DiscreteDistributionData const &, bool  >(""));
    class_699cdb7e19d058c5bb8d30066fdc55ac.def("extract_model", method_pointer_09f4b4034c225520bfc427bc9432f05a, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_699cdb7e19d058c5bb8d30066fdc55ac.def("ascii_read", method_pointer_59d66a506b4f5e198e7598f67a1365fd, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_699cdb7e19d058c5bb8d30066fdc55ac.def("get_distribution", method_pointer_3dc769f178ff5cacaf81003605ec49dc, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_699cdb7e19d058c5bb8d30066fdc55ac.staticmethod("ascii_read");

    if(autowig::Held< class ::stat_tool::DiscreteDistributionData >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteDistributionData >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteDistributionData >::Type, autowig::Held< class ::stat_tool::FrequencyDistribution >::Type >();
    }

}