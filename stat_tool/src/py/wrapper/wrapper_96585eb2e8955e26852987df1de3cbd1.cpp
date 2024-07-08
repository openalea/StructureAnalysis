#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::DiscreteParametricModel const volatile * get_pointer<class ::stat_tool::DiscreteParametricModel const volatile >(class ::stat_tool::DiscreteParametricModel const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_96585eb2e8955e26852987df1de3cbd1()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::DiscreteParametricModel::*method_pointer_32615e4dbd38500b914c96f21f414665)(class ::stat_tool::StatError &) const = &::stat_tool::DiscreteParametricModel::extract_data;
    class ::stat_tool::DiscreteParametricModel * (*method_pointer_6b4a5891c54d5b0cb3a860a4b25d4a1a)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, double ) = ::stat_tool::DiscreteParametricModel::ascii_read;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::DiscreteParametricModel::*method_pointer_10571acabe925b6aaff7cc61716a1496)(class ::stat_tool::StatError &, int ) const = &::stat_tool::DiscreteParametricModel::simulation;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::DiscreteParametricModel::*method_pointer_4c807bb6b90358228851305f7e6bc996)() const = &::stat_tool::DiscreteParametricModel::get_frequency_distribution;
    boost::python::class_< class ::stat_tool::DiscreteParametricModel, autowig::Held< class ::stat_tool::DiscreteParametricModel >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::DiscreteParametric > > class_96585eb2e8955e26852987df1de3cbd1("DiscreteParametricModel", "Discrete parametric distribution\n\n", boost::python::no_init);
    class_96585eb2e8955e26852987df1de3cbd1.def(boost::python::init< int , enum ::stat_tool::discrete_parametric , int , int , double , double  >(""));
    class_96585eb2e8955e26852987df1de3cbd1.def(boost::python::init< enum ::stat_tool::discrete_parametric , int , int , double , double , double  >(""));
    class_96585eb2e8955e26852987df1de3cbd1.def(boost::python::init< class ::stat_tool::FrequencyDistribution const & >(""));
    class_96585eb2e8955e26852987df1de3cbd1.def(boost::python::init< class ::stat_tool::Distribution const & >(""));
    class_96585eb2e8955e26852987df1de3cbd1.def(boost::python::init< class ::stat_tool::DiscreteParametric const & >(""));
    class_96585eb2e8955e26852987df1de3cbd1.def(boost::python::init< class ::stat_tool::Distribution const &, class ::stat_tool::FrequencyDistribution const * >(""));
    class_96585eb2e8955e26852987df1de3cbd1.def(boost::python::init< class ::stat_tool::DiscreteParametric const &, class ::stat_tool::FrequencyDistribution const * >(""));
    class_96585eb2e8955e26852987df1de3cbd1.def(boost::python::init< class ::stat_tool::DiscreteParametricModel const &, bool  >(""));
    class_96585eb2e8955e26852987df1de3cbd1.def("extract_data", method_pointer_32615e4dbd38500b914c96f21f414665, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_96585eb2e8955e26852987df1de3cbd1.def("ascii_read", method_pointer_6b4a5891c54d5b0cb3a860a4b25d4a1a, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_96585eb2e8955e26852987df1de3cbd1.def("simulation", method_pointer_10571acabe925b6aaff7cc61716a1496, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_96585eb2e8955e26852987df1de3cbd1.def("get_frequency_distribution", method_pointer_4c807bb6b90358228851305f7e6bc996, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_96585eb2e8955e26852987df1de3cbd1.staticmethod("ascii_read");

    if(autowig::Held< class ::stat_tool::DiscreteParametricModel >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteParametricModel >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::DiscreteParametricModel >::Type, autowig::Held< class ::stat_tool::DiscreteParametric >::Type >();
    }

}