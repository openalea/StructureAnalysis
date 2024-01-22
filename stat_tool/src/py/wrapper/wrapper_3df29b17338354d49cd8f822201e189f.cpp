#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::CompoundData const volatile * get_pointer<class ::stat_tool::CompoundData const volatile >(class ::stat_tool::CompoundData const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_3df29b17338354d49cd8f822201e189f()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::stat_tool::DiscreteDistributionData * (::stat_tool::CompoundData::*method_pointer_f030d5b9a4fd569c81737fcc1a40e44d)(class ::stat_tool::StatError &, enum ::stat_tool::compound_distribution ) const = &::stat_tool::CompoundData::extract;
    class ::stat_tool::Compound * (::stat_tool::CompoundData::*method_pointer_b4247691758657039c75f6e8a5005fc5)() const = &::stat_tool::CompoundData::get_compound;
    class ::stat_tool::FrequencyDistribution * (::stat_tool::CompoundData::*method_pointer_71f25a4d9a3d5dd48e19ef196515a723)() const = &::stat_tool::CompoundData::get_sum_frequency_distribution;
    class ::stat_tool::FrequencyDistribution * (::stat_tool::CompoundData::*method_pointer_7eb4fdc1dfe05a79952798f599446ac0)() const = &::stat_tool::CompoundData::get_frequency_distribution;
    boost::python::class_< class ::stat_tool::CompoundData, autowig::Held< class ::stat_tool::CompoundData >::Type, boost::python::bases< class ::stat_tool::StatInterface, class ::stat_tool::FrequencyDistribution > > class_3df29b17338354d49cd8f822201e189f("CompoundData", "Data structure corresponding to a compound distribution\n\n", boost::python::no_init);
    class_3df29b17338354d49cd8f822201e189f.def(boost::python::init<  >(""));
    class_3df29b17338354d49cd8f822201e189f.def(boost::python::init< class ::stat_tool::FrequencyDistribution const &, class ::stat_tool::Compound const & >(""));
    class_3df29b17338354d49cd8f822201e189f.def(boost::python::init< class ::stat_tool::Compound const & >(""));
    class_3df29b17338354d49cd8f822201e189f.def(boost::python::init< class ::stat_tool::CompoundData const &, bool  >(""));
    class_3df29b17338354d49cd8f822201e189f.def("extract", method_pointer_f030d5b9a4fd569c81737fcc1a40e44d, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_3df29b17338354d49cd8f822201e189f.def("get_compound", method_pointer_b4247691758657039c75f6e8a5005fc5, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_3df29b17338354d49cd8f822201e189f.def("get_sum_frequency_distribution", method_pointer_71f25a4d9a3d5dd48e19ef196515a723, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_3df29b17338354d49cd8f822201e189f.def("get_frequency_distribution", method_pointer_7eb4fdc1dfe05a79952798f599446ac0, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");

    if(autowig::Held< class ::stat_tool::CompoundData >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::CompoundData >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::implicitly_convertible< autowig::Held< class ::stat_tool::CompoundData >::Type, autowig::Held< class ::stat_tool::FrequencyDistribution >::Type >();
    }

}