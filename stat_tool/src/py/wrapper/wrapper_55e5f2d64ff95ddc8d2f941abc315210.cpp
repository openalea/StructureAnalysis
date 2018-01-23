#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::initializer_list< class ::stat_tool::FrequencyDistribution > const volatile * get_pointer<class ::std::initializer_list< class ::stat_tool::FrequencyDistribution > const volatile >(class ::std::initializer_list< class ::stat_tool::FrequencyDistribution > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_55e5f2d64ff95ddc8d2f941abc315210()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    ::std::initializer_list< class ::stat_tool::FrequencyDistribution >::size_type  (::std::initializer_list< ::stat_tool::FrequencyDistribution >::*method_pointer_ea77dfc0152e58c0aae730de5648c6df)() const = &::std::initializer_list< class ::stat_tool::FrequencyDistribution >::size;
    ::std::initializer_list< class ::stat_tool::FrequencyDistribution >::const_iterator  (::std::initializer_list< ::stat_tool::FrequencyDistribution >::*method_pointer_c5a96a56b49952ef8f2b8a6084121a2f)() const = &::std::initializer_list< class ::stat_tool::FrequencyDistribution >::begin;
    ::std::initializer_list< class ::stat_tool::FrequencyDistribution >::const_iterator  (::std::initializer_list< ::stat_tool::FrequencyDistribution >::*method_pointer_a2a432d41290575dad0409071747f892)() const = &::std::initializer_list< class ::stat_tool::FrequencyDistribution >::end;
    boost::python::class_< class ::std::initializer_list< class ::stat_tool::FrequencyDistribution >, autowig::Held< class ::std::initializer_list< class ::stat_tool::FrequencyDistribution > >::Type > class_55e5f2d64ff95ddc8d2f941abc315210("_InitializerList_55e5f2d64ff95ddc8d2f941abc315210", "", boost::python::no_init);
    class_55e5f2d64ff95ddc8d2f941abc315210.def("__len__", method_pointer_ea77dfc0152e58c0aae730de5648c6df, "");
    class_55e5f2d64ff95ddc8d2f941abc315210.def("begin", method_pointer_c5a96a56b49952ef8f2b8a6084121a2f, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_55e5f2d64ff95ddc8d2f941abc315210.def("end", method_pointer_a2a432d41290575dad0409071747f892, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");

}