#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::allocator< class ::stat_tool::FrequencyDistribution > const volatile * get_pointer<class ::std::allocator< class ::stat_tool::FrequencyDistribution > const volatile >(class ::std::allocator< class ::stat_tool::FrequencyDistribution > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_f2215828ca3a504b8e3d080988e9947f()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< class ::std::allocator< class ::stat_tool::FrequencyDistribution >, autowig::Held< class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::Type > class_f2215828ca3a504b8e3d080988e9947f("_Allocator_f2215828ca3a504b8e3d080988e9947f", "", boost::python::no_init);
    class_f2215828ca3a504b8e3d080988e9947f.def(boost::python::init<  >(""));
    class_f2215828ca3a504b8e3d080988e9947f.def(boost::python::init< class ::std::allocator< class ::stat_tool::FrequencyDistribution > const & >(""));

}