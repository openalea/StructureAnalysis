#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::allocator< enum ::stat_tool::discrete_parametric > const volatile * get_pointer<class ::std::allocator< enum ::stat_tool::discrete_parametric > const volatile >(class ::std::allocator< enum ::stat_tool::discrete_parametric > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_a8a9526a25e8542295dd2df108ca409c()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< class ::std::allocator< enum ::stat_tool::discrete_parametric >, autowig::Held< class ::std::allocator< enum ::stat_tool::discrete_parametric > >::Type > class_a8a9526a25e8542295dd2df108ca409c("_Allocator_a8a9526a25e8542295dd2df108ca409c", "", boost::python::no_init);
    class_a8a9526a25e8542295dd2df108ca409c.def(boost::python::init<  >(""));
    class_a8a9526a25e8542295dd2df108ca409c.def(boost::python::init< class ::std::allocator< enum ::stat_tool::discrete_parametric > const & >(""));

}