#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::allocator< class ::stat_tool::DiscreteParametric > const volatile * get_pointer<class ::std::allocator< class ::stat_tool::DiscreteParametric > const volatile >(class ::std::allocator< class ::stat_tool::DiscreteParametric > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_1f29b58eee9156678af1ee8254be47f4()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< class ::std::allocator< class ::stat_tool::DiscreteParametric >, autowig::Held< class ::std::allocator< class ::stat_tool::DiscreteParametric > >::Type > class_1f29b58eee9156678af1ee8254be47f4("_Allocator_1f29b58eee9156678af1ee8254be47f4", "", boost::python::no_init);
    class_1f29b58eee9156678af1ee8254be47f4.def(boost::python::init<  >(""));
    class_1f29b58eee9156678af1ee8254be47f4.def(boost::python::init< class ::std::allocator< class ::stat_tool::DiscreteParametric > const & >(""));

}