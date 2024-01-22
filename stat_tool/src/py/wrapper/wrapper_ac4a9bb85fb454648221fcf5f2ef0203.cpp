#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> struct ::std::allocator_traits< class ::std::allocator< char > > const volatile * get_pointer<struct ::std::allocator_traits< class ::std::allocator< char > > const volatile >(struct ::std::allocator_traits< class ::std::allocator< char > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_ac4a9bb85fb454648221fcf5f2ef0203()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< struct ::std::allocator_traits< class ::std::allocator< char > >, autowig::Held< struct ::std::allocator_traits< class ::std::allocator< char > > >::Type > class_ac4a9bb85fb454648221fcf5f2ef0203("_AllocatorTraits_ac4a9bb85fb454648221fcf5f2ef0203", "", boost::python::no_init);

}