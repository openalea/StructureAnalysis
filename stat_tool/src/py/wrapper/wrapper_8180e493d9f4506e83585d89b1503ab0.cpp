#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::allocator< enum ::stat_tool::process_distribution > const volatile * get_pointer<class ::std::allocator< enum ::stat_tool::process_distribution > const volatile >(class ::std::allocator< enum ::stat_tool::process_distribution > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_8180e493d9f4506e83585d89b1503ab0()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< class ::std::allocator< enum ::stat_tool::process_distribution >, autowig::Held< class ::std::allocator< enum ::stat_tool::process_distribution > >::Type > class_8180e493d9f4506e83585d89b1503ab0("_Allocator_8180e493d9f4506e83585d89b1503ab0", "", boost::python::no_init);
    class_8180e493d9f4506e83585d89b1503ab0.def(boost::python::init<  >(""));
    class_8180e493d9f4506e83585d89b1503ab0.def(boost::python::init< class ::std::allocator< enum ::stat_tool::process_distribution > const & >(""));

}