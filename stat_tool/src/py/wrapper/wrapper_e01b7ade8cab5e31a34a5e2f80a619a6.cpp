#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::allocator< int > const volatile * get_pointer<class ::std::allocator< int > const volatile >(class ::std::allocator< int > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_e01b7ade8cab5e31a34a5e2f80a619a6()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::class_< class ::std::allocator< int >, autowig::Held< class ::std::allocator< int > >::Type > class_e01b7ade8cab5e31a34a5e2f80a619a6("_Allocator_e01b7ade8cab5e31a34a5e2f80a619a6", "", boost::python::no_init);
    class_e01b7ade8cab5e31a34a5e2f80a619a6.def(boost::python::init<  >(""));
    class_e01b7ade8cab5e31a34a5e2f80a619a6.def(boost::python::init< class ::std::allocator< int > const & >(""));

}