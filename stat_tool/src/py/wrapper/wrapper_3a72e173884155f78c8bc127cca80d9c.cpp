#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::ios_base::Init const volatile * get_pointer<class ::std::ios_base::Init const volatile >(class ::std::ios_base::Init const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_3a72e173884155f78c8bc127cca80d9c()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    std::string name_5647113ef4105dfab0588ffcaf6c479b = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + "._ios_base");
    boost::python::object module_5647113ef4105dfab0588ffcaf6c479b(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_5647113ef4105dfab0588ffcaf6c479b.c_str()))));
    boost::python::scope().attr("_ios_base") = module_5647113ef4105dfab0588ffcaf6c479b;
    boost::python::scope scope_5647113ef4105dfab0588ffcaf6c479b = module_5647113ef4105dfab0588ffcaf6c479b;
    boost::python::class_< class ::std::ios_base::Init, autowig::Held< class ::std::ios_base::Init >::Type > class_3a72e173884155f78c8bc127cca80d9c("Init", "", boost::python::no_init);
    class_3a72e173884155f78c8bc127cca80d9c.def(boost::python::init<  >(""));

}