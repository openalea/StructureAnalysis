#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::basic_ostream< char, struct ::std::char_traits< char > >::sentry const volatile * get_pointer<class ::std::basic_ostream< char, struct ::std::char_traits< char > >::sentry const volatile >(class ::std::basic_ostream< char, struct ::std::char_traits< char > >::sentry const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_79c03425b8505668b16ffdd958127107()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    std::string name_e1391944268253558f04b6f996bb5a8b = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".__basic_ostream_e1391944268253558f04b6f996bb5a8b");
    boost::python::object module_e1391944268253558f04b6f996bb5a8b(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_e1391944268253558f04b6f996bb5a8b.c_str()))));
    boost::python::scope().attr("__basic_ostream_e1391944268253558f04b6f996bb5a8b") = module_e1391944268253558f04b6f996bb5a8b;
    boost::python::scope scope_e1391944268253558f04b6f996bb5a8b = module_e1391944268253558f04b6f996bb5a8b;
    boost::python::class_< class ::std::basic_ostream< char, struct ::std::char_traits< char > >::sentry, autowig::Held< class ::std::basic_ostream< char, struct ::std::char_traits< char > >::sentry >::Type, boost::noncopyable > class_79c03425b8505668b16ffdd958127107("Sentry", "", boost::python::no_init);
    class_79c03425b8505668b16ffdd958127107.def(boost::python::init< class ::std::basic_ostream< char, struct ::std::char_traits< char > > & >(""));

}