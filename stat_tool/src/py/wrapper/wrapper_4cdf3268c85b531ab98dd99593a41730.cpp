#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::StatError const volatile * get_pointer<class ::stat_tool::StatError const volatile >(class ::stat_tool::StatError const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_4cdf3268c85b531ab98dd99593a41730()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > >  (::stat_tool::StatError::*method_pointer_3d5f3866720253f79f91599400394e9e)(enum ::stat_tool::error_type ) const = &::stat_tool::StatError::ascii_write;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::StatError::*method_pointer_1b338da12c2f5a6a8d94f1d3bc88e918)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, enum ::stat_tool::error_type ) const = &::stat_tool::StatError::ascii_write;
    void  (::stat_tool::StatError::*method_pointer_6aca1e4328215b179690c13e4babbc77)() = &::stat_tool::StatError::init;
    int  (::stat_tool::StatError::*method_pointer_428a45e45b335e5ca2005b2914d07101)() const = &::stat_tool::StatError::get_nb_error;
    int  (::stat_tool::StatError::*method_pointer_9007bd748c2c5e138d772eb3fd7ec146)() const = &::stat_tool::StatError::get_max_nb_error;
    boost::python::class_< class ::stat_tool::StatError, autowig::Held< class ::stat_tool::StatError >::Type > class_4cdf3268c85b531ab98dd99593a41730("StatError", "Class for error management\n\n", boost::python::no_init);
    class_4cdf3268c85b531ab98dd99593a41730.def(boost::python::init< int  >(""));
    class_4cdf3268c85b531ab98dd99593a41730.def("ascii_write", method_pointer_3d5f3866720253f79f91599400394e9e, "");
    class_4cdf3268c85b531ab98dd99593a41730.def("ascii_write", method_pointer_1b338da12c2f5a6a8d94f1d3bc88e918, boost::python::return_internal_reference<>(), "");
    class_4cdf3268c85b531ab98dd99593a41730.def("init", method_pointer_6aca1e4328215b179690c13e4babbc77, "");
    class_4cdf3268c85b531ab98dd99593a41730.def("get_nb_error", method_pointer_428a45e45b335e5ca2005b2914d07101, "");
    class_4cdf3268c85b531ab98dd99593a41730.def("get_max_nb_error", method_pointer_9007bd748c2c5e138d772eb3fd7ec146, "");

}