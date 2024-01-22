#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::Test const volatile * get_pointer<class ::stat_tool::Test const volatile >(class ::stat_tool::Test const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_6f8b5760ad2a51adaa337840a6a109d1()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::Test::*method_pointer_5d1e82f49af15250bdfc9ed3c007426f)(class ::stat_tool::Test const &) = &::stat_tool::Test::copy;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Test::*method_pointer_9566c53b184955ba904b8b0abf1007bb)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool , bool ) const = &::stat_tool::Test::ascii_print;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::Test::*method_pointer_d451faba26305b43b033e92605e37fed)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::Test::spreadsheet_print;
    void  (::stat_tool::Test::*method_pointer_fe1275e4546a5da09641c670d5b6177d)() = &::stat_tool::Test::standard_normal_critical_probability_computation;
    void  (::stat_tool::Test::*method_pointer_1828d9c2d2415ae68b778f6d26c05036)() = &::stat_tool::Test::standard_normal_value_computation;
    void  (::stat_tool::Test::*method_pointer_a2f30cdc9f8656f6a5b9b12990f0187c)() = &::stat_tool::Test::chi2_critical_probability_computation;
    void  (::stat_tool::Test::*method_pointer_e2b31f57c62459d2b55949f85f609ab5)() = &::stat_tool::Test::chi2_value_computation;
    void  (::stat_tool::Test::*method_pointer_0496cfc256345921a65582f307882888)() = &::stat_tool::Test::F_critical_probability_computation;
    void  (::stat_tool::Test::*method_pointer_395c5b3a9b045414a12c16407e5beb09)() = &::stat_tool::Test::F_value_computation;
    void  (::stat_tool::Test::*method_pointer_56006b70c92b5107b03c64bf8d3a993a)() = &::stat_tool::Test::t_critical_probability_computation;
    void  (::stat_tool::Test::*method_pointer_6dc9bf3a118e5863904b7e8efab67296)() = &::stat_tool::Test::t_value_computation;
    boost::python::class_< class ::stat_tool::Test, autowig::Held< class ::stat_tool::Test >::Type > class_6f8b5760ad2a51adaa337840a6a109d1("Test", "Test of hypothesis\n\n", boost::python::no_init);
    class_6f8b5760ad2a51adaa337840a6a109d1.def("copy", method_pointer_5d1e82f49af15250bdfc9ed3c007426f, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("ascii_print", method_pointer_9566c53b184955ba904b8b0abf1007bb, boost::python::return_internal_reference<>(), "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("spreadsheet_print", method_pointer_d451faba26305b43b033e92605e37fed, boost::python::return_internal_reference<>(), "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("standard_normal_critical_probability_computation", method_pointer_fe1275e4546a5da09641c670d5b6177d, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("standard_normal_value_computation", method_pointer_1828d9c2d2415ae68b778f6d26c05036, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("chi_2__critical_probability_computation", method_pointer_a2f30cdc9f8656f6a5b9b12990f0187c, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("chi_2__value_computation", method_pointer_e2b31f57c62459d2b55949f85f609ab5, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("f__critical_probability_computation", method_pointer_0496cfc256345921a65582f307882888, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("f__value_computation", method_pointer_395c5b3a9b045414a12c16407e5beb09, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("t_critical_probability_computation", method_pointer_56006b70c92b5107b03c64bf8d3a993a, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def("t_value_computation", method_pointer_6dc9bf3a118e5863904b7e8efab67296, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def_readwrite("ident", &::stat_tool::Test::ident, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def_readwrite("one_side", &::stat_tool::Test::one_side, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def_readwrite("df_1", &::stat_tool::Test::df1, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def_readwrite("df_2", &::stat_tool::Test::df2, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def_readwrite("value", &::stat_tool::Test::value, "");
    class_6f8b5760ad2a51adaa337840a6a109d1.def_readwrite("critical_probability", &::stat_tool::Test::critical_probability, "");

}