#include "_stat_tool.h"



namespace autowig
{

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::stat_tool::SinglePlot const volatile * get_pointer<class ::stat_tool::SinglePlot const volatile >(class ::stat_tool::SinglePlot const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_bd64fe19cf7b5190a3344dcce6986fc3()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    void  (::stat_tool::SinglePlot::*method_pointer_6a7e6a436fdf5223936bab1c87e3a925)(float , float ) = &::stat_tool::SinglePlot::add_point;
    void  (::stat_tool::SinglePlot::*method_pointer_263879913028539985e7fac321db4e3b)(::stat_tool::PlotPoint const &) = &::stat_tool::SinglePlot::add_point;
    void  (::stat_tool::SinglePlot::*method_pointer_03589cd205885ff4909a4223cd2a2476)(float , float , class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const &) = &::stat_tool::SinglePlot::add_text;
    float  (::stat_tool::SinglePlot::*method_pointer_0284f38ae21d5fbe827f5eab44dce8cb)(int ) = &::stat_tool::SinglePlot::get_x;
    float  (::stat_tool::SinglePlot::*method_pointer_4dac2808834d5d63ad5a8a2212172f8e)(int ) = &::stat_tool::SinglePlot::get_y;
    class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > >  (::stat_tool::SinglePlot::*method_pointer_2f87b39079bc5033b75d79a6ff4ccc2e)(int ) = &::stat_tool::SinglePlot::get_label;
    int  (::stat_tool::SinglePlot::*method_pointer_2ec6352b78ff512f84da5d6d1caed6a9)() = &::stat_tool::SinglePlot::get_size;
    int  (::stat_tool::SinglePlot::*method_pointer_9bfcc702128d5ec79deb061830e16ef4)() = &::stat_tool::SinglePlot::size;
    boost::python::class_< class ::stat_tool::SinglePlot, autowig::Held< class ::stat_tool::SinglePlot >::Type > class_bd64fe19cf7b5190a3344dcce6986fc3("SinglePlot", "", boost::python::no_init);
    class_bd64fe19cf7b5190a3344dcce6986fc3.def(boost::python::init<  >(""));
    class_bd64fe19cf7b5190a3344dcce6986fc3.def(boost::python::init< class ::stat_tool::SinglePlot const & >(""));
    class_bd64fe19cf7b5190a3344dcce6986fc3.def("add_point", method_pointer_6a7e6a436fdf5223936bab1c87e3a925, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def("add_point", method_pointer_263879913028539985e7fac321db4e3b, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def("add_text", method_pointer_03589cd205885ff4909a4223cd2a2476, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def("get_x", method_pointer_0284f38ae21d5fbe827f5eab44dce8cb, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def("get_y", method_pointer_4dac2808834d5d63ad5a8a2212172f8e, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def("get_label", method_pointer_2f87b39079bc5033b75d79a6ff4ccc2e, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def("get_size", method_pointer_2ec6352b78ff512f84da5d6d1caed6a9, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def("__len__", method_pointer_9bfcc702128d5ec79deb061830e16ef4, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def_readwrite("data", &::stat_tool::SinglePlot::data, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def_readwrite("data_and_text", &::stat_tool::SinglePlot::data_and_text, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def_readwrite("legend", &::stat_tool::SinglePlot::legend, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def_readwrite("style", &::stat_tool::SinglePlot::style, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def_readwrite("color", &::stat_tool::SinglePlot::color, "");
    class_bd64fe19cf7b5190a3344dcce6986fc3.def_readwrite("label", &::stat_tool::SinglePlot::label, "");

}