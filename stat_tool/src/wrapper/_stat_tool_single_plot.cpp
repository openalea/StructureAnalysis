#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_single_plot()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        void (::stat_tool::SinglePlot::*method_pointer_6a7e6a436fdf5223936bab1c87e3a925)(float, float) = &::stat_tool::SinglePlot::add_point;
        void (::stat_tool::SinglePlot::*method_pointer_844b03359db15d14a628d52100fd1a25)(struct ::std::pair<float, float> const &) = &::stat_tool::SinglePlot::add_point;
        void (::stat_tool::SinglePlot::*method_pointer_303fb6128f8e56a9bcbeda47ac209905)(float, float, class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > const &) = &::stat_tool::SinglePlot::add_text;
        float (::stat_tool::SinglePlot::*method_pointer_0284f38ae21d5fbe827f5eab44dce8cb)(int) = &::stat_tool::SinglePlot::get_x;
        float (::stat_tool::SinglePlot::*method_pointer_4dac2808834d5d63ad5a8a2212172f8e)(int) = &::stat_tool::SinglePlot::get_y;
        class ::std::basic_string<char, std::char_traits<char>, std::allocator<char> > (::stat_tool::SinglePlot::*method_pointer_2f87b39079bc5033b75d79a6ff4ccc2e)(int) = &::stat_tool::SinglePlot::get_label;
        int (::stat_tool::SinglePlot::*method_pointer_2ec6352b78ff512f84da5d6d1caed6a9)() = &::stat_tool::SinglePlot::get_size;
        int (::stat_tool::SinglePlot::*method_pointer_9bfcc702128d5ec79deb061830e16ef4)() = &::stat_tool::SinglePlot::size;
        struct ::std::_List_const_iterator<std::pair<float, float> > (::stat_tool::SinglePlot::*method_pointer_48bd08842095501da93a8d4ef039932e)() = &::stat_tool::SinglePlot::begin;
        struct ::std::_List_const_iterator<std::pair<float, float> > (::stat_tool::SinglePlot::*method_pointer_73c1e2ab42fe5cfda0f730a9dcb393cd)() = &::stat_tool::SinglePlot::end;
        boost::python::class_< class ::stat_tool::SinglePlot, std::shared_ptr< class ::stat_tool::SinglePlot > >("SinglePlot", boost::python::no_init)
            .def(boost::python::init<  >())
            .def(boost::python::init< class ::stat_tool::SinglePlot const & >())
            .def("add_point", method_pointer_6a7e6a436fdf5223936bab1c87e3a925)
            .def("add_point", method_pointer_844b03359db15d14a628d52100fd1a25)
            .def("add_text", method_pointer_303fb6128f8e56a9bcbeda47ac209905)
            .def("get_x", method_pointer_0284f38ae21d5fbe827f5eab44dce8cb)
            .def("get_y", method_pointer_4dac2808834d5d63ad5a8a2212172f8e)
            .def("get_label", method_pointer_2f87b39079bc5033b75d79a6ff4ccc2e)
            .def("get_size", method_pointer_2ec6352b78ff512f84da5d6d1caed6a9)
            .def("size", method_pointer_9bfcc702128d5ec79deb061830e16ef4)
            .def("begin", method_pointer_48bd08842095501da93a8d4ef039932e)
            .def("end", method_pointer_73c1e2ab42fe5cfda0f730a9dcb393cd)
            .def_readwrite("data", &::stat_tool::SinglePlot::data)
            .def_readwrite("data_and_text", &::stat_tool::SinglePlot::data_and_text)
            .def_readwrite("legend", &::stat_tool::SinglePlot::legend)
            .def_readwrite("style", &::stat_tool::SinglePlot::style)
            .def_readwrite("color", &::stat_tool::SinglePlot::color)
            .def_readwrite("label", &::stat_tool::SinglePlot::label);
}