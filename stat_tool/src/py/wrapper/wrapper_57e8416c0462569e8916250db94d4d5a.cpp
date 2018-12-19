#include "_stat_tool.h"



namespace autowig
{
    class Wrap_57e8416c0462569e8916250db94d4d5a : public ::stat_tool::StatInterface, public boost::python::wrapper< class ::stat_tool::StatInterface >
    {
        public:
            
            virtual bool  plot_write(class ::stat_tool::StatError & param_0, char const * param_1, char const * param_2) const
            { return this->get_override("plot_write")(param_0, param_1, param_2); }
                        
            virtual bool  spreadsheet_write(class ::stat_tool::StatError & param_0, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const param_1) const
            { return this->get_override("spreadsheet_write")(param_0, param_1); }
                        
            virtual bool  ascii_write(class ::stat_tool::StatError & param_0, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const param_1, bool  param_2) const
            { return this->get_override("ascii_write")(param_0, param_1, param_2); }
                        
            virtual class ::std::basic_ostream< char, struct ::std::char_traits< char > > & ascii_write(class ::std::basic_ostream< char, struct ::std::char_traits< char > > & param_0, bool  param_1) const
            {
                 ::std::basic_ostream< char, ::std::char_traits< char > >* result = this->get_override("ascii_write")(param_0, param_1);
                 return *result;
            }                 
                        
            virtual class ::std::basic_ostream< char, struct ::std::char_traits< char > > & line_write(class ::std::basic_ostream< char, struct ::std::char_traits< char > > & param_0) const
            {
                 ::std::basic_ostream< char, ::std::char_traits< char > >* result = this->get_override("line_write")(param_0);
                 return *result;
            }                 
                        

        protected:
            

        private:
            

    };

}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> autowig::Wrap_57e8416c0462569e8916250db94d4d5a const volatile * get_pointer<autowig::Wrap_57e8416c0462569e8916250db94d4d5a const volatile >(autowig::Wrap_57e8416c0462569e8916250db94d4d5a const volatile *c) { return c; }
    template <> class ::stat_tool::StatInterface const volatile * get_pointer<class ::stat_tool::StatInterface const volatile >(class ::stat_tool::StatInterface const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_57e8416c0462569e8916250db94d4d5a()
{

    std::string name_0cdd446515295e8e8373e99f328c3748 = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
    boost::python::object module_0cdd446515295e8e8373e99f328c3748(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_0cdd446515295e8e8373e99f328c3748.c_str()))));
    boost::python::scope().attr("stat_tool") = module_0cdd446515295e8e8373e99f328c3748;
    boost::python::scope scope_0cdd446515295e8e8373e99f328c3748 = module_0cdd446515295e8e8373e99f328c3748;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::StatInterface::*method_pointer_109e8f598f3c51bca3b606de4d27c07e)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &) const = &::stat_tool::StatInterface::line_write;
    class ::std::basic_ostream< char, struct ::std::char_traits< char > > & (::stat_tool::StatInterface::*method_pointer_5f38643402d45a93b04e0b6016f90e73)(class ::std::basic_ostream< char, struct ::std::char_traits< char > > &, bool ) const = &::stat_tool::StatInterface::ascii_write;
    class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > >  (::stat_tool::StatInterface::*method_pointer_a928a20806fb5a658f4b4fa2eb451e9d)(bool ) const = &::stat_tool::StatInterface::ascii_write;
    bool  (::stat_tool::StatInterface::*method_pointer_68c734f08f13501ba14d57e85abbc5e4)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const, bool ) const = &::stat_tool::StatInterface::ascii_write;
    bool  (::stat_tool::StatInterface::*method_pointer_358914f49a5c5b7fabe0adf15f780079)(class ::stat_tool::StatError &, class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > > const) const = &::stat_tool::StatInterface::spreadsheet_write;
    ::stat_tool::MultiPlotSet * (::stat_tool::StatInterface::*method_pointer_dcd5f5dbeaed52b1a46267335ba3f07b)() const = &::stat_tool::StatInterface::get_plotable;
    boost::python::class_< autowig::Wrap_57e8416c0462569e8916250db94d4d5a, autowig::Held< autowig::Wrap_57e8416c0462569e8916250db94d4d5a >::Type, boost::noncopyable > class_57e8416c0462569e8916250db94d4d5a("StatInterface", "Abstract class defining a common interface\n\n", boost::python::no_init);
    class_57e8416c0462569e8916250db94d4d5a.def("line_write", boost::python::pure_virtual(method_pointer_109e8f598f3c51bca3b606de4d27c07e), boost::python::return_internal_reference<>(), "");
    class_57e8416c0462569e8916250db94d4d5a.def("ascii_write", boost::python::pure_virtual(method_pointer_5f38643402d45a93b04e0b6016f90e73), boost::python::return_internal_reference<>(), "");
    class_57e8416c0462569e8916250db94d4d5a.def("ascii_write", method_pointer_a928a20806fb5a658f4b4fa2eb451e9d, "");
    class_57e8416c0462569e8916250db94d4d5a.def("ascii_write", boost::python::pure_virtual(method_pointer_68c734f08f13501ba14d57e85abbc5e4), "");
    class_57e8416c0462569e8916250db94d4d5a.def("spreadsheet_write", boost::python::pure_virtual(method_pointer_358914f49a5c5b7fabe0adf15f780079), "");
    class_57e8416c0462569e8916250db94d4d5a.def("get_plotable", method_pointer_dcd5f5dbeaed52b1a46267335ba3f07b, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    if(autowig::Held< class ::stat_tool::StatInterface >::is_class)
    {
        boost::python::implicitly_convertible< autowig::Held< autowig::Wrap_57e8416c0462569e8916250db94d4d5a >::Type, autowig::Held< class ::stat_tool::StatInterface >::Type >();
        boost::python::register_ptr_to_python< autowig::Held< class ::stat_tool::StatInterface >::Type >();
    }    

}