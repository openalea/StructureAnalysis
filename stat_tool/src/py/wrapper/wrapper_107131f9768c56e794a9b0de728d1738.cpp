#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_c38bd1a9de165b4a86de6b417f764214(class ::std::vector< double, class ::std::allocator< double > > & instance, ::std::vector< double, class ::std::allocator< double > >::size_type  param_in_0, double param_out) { instance.operator[](param_in_0) = param_out; }
    void method_decorator_beaa574a87dd59068311360e573c777c(class ::std::vector< double, class ::std::allocator< double > > & instance, ::std::vector< double, class ::std::allocator< double > >::size_type  param_in_0, double param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_22f52ccf4b785404ada1e9bcb9fa01aa(class ::std::vector< double, class ::std::allocator< double > > & instance, double param_out) { instance.front() = param_out; }
    void method_decorator_cae6ded0197f5723954a48c02008c60c(class ::std::vector< double, class ::std::allocator< double > > & instance, double param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< double, class ::std::allocator< double > > const volatile * get_pointer<class ::std::vector< double, class ::std::allocator< double > > const volatile >(class ::std::vector< double, class ::std::allocator< double > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_107131f9768c56e794a9b0de728d1738()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_3fcbcdd89dc951b39073e310b61a04e9)(::std::vector< double, class ::std::allocator< double > >::size_type , ::std::vector< double, class ::std::allocator< double > >::value_type const &) = &::std::vector< double, class ::std::allocator< double > >::assign;
    ::std::vector< double, class ::std::allocator< double > >::size_type  (::std::vector< double, ::std::allocator< double > >::*method_pointer_2762c9c3801851f798933620de2a38c5)() const = &::std::vector< double, class ::std::allocator< double > >::size;
    ::std::vector< double, class ::std::allocator< double > >::size_type  (::std::vector< double, ::std::allocator< double > >::*method_pointer_441a8fd7a6105f52ab269c45a59d8985)() const = &::std::vector< double, class ::std::allocator< double > >::max_size;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_cc6b2e20d84b5d1fb99f0ecb401c5dbe)(::std::vector< double, class ::std::allocator< double > >::size_type ) = &::std::vector< double, class ::std::allocator< double > >::resize;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_d9f1c936f872578eb9f098967a8c1d30)(::std::vector< double, class ::std::allocator< double > >::size_type , ::std::vector< double, class ::std::allocator< double > >::value_type const &) = &::std::vector< double, class ::std::allocator< double > >::resize;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_974010b191b3526f93b53ba098c213a9)() = &::std::vector< double, class ::std::allocator< double > >::shrink_to_fit;
    ::std::vector< double, class ::std::allocator< double > >::size_type  (::std::vector< double, ::std::allocator< double > >::*method_pointer_5c9706119c135c8ca0dbcadbab171935)() const = &::std::vector< double, class ::std::allocator< double > >::capacity;
    bool  (::std::vector< double, ::std::allocator< double > >::*method_pointer_fbb341f4fc855b39aa2eeb3c29aefdd1)() const = &::std::vector< double, class ::std::allocator< double > >::empty;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_56d1f2045f5952b983dea94ad5c12052)(::std::vector< double, class ::std::allocator< double > >::size_type ) = &::std::vector< double, class ::std::allocator< double > >::reserve;
    ::std::vector< double, class ::std::allocator< double > >::reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_c38bd1a9de165b4a86de6b417f764214)(::std::vector< double, class ::std::allocator< double > >::size_type ) = &::std::vector< double, class ::std::allocator< double > >::operator[];
    ::std::vector< double, class ::std::allocator< double > >::const_reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_6ed16026f8925c70aefab0853e8876a8)(::std::vector< double, class ::std::allocator< double > >::size_type ) const = &::std::vector< double, class ::std::allocator< double > >::operator[];
    ::std::vector< double, class ::std::allocator< double > >::reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_beaa574a87dd59068311360e573c777c)(::std::vector< double, class ::std::allocator< double > >::size_type ) = &::std::vector< double, class ::std::allocator< double > >::at;
    ::std::vector< double, class ::std::allocator< double > >::const_reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_25da210bb3c75c269af3895505a003fd)(::std::vector< double, class ::std::allocator< double > >::size_type ) const = &::std::vector< double, class ::std::allocator< double > >::at;
    ::std::vector< double, class ::std::allocator< double > >::reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_22f52ccf4b785404ada1e9bcb9fa01aa)() = &::std::vector< double, class ::std::allocator< double > >::front;
    ::std::vector< double, class ::std::allocator< double > >::const_reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_fb2421b769975570ac8387adada3a0b2)() const = &::std::vector< double, class ::std::allocator< double > >::front;
    ::std::vector< double, class ::std::allocator< double > >::reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_cae6ded0197f5723954a48c02008c60c)() = &::std::vector< double, class ::std::allocator< double > >::back;
    ::std::vector< double, class ::std::allocator< double > >::const_reference  (::std::vector< double, ::std::allocator< double > >::*method_pointer_e34e38388a585de3aaa97622256cd858)() const = &::std::vector< double, class ::std::allocator< double > >::back;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_f6adc78c87285d669260e4289dd02852)(::std::vector< double, class ::std::allocator< double > >::value_type const &) = &::std::vector< double, class ::std::allocator< double > >::push_back;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_0d34433fdbe7506180094ce87fae6cab)() = &::std::vector< double, class ::std::allocator< double > >::pop_back;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_331520d5c9a25683a40120b524320d20)(class ::std::vector< double, class ::std::allocator< double > > &) = &::std::vector< double, class ::std::allocator< double > >::swap;
    void  (::std::vector< double, ::std::allocator< double > >::*method_pointer_a4208670622158c3a2376d15176ceee2)() = &::std::vector< double, class ::std::allocator< double > >::clear;
    boost::python::class_< class ::std::vector< double, class ::std::allocator< double > >, autowig::Held< class ::std::vector< double, class ::std::allocator< double > > >::Type > class_107131f9768c56e794a9b0de728d1738("_Vector_107131f9768c56e794a9b0de728d1738", "", boost::python::no_init);
    class_107131f9768c56e794a9b0de728d1738.def(boost::python::init<  >(""));
    class_107131f9768c56e794a9b0de728d1738.def(boost::python::init< ::std::vector< double, class ::std::allocator< double > >::allocator_type const & >(""));
    class_107131f9768c56e794a9b0de728d1738.def(boost::python::init< ::std::vector< double, class ::std::allocator< double > >::size_type , ::std::vector< double, class ::std::allocator< double > >::allocator_type const & >(""));
    class_107131f9768c56e794a9b0de728d1738.def(boost::python::init< ::std::vector< double, class ::std::allocator< double > >::size_type , ::std::vector< double, class ::std::allocator< double > >::value_type const &, ::std::vector< double, class ::std::allocator< double > >::allocator_type const & >(""));
    class_107131f9768c56e794a9b0de728d1738.def(boost::python::init< class ::std::vector< double, class ::std::allocator< double > > const & >(""));
    class_107131f9768c56e794a9b0de728d1738.def(boost::python::init< class ::std::vector< double, class ::std::allocator< double > > const &, ::std::vector< double, class ::std::allocator< double > >::allocator_type const & >(""));
    class_107131f9768c56e794a9b0de728d1738.def("assign", method_pointer_3fcbcdd89dc951b39073e310b61a04e9, "");
    class_107131f9768c56e794a9b0de728d1738.def("__len__", method_pointer_2762c9c3801851f798933620de2a38c5, "");
    class_107131f9768c56e794a9b0de728d1738.def("max_size", method_pointer_441a8fd7a6105f52ab269c45a59d8985, "");
    class_107131f9768c56e794a9b0de728d1738.def("resize", method_pointer_cc6b2e20d84b5d1fb99f0ecb401c5dbe, "");
    class_107131f9768c56e794a9b0de728d1738.def("resize", method_pointer_d9f1c936f872578eb9f098967a8c1d30, "");
    class_107131f9768c56e794a9b0de728d1738.def("shrink_to_fit", method_pointer_974010b191b3526f93b53ba098c213a9, "");
    class_107131f9768c56e794a9b0de728d1738.def("capacity", method_pointer_5c9706119c135c8ca0dbcadbab171935, "");
    class_107131f9768c56e794a9b0de728d1738.def("empty", method_pointer_fbb341f4fc855b39aa2eeb3c29aefdd1, "");
    class_107131f9768c56e794a9b0de728d1738.def("reserve", method_pointer_56d1f2045f5952b983dea94ad5c12052, "");
    class_107131f9768c56e794a9b0de728d1738.def("__getitem__", method_pointer_c38bd1a9de165b4a86de6b417f764214, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("__getitem__", autowig::method_decorator_c38bd1a9de165b4a86de6b417f764214);
    class_107131f9768c56e794a9b0de728d1738.def("__getitem__", method_pointer_6ed16026f8925c70aefab0853e8876a8, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("at", method_pointer_beaa574a87dd59068311360e573c777c, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("at", autowig::method_decorator_beaa574a87dd59068311360e573c777c);
    class_107131f9768c56e794a9b0de728d1738.def("at", method_pointer_25da210bb3c75c269af3895505a003fd, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("front", method_pointer_22f52ccf4b785404ada1e9bcb9fa01aa, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("front", autowig::method_decorator_22f52ccf4b785404ada1e9bcb9fa01aa);
    class_107131f9768c56e794a9b0de728d1738.def("front", method_pointer_fb2421b769975570ac8387adada3a0b2, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("back", method_pointer_cae6ded0197f5723954a48c02008c60c, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("back", autowig::method_decorator_cae6ded0197f5723954a48c02008c60c);
    class_107131f9768c56e794a9b0de728d1738.def("back", method_pointer_e34e38388a585de3aaa97622256cd858, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_107131f9768c56e794a9b0de728d1738.def("push_back", method_pointer_f6adc78c87285d669260e4289dd02852, "");
    class_107131f9768c56e794a9b0de728d1738.def("pop_back", method_pointer_0d34433fdbe7506180094ce87fae6cab, "");
    class_107131f9768c56e794a9b0de728d1738.def("swap", method_pointer_331520d5c9a25683a40120b524320d20, "");
    class_107131f9768c56e794a9b0de728d1738.def("clear", method_pointer_a4208670622158c3a2376d15176ceee2, "");

    struct vector_107131f9768c56e794a9b0de728d1738_from_python
    {
        vector_107131f9768c56e794a9b0de728d1738_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< double, class ::std::allocator< double > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< double, class ::std::allocator< double > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< double, class ::std::allocator< double > >();
            data->convertible = storage;
            class ::std::vector< double, class ::std::allocator< double > >& result = *((class ::std::vector< double, class ::std::allocator< double > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< double  >(py_elem_obj));
            }
        }
    };

    vector_107131f9768c56e794a9b0de728d1738_from_python();
}