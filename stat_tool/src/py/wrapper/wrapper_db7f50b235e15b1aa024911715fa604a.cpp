#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_f0d4c383754f5a8bbec28ab41c2917cd(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > & instance, ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  param_out) { instance.operator[](param_in_0) = param_out; }
    void method_decorator_f86e78795f185656ac6bee78a061be87(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > & instance, ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_ca5b43b1a06d5433af0db92d4696d989(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > & instance, const ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_638cfa59ea435ba1809aa29a95560e18(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > & instance, const ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > const volatile * get_pointer<class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > const volatile >(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_db7f50b235e15b1aa024911715fa604a()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_98cad235f33b5012ada8f95493b92544)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type , ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::value_type const &) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::assign;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_c632a5e3dd8c5dda9749ea255f35edc2)(class ::std::initializer_list< class ::stat_tool::Vectors > ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::assign;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_f67ce359a30f58d28fc366cdb4586cb5)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_f5e1d945560e5c6f88beb5ff37891076)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::max_size;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_663ee50fa32752d6a0e5fd0d42c3ef27)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::resize;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_d7fa969c8ac55156ad585b75052d6813)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type , ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::value_type const &) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::resize;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_234f087679ed5e8f893a55d166cfb7e4)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::shrink_to_fit;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_240c6bd6df49518abe5eb58a8a9f2b69)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::capacity;
    bool  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_e40e6bfdea1f5d8f956e09439520ccff)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::empty;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_6b4fc5e0cbf9532ba596ea620e8a320d)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reserve;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_f0d4c383754f5a8bbec28ab41c2917cd)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::operator[];
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_4a294f10be4e590bbf39c732f582c50d)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::operator[];
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_f86e78795f185656ac6bee78a061be87)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::at;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_bb66332ae2af57cfbeee55facb131d36)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::at;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_ca5b43b1a06d5433af0db92d4696d989)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::front;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_4ba11bae9fe355c59d868c0fb5e5c153)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::front;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_638cfa59ea435ba1809aa29a95560e18)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::back;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_17cf274358a750c7b74c9180e1c08722)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::back;
    class ::stat_tool::Vectors * (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_06a9b89cac2253798b62400a3d3cd2a7)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::data;
    class ::stat_tool::Vectors const * (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_d249b5a65b1c56358bd76c81a737ca70)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::data;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_6c5f935b824250eb952d96814b487b52)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::value_type const &) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::push_back;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_c00e309f8e635de2a9e0758330f081c5)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::pop_back;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_63a6c8fedcf15f55afbd0706c0d0cdcb)(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > &) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::swap;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_956846ac8e1557aba3fab6a2431160ef)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::clear;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_62f66319367f56c49ff63b245b9a0486)(class ::stat_tool::Vectors const *, class ::stat_tool::Vectors const *) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::assign;
    boost::python::class_< class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >, autowig::Held< class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > >::Type > class_db7f50b235e15b1aa024911715fa604a("_Vector_db7f50b235e15b1aa024911715fa604a", "", boost::python::no_init);
    class_db7f50b235e15b1aa024911715fa604a.def(boost::python::init<  >(""));
    class_db7f50b235e15b1aa024911715fa604a.def(boost::python::init< ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::allocator_type const & >(""));
    class_db7f50b235e15b1aa024911715fa604a.def(boost::python::init< ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type , ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::allocator_type const & >(""));
    class_db7f50b235e15b1aa024911715fa604a.def(boost::python::init< ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type , ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::value_type const &, ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::allocator_type const & >(""));
    class_db7f50b235e15b1aa024911715fa604a.def(boost::python::init< class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > const & >(""));
    class_db7f50b235e15b1aa024911715fa604a.def(boost::python::init< class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > const &, ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::allocator_type const & >(""));
    class_db7f50b235e15b1aa024911715fa604a.def(boost::python::init< class ::std::initializer_list< class ::stat_tool::Vectors > , ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::allocator_type const & >(""));
    class_db7f50b235e15b1aa024911715fa604a.def("assign", method_pointer_98cad235f33b5012ada8f95493b92544, "");
    class_db7f50b235e15b1aa024911715fa604a.def("assign", method_pointer_c632a5e3dd8c5dda9749ea255f35edc2, "");
    class_db7f50b235e15b1aa024911715fa604a.def("__len__", method_pointer_f67ce359a30f58d28fc366cdb4586cb5, "");
    class_db7f50b235e15b1aa024911715fa604a.def("max_size", method_pointer_f5e1d945560e5c6f88beb5ff37891076, "");
    class_db7f50b235e15b1aa024911715fa604a.def("resize", method_pointer_663ee50fa32752d6a0e5fd0d42c3ef27, "");
    class_db7f50b235e15b1aa024911715fa604a.def("resize", method_pointer_d7fa969c8ac55156ad585b75052d6813, "");
    class_db7f50b235e15b1aa024911715fa604a.def("shrink_to_fit", method_pointer_234f087679ed5e8f893a55d166cfb7e4, "");
    class_db7f50b235e15b1aa024911715fa604a.def("capacity", method_pointer_240c6bd6df49518abe5eb58a8a9f2b69, "");
    class_db7f50b235e15b1aa024911715fa604a.def("empty", method_pointer_e40e6bfdea1f5d8f956e09439520ccff, "");
    class_db7f50b235e15b1aa024911715fa604a.def("reserve", method_pointer_6b4fc5e0cbf9532ba596ea620e8a320d, "");
    class_db7f50b235e15b1aa024911715fa604a.def("__getitem__", method_pointer_f0d4c383754f5a8bbec28ab41c2917cd, boost::python::return_internal_reference<>(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("__getitem__", autowig::method_decorator_f0d4c383754f5a8bbec28ab41c2917cd);
    class_db7f50b235e15b1aa024911715fa604a.def("__getitem__", method_pointer_4a294f10be4e590bbf39c732f582c50d, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("at", method_pointer_f86e78795f185656ac6bee78a061be87, boost::python::return_internal_reference<>(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("at", autowig::method_decorator_f86e78795f185656ac6bee78a061be87);
    class_db7f50b235e15b1aa024911715fa604a.def("at", method_pointer_bb66332ae2af57cfbeee55facb131d36, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("front", method_pointer_ca5b43b1a06d5433af0db92d4696d989, boost::python::return_internal_reference<>(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("front", autowig::method_decorator_ca5b43b1a06d5433af0db92d4696d989);
    class_db7f50b235e15b1aa024911715fa604a.def("front", method_pointer_4ba11bae9fe355c59d868c0fb5e5c153, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("back", method_pointer_638cfa59ea435ba1809aa29a95560e18, boost::python::return_internal_reference<>(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("back", autowig::method_decorator_638cfa59ea435ba1809aa29a95560e18);
    class_db7f50b235e15b1aa024911715fa604a.def("back", method_pointer_17cf274358a750c7b74c9180e1c08722, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("data", method_pointer_06a9b89cac2253798b62400a3d3cd2a7, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("data", method_pointer_d249b5a65b1c56358bd76c81a737ca70, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("push_back", method_pointer_6c5f935b824250eb952d96814b487b52, "");
    class_db7f50b235e15b1aa024911715fa604a.def("pop_back", method_pointer_c00e309f8e635de2a9e0758330f081c5, "");
    class_db7f50b235e15b1aa024911715fa604a.def("swap", method_pointer_63a6c8fedcf15f55afbd0706c0d0cdcb, "");
    class_db7f50b235e15b1aa024911715fa604a.def("clear", method_pointer_956846ac8e1557aba3fab6a2431160ef, "");
    class_db7f50b235e15b1aa024911715fa604a.def("assign", method_pointer_62f66319367f56c49ff63b245b9a0486, "");

    struct vector_db7f50b235e15b1aa024911715fa604a_from_python
    {
        vector_db7f50b235e15b1aa024911715fa604a_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >();
            data->convertible = storage;
            class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >& result = *((class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< class ::stat_tool::Vectors  >(py_elem_obj));
            }
        }
    };

    vector_db7f50b235e15b1aa024911715fa604a_from_python();
}