#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_65d73148dc4a508b9db1cda860b09e83(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > & instance, ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_85254a212bfb5de2a8d1e892bebff8e4(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > & instance, const ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_9d49cf5d7771526a92d0faa0353ddae6(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > & instance, const ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  param_out) { instance.back() = param_out; }
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
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_83baa9fc9c1157d9b24615b6b82790b7)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type , ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::assign;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_8c95d46dc6e35c1a92ed5968642ca1ad)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_59c91671f7b5536589a5d44b373d839d)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::capacity;
    bool  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_00a442cf8a5f53f38bfff343e8f23f7b)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::empty;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_0b5530b8e20050d2950f7adc461f5da8)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::max_size;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_a79b8859380452a392bdc2adde215aed)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reserve;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_65d73148dc4a508b9db1cda860b09e83)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::at;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_bf336e01fc2458ccafb095c757d43b9f)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::size_type ) const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::at;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_85254a212bfb5de2a8d1e892bebff8e4)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::front;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_6071fe18e9b253589ef7c98a5c88134b)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::front;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_9d49cf5d7771526a92d0faa0353ddae6)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::back;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_b5e6f2029fa05132a3d3bfccab670a56)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::back;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::value_type * (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_c426cdc96532577990ba57349cd35833)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::data;
    ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::value_type const * (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_617fba872ea3521a852a8ee834d0335c)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::data;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_0b22f529f5435db6af82ef86fe3e03a0)(::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::const_reference ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::push_back;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_97d661725bf454f68d5a0ae28aca235d)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::pop_back;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_3f53e626eb205b5393375291315b1cb1)() = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::clear;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_1ec767cb94f8580f83d9f9cd8499d849)(class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > &) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::swap;
    bool  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_ba51b1ea970f5b49a87e7b6ab5fe8c58)() const = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::__invariants;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_9c26dc55b84a5ccdb4993df36c9aa4c4)(class ::stat_tool::Vectors *, class ::stat_tool::Vectors *) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::assign;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_2f262c3143735ccc9f9d861d3de48ed6)(class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > , class ::std::move_iterator< class ::std::__wrap_iter< class ::stat_tool::Vectors * > > ) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::assign;
    void  (::std::vector< ::stat_tool::Vectors, ::std::allocator< ::stat_tool::Vectors > >::*method_pointer_a7cbd11ddd77559083f8f6b6f1d70e5d)(class ::stat_tool::Vectors const *, class ::stat_tool::Vectors const *) = &::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >::assign;
    boost::python::class_< class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > >, autowig::Held< class ::std::vector< class ::stat_tool::Vectors, class ::std::allocator< class ::stat_tool::Vectors > > >::Type > class_db7f50b235e15b1aa024911715fa604a("_Vector_db7f50b235e15b1aa024911715fa604a", "", boost::python::no_init);
    class_db7f50b235e15b1aa024911715fa604a.def("assign", method_pointer_83baa9fc9c1157d9b24615b6b82790b7, "");
    class_db7f50b235e15b1aa024911715fa604a.def("__len__", method_pointer_8c95d46dc6e35c1a92ed5968642ca1ad, "");
    class_db7f50b235e15b1aa024911715fa604a.def("capacity", method_pointer_59c91671f7b5536589a5d44b373d839d, "");
    class_db7f50b235e15b1aa024911715fa604a.def("empty", method_pointer_00a442cf8a5f53f38bfff343e8f23f7b, "");
    class_db7f50b235e15b1aa024911715fa604a.def("max_size", method_pointer_0b5530b8e20050d2950f7adc461f5da8, "");
    class_db7f50b235e15b1aa024911715fa604a.def("reserve", method_pointer_a79b8859380452a392bdc2adde215aed, "");
    class_db7f50b235e15b1aa024911715fa604a.def("at", method_pointer_65d73148dc4a508b9db1cda860b09e83, boost::python::return_internal_reference<>(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("at", autowig::method_decorator_65d73148dc4a508b9db1cda860b09e83);
    class_db7f50b235e15b1aa024911715fa604a.def("at", method_pointer_bf336e01fc2458ccafb095c757d43b9f, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("front", method_pointer_85254a212bfb5de2a8d1e892bebff8e4, boost::python::return_internal_reference<>(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("front", autowig::method_decorator_85254a212bfb5de2a8d1e892bebff8e4);
    class_db7f50b235e15b1aa024911715fa604a.def("front", method_pointer_6071fe18e9b253589ef7c98a5c88134b, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("back", method_pointer_9d49cf5d7771526a92d0faa0353ddae6, boost::python::return_internal_reference<>(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("back", autowig::method_decorator_9d49cf5d7771526a92d0faa0353ddae6);
    class_db7f50b235e15b1aa024911715fa604a.def("back", method_pointer_b5e6f2029fa05132a3d3bfccab670a56, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("data", method_pointer_c426cdc96532577990ba57349cd35833, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("data", method_pointer_617fba872ea3521a852a8ee834d0335c, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_db7f50b235e15b1aa024911715fa604a.def("push_back", method_pointer_0b22f529f5435db6af82ef86fe3e03a0, "");
    class_db7f50b235e15b1aa024911715fa604a.def("pop_back", method_pointer_97d661725bf454f68d5a0ae28aca235d, "");
    class_db7f50b235e15b1aa024911715fa604a.def("clear", method_pointer_3f53e626eb205b5393375291315b1cb1, "");
    class_db7f50b235e15b1aa024911715fa604a.def("swap", method_pointer_1ec767cb94f8580f83d9f9cd8499d849, "");
    class_db7f50b235e15b1aa024911715fa604a.def("invariants", method_pointer_ba51b1ea970f5b49a87e7b6ab5fe8c58, "");
    class_db7f50b235e15b1aa024911715fa604a.def("assign", method_pointer_9c26dc55b84a5ccdb4993df36c9aa4c4, "");
    class_db7f50b235e15b1aa024911715fa604a.def("assign", method_pointer_2f262c3143735ccc9f9d861d3de48ed6, "");
    class_db7f50b235e15b1aa024911715fa604a.def("assign", method_pointer_a7cbd11ddd77559083f8f6b6f1d70e5d, "");

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