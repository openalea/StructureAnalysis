#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_e8049c6a5d7b50d2b9d74d7dcc0a38eb(class ::std::vector< int, class ::std::allocator< int > > & instance, ::std::vector< int, class ::std::allocator< int > >::size_type  param_in_0, int param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_56536485a2435a55b80399050e485a93(class ::std::vector< int, class ::std::allocator< int > > & instance, int param_out) { instance.front() = param_out; }
    void method_decorator_d98ee2218c6e582799423a1894af6aaf(class ::std::vector< int, class ::std::allocator< int > > & instance, int param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< int, class ::std::allocator< int > > const volatile * get_pointer<class ::std::vector< int, class ::std::allocator< int > > const volatile >(class ::std::vector< int, class ::std::allocator< int > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_6b9ae5eac40858c9a0f5e6e21c15d1d3()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_fec8568f38ea57d29eda0d0e7ef7bb61)(::std::vector< int, class ::std::allocator< int > >::size_type , ::std::vector< int, class ::std::allocator< int > >::const_reference ) = &::std::vector< int, class ::std::allocator< int > >::assign;
    ::std::vector< int, class ::std::allocator< int > >::size_type  (::std::vector< int, ::std::allocator< int > >::*method_pointer_77f67c3020965bcdbbe36f6d3c79297d)() const = &::std::vector< int, class ::std::allocator< int > >::size;
    ::std::vector< int, class ::std::allocator< int > >::size_type  (::std::vector< int, ::std::allocator< int > >::*method_pointer_f0284fcb229051f984da13c1d83c94ff)() const = &::std::vector< int, class ::std::allocator< int > >::capacity;
    bool  (::std::vector< int, ::std::allocator< int > >::*method_pointer_41ea14132fb654cb855fbf3af42e7228)() const = &::std::vector< int, class ::std::allocator< int > >::empty;
    ::std::vector< int, class ::std::allocator< int > >::size_type  (::std::vector< int, ::std::allocator< int > >::*method_pointer_a7ab2e475ccc578db3d543eaa3bc4544)() const = &::std::vector< int, class ::std::allocator< int > >::max_size;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_8a73c7d270b05cbbb64fba92e80dd145)(::std::vector< int, class ::std::allocator< int > >::size_type ) = &::std::vector< int, class ::std::allocator< int > >::reserve;
    ::std::vector< int, class ::std::allocator< int > >::reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_e8049c6a5d7b50d2b9d74d7dcc0a38eb)(::std::vector< int, class ::std::allocator< int > >::size_type ) = &::std::vector< int, class ::std::allocator< int > >::at;
    ::std::vector< int, class ::std::allocator< int > >::const_reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_a5827e4f1131537eb23a38bf28b92da6)(::std::vector< int, class ::std::allocator< int > >::size_type ) const = &::std::vector< int, class ::std::allocator< int > >::at;
    ::std::vector< int, class ::std::allocator< int > >::reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_56536485a2435a55b80399050e485a93)() = &::std::vector< int, class ::std::allocator< int > >::front;
    ::std::vector< int, class ::std::allocator< int > >::const_reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_e6ec0047097a5395a29ff08fc1332603)() const = &::std::vector< int, class ::std::allocator< int > >::front;
    ::std::vector< int, class ::std::allocator< int > >::reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_d98ee2218c6e582799423a1894af6aaf)() = &::std::vector< int, class ::std::allocator< int > >::back;
    ::std::vector< int, class ::std::allocator< int > >::const_reference  (::std::vector< int, ::std::allocator< int > >::*method_pointer_bd02b15c21da5e05bc2d5d2273777ed8)() const = &::std::vector< int, class ::std::allocator< int > >::back;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_059a7a64cad45bb883f68794a46f0a1a)(::std::vector< int, class ::std::allocator< int > >::const_reference ) = &::std::vector< int, class ::std::allocator< int > >::push_back;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_b4a17f84c3ae5d36876dc92278aabd38)() = &::std::vector< int, class ::std::allocator< int > >::pop_back;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_be385133d5415036ba8aa12cd2b99c76)() = &::std::vector< int, class ::std::allocator< int > >::clear;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_85d8e257ec345360aaa8fea08c1324eb)(class ::std::vector< int, class ::std::allocator< int > > &) = &::std::vector< int, class ::std::allocator< int > >::swap;
    bool  (::std::vector< int, ::std::allocator< int > >::*method_pointer_63d4eabca8b55be8a1aadb729c70772b)() const = &::std::vector< int, class ::std::allocator< int > >::__invariants;
    void  (::std::vector< int, ::std::allocator< int > >::*method_pointer_91864741fe0a5094b7f9e99332558507)(class ::std::move_iterator< class ::std::__wrap_iter< int * > > , class ::std::move_iterator< class ::std::__wrap_iter< int * > > ) = &::std::vector< int, class ::std::allocator< int > >::assign;
    boost::python::class_< class ::std::vector< int, class ::std::allocator< int > >, autowig::Held< class ::std::vector< int, class ::std::allocator< int > > >::Type > class_6b9ae5eac40858c9a0f5e6e21c15d1d3("_Vector_6b9ae5eac40858c9a0f5e6e21c15d1d3", "", boost::python::no_init);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("assign", method_pointer_fec8568f38ea57d29eda0d0e7ef7bb61, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("__len__", method_pointer_77f67c3020965bcdbbe36f6d3c79297d, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("capacity", method_pointer_f0284fcb229051f984da13c1d83c94ff, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("empty", method_pointer_41ea14132fb654cb855fbf3af42e7228, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("max_size", method_pointer_a7ab2e475ccc578db3d543eaa3bc4544, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("reserve", method_pointer_8a73c7d270b05cbbb64fba92e80dd145, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("at", method_pointer_e8049c6a5d7b50d2b9d74d7dcc0a38eb, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("at", autowig::method_decorator_e8049c6a5d7b50d2b9d74d7dcc0a38eb);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("at", method_pointer_a5827e4f1131537eb23a38bf28b92da6, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("front", method_pointer_56536485a2435a55b80399050e485a93, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("front", autowig::method_decorator_56536485a2435a55b80399050e485a93);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("front", method_pointer_e6ec0047097a5395a29ff08fc1332603, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("back", method_pointer_d98ee2218c6e582799423a1894af6aaf, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("back", autowig::method_decorator_d98ee2218c6e582799423a1894af6aaf);
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("back", method_pointer_bd02b15c21da5e05bc2d5d2273777ed8, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("push_back", method_pointer_059a7a64cad45bb883f68794a46f0a1a, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("pop_back", method_pointer_b4a17f84c3ae5d36876dc92278aabd38, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("clear", method_pointer_be385133d5415036ba8aa12cd2b99c76, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("swap", method_pointer_85d8e257ec345360aaa8fea08c1324eb, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("invariants", method_pointer_63d4eabca8b55be8a1aadb729c70772b, "");
    class_6b9ae5eac40858c9a0f5e6e21c15d1d3.def("assign", method_pointer_91864741fe0a5094b7f9e99332558507, "");

    struct vector_6b9ae5eac40858c9a0f5e6e21c15d1d3_from_python
    {
        vector_6b9ae5eac40858c9a0f5e6e21c15d1d3_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< int, class ::std::allocator< int > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< int, class ::std::allocator< int > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< int, class ::std::allocator< int > >();
            data->convertible = storage;
            class ::std::vector< int, class ::std::allocator< int > >& result = *((class ::std::vector< int, class ::std::allocator< int > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< int  >(py_elem_obj));
            }
        }
    };

    vector_6b9ae5eac40858c9a0f5e6e21c15d1d3_from_python();
}