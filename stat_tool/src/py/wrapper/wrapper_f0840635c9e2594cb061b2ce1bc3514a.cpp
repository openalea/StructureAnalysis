#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_1470d2a590d75f1690eef31561d13f47(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > & instance, ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  param_out) { instance.operator[](param_in_0) = param_out; }
    void method_decorator_f4d77e53ae055daf8be34daf5c4f9b3e(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > & instance, ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_340c3b2eff2a5503a6698a72f0dbaa4d(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > & instance, const ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_1e1610524a6053239b15b5f3478e55db(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > & instance, const ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > const volatile * get_pointer<class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > const volatile >(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_f0840635c9e2594cb061b2ce1bc3514a()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_ed4258ed4d695402b3e64bc4a38ad9f6)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type , ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::value_type const &) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::assign;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_f743285ed70f5cee9a26ea28455ed82e)(class ::std::initializer_list< class ::stat_tool::FrequencyDistribution > ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::assign;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_11a81fbbadfd5463bf60c597c5a11e8f)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_b012012fe0415bbaaa8be0cda73464d3)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::max_size;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_2d9664c3d29a5445871432890353093d)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::resize;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_6207efc0b9645e4591f8efda397fd39a)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type , ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::value_type const &) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::resize;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_651df6d16b9a5144a998c6dd55e2c65b)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::shrink_to_fit;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_4c43f568dab15a63bce6cc29d444151f)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::capacity;
    bool  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_0dfdad518c7d5efcae8c5dcc8731708a)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::empty;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_acfebb82196052f38d8eb36f40e28cea)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reserve;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_1470d2a590d75f1690eef31561d13f47)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::operator[];
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_92083c64a1ae5f6a951ec2e233f522df)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::operator[];
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_f4d77e53ae055daf8be34daf5c4f9b3e)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::at;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_b5af4810e3ba58cba76d4bb466f4dc1c)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type ) const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::at;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_340c3b2eff2a5503a6698a72f0dbaa4d)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::front;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_bd575619154f5b7a83d69551a0972f8f)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::front;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_1e1610524a6053239b15b5f3478e55db)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::back;
    ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::const_reference  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_7d2ce74aa4535e1898aa54307638b950)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::back;
    class ::stat_tool::FrequencyDistribution * (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_6f323f454d5b5cc78e96fb3b3206c5ae)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::data;
    class ::stat_tool::FrequencyDistribution const * (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_87673768c62f5c4a96987fa8a8f33c73)() const = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::data;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_8dc8ae7e87705ee88daf63470ca2c95f)(::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::value_type const &) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::push_back;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_2723d21b3ef656e5b39d56db63779b7d)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::pop_back;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_8c3ab5be63135191b7cf3beb58c9fdf4)(class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > &) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::swap;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_07472d7e4ed95d44b89757ab17559373)() = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::clear;
    void  (::std::vector< ::stat_tool::FrequencyDistribution, ::std::allocator< ::stat_tool::FrequencyDistribution > >::*method_pointer_b278ac5e51915416b63e2bac5f342106)(class ::stat_tool::FrequencyDistribution const *, class ::stat_tool::FrequencyDistribution const *) = &::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::assign;
    boost::python::class_< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >, autowig::Held< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > >::Type > class_f0840635c9e2594cb061b2ce1bc3514a("_Vector_f0840635c9e2594cb061b2ce1bc3514a", "", boost::python::no_init);
    class_f0840635c9e2594cb061b2ce1bc3514a.def(boost::python::init<  >(""));
    class_f0840635c9e2594cb061b2ce1bc3514a.def(boost::python::init< ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::allocator_type const & >(""));
    class_f0840635c9e2594cb061b2ce1bc3514a.def(boost::python::init< ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type , ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::allocator_type const & >(""));
    class_f0840635c9e2594cb061b2ce1bc3514a.def(boost::python::init< ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::size_type , ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::value_type const &, ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::allocator_type const & >(""));
    class_f0840635c9e2594cb061b2ce1bc3514a.def(boost::python::init< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > const & >(""));
    class_f0840635c9e2594cb061b2ce1bc3514a.def(boost::python::init< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > const &, ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::allocator_type const & >(""));
    class_f0840635c9e2594cb061b2ce1bc3514a.def(boost::python::init< class ::std::initializer_list< class ::stat_tool::FrequencyDistribution > , ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >::allocator_type const & >(""));
    class_f0840635c9e2594cb061b2ce1bc3514a.def("assign", method_pointer_ed4258ed4d695402b3e64bc4a38ad9f6, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("assign", method_pointer_f743285ed70f5cee9a26ea28455ed82e, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("__len__", method_pointer_11a81fbbadfd5463bf60c597c5a11e8f, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("max_size", method_pointer_b012012fe0415bbaaa8be0cda73464d3, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("resize", method_pointer_2d9664c3d29a5445871432890353093d, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("resize", method_pointer_6207efc0b9645e4591f8efda397fd39a, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("shrink_to_fit", method_pointer_651df6d16b9a5144a998c6dd55e2c65b, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("capacity", method_pointer_4c43f568dab15a63bce6cc29d444151f, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("empty", method_pointer_0dfdad518c7d5efcae8c5dcc8731708a, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("reserve", method_pointer_acfebb82196052f38d8eb36f40e28cea, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("__getitem__", method_pointer_1470d2a590d75f1690eef31561d13f47, boost::python::return_internal_reference<>(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("__getitem__", autowig::method_decorator_1470d2a590d75f1690eef31561d13f47);
    class_f0840635c9e2594cb061b2ce1bc3514a.def("__getitem__", method_pointer_92083c64a1ae5f6a951ec2e233f522df, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("at", method_pointer_f4d77e53ae055daf8be34daf5c4f9b3e, boost::python::return_internal_reference<>(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("at", autowig::method_decorator_f4d77e53ae055daf8be34daf5c4f9b3e);
    class_f0840635c9e2594cb061b2ce1bc3514a.def("at", method_pointer_b5af4810e3ba58cba76d4bb466f4dc1c, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("front", method_pointer_340c3b2eff2a5503a6698a72f0dbaa4d, boost::python::return_internal_reference<>(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("front", autowig::method_decorator_340c3b2eff2a5503a6698a72f0dbaa4d);
    class_f0840635c9e2594cb061b2ce1bc3514a.def("front", method_pointer_bd575619154f5b7a83d69551a0972f8f, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("back", method_pointer_1e1610524a6053239b15b5f3478e55db, boost::python::return_internal_reference<>(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("back", autowig::method_decorator_1e1610524a6053239b15b5f3478e55db);
    class_f0840635c9e2594cb061b2ce1bc3514a.def("back", method_pointer_7d2ce74aa4535e1898aa54307638b950, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("data", method_pointer_6f323f454d5b5cc78e96fb3b3206c5ae, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("data", method_pointer_87673768c62f5c4a96987fa8a8f33c73, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("push_back", method_pointer_8dc8ae7e87705ee88daf63470ca2c95f, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("pop_back", method_pointer_2723d21b3ef656e5b39d56db63779b7d, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("swap", method_pointer_8c3ab5be63135191b7cf3beb58c9fdf4, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("clear", method_pointer_07472d7e4ed95d44b89757ab17559373, "");
    class_f0840635c9e2594cb061b2ce1bc3514a.def("assign", method_pointer_b278ac5e51915416b63e2bac5f342106, "");

    struct vector_f0840635c9e2594cb061b2ce1bc3514a_from_python
    {
        vector_f0840635c9e2594cb061b2ce1bc3514a_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >();
            data->convertible = storage;
            class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >& result = *((class ::std::vector< class ::stat_tool::FrequencyDistribution, class ::std::allocator< class ::stat_tool::FrequencyDistribution > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< class ::stat_tool::FrequencyDistribution  >(py_elem_obj));
            }
        }
    };

    vector_f0840635c9e2594cb061b2ce1bc3514a_from_python();
}