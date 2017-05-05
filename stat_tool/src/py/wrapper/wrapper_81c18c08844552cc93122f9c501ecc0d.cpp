#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_f2609e8dda33526cb8e0f6bca07d716a(class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > & instance, ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size_type  param_in_0, ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::reference  param_out)     { instance.at(param_in_0) = param_out; }
    void method_decorator_75b68accdb275f9aa5cefe8388efe272(class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > & instance, ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::reference  param_out)     { instance.front() = param_out; }
    void method_decorator_1325714ee0e355b083afd715c0f1e71e(class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > & instance, ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::reference  param_out)     { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > const volatile * get_pointer<class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > const volatile >(class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_81c18c08844552cc93122f9c501ecc0d()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_31334cf4600e562a858f05cd286b1918)(::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size_type , ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::value_type const &) = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::assign;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size_type  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_e6138d5e4fd85aa9ae13851f2e356e99)() const = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size_type  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_b426f501914b58dba5503a2c78d8c27e)() const = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::max_size;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size_type  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_b36744b435fd53c4a506fa47082786a4)() const = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::capacity;
    bool  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_4cd9c603c30f5064997ced23a140d43d)() const = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::empty;
    void  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_dd72e1edaf885b97a750bf45b76925a3)(::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size_type ) = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::reserve;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::reference  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_f2609e8dda33526cb8e0f6bca07d716a)(::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size_type ) = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::at;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::const_reference  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_f220f529196f522cb784fc409db55427)(::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::size_type ) const = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::at;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::reference  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_75b68accdb275f9aa5cefe8388efe272)() = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::front;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::const_reference  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_838f7ed9e7465a5b9900e94eca4e77b9)() const = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::front;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::reference  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_1325714ee0e355b083afd715c0f1e71e)() = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::back;
    ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::const_reference  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_5d0c908447775f05b574af4368a59ab5)() const = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::back;
    void  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_cb59d2f9dddf56989045959fcfff443b)(::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::value_type const &) = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::push_back;
    void  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_6e0f1e7fedc1536c9b7315248b0d1a9a)() = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::pop_back;
    void  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_cddc307e784756dbb7bc1cebd4104c27)(class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > &) = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::swap;
    void  (::std::vector< enum ::stat_tool::process_distribution, ::std::allocator< enum ::stat_tool::process_distribution > >::*method_pointer_ad1aad8cd05c59b4a9d332d5beaca2d4)() = &::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >::clear;
    boost::python::class_< class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >, autowig::Held< class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > >::Type > class_81c18c08844552cc93122f9c501ecc0d("_Vector_81c18c08844552cc93122f9c501ecc0d", "", boost::python::no_init);
    class_81c18c08844552cc93122f9c501ecc0d.def(boost::python::init<  >(""));
    class_81c18c08844552cc93122f9c501ecc0d.def(boost::python::init< class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > const & >(""));
    class_81c18c08844552cc93122f9c501ecc0d.def("assign", method_pointer_31334cf4600e562a858f05cd286b1918, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("__len__", method_pointer_e6138d5e4fd85aa9ae13851f2e356e99, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("max_size", method_pointer_b426f501914b58dba5503a2c78d8c27e, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("capacity", method_pointer_b36744b435fd53c4a506fa47082786a4, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("empty", method_pointer_4cd9c603c30f5064997ced23a140d43d, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("reserve", method_pointer_dd72e1edaf885b97a750bf45b76925a3, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("at", method_pointer_f2609e8dda33526cb8e0f6bca07d716a, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_81c18c08844552cc93122f9c501ecc0d.def("at", autowig::method_decorator_f2609e8dda33526cb8e0f6bca07d716a);
    class_81c18c08844552cc93122f9c501ecc0d.def("at", method_pointer_f220f529196f522cb784fc409db55427, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_81c18c08844552cc93122f9c501ecc0d.def("front", method_pointer_75b68accdb275f9aa5cefe8388efe272, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_81c18c08844552cc93122f9c501ecc0d.def("front", autowig::method_decorator_75b68accdb275f9aa5cefe8388efe272);
    class_81c18c08844552cc93122f9c501ecc0d.def("front", method_pointer_838f7ed9e7465a5b9900e94eca4e77b9, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_81c18c08844552cc93122f9c501ecc0d.def("back", method_pointer_1325714ee0e355b083afd715c0f1e71e, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_81c18c08844552cc93122f9c501ecc0d.def("back", autowig::method_decorator_1325714ee0e355b083afd715c0f1e71e);
    class_81c18c08844552cc93122f9c501ecc0d.def("back", method_pointer_5d0c908447775f05b574af4368a59ab5, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_81c18c08844552cc93122f9c501ecc0d.def("push_back", method_pointer_cb59d2f9dddf56989045959fcfff443b, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("pop_back", method_pointer_6e0f1e7fedc1536c9b7315248b0d1a9a, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("swap", method_pointer_cddc307e784756dbb7bc1cebd4104c27, "");
    class_81c18c08844552cc93122f9c501ecc0d.def("clear", method_pointer_ad1aad8cd05c59b4a9d332d5beaca2d4, "");

    struct vector_81c18c08844552cc93122f9c501ecc0d_from_python
    {
        vector_81c18c08844552cc93122f9c501ecc0d_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >();
            data->convertible = storage;
            class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >& result = *((class ::std::vector< enum ::stat_tool::process_distribution, class ::std::allocator< enum ::stat_tool::process_distribution > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< enum ::stat_tool::process_distribution  >(py_elem_obj));
            }
        }
    };

    vector_81c18c08844552cc93122f9c501ecc0d_from_python();
}