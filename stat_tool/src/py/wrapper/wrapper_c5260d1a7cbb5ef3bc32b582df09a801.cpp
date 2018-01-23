#include "_stat_tool.h"



namespace autowig
{

    void method_decorator_b2dc7c9864fc5229847fac6a174f3f0d(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > & instance, ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  param_out) { instance.operator[](param_in_0) = param_out; }
    void method_decorator_d6099c3104af5a46bbaddcfe1bee30f1(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > & instance, ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  param_in_0, const ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  param_out) { instance.at(param_in_0) = param_out; }
    void method_decorator_4c47ded6a4d75c618d527e5817501594(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > & instance, const ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  param_out) { instance.front() = param_out; }
    void method_decorator_45c9d8576617591193097d21f86d935c(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > & instance, const ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  param_out) { instance.back() = param_out; }
}

#if defined(_MSC_VER)
    #if (_MSC_VER == 1900)
namespace boost
{
    template <> class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const volatile * get_pointer<class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const volatile >(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const volatile *c) { return c; }
}
    #endif
#endif



void wrapper_c5260d1a7cbb5ef3bc32b582df09a801()
{

    std::string name_a5e4e9231d6351ccb0e06756b389f0af = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".std");
    boost::python::object module_a5e4e9231d6351ccb0e06756b389f0af(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(name_a5e4e9231d6351ccb0e06756b389f0af.c_str()))));
    boost::python::scope().attr("std") = module_a5e4e9231d6351ccb0e06756b389f0af;
    boost::python::scope scope_a5e4e9231d6351ccb0e06756b389f0af = module_a5e4e9231d6351ccb0e06756b389f0af;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_5774436e01fe5afe947d5259394933f3)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type , ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::value_type const &) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::assign;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_00068b27e16956399acdaecd87bdfeb2)(class ::std::initializer_list< class ::stat_tool::DiscreteParametric > ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::assign;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_f64ccbdd7dca5e55968d668ad56b50b4)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_bb905999b8e052d3a9fcfc89aa73ebe9)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::max_size;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_afce787ebde35198ace57559280f3ccd)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::resize;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_eba2c57e4afb5f3f8d4328bc97545974)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type , ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::value_type const &) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::resize;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_3f93bd0f08f45c6599b15767916806e5)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::shrink_to_fit;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_2e044c9a2cbf5c1884e6d3e1fcaa15a8)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::capacity;
    bool  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_c43e6c1cbaae57cab981ba1a5cd94ae0)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::empty;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_30562a61977953f19efb9dd6f738bfd6)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reserve;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_b2dc7c9864fc5229847fac6a174f3f0d)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::operator[];
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_826e60df90a558088fb0873fb85d93df)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::operator[];
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_d6099c3104af5a46bbaddcfe1bee30f1)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::at;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_d9251c6ecb615be9bf74204f77e1363c)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type ) const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::at;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_4c47ded6a4d75c618d527e5817501594)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::front;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_653eeb84630c5e9589e00743b13f1890)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::front;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_45c9d8576617591193097d21f86d935c)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::back;
    ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::const_reference  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_1fa248dc911d5f3893d587f5d33a9438)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::back;
    class ::stat_tool::DiscreteParametric * (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_f1f8d572818c56eb86cd1b9011d7e6cd)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::data;
    class ::stat_tool::DiscreteParametric const * (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_302ff38f4cdb5748bf66b586649c6f16)() const = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::data;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_9133b69f693a50e5940f2a06dc76055c)(::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::value_type const &) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::push_back;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_193f9807b02d5391aaedec88a7a7973c)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::pop_back;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_d5bf4a0615a55781aa02a1c8feb03ded)(class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > &) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::swap;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_18c7b7b548a65b37ab7accafad2b0228)() = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::clear;
    void  (::std::vector< ::stat_tool::DiscreteParametric, ::std::allocator< ::stat_tool::DiscreteParametric > >::*method_pointer_5af19a1366e8542fb85129c0b4da73ba)(class ::stat_tool::DiscreteParametric const *, class ::stat_tool::DiscreteParametric const *) = &::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::assign;
    boost::python::class_< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >, autowig::Held< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > >::Type > class_c5260d1a7cbb5ef3bc32b582df09a801("_Vector_c5260d1a7cbb5ef3bc32b582df09a801", "", boost::python::no_init);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def(boost::python::init<  >(""));
    class_c5260d1a7cbb5ef3bc32b582df09a801.def(boost::python::init< ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::allocator_type const & >(""));
    class_c5260d1a7cbb5ef3bc32b582df09a801.def(boost::python::init< ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type , ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::allocator_type const & >(""));
    class_c5260d1a7cbb5ef3bc32b582df09a801.def(boost::python::init< ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::size_type , ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::value_type const &, ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::allocator_type const & >(""));
    class_c5260d1a7cbb5ef3bc32b582df09a801.def(boost::python::init< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const & >(""));
    class_c5260d1a7cbb5ef3bc32b582df09a801.def(boost::python::init< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > const &, ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::allocator_type const & >(""));
    class_c5260d1a7cbb5ef3bc32b582df09a801.def(boost::python::init< class ::std::initializer_list< class ::stat_tool::DiscreteParametric > , ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >::allocator_type const & >(""));
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("assign", method_pointer_5774436e01fe5afe947d5259394933f3, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("assign", method_pointer_00068b27e16956399acdaecd87bdfeb2, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("__len__", method_pointer_f64ccbdd7dca5e55968d668ad56b50b4, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("max_size", method_pointer_bb905999b8e052d3a9fcfc89aa73ebe9, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("resize", method_pointer_afce787ebde35198ace57559280f3ccd, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("resize", method_pointer_eba2c57e4afb5f3f8d4328bc97545974, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("shrink_to_fit", method_pointer_3f93bd0f08f45c6599b15767916806e5, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("capacity", method_pointer_2e044c9a2cbf5c1884e6d3e1fcaa15a8, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("empty", method_pointer_c43e6c1cbaae57cab981ba1a5cd94ae0, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("reserve", method_pointer_30562a61977953f19efb9dd6f738bfd6, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("__getitem__", method_pointer_b2dc7c9864fc5229847fac6a174f3f0d, boost::python::return_internal_reference<>(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("__getitem__", autowig::method_decorator_b2dc7c9864fc5229847fac6a174f3f0d);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("__getitem__", method_pointer_826e60df90a558088fb0873fb85d93df, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("at", method_pointer_d6099c3104af5a46bbaddcfe1bee30f1, boost::python::return_internal_reference<>(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("at", autowig::method_decorator_d6099c3104af5a46bbaddcfe1bee30f1);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("at", method_pointer_d9251c6ecb615be9bf74204f77e1363c, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("front", method_pointer_4c47ded6a4d75c618d527e5817501594, boost::python::return_internal_reference<>(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("front", autowig::method_decorator_4c47ded6a4d75c618d527e5817501594);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("front", method_pointer_653eeb84630c5e9589e00743b13f1890, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("back", method_pointer_45c9d8576617591193097d21f86d935c, boost::python::return_internal_reference<>(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("back", autowig::method_decorator_45c9d8576617591193097d21f86d935c);
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("back", method_pointer_1fa248dc911d5f3893d587f5d33a9438, boost::python::return_value_policy< boost::python::return_by_value >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("data", method_pointer_f1f8d572818c56eb86cd1b9011d7e6cd, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("data", method_pointer_302ff38f4cdb5748bf66b586649c6f16, boost::python::return_value_policy< boost::python::reference_existing_object >(), "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("push_back", method_pointer_9133b69f693a50e5940f2a06dc76055c, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("pop_back", method_pointer_193f9807b02d5391aaedec88a7a7973c, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("swap", method_pointer_d5bf4a0615a55781aa02a1c8feb03ded, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("clear", method_pointer_18c7b7b548a65b37ab7accafad2b0228, "");
    class_c5260d1a7cbb5ef3bc32b582df09a801.def("assign", method_pointer_5af19a1366e8542fb85129c0b4da73ba, "");

    struct vector_c5260d1a7cbb5ef3bc32b582df09a801_from_python
    {
        vector_c5260d1a7cbb5ef3bc32b582df09a801_from_python()
        {
            boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > >());
        }

        static void* convertible(PyObject* obj_ptr)
        { return obj_ptr; }

        static void construct(PyObject* obj_ptr, boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
            void* storage = ((boost::python::converter::rvalue_from_python_storage< class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > > >*)data)->storage.bytes;
            new (storage) class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >();
            data->convertible = storage;
            class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >& result = *((class ::std::vector< class ::stat_tool::DiscreteParametric, class ::std::allocator< class ::stat_tool::DiscreteParametric > >*)storage);
            unsigned int i = 0;
            for(;; i++)
            {
                boost::python::handle<> py_elem_hdl(boost::python::allow_null(PyIter_Next(obj_iter.get())));
                if(PyErr_Occurred())
                { boost::python::throw_error_already_set(); }
                if(!py_elem_hdl.get())
                { break; }
                boost::python::object py_elem_obj(py_elem_hdl);
                result.push_back(boost::python::extract< class ::stat_tool::DiscreteParametric  >(py_elem_obj));
            }
        }
    };

    vector_c5260d1a7cbb5ef3bc32b582df09a801_from_python();
}