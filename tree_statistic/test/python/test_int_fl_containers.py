# a test for the class int_fl_containers: estimation - syntax
import sys, os
import openalea.tree_statistic.int_fl_containers as int_fl_containers
w=int_fl_containers.Int_fl_container(1,0)
print w.NbInt()
print w.NbFloat()
print w
w=int_fl_containers.Int_fl_container(1,2)
print w.NbInt()
print w.NbFloat()
print w
print w.Int(0)
# print w.Double(0)
# w.Double(0)=1.5;
val=2
print w