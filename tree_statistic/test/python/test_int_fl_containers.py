# a test for the class int_fl_containers: estimation - syntax
import sys, os
AMLLIBDIR=os.getenv("HOME")+"/devlp/AMAPmod/build-linux/lib"
sys.path.append(AMLLIBDIR)
map(lambda x:sys.path.append(os.getenv("AMAPDIR")+"/STAT_TREES/Python/"+x),
    ["Int_fl_containers/"])
sys.path.append(os.getenv("AMAPDIR")+"/STAT/Python/")
import int_fl_containers
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