"""Abstract base class used by test_mixture, test_compound, etc"""

from openalea.stat_tool import *
from openalea.stat_tool.plot import DISABLE_PLOT
import os
from openalea.stat_tool.output import Display, Save


__revision__ = "$Id: tools.py 9000 2010-05-26 09:21:23Z cokelaer $"

def runTestClass(myclass):
    functions = [x for x in dir(myclass) if x.startswith('test')]
    for function in functions:
        getattr(myclass, function)()

class interface():
    """Interface to be used by test file that perform tests on the following
    data structure: compound, convolution, mixture, histogram, vector

    :param data: a data that will be filled using the build_data structure
    :para filename: a filename to a file containing the relevant data structure
    :param structure: reference to a data structure Class that is not instantiated.

    :Usage:
    In you test file, add ::

        >>> from tools import interface

    Then, if we consider the Compound class case, create a class as follows::

        class Test(interface):
            def __init__(self):
                self.data = self.build_data()
                self.filename = "compound1.cd"
                self.structure = Compound

            def build_data(self):
                d1 = Binomial(2, 5, 0.5)
                d2 = NegativeBinomial(0, 2, 0.5)
                data = Compound(d1, d2)
                return data

            def test_empty(self):
                self.empty()

    """
    def __init__(self, data=None, filename=None, Structure=None):
        self.data = data
        self.filename = filename
        self.structure = Structure

    def build_data(self):
        raise NotImplementedError()

    def empty(self):
        """Test that empty constructor fails"""
        try:
            _m = self.structure()
            assert False
        except TypeError:
            assert True

    def constructor_from_file(self):
        """Test constructor from file"""
        if self.filename == None:
            return None
        else:
            c = self.structure(self.filename)
            assert c
            return c

    def constructor_from_file_failure(self):
        """run constructor with filename argument"""
        try:
            _h = self.structure("whatever_wrong_filename.txt")
            assert False
        except IOError:
            assert True

    def print_data(self):
        """test that print command exists"""
        print self.data

    def display(self):
        """check that .display/Display is callable"""
        data = self.data
        data.display()
        Display(data)
        assert data.display()==Display(data)

    def display_versus_ascii_write(self):
        """check that display is equivalent to ascii_write"""
        assert Display(self.data) == self.data.ascii_write(False)

    def display_versus_str(self):
        """check that display and str are equivalent"""
        data = self.data
        s = str(data)
        assert Display(data) == s

    def plot(self):
        """run plotting routines """
        if DISABLE_PLOT == False:
            self.data.plot()

    def save(self, Format=None, skip_reading=False):
        """In the Vector case, Format should be Data.
        :param skip_reading: some class do not have Filename Constructor;
            skip_reading can be set to False to prevent code to be run.

        .. todo:: This is surely a bug. to be checked"""


        c1 = self.data

        try:
            os.remove('test1.dat')
        except:
            pass
        try:
            os.remove('test2.dat')
        except:
            pass

        if Format is None:
            c1.save('test1.dat')
            Save(c1, 'test2.dat')
        else:
            c1.save('test1.dat', Format="Data")
            Save(c1, 'test2.dat', Format="Data")

        if skip_reading:
            pass
        else:
            c1_read = self.structure('test1.dat')
            c2_read = self.structure('test2.dat')

            print c1_read

            assert c1 and c1_read and c2_read
            assert str(c1_read) == str(c2_read)

        #os.remove('test1.dat')
        #os.remove('test2.dat')

    def plot_write(self):
        h = self.data
        h.plot_write('test', 'title')

    def file_ascii_write(self):
        h = self.data
        h.file_ascii_write('test.dat', True)
        os.remove('test.dat')

    def file_ascii_data_write(self):
        h = self.data
        h.file_ascii_data_write('test.dat', True)
        os.remove('test.dat')

    def spreadsheet_write(self):
        h = self.data
        h.spreadsheet_write('test.dat')
        os.remove('test.dat')

    def survival_ascii_write(self):
        d = self.data
        d.survival_ascii_write()

    def survival_plot_write(self):
        d = self.data
        d.survival_plot_write('test','test')

    def survival_file_ascii_write(self):
        d = self.data
        d.survival_file_ascii_write()

    def survival_spreadsheet_write(self):
        d = self.data
        d.survival_spreadsheet_write('test.xsl')
        os.remove('test.xsl')

    def simulate(self):
        """Test the simulate method"""
        m = self.data
        s = m.simulate(1000)
        s2 = Simulate(m, 1000)
        assert len(s) == 1000
        assert len(s2) == 1000
        assert str(s2)
        assert str(s)
        return s

    def test_len(self):
        raise NotImplementedError()

    def test_extract(self):
        raise NotImplementedError()

    def test_extract_data(self):
        raise NotImplementedError()


