# -*- coding: utf-8 -*-
"""Abstract base class used by test_mixture, test_compound, etc

author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__version__ = "$Id: tools.py 11452 2011-11-25 15:49:11Z jbdurand $"


from openalea.stat_tool import Simulate
from openalea.stat_tool.plot import DISABLE_PLOT
import os
from openalea.stat_tool.output import Display, Save


def runTestClass(myclass):
    functions = [x for x in dir(myclass) if x.startswith('test')]
    for function in functions:
        getattr(myclass, function)()


def _remove_file(filename):
    """alias to remove a file"""
    try:
        os.remove(filename)
    except:
        pass


class interface():
    """Interface to be used by test files related to data structure such as
    compound, convolution, mixture, histogram, vectors.

    :param data: a data that will be filled using the build_data structure
    :param filename: a filename to a file containing the relevant data structure
    :param structure: reference to a data structure Class that is not instantiated.

    :Usage:
    In you test file, add::

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

    """
    def __init__(self, data=None, filename=None, Structure=None):

        if data is None:
            raise AttributeError("data must be provided")
        if Structure is None:
            raise AttributeError("Structure  must be provided")


        self.data = data
        self.filename = filename
        self.structure = Structure
        self.N = 1000

    def build_data(self):
        raise NotImplementedError()

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
        except Exception:
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
        """In the Vector case, Format should be set to Data.
        :param skip_reading: some class do not have Filename Constructor;
            skip_reading can be set to False to prevent code to be run.

        .. todo:: This is surely a bug. to be checked"""


        c1 = self.data
        #_remove_file('test1.dat')
        #_remove_file('test2.dat')

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

            assert c1 and c1_read and c2_read
            assert str(c1_read) == str(c2_read)

        _remove_file('test1.dat')
        _remove_file('test2.dat')

    def plot_write(self):
        h = self.data
        h.plot_write('test', 'title')
        _remove_file('test.print')
        _remove_file('test.plot')
        _remove_file('test0.dat')

    def file_ascii_write(self):
        h = self.data
        h.file_ascii_write('test.dat', True)
        _remove_file('test.dat')

    def file_ascii_data_write(self):
        h = self.data
        h.file_ascii_data_write('test.dat', True)
        _remove_file('test.dat')

    def spreadsheet_write(self):
        h = self.data
        h.spreadsheet_write('test.dat')
        _remove_file('test.dat')

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
        _remove_file('test.xsl')

    def simulate(self, N=-1):
        """Test the simulate method"""
        if N == -1:
            N = self.N
        m = self.data
        s = m.simulate(N)
        s2 = Simulate(m, N)
        assert len(s) == N
        assert len(s2) == N
        assert str(s2)
        assert str(s)
        return s

    def test_len(self):
        raise NotImplementedError()

    def test_extract(self):
        raise NotImplementedError()

    def test_extract_data(self):
        raise NotImplementedError()

