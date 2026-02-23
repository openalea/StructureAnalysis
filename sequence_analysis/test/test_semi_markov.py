"""Semi markov tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""

__revision__ = "$Id$"

import pytest

from openalea.sequence_analysis import SemiMarkov, Simulate


from .tools import interface
from .tools import robust_path as get_shared_data


@pytest.fixture
def interface_instance():
    filename = str(get_shared_data("test_semi_markov.dat"))
    return interface(
        data=SemiMarkov(filename),
        filename=filename,
        Structure=SemiMarkov,
    )


class TestSemiMarkov:
    """a simple unittest class"""

    def test_constructor_from_file(self, interface_instance):
        interface_instance.constructor_from_file()

    def build_data(self):
        """todo: check identifier output. should be a list"""
        # build a list of 2 sequences with a variable that should be identical
        # to sequences1.seq
        return SemiMarkov(str(get_shared_data("test_semi_markov.dat")))

    def test_simulate(self, interface_instance):
        Simulate(interface_instance.data, 1, 1000, True)
        pass

    def test_empty(self, interface_instance):
        interface_instance.empty()

    def test_constructor_from_file_failure(self, interface_instance):
        interface_instance.constructor_from_file_failure()

    def test_print(self, interface_instance):
        interface_instance.print_data()

    def test_display(self, interface_instance):
        interface_instance.display()
        interface_instance.display_versus_ascii_write()
        interface_instance.display_versus_str()

    def test_len(self, interface_instance):
        seq = interface_instance.data
        pass

    def test_plot(self, interface_instance):
        interface_instance.plot()

    def test_save(self, interface_instance):
        interface_instance.save(skip_reading=True)

    def test_plot_write(self, interface_instance):
        interface_instance.plot_write()

    def test_file_ascii_write(self, interface_instance):
        interface_instance.file_ascii_write()

    def test_spreadsheet_write(self, interface_instance):
        interface_instance.spreadsheet_write()

    def test_thresholding(self, interface_instance):
        interface_instance.data.thresholding(1)

    def test_extract(self, interface_instance):
        pass
        # interface_instance.data.extract(0,1,1)

    def test_extract_data(self, interface_instance):
        interface_instance.data.extract_data()
