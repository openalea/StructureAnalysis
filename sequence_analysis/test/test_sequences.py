"""Sequences tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id: test_sequences.py 9885 2010-11-06 18:19:34Z cokelaer $"


from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.sequences import Sequences, Split, SaveMTG
from openalea.sequence_analysis.data_transform import Cumulate, Difference, \
    ExtractVectors, IndexParameterExtract, IndexParameterSelect, RecurrenceTimeSequences

from openalea.stat_tool.data_transform import *
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.cluster import Transcode, Cluster

from tools import interface
from tools import runTestClass

from openalea.sequence_analysis.sequences import Sequences, IndexParameterType
from openalea.sequence_analysis import get_shared_data


class Test(interface):
    """a simple unittest class


    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           get_shared_data("sequences1.seq"),
                           Sequences)
        self.seqn = self.build_seqn()
        self.seq1 = self.build_seq1()

    def build_data(self):
        """todo: check identifier output. should be a list """
        # build a list of 2 sequences with a variable that should be identical
        # to sequences1.seq
        data = Sequences([
                    [1, 0, 0, 0, 1, 1, 2, 0, 2, 2, 2, 1, 1, 0, 1, 0,
                     1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 2, 2, 2, 1],
                    [0, 0, 0, 1, 1, 0, 2, 0, 2, 2 ,2 ,1 ,1 ,1 ,1 ,0 ,1
                     ,0 ,0 ,0 ,0 ,0]])
        assert data

        assert data.nb_sequence == 2
        assert data.nb_variable == 1
        assert data.cumul_length == 52
        assert data.max_length == 30

        assert [0, 1] == data.get_identifiers()

        return data

    def build_seqn(self):
        s = Sequences([[[1,1,1],[12,12,12]],[[2,2,2],[22,23,24]]])
        return s

    def build_seq1(self):
        s = Sequences([[1,1,1],[2,2,2]])
        return s

    def build_seq_wrong_identifiers(self):
        try:
            s = Sequences([[1,1,1],[2,2,2]], Identifiers=[-1,1])
            assert False
        except:
            assert True

    def _test_empty(self):
        self.empty()

    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def test_print(self):
        self.print_data()

    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()

    def test_len(self):
        seq = self.data
        assert len(seq) == 2
        assert len(seq) == seq.nb_sequence

    def test_plot(self):
        self.plot()

    def test_save(self):
        self.save(skip_reading=True)

    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()

    def test_spreadsheet_write(self):
        self.spreadsheet_write()

    def test_extract(self):
        #todo
        seqn = self.seqn
        assert seqn.extract_value(1)
        assert seqn.extract_value(2)

    def test_extract_data(self):
        pass

    def test_index_parameter_type(self):

        seq1 = Sequences([[1.,1,1],[2.,2,2.]])
        assert IndexParameterType(seq1)=='IMPLICIT_TYPE'
        seq1 = Sequences([[1.,1,1],[2.,2,2.]], IndexParameterType="TIME")
        assert IndexParameterType(seq1)=='TIME'
        seq1 = Sequences([[1.,1,1],[2.,2,2.]],IndexParameterType="POSITION" )
        assert IndexParameterType(seq1)=='POSITION'


    def test_constructors(self):
        # heterogeneous or homogeneous type
        seq1 = Sequences([1, 2, 3, 4])
        seq1 = Sequences([1, 2, 3, 4.])
        assert seq1.nb_sequence == 1
        assert seq1.nb_variable == 1

        # single sequence multivariate
        seq2 = Sequences([[1,2],[3,4], [5,6]])
        assert seq2.nb_sequence==1
        assert seq2.nb_variable==2
        #ambiguous case (length>5) 
        seq2 = Sequences([[1,2,3,4,5,6],[3,4,3,4,5,6], [5,6,4,5,6,7]])
        assert seq2.nb_sequence==3
        assert seq2.nb_variable==1


        # univariates sequences
        seq3 = Sequences([[1,2],[3,4], [5,6,7]])
        assert seq3.nb_sequence==3
        assert seq3.nb_variable==1

        # general case
        seq4 = Sequences([ [[1,2],[3,4]], [[21,22],[23,24]], [[31,32],[33,34], [35,36] ]])
        assert seq4.nb_sequence==3
        assert seq4.nb_variable==2

        seq4 = Sequences([ [[1,2],[3,4]], [[21,22],[23,24]], [[31,32],[33,34], [35,36] ]], 
            VertexIdentifiers=[[1,2],[3,4],[5,6,7]], Identifiers=[1,2,3])

        seq4 = Sequences([ [[1,2],[3,4]], [[21,22],[23,24]], [[31,32],[33,34], [35,36] ]], 
            IndexParameterType="POSITION", IndexParameter=[[0,1,10], [2,3,11], [4,5,6,12]])

        seq4 = Sequences([ [[1,2],[3,4]], [[21,22],[23,24]], [[31,32],[33,34], [35,36] ]], 
            IndexParameterType="TIME", IndexParameter=[[0,1], [2,3], [4,5,6]])




    def test_constructor_one_sequence(self):
        # two sequences with 2 variables (int)
        s = Sequences([[1,1,1],[2,2,2]])
        assert s
        # two sequences with 2 variables (real)
        s = Sequences([[1.,1.,1.],[2.,2.,2.]])
        assert s
        # two sequences with 2 variables mix of (real and int)
        # works because the first one is float so all others are assume to be float as well
        s = Sequences([[1.,1.,1.],[2.,2.,2]])
        assert s
        # here it fails because the first number is int but others may be float
        try:
            s = Sequences([[1,1.,1.],[2.,2.,2]])
            assert False
        except:
            assert True

    def test_constructor_two_sequences(self):
        # two sequences with 2 variables (int)
        s = Sequences([[[1,1,1],[12,12,12]],[[2,2,2],[22,23,24]]])
        assert s
        s = Sequences([[[1.,1.,1.],[1.,1.,1.]],[[2.,2.,2.],[2.,2.,2.]]])
        assert s

    def test_container(self):
        # first index is the sequence and second index is the variable
        s = self.seqn
        assert s[0,0] == [1,1,1]
        assert s[0,1] == [12,12,12]
        assert s[1,0] == [2,2,2]
        assert s[1,1] == [22,23,24]
        assert s[1,1][1] == 23
        assert len(s) == 2

        s = self.seq1
        assert s[0,0] == [1,1,1]
        assert s[0,1] == [2,2,2]
        assert s[0,0][0] == 1
        assert len(s) == 1



    def test_value_select(self):
        "test_value_select implemented but need to be checked"
        seqn = self.seqn
        a = seqn.value_select(1, 1, 2,True)
        assert a
        assert str(ValueSelect(seqn, 1, 1, 2)) == str(seqn.value_select(1,1,2, True))

    def test_select_variable(self):
        "test_select_variable implemented but need to be checked (index issue)"
        # !!!!!!!! NEED to CHECK THE INDEX 0, 1 , ... or 1,2,....
        # what about identifiers ?
        # Variable seems to start at 1 not 0
        s = self.seqn
        #select variable 1
        select = s.select_variable([1], keep=True)
        assert select[0,0] == [1]
        assert select[1,0] == [2]


    def test_select_individual(self):
        #select one or several sequences
        s = self.seqn

        # select all
        select = s.select_individual([0,1], keep=True)
        assert s.display() == select.markovian_sequences().display()

        #select first sequence only
        select = s.select_individual([0], keep=True)
        assert select[0,0] == [1, 1, 1]
        assert select[0,1] == [12, 12, 12]
        try:
            select[1,0]
            assert False
        except:
            assert True

        #select second sequence only
        select = s.select_individual([1], keep=True)
        assert select[0,0] == [2, 2, 2]
        assert select[0,1] == [22, 23, 24]
        try:
            select[1,0]
            assert False
        except:
            assert True

    def test_shift_seqn(self):
        s = self.seqn
        shifted = s.shift(1,2)
        assert shifted[0,0] == [3,1,1]
        assert shifted[0,1] == [14,12,12]
        assert shifted[1,0] == [4,2,2]
        assert shifted[1,1] == [24,23,24]


    def test_shift_seq1(self):
        s = self.seq1
        shifted = s.shift(1,2)
        assert shifted[0,0] == [3, 1, 1]
        assert shifted[0,1] == [4, 2, 2]

    def test_merge (self):
        s1 = self.seqn
        s2 = self.seqn
        s3 = self.seq1

        sall = s1.merge([s2])
        assert sall.nb_sequence == 4
        assert sall.nb_variable == 3

        sall = s1.merge([s3])
        assert sall.nb_sequence == 3
        assert sall.nb_variable == 3

    def test_merge_and_Merge(self):
        s1 = self.seqn
        s2 = self.seqn

        a = s1.merge([s2])
        b = s2.merge([s1])
        v = Merge(s1, s2)

        assert str(a) == str(b)
        assert str(a) == str(v)

    def test_imcompatible_merge(self):
        s1 = self.seqn
        s3 = self.seq1
        try:
            s1.merge(s3)
            assert False
        except:
            assert True

    def test_merge_variable(self):
        import copy
        s1 = self.seqn
        s2 = self.seqn
        s3 = self.seq1

        sall = s1.merge_variable([s2],1) # why 1 ? same result with 2 !
        assert sall.nb_sequence == 2
        assert sall.nb_variable == 6

        # sall =  s1.merge_variable([s3],1)
        #assert sall.nb_sequence == 2
        #assert sall.nb_variable == 3

    def test_merge_variable_and_MergeVariable(self):
        s1 = self.seqn
        s2 = self.seqn
        s3 = self.seq1

        a = s1.merge_variable([s2],1)
        b = s2.merge_variable([s1],1)
        v = MergeVariable(s1,s2)

        assert str(a) == str(b)
        assert str(a) == str(v)


    def test_cluster_step(self):
        seq1 = Sequences([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        assert str(Cluster(seq1, "Step", 1, 2)) == str(seq1.cluster_step(1, 2, False))
        seqn = Sequences([[[1, 2, 3], [1, 3, 1]], [[4, 5, 6], [7,8,9]]])
        assert str(Cluster(seqn, "Step", 1, 2)) == str(seqn.cluster_step(1, 2, False))

    def test_cluster_limit(self):
        seq1 = Sequences([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        assert str(Cluster(seq1, "Limit", 1, [2])) == \
            str(seq1.cluster_limit(1, [2], False))
        seqn = Sequences([[[1, 2, 3], [1, 3, 1]], [[4, 5, 6], [7,8,9]]])
        assert str(Cluster(seqn, "Limit", 1, [2, 4, 6])) == \
            str(seqn.cluster_limit(1, [2, 4 ,6], False))

    def test_transcode(self):
        """This functionality need to be checked.

        See also the vector case!"""
        seq = self.seqn
        assert  str(seq.transcode(1, [0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,01,1,1,1,1,0,0], False))==\
            str(Transcode(seq, 1, [0,1,1,1,1,1,1,1,1,1,1,1,0,0,0,01,1,1,1,1,0,0]))

    def test_reverse(self):
        """reverse to be checked. seems to give same output as input"""
        s = self.seqn
        s.reverse()


    def test_max_length(self):
        s = self.data
        assert s.max_length==30

    def test_get_max_value(self):
        s = self.data
        assert s.get_max_value(0)==2

    def test_get_min_value(self):
        s = self.data
        assert s.get_min_value(0)==0

    def test_get_length(self):
        s = self.data
        assert s.get_length(0)==30
        assert s.get_length(1)==22

    def test_difference(self):
        data = self.data
        assert str(Difference(data, 1)) == str(data.difference(1, False))
        res = Difference(data, 1)
        assert res.cumul_length == 50

    def test_cumulate(self):
        #see also test_cumulate for more tests
        s = self.data
        res = Cumulate(s)
        assert res.cumul_length == 52


    def test_extract_vectors(self):
        """see test_extract_vectors"""
        ExtractVectors(self.data, "Length")

    def _test_index_parameter_extract(self):
        """fixme: markovian_sequences should be in wrapper ? """
        aml = self.data.index_parameter_extract(0, 29).markovian_sequences()
        mod = IndexParameterExtract(self.data, 0, MaxIndex=29)
        assert str(aml) == str(mod)


    def test_index_parameter_select(self):
        """test to be done"""
        pass

    def test_recurrence_time_sequences(self):
        aml = self.data.recurrence_time_sequences(1, 1)
        mod = RecurrenceTimeSequences(self.data, 1, 1)
        assert str(aml.markovian_sequences()) == str(mod)

    def test_remove_run(self):
        """test to be done"""
        pass
    def test_transform_position(self):
        """test to be done"""
        pass
    def test_segmentation_extract(self):
        """test to be done"""
        pass
    def test_variable_scaling(self):
        """test to be done"""
        pass
    def test_remove_index_parameter(self):
        """test to be done"""
        self.data


    def test_write_mtg(self):
        import os
        self.data.mtg_write('test.mtg', [1,2])
        os.remove('test.mtg')

    def test_SaveMTG(self):
        SaveMTG(self.data,Filename='test.mtg', Type=['N'])

    def test_split(self):
        #markovian sequences
        data = Sequences(get_shared_data('vanille_m.seq'))
        Split(data, 2)

    def test_initial_run(self):
        from openalea.sequence_analysis import ComputeInitialRun
        #markovian sequences
        data = Sequences(get_shared_data('vanille_m.seq'))
        ComputeInitialRun(data)

"""

seq.segmentation_extract
seq.pointwise_average
seq.cross
seq.get_index_parameter_type
seq.sojourn_time_sequences
seq.remove_index_parameter
seq.round
"""

if __name__ == "__main__":
    runTestClass(Test())
