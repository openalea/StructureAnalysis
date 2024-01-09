"""Build an HTML report given a model."""
from .config import OUTPUT_PATH
from .transition_graph_display import TransitionGraphDisplay
import visuscanpath
import Gnuplot
import logging
import matplotlib
import numpy as np
import os
import seaborn as sns
import unidecode
import webbrowser


class Htmlreport(object):
    """Build an HTML report given a model."""
    def __init__(self, model, png_path, graph_probability_threshold=0.01,
                 graph_decimals_displayed=3,
                 number_of_text_restorations_to_display=100, 
                 subj=None, text=None, phase=False,
                 output_path=None):

        if output_path is None:
            self._output_path = OUTPUT_PATH
        else:
            self._output_path = output_path
        self._model = model
        self._png_path = png_path
        self._graph_probability_threshold = graph_probability_threshold
        self._graph_decimals_displayed = graph_decimals_displayed
        self._number_of_text_restorations_to_display = number_of_text_restorations_to_display
        self._transition_graph_file_path = os.path.join(self._output_path, self._model.model_id + '-transition-graph')
        self._html_report_file_path = os.path.join(self._output_path, self._model.model_id + '-report' + '.html')
        self._model.eye_movement_data.restored_data_file_path = os.path.join(self._model._tmp_path,
                                                                             self._model.model_id + '-restored-data' +
                                                                             '.xls')
        self._fdur_file_path = os.path.join(self._output_path, self._model.model_id + '-boxplot-fdur.png')
        self._sacamp_file_path = os.path.join(self._output_path, self._model.model_id + '-boxplot-sacamp.png')
        self._winc_file_path = os.path.join(self._output_path, self._model.model_id + '-boxplot-winc.png')
        self._cinc_file_path = os.path.join(self._output_path, self._model.model_id + '-boxplot-cinc.png')
        self._cosinst_file_path = os.path.join(self._output_path, self._model.model_id + '-boxplot-cosinst.png')
        self._coscum_file_path = os.path.join(self._output_path, self._model.model_id + '-boxplot-coscum.png')
        self._reading_speed_file_path = os.path.join(self._output_path, self._model.model_id + '-boxplot-readingspeed.png')
        self._phase = phase
        self._nb_states = self._model.k
        if phase:
            self._nb_states -= 1

        self._model.update_restored_data()
        #self._model.eye_movement_data.restored_data.to_excel(self._model.eye_movement_data.restored_data_file_path,
        #                                                   index=False)
        self._generate_transition_graph()
        self._generate_sojourn_distribution_plots()
        self._generate_emission_distribution_plots()
        # self._add_reading_speed_per_strategy_column()
        # self._generate_indicators_plots()
        self._generate_text_restoration_plots(subj, text)
        self._prepare_html_code()

    def _generate_transition_graph(self):
        TransitionGraphDisplay(self._model.hsmc_file, self._transition_graph_file_path,
                               probability_threshold=self._graph_probability_threshold,
                               decimals_displayed=self._graph_decimals_displayed, display=False)

    def _generate_sojourn_distribution_plots(self):
        for i in xrange(0, self._model.k):
            try:
                plot_file = os.path.join(self._output_path, self._model.model_id +
                                         '-sojourn-distribution-state-' + str(i))
                self._model.hsmm.extract(5, 0, i).plot_write(plot_file, "")
                g = Gnuplot.Gnuplot()
                g("set terminal png size 1024,768")
                g("set output '" + plot_file + ".png'")
                g("load '" + plot_file + ".plot'")
            except Exception as e:
                logging.warning(e)

    def _generate_emission_distribution_plots(self):
        for i in xrange(0, len(self._model.parameters.output_processes)):
            for j in xrange(0, self._model.k):
                    try:
                        plot_file = os.path.join(self._output_path, self._model.model_id + '-output-process-' +
                                                 str(i + 1) + '-state-' + str(j))
                        self._model.hsmm.extract(1, i + 1, j).plot_write(
                            plot_file, "")
                        g = Gnuplot.Gnuplot()
                        g("set terminal png size 1024,768")
                        g("set output '" + plot_file + ".png'")
                        g("load '" + plot_file + ".plot'")
                    except Exception, e:
                        logging.warning(e)

    def _generate_indicators_plots(self):
        phase_or_state = 'STATES'
        if self._phase:
            phase_or_state = 'PHASE'

        matplotlib.rcParams['figure.figsize'] = (20.0, 10.0)
        df = self._model.eye_movement_data.restored_data

        plot = sns.violinplot(x="PHASE", y="FDUR", data=df)
        plot.get_figure().savefig(self._fdur_file_path)

        plot = sns.violinplot(x="PHASE", y="SACAMP", data=df)
        plot.get_figure().savefig(self._sacamp_file_path)

        df.WINC = df.WINC.astype(float)
        plot = sns.violinplot(x="PHASE", y="WINC", data=df)
        plot.get_figure().savefig(self._winc_file_path)

        if "CINC" in list(df.keys()):
            df.CINC = df.CINC.astype(float)
            plot = sns.violinplot(x="PHASE", y="CINC", data=df)
            plot.get_figure().savefig(self._cinc_file_path)

        if "COSINST" in list(df.keys()):
            import pandas
            df['COS_INST'] = pandas.Series(df.COSINST.astype(float), index=df.index)
        elif "COS_INST" in list(df.keys()):
            df.COS_INST = df.COS_INST.astype(float)
        plot = sns.violinplot(x="PHASE", y="COS_INST", data=df)
        plot.get_figure().savefig(self._cosinst_file_path)

        if "COSCUM" in list(df.keys()):
            df['COS_CUM'] = pandas.Series(df.COSCUM.astype(float), index=df.index)
        elif "COS_CUM" in list(df.keys()):
            df.COS_CUM = df.COS_CUM.astype(float)
        plot = sns.violinplot(x="PHASE", y="COS_CUM", data=df)
        plot.get_figure().savefig(self._coscum_file_path)

        if "CINC" in list(df.keys()):
            df.CINC = df.CINC.astype(float)
            plot = sns.violinplot(x="PHASE", y="CINC", data=df)
            plot.get_figure().savefig(self._cinc_file_path)

        plot = sns.violinplot(x="PHASE", y="READINGSPEED", data=df)
        plot.get_figure().savefig(self._reading_speed_file_path)

        """
        df['WINCTMP'] = df['WINC'] > 0
        logging.info(df[[phase_or_state, 'WINCTMP']].groupby(["PHASE"]).sum())
        logging.info(df[[phase_or_state, 'WINCTMP']].groupby(["PHASE"]).count())
        logging.info(df[[phase_or_state, 'FDUR']].groupby(["PHASE"]).sum())

        # proportions of saccade directions per reading strategies
        sacdir = df[['SACDIR', phase_or_state]].groupby([phase_or_state, 'SACDIR']).size()
        logging.info(sacdir)
        """

    def _generate_text_restoration_plots(self, subj=None, text=None):
        import random
        random.seed(1)
        np.random.seed(1)
        visuscanpath.config.X_COL = 'X'
        visuscanpath.config.Y_COL = 'Y'
        visuscanpath.config.FIXATION_DURATION_COL = 'FDUR'
        visuscanpath.config.SUBJECT_COL = 'SUBJ_NAME'
        visuscanpath.config.TEXT_NAME_COL = 'TEXT'
        rdf = self._model.eye_movement_data.restored_data
        """ favorite F texts to display
        self.text_restorations_display_list = [[1, 1], [1, 3], [1, 5], [1, 6],
                                                 [1, 9], [1, 54], [2, 5],
                                                 [2, 6], [2, 9], [2, 54],
                                                 [4, 3], [4, 5], [4, 6],
                                                 [4, 9], [4, 54], [6, 9]]
        """
        rdf['READMODETMP'] = rdf['READMODE']
        """
        rdf.at[rdf.READMODETMP == 0, 'READMODETMP'] = 'LReg'
        rdf.at[rdf.READMODETMP == 1, 'READMODETMP'] = 'Reg'
        rdf.at[rdf.READMODETMP == 2, 'READMODETMP'] = 'Ref'
        rdf.at[rdf.READMODETMP == 3, 'READMODETMP'] = 'Pr'
        rdf.at[rdf.READMODETMP == 4, 'READMODETMP'] = 'LPr'
        """
        phase_or_state = 'STATES'
        if self._phase:
            phase_or_state = 'PHASE'
        
        def subj_name(s):
            res = str(s)
            if s < 10:
                res = "s0" + res
            else:        
                res = "s" + res
            return res
            
        if not("SUBJ_NAME" in list(rdf.keys())):
            import pandas
            rdf['SUBJ_NAME'] = pandas.Series(list(map(subj_name, rdf.SUBJ)), index=rdf.index)            

        self.text_restorations_display_list = []
        subj_text_list = rdf[['SUBJ_NAME', 'TEXT']].groupby(['SUBJ_NAME', 'TEXT']).count().reset_index()
        if (subj is None):
            rdm_index = np.random.choice(subj_text_list.index, self._number_of_text_restorations_to_display, replace=True)
        else:
            assert len(subj) == len(text), "Lengths of nb texts and nb sujb do not match"
            text_list = sorted(list(set(subj_text_list.TEXT)))
            rdm_index = []
            for i in subj_text_list.index:
                for j in range(len(subj)):
                    if ((subj_text_list.loc[i][0] == subj_name(subj[j])) and
                        (subj_text_list.loc[i][1] == text_list[text[j]])):
                        rdm_index += [i]
        for i in rdm_index:
            self.text_restorations_display_list.append([subj_text_list.loc[i][0], subj_text_list.loc[i][1]])
            visuscanpath.plot_scanpath(dataframe=self._model.eye_movement_data.restored_data, display=False,
                                       img_path=os.path.join(self._png_path,
                                                             unidecode.unidecode(subj_text_list.loc[i][1]) + ".png"),
                                       subject=subj_text_list.loc[i][0],
                                       text_name=subj_text_list.loc[i][1],
                                       #print_col='READMODETMP',
                                       hue=phase_or_state,
                                       output_file=os.path.join(self._output_path, self._model.model_id + "-subj-" +
                                                                subj_text_list.loc[i][0] +
                                                                "-text-" +
                                                                unidecode.unidecode(subj_text_list.loc[i][1]) + ".png"),

                                       )

    def _prepare_html_code(self):
        phase_or_state = 'STATES'
        if self._phase:
            phase_or_state = 'PHASE'

        self._html_code = '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0'
        self._html_code += 'Transitional//EN"'
        self._html_code += '"http://www.w3.org/TR/xhtml1/DTD/xhtml1-'
        self._html_code += 'transitional.dtd">\n'
        self._html_code += '<html xmlns="http://www.w3.org/1999/xhtml">\n'
        self._html_code += '<head>\n'
        self._html_code += '<meta http-equiv="Content-Type"'
        self._html_code += 'content="text/html;charset=utf-8"/>\n'
        self._html_code += '<meta http-equiv="Content-Style-Type"'
        self._html_code += 'content="text/css" />\n'
        self._html_code += '<meta name="generator" content="pandoc" />\n'
        self._html_code += '<title>EMA report '
        self._html_code += self._model.model_id
        self._html_code += '</title>\n'
        self._html_code += '<style type="text/css">code{white-space: pre;}'
        self._html_code += '</style>\n'
        self._html_code += '<link rel="stylesheet" type="text/css"'
        self._html_code += 'href="'
        self._html_code += os.path.join(self._output_path, 'style.css')
        self._html_code += '">\n'
        self._html_code += '</head>\n'
        self._html_code += '<body>\n'
        self._html_code += '<div id="transition_graph_wrap">\n'
        self._html_code += '<div id ="transition_graph_img">\n'
        self._html_code += '<h1>Initialization and transition probabilities '
        self._html_code += 'graph</h1><p>(edge display threshold : '
        self._html_code += str(self._graph_probability_threshold)
        self._html_code += ')</p>\n'
        self._html_code += '<img src="'
        self._html_code += os.path.join(self._output_path, self._model.model_id + '-transition-graph.png')
        self._html_code += '" alt="N/D">\n'
        self._html_code += '</div>\n'
        self._html_code += '<div id="model_criterion">'
        self._html_code += '<h1>Model criterion</h1>\n'
        if bool(self._model.criterion):
            for key, value in self._model.criterion.iteritems():
                self._html_code += '<p>'
                self._html_code += key + ' : ' + str(value)
                self._html_code += '</p>\n'
        self._html_code += '</div>\n'
        self._html_code += '</div>\n'
        self._html_code += '<div id="sojourn_emission_plot_wrap">\n'
        self._html_code += '<div><h1>Sojourn distribution bar charts</h1></div>\n'
        for i in xrange(0, len(self._model.parameters.output_processes)):
            self._html_code += '<div><h1>Emission distribution bar charts : Output process ' + str(i) + '</h1></div>\n'
        # sojourn
        for i in xrange(0, self._model.k):
            self._html_code += '<div><img src="'
            self._html_code += os.path.join(self._output_path, self._model.model_id + '-sojourn-distribution-state-' +
                                            str(i) + '.png')
            self._html_code += '" alt="N/D"></div>\n'
            # emission
            for j in xrange(0, len(self._model.parameters.output_processes)):
                self._html_code += '<div><img src="'
                self._html_code += os.path.join(self._output_path, self._model.model_id + '-output-process-' +
                                                str(j + 1) + '-state-' + str(i) + '.png')
                self._html_code += '" alt="N/D"></div>\n'
        self._html_code += '</div>\n'
        """
        self._html_code += '<div id="boxplots">\n'
        self._html_code += '<h1>Indicators\' boxplots</h1>'
        # indicators
        self._html_code += '<div id="fdur_wrap">\n'
        self._html_code += '<div><img src="'
        self._html_code += os.path.join(
            self._output_path, self._model.model_id + '-boxplot-fdur.png')
        self._html_code += '" alt="whoops"></div>\n'
        self._html_code += '</div>\n'
        self._html_code += '<div id="sacamp_wrap">\n'
        self._html_code += '<div><img src="'
        self._html_code += os.path.join(
            self._output_path, self._model.model_id + '-boxplot-sacamp.png')
        self._html_code += '" alt="whoops"></div>\n'
        self._html_code += '</div>\n'
        self._html_code += '<div id="winc_wrap">\n'
        self._html_code += '<div><img src="'
        self._html_code += os.path.join(
            self._output_path, self._model.model_id + '-boxplot-winc.png')
        self._html_code += '" alt="whoops"></div>\n'
        self._html_code += '</div>\n'
        self._html_code += '<div id="cinc_wrap">\n'
        self._html_code += '<div><img src="'
        self._html_code += os.path.join(
            self._output_path, self._model.model_id + '-boxplot-cinc.png')
        self._html_code += '" alt="whoops"></div>\n'
        self._html_code += '</div>\n'
        self._html_code += '<div id="cosinst_wrap">\n'
        self._html_code += '<div><img src="'
        self._html_code += os.path.join(
            self._output_path, self._model.model_id + '-boxplot-cosinst.png')
        self._html_code += '" alt="whoops"></div>\n'
        self._html_code += '</div>\n'
        self._html_code += '<div id="fdur_wrap">\n'
        self._html_code += '<div><img src="'
        self._html_code += os.path.join(
            self._output_path, self._model.model_id + '-boxplot-coscum.png')
        self._html_code += '" alt="whoops"></div>\n'
        self._html_code += '</div>\n'
        self._html_code += '<div id="readingspeed_wrap">\n'
        self._html_code += '<div><img src="'
        self._html_code += os.path.join(
            self._output_path, self._model.model_id + '-boxplot-readingspeed.png')
        self._html_code += '" alt="whoops"></div>\n'
        self._html_code += '</div>\n'
        """
        # indicator table
        self._html_code += '<h1> Indicators :</h1>\n'
        self._html_code += '<div id="indicator_wrap">\n'
        self._html_code += '<table style="width:100%">\n'
        self._html_code += '<tr>\n'
        self._html_code += '<th> Indicators </th>\n'
        for i in xrange(0, self._nb_states):
            self._html_code += '<th>strategy ' + str(i) + '</th>'
        self._html_code += '</tr>\n'
        # fixation duration
        self._html_code += '<tr>\n'
        self._html_code += '<td>fixation duration (ms)</td>\n'
        fdur_mean = self._model.eye_movement_data.restored_data[[phase_or_state, 'FDUR']].groupby(phase_or_state).mean()
        fdur_std = self._model.eye_movement_data.restored_data[[phase_or_state, 'FDUR']].groupby(phase_or_state).std()
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(fdur_mean.loc[i][0]) + ' &plusmn; '\
                                   + "{:.2f}".format(fdur_std.loc[i][0]) + '</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'
        # saccade amplitude
        self._html_code += '<tr>\n'
        self._html_code += '<td>saccade amplitude (px)</td>\n'
        sacamp_mean = self._model.eye_movement_data.restored_data[[phase_or_state, 'SACAMP']].groupby(phase_or_state).mean()
        sacamp_std = self._model.eye_movement_data.restored_data[[phase_or_state, 'SACAMP']].groupby(phase_or_state).std()
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(sacamp_mean.loc[i][0]) + ' &plusmn; ' \
                                   + "{:.2f}".format(sacamp_std.loc[i][0]) + '</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'
        # reading speed
        df = self._model.eye_movement_data.restored_data
        if "COSINST" in list(df.keys()):
            import pandas
            df['COS_INST'] = pandas.Series(df.COSINST.astype(float), index=df.index)
        if "COSCUM" in list(df.keys()):
            import pandas
            df['COS_CUM'] = pandas.Series(df.COSCUM.astype(float), index=df.index)

        if "NEW_READ_WORDS" in list(df.keys()):
            rs = self._model.eye_movement_data.restored_data[[phase_or_state, 'NEW_READ_WORDS']].groupby(phase_or_state).sum().reset_index()['NEW_READ_WORDS'] /\
                (self._model.eye_movement_data.restored_data[[phase_or_state, 'FDUR']].groupby(phase_or_state).sum().reset_index()['FDUR'] / 1000 / 60)
            self._html_code += '<tr>\n'
            self._html_code += '<td>reading speed (wpm)</td>\n'
            for i in xrange(0, self._nb_states):
                try:
                    self._html_code += '<td>' + "{:.2f}".format(rs[i]) + '</td>\n'
                except Exception, e:
                    self._html_code += '<td>0</td>\n'
            self._html_code += '</tr>\n'
        # word increment
        self._html_code += '<tr>\n'
        self._html_code += '<td>word increment</td>\n'
        winc_mean = self._model.eye_movement_data.restored_data[[phase_or_state, 'WINC']].groupby(phase_or_state).mean()
        winc_std = self._model.eye_movement_data.restored_data[[phase_or_state, 'WINC']].groupby(phase_or_state).std()
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(winc_mean.loc[i][0]) + ' &plusmn; '\
                                   + "{:.2f}".format(winc_std.loc[i][0]) + '</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'
        # character increment
        if "CINC" in list(df.keys()):
            self._html_code += '<tr>\n'
            self._html_code += '<td>character increment</td>\n'
            cinc_mean = self._model.eye_movement_data.restored_data[[phase_or_state, 'CINC']].groupby(phase_or_state).mean()
            cinc_std = self._model.eye_movement_data.restored_data[[phase_or_state, 'CINC']].groupby(phase_or_state).std()
            for i in xrange(0, self._nb_states):
                try:
                    self._html_code += '<td>' + "{:.2f}".format(cinc_mean.loc[i][0]) + ' &plusmn; '\
                                    + "{:.2f}".format(cinc_std.loc[i][0]) + '</td>\n'
                except Exception, e:
                    self._html_code += '<td>0</td>\n'
            self._html_code += '</tr>\n'
        # cos inst
        self._html_code += '<tr>\n'
        self._html_code += '<td>instant cosine</td>\n'
        cos_inst_mean = self._model.eye_movement_data.restored_data[[phase_or_state, 'COS_INST']].groupby(phase_or_state).mean()
        cos_inst_std = self._model.eye_movement_data.restored_data[[phase_or_state, 'COS_INST']].groupby(phase_or_state).std()
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(cos_inst_mean.loc[i][0]) + ' &plusmn; '\
                                   + "{:.2f}".format(cos_inst_std.loc[i][0]) + '</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'
        # cos cum
        self._html_code += '<tr>\n'
        self._html_code += '<td>cumulated cosine</td>\n'
        cos_cum_mean = self._model.eye_movement_data.restored_data[[phase_or_state, 'COS_CUM']].groupby(phase_or_state).mean()
        cos_cum_std = self._model.eye_movement_data.restored_data[[phase_or_state, 'COS_CUM']].groupby(phase_or_state).std()
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(cos_cum_mean.loc[i][0]) + ' &plusmn; '\
                                   + "{:.2f}".format(cos_cum_std.loc[i][0]) + '</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'
        self._html_code += '</table>\n'
        self._html_code += '</div>\n'

        # saccade direction table
        self._html_code += '<h1>Saccade directions :</h1>\n'
        self._html_code += '<div id="saccade_dir_wrap">\n'
        sacdir = self._model.eye_movement_data.restored_data[['SACDIR', phase_or_state]].groupby([phase_or_state, 'SACDIR']).size()
        freq = []
        for i in xrange(0, 5):
            freq.append(sacdir[i] / sacdir[i].sum())
        self._html_code += '<table style="width:100%">\n'
        self._html_code += '<tr>\n'
        self._html_code += '<th>Saccade Direction</th>\n'
        for i in xrange(0, self._nb_states):
            self._html_code += '<th>strategy ' + str(i) + ' - freq (cumul)</th>'
        self._html_code += '</tr>\n'
        self._html_code += '<tr>\n'
        self._html_code += '<td>backward</td>\n'
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(freq[i][0]) + '% (' + str(sacdir[i][0]) + ')</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'

        self._html_code += '<tr>\n'
        self._html_code += '<td>downward</td>\n'
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(freq[i][1]) + '% (' + str(sacdir[i][1]) + ')</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'

        self._html_code += '<tr>\n'
        self._html_code += '<td>forward</td>\n'
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(freq[i][2]) + '% (' + str(sacdir[i][2]) + ')</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'

        self._html_code += '<tr>\n'
        self._html_code += '<td>upward</td>\n'
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(freq[i][4]) + '% (' + str(sacdir[i][4]) + ')</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'

        self._html_code += '<tr>\n'
        self._html_code += '<td>last</td>\n'
        for i in xrange(0, self._nb_states):
            try:
                self._html_code += '<td>' + "{:.2f}".format(freq[i][3]) + '% (' + str(sacdir[i][3]) + ')</td>\n'
            except Exception, e:
                self._html_code += '<td>0</td>\n'
        self._html_code += '</tr>\n'
        self._html_code += '</table>\n'
        self._html_code += '</div>\n'
        # scanpath
        self._html_code += '<div><h1>Text restorations</h1></div>\n'
        self._html_code += '<div id="restauration_wrap">\n'
        for i in xrange(0, len(self.text_restorations_display_list)):
            self._html_code += '<div>'
            self._html_code += '<img src="'
            self._html_code += os.path.join(
                self._output_path, self._model.model_id + '-subj-'
                + str(self.text_restorations_display_list[i][0]) + '-text-'
                + unidecode.unidecode(self.text_restorations_display_list[i][1]) + '.png')
            self._html_code += '" alt="N/D"></div>\n'
        self._html_code += '</div>\n'
        self._html_code += '</body>\n'
        self._html_code += '</html>\n'

    @property
    def html_code(self):
        return self._html_code

    def make_css(self):
        with open(os.path.join(self._output_path, 'style.css'), 'w') as f:
            f.write('table, th, td {\nborder: 1px solid black;\nborder-collapse: collapse;\n}th, td {\npadding: 5px;\ntext-align: left;}\n')
            f.write('div\n')
            f.write('{\n')
            f.write('float: left;\n')
            f.write('width: 100%;\n')
            f.write('text-align: center;\n')
            f.write('}\n')
            f.write('#transition_graph_wrap > div\n')
            f.write('{\n')
            f.write('width : 80%;\n')
            f.write('}\n')
            f.write('div#model_criterion\n')
            f.write('{\n')
            f.write('width: 20%;\n')
            f.write('}\n')
            f.write('#sojourn_emission_plot_wrap\n')
            f.write('{\n')
            f.write('width: 100%;\n')
            f.write('}\n')
            f.write('#sojourn_emission_plot_wrap > div\n')
            f.write('{\n')
            f.write('width: ' + str(100. / (len(self._model.parameters.output_processes) + 1) - 1) + '%;\n')
            f.write('padding-bottom: 0.5%;\n')
            f.write('margin: 0.5%;\n')
            f.write('}\n')
            f.write('#restauration_wrap > div\n')
            f.write('{\n')
            f.write('width: 49%;\n')
            f.write('padding-bottom: 0.5%;\n')
            f.write('margin: 0.5%;\n')
            f.write('}\n')
            f.write('div > img\n')
            f.write('{\n')
            f.write('max-height: 100%;\n')
            f.write('max-width: 100%;\n')
            f.write('display: block;\n')
            f.write('margin-left: auto;\n')
            f.write('margin-right: auto;\n')
            f.write('}\n')

    def make_html(self, open_in_web_browser=False):
        self.make_css()
        with open(self._html_report_file_path, 'w') as f:
            f.write(self._html_code)
        if open_in_web_browser:
            webbrowser.open(url=self._html_report_file_path, new=2)

    def make_latex(self):
        logging.info('Not implemented yet.')
