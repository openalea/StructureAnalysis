"""Build an HTML report given a model."""

from .config import OUTPUT_PATH
from .transition_graph_display import TransitionGraphDisplay

import Gnuplot
import logging
import matplotlib
import numpy as np
import os
import seaborn as sns
import unidecode
import webbrowser
import shutil

class Htmlreport(object):
    """Build an HTML report given a model."""
    def __init__(self, model, graph_probability_threshold=0.01,
                 graph_decimals_displayed=3,
                 output_path=None):

        if output_path is None:
            self._output_path = OUTPUT_PATH
        else:
            self._output_path = output_path
        self._model = model
        self._graph_probability_threshold = graph_probability_threshold
        self._graph_decimals_displayed = graph_decimals_displayed
        self._transition_graph_file_path = os.path.join(self._output_path, self._model.model_id + '-transition-graph')
        self._html_report_file_path = os.path.join(self._output_path, self._model.model_id + '-report' + '.html')

        self._nb_states = self._model.k

        self._model.update_restored_data()
        #self._model.eye_movement_data.restored_data.to_excel(self._model.eye_movement_data.restored_data_file_path,
        #                                                   index=False)
        self._generate_transition_graph()
        self._generate_sojourn_distribution_plots()
        self._generate_emission_distribution_plots()
        # self._add_reading_speed_per_strategy_column()
        # self._generate_indicators_plots()
        self._prepare_html_code()

    def _generate_transition_graph(self):
        TransitionGraphDisplay(self._model.hsmc_file, self._transition_graph_file_path,
                               probability_threshold=self._graph_probability_threshold,
                               decimals_displayed=self._graph_decimals_displayed, display=False)

    def _generate_sojourn_distribution_plots(self):
        for i in xrange(0, self._model.k):
            try:
                plot_file = os.path.join(".", self._model.model_id +
                                         '-sojourn-distribution-state-' + str(i))
                self._model.hsmm.extract(5, 0, i).plot_write(plot_file, "")
                g = Gnuplot.Gnuplot()
                g("set terminal png size 1024,768")
                g("set output '" + plot_file + ".png'")
                g("load '" + plot_file + ".plot'")
                g.close()
                for f in [plot_file + ".png", plot_file + ".plot", plot_file + ".print"]:
                    try:
                        shutil.move(f, self._output_path)
                    except:
                        os.remove(os.path.join(self._output_path, f))
                        shutil.move(f, self._output_path)
                for f in [plot_file + "0.dat", plot_file + "1.dat"]:
                    try:
                        shutil.move(f, self._output_path)
                    except:
                        os.remove(os.path.join(self._output_path, f))
                        shutil.move(f, self._output_path)
            except Exception as e:
                logging.warning(e)

    def _generate_emission_distribution_plots(self):
        for i in xrange(0, len(self._model.parameters.output_processes)):
            for j in xrange(0, self._model.k):
                    try:
                        plot_file = os.path.join(".", self._model.model_id + '-output-process-' +
                                                 str(i + 1) + '-state-' + str(j))
                        self._model.hsmm.extract(1, i + 1, j).plot_write(
                            plot_file, "")
                        g = Gnuplot.Gnuplot()
                        g("set terminal png size 1024,768")
                        g("set output '" + plot_file + ".png'")
                        g("load '" + plot_file + ".plot'")
                        g.close()
                        for f in [plot_file + ".png", plot_file + ".plot", plot_file + ".print"]:
                            try:
                                shutil.move(f, self._output_path)
                            except:
                                os.remove(os.path.join(self._output_path, f))
                                shutil.move(f, self._output_path)
                        for f in [plot_file + "0.dat", plot_file + "1.dat"]:
                            try:
                                shutil.move(f, self._output_path)
                            except:
                                os.remove(os.path.join(self._output_path, f))
                                shutil.move(f, self._output_path)
                    except Exception, e:
                        logging.warning(e)


    def _prepare_html_code(self):
        import pandas as pd
        phase_or_state = 'STATES'

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
        self._html_code += os.path.join(".", 'style.css')
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
        self._html_code += os.path.join(".", self._model.model_id + '-transition-graph.png')
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
            self._html_code += os.path.join(".", self._model.model_id + '-sojourn-distribution-state-' +
                                            str(i) + '.png')
            self._html_code += '" alt="N/D"></div>\n'
            # emission
            for j in xrange(0, len(self._model.parameters.output_processes)):
                self._html_code += '<div><img src="'
                self._html_code += os.path.join(".", self._model.model_id + '-output-process-' +
                                                str(j + 1) + '-state-' + str(i) + '.png')
                self._html_code += '" alt="N/D"></div>\n'
        self._html_code += '</div>\n'
        # indicator table
        self._html_code += '</body>\n'
        self._html_code += '</html>\n'

    @property
    def html_code(self):
        return self._html_code

    def make_css(self):
        with open(os.path.join(".", 'style.css'), 'w') as f:
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
        current_pos_file = "." + os.sep + 'style.css'
        try:
            shutil.move(current_pos_file, self._output_path)
        except:
            os.remove(os.path.join(self._output_path, "style.css"))
            shutil.move(current_pos_file, self._output_path)            

    def make_html(self, open_in_web_browser=False):
        self.make_css()
        with open(self._html_report_file_path, 'w') as f:
            f.write(self._html_code)
        if open_in_web_browser:
            webbrowser.open(url=self._html_report_file_path, new=2)

    def make_latex(self):
        logging.info('Not implemented yet.')
