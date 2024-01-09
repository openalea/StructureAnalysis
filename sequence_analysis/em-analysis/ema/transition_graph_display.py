from .config import OUTPUT_PATH
import functools
import graphviz as gv
import numpy as np
import os


class TransitionGraphDisplay:

    def __init__(self, hsmc_file, output_file=os.path.join(
                 OUTPUT_PATH, 'transition_graph'),
                 probability_threshold=0.01, decimals_displayed=4, display=True):
        self._hsmc_file = hsmc_file
        self.output_file = output_file
        self.decimals_displayed = decimals_displayed
        self.probability_threshold = probability_threshold
        self.graph = self.build_graph()
        self.save_graph()
        if display:
            self.display_graph()

    @staticmethod
    def add_nodes(graph, nodes):
        for n in nodes:
            if isinstance(n, tuple):
                graph.node(n[0], **n[1])
            else:
                graph.node(n)
        return graph

    @staticmethod
    def add_edges(graph, edges):
        for e in edges:
            if isinstance(e[0], tuple):
                graph.edge(*e[0], **e[1])
            else:
                graph.edge(*e)
        return graph

    @staticmethod
    def apply_styles(graph, styles):
        graph.graph_attr.update(
            ('graph' in styles and styles['graph']) or {}
        )
        graph.node_attr.update(
            ('nodes' in styles and styles['nodes']) or {}
        )
        graph.edge_attr.update(
            ('edges' in styles and styles['edges']) or {}
        )
        return graph

    def get_transition_states_style(self):
        style = {
            'graph': {
                'label': '',
                'fontsize': '16',
                #'fontcolor': 'white',
                #'bgcolor': '#333333',
                'rankdir': 'BT',
            },
            'nodes': {
                'fontname': 'Helvetica',
                'shape': 'circle',
                'fontcolor': 'white',
                'color': 'white',
                'style': 'filled',  # invis
                'fillcolor': '#006699',
            },
            'edges': {
                #'style': 'dashed',
                #'color': 'white',
                'arrowhead': 'open',
                'fontname': 'Courier',
                'fontsize': '12',
                #'fontcolor': 'white',
            }
        }
        return style

    @staticmethod
    def get_initialization_states_style():
        style = {
            'nodes': {
                'style': 'invis',
            }
        }
        return style

    def build_graph(self):
        initialization_probabilities = []
        transition_probabilities = []
        k = -1

        with open(self._hsmc_file) as f:
            l = f.readline()
            while(bool(l)):
                if '#' in l:
                    pass
                if 'STATES' in l:
                    k = int(l.split()[0])
                if 'INITIAL_PROBABILITIES' in l:
                    l = f.readline()
                    initialization_probabilities = np.array(
                        map(float, l.split()))
                    assert(len(initialization_probabilities) == k)
                if 'TRANSITION_PROBABILITIES' in l:
                    transition_probabilities = np.zeros(
                        shape=(k, k))
                    for i in xrange(0, k):
                        l = f.readline()
                        row = map(float, l.split())
                        transition_probabilities[i, :] = row
                    break
                l = f.readline()

        initialization_probabilities = np.around(
            initialization_probabilities, decimals=self.decimals_displayed)
        transition_probabilities = np.around(
            transition_probabilities, decimals=self.decimals_displayed)

        digraph = functools.partial(gv.Digraph, format='png')
        g = digraph()
        nodes = [str(i) for i in xrange(0, k)]
        edges = []
        for i in xrange(0, k):
            for j in xrange(0, k):
                if transition_probabilities[i, j] > max(0, self.probability_threshold):
                    # grey lvl in rgb(x, x, x)
                    edge_color = (1 - transition_probabilities[i, j]) * 255
                    # convertion in hexadecimal
                    edge_color = '#%02x%02x%02x' % (edge_color, edge_color, edge_color)
                    edges.append(((str(i), str(j)), {'label': str(transition_probabilities[i, j]), 'color': edge_color}))
                else:
                    edge_color = '#ffffff'
                    edges.append(((str(i), str(j)), {'color': edge_color}))

        invisible_nodes = [str(i) for i in xrange(k, k*2)]

        g_invisible = digraph()
        g_invisible = self.add_nodes(g_invisible, invisible_nodes)
        g_invisible = self.apply_styles(g_invisible, self.get_initialization_states_style())

        g = self.add_nodes(g, nodes)
        g = self.add_edges(g, edges)
        g = self.apply_styles(g, self.get_transition_states_style())

        g.subgraph(g_invisible)

        for i in xrange(0, k):
            parent_id = i+k
            if initialization_probabilities[i] > self.probability_threshold:
                # grey lvl in rgb(x, x, x)
                edge_color = (1 - initialization_probabilities[i] + 0.2) * 255
                # convertion in hexadecimal
                edge_color = '#%02x%02x%02x' % (255, edge_color, edge_color)
                g.edge(str(parent_id), str(i), weight='2', color=edge_color, label=str(initialization_probabilities[i]))
            else:
                edge_color = '#ffffff'
                g.edge(str(parent_id), str(i), weight='2', color=edge_color)
        return g

    def save_graph(self):
        self.graph.render(self.output_file)
