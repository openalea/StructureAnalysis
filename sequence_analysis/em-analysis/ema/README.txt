Detail of the src/ema/ folder
-----------------------------
__init__.py: define module ema, logger and global vars s.a. PROJECT_PATH.
eye_movement_data.py: handle eye movement data with a pandas dataframe. Contains several useful functions to perform analysis
html_report.py: run reporting tools on a model, display and plot model params, indicators, criterion
indicator.py: run by html_report to build indicators s.a. fixations per state
initial_probabilities.py: stores parameter pi
model.py: handle parameter interaction between C++ and python wrappers by storing and parsing parameter files (*.hsmc). It also has several useful functions helping perform analysis
model_parameters.py: stores all model parameters
observation_distribution.py: stores observation distributions
occupancy_distribution.py: stores occupancy distributions
output_process.py: stores output processes (and therefore observation_distributions)
png_plot.py: provide functions for plotting scanpaths
sequential_data.py: superclass of eye_movement_data when data provided is sequential and not a dataframe
transition_graph_display.py: provide functions for plotting transition probabilities as a graph
xl_io.py: provide functions to write data into a xl file