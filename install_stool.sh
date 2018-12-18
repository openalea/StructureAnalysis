CONDA_PY=27 conda build --py 2.7 bin/conda/libstat_tool -c statiskit -c statiskit/label/temp -c defaults --override-channels
CONDA_PY=27 conda install --py 2.7 -y libstat_tool --use-local -c statiskit -c statiskit/label/temp -c defaults --override-channels

python stwig.py

CONDA_PY=27 conda build --python=2.7 stat_tool/bin/conda/python-stat_tool -c statiskit -c statiskit/label/temp -c defaults --override-channels
CONDA_PY=27 conda install -y python-stat_tool --use-local -c statiskit  -c statiskit/label/temp -c defaults --override-channels
