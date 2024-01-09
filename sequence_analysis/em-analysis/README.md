# vplants.ema

## Description
This package aims at analyzing oculometric sequences.
It uses the original vplants libraries which are necessary for analyzing sequence (tools, stat_tools, sequence_analysis)
(see lib/).

## What's new ?
[Change log](ChangeLog.txt)

## Compilation Tools
python >= 2.4

scons >= 0.97

g++ >= 3.0 (on linux)

mingw (on windows)

openalea.config >= 0.2

openalea.sconsx >= 0.4

openalea.distx >= 0.3

qt >= 4.2 (on windows)

sphinx >= 1.5.2

## Pip requirements
[pip requirements](requirements.txt)


## Installing ema on ArchLinux

### OS specific dependencies
```bash
sudo pacman -S pip2 scons qt4 boost tk graphviz gnuplot yaourt git
```
```bash
yaourt -S termcap
```

### Python dependencies (installable via pip2)
```bash
sudo -H pip2 install numpy pandas unidecode ggplot xlrd xlwt graphviz pillow nose nbdime
```


### Installing and setting up ema
Clone ema repo
```bash
git clone https://github.com/PyENE/em-analysis.git
```

Gnuplot
```bash
python2.7 ema/lib/gnuplot-py-1.8/setup.py install --user
```

Openalea
```bash
python2.7 ema/lib/openalea/multisetup.py develop --user
```

Vplants.tools
```bash
python2.7 ema/lib/vplants/tools/setup.py  develop --user
```

Vplants.stat_tools
```bash
python2.7 ema/lib/vplants/stat_tools/setup.py  develop --user
```

Vplants.sequence_analysis
```bash
python2.7 ema/lib/vplants/sequence_analysis/setup.py  develop --user
```

Vplants.ema
```bash
python2.7 ema/setup.py develop --user
```

Export libs to LD_LIBRARY_PATH
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/$HOME/.local/lib/python2.7/site-packages/
```

