
this package contains software used by the Sonification of Solar 
Harmonics (SoSH) Project for visualizations and for retrieving solar 
data from the web.  to run it, you will need to install python.  the 
distribution we recommend is anaconda, available from 
https://www.anaconda.com/download .  alternatively, for a more minimal 
installation, one may download conda (anaconda's package manager) from 
https://conda.io/en/latest/miniconda.html and follow the instructions 
given here to install only the modules needed for this project.

three files are provided that specify the environments needed to 
use this package:
soshdraw.env.yml - installs modules needed to generate graphics.
soshdata.env.yml - installs modules needed to download and convert data.
soshfull.env.yml - installs both.

to create and use an environment, you can run, for example,

conda env create -f soshfull.env.yml
conda activate sosh

alternatively, if you have installed the full anaconda, you will 
already have all the modules needed for graphics.  to install the 
two additional modules needed to download data and convert it to 
wav format you can run

conda config --add channels conda-forge
conda install drms
conda install pysoundfile

NOTE: we are currently unable to install pysoundfile for windows 
using conda.  however, you may install it from the python 
package index using pip.  for example, at your windows command prompt 
you can run

type soshdata.env.ylm | find /v "pysoundfile" > winenv.yml
conda env create -f winenv.yml
conda activate soshdata
pip install pysoundfile

for further instructions, see quickstart_visual.txt or 
quickstart_data.txt .  the file instructions_visual.txt contains 
the full documentation for the visualizations.

all included files:
addharmonics.py   instructions_visual.txt  soshdata.env.yml
addradial.py      model.surface.modes      soshdata.py
daynumbers.txt    quickstart_data.txt      soshdraw.env.yml
drawharmonics.py  quickstart_visual.txt    soshfull.env.yml
drawradial.py     README_visual.txt        sosh.py

