# quant-condensate
Cluster analysis for Courchaine et al, 2020


## Prereq's:
1. Install Python 3.6 using miniconda or anaconda
2. Install python-microscopy environment (instructions at http://python-microscopy.org/doc/Installation/InstallationWithAnaconda.html, though the manual route may be necessary at the moment if you want to load in two-f tiffs as this requires changing a line in pyme (andrew will change pyme to fix this when he can, though reminders never hurt)
3. Make sure you have scipy>=1.3 installed.
4. Clone this repository and install via setuptools (run `python setup.py develop` or `python setup.py install` in the top directory)
5. Register the plugin modules with PYME by running `"`python quant_condensate/install_plugin.py", again from the top directory
