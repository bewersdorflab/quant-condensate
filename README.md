# pyme-plugin
Template for making plug-ins for the Python-Microscopy Environment (PYME)

## Instructions
1. Clone the repository, add your recipe module files, dsviewer module files, and visgui module files to the respective
folders in `pyme-plugin/package_name`. 
2. Change the name of the ```pyme-plugin/package_name``` directory.
3. In the `_etc/PYME/plugins/` folder, go into the `dsviewer`, `recipes`, and `visgui` folders and edit
`package_name.txt` to include your actual package name, and to point to each of your module files with a new line. For example, two dsviewer module files could have a `package_name.txt` file of:
```
mynewpackage.dsviewer_modules.my_core_functions
mynewpackage.dsviewer_modules.my_other_dsview_modules
```
4. Change the filenames of the `package_name.txt` files to match your package.
5. Run `python install_plugin.py` to make PYME aware of your modules.
6. Install your package however you need to install it e.g. distutils, etc.
 
