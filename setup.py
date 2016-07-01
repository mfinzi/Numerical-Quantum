# -*- coding: utf-8 -*-

# Small setup file to freeze the py's and dependent libraries (numpy, scipy, matplotlib, etc)
# into a self contained package using cx_Freeze. Call python setup.py bdist_msi to make the windows
# installer
import scipy.special
import sys, os
import scipy
from cx_Freeze import setup, Executable

base = None
if sys.platform == 'win32':
    base = 'Win32GUI'

# Scipy libraries fail without these two dll's that come from numpy, but for some reason,
# they aren't pulled in by the regular include
includeList =["libifcoremd.dll","libmmd.dll"]
# Also, scipy's libraries use the __file__ command which differs in a zipped folder,
# so when CX_Freeze zips them, it breaks the library. Instead we copy over all the files from scipy
# assumes you have scipy installed
scipy_path = os.path.dirname(scipy.__file__)
includeList.append(scipy_path)


executables = [
    Executable('fullApp.py', base=base,
    	shortcutName="QM Sim",
    	shortcutDir="DesktopFolder",
    	icon = "ih.ico")
]



build_exe_options = {
	'excludes' : ['collections.sys',
				  'collections._weakref'],
	'include_files' : includeList,
	"packages":[],
	"icon":["ih.ico"]
}

shortcutTable = [
	("DesktopShortcut",        # Shortcut
     "DesktopFolder",          # Directory_
     "QM Sim",           # Name
     "TARGETDIR",              # Component_
     "[TARGETDIR]fullApp.exe",# Target
     None,                     # Arguments
     None,                     # Description
     None,                     # Hotkey
     None,#"ih.ico",      # Icon
     None,                     # IconIndex
     None,                     # ShowCmd
     'TARGETDIR'               # WkDir
     )
    ]

setup(name='wave_sim',
      version='0.45',
      description='QM Sim',
      author='Marc Finzi',
      author_email='mfinzi@hmc.edu',
      executables=executables,
      options = {"build_exe": build_exe_options,
      			 "bdist_msi": {'data': {"Shortcut": shortcutTable}}
      			 }
      )
