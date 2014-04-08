#!/usr/bin/env python

from distutils.core import setup

class Dummy:
    pass
pkginfo = Dummy()
execfile('SFit/__pkginfo__.py', pkginfo.__dict__)                               

setup (name = "ScrewFit",
       version = "1.0.4",
       description = "Secondary Structure Determination and Characterization",
       long_description =
"""ScrewFit is an Open Source program for the determination and
characterization of secondary structure of proteins. It provides an
efficient description of structural motifs based on simple geometrical
assumptions.
""",
       author = "Paolo Calligari",
       author_email = "paolo.calligari@sissa.it",
       url = "http://dirac.cnrs-orleans.fr/plone/software/screwfit",
       license = "CeCILL",
       packages = ['SFit'],
       scripts = ['screwfit'],
       )

