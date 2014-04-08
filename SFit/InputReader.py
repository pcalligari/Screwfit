#!/usr/local/bin//python

#
# Written by Paolo Calligari (paolo.calligari@ens.fr)
#
#

from MMTK import *
from MMTK.Trajectory import Trajectory 
from MMTK.PDB import PDBConfiguration
from MMTK.PDBMoleculeFactory import PDBMoleculeFactory
from Scientific.IO.ArrayIO import *
import sys
import string
import os
from getopt import *

from SFit.Utils import printHelp


# DEBUG
import pdb

class Input_data:
    
    def __init__(self):
        self.models = None
        self.input_file = None
        self.title = None
        self.output = None
        self.plot = None
        self.database = None
        self.pdb = None
        self.sequence = None 
        self.variation = None
        self.db_name = None
        self.universe = None
        self.trajectory = None
        self.traj_skip = None
        
    def inputFileRead(self,filename):
        """
        read and process an specification file
        """
        keywords = ['models','input_file','title','plot','output','traj_skip','sequence']
        newvars = {}
        file_text = Utility.readURL(filename)
        exec file_text in vars(sys.modules['__builtin__']), newvars
        for key in keywords:
            if newvars.has_key(key): setattr(self,key,newvars[key])
            else: setattr(self,key,None)
            
        if type(self.models) == str:
            listmodels=map(lambda x : int(x), string.split(self.models,":"))
            if len(listmodels) == 3:
                self.models = range(listmodels[0], listmodels[1]+1, listmodels[2])
            elif len(listmodels) == 2:
                self.models = range(listmodels[0],listmodels[1]+1)
            else:
                raise TypeError("Please verify your selection of models")
        else:
            self.models = self.models
            
        if string.split(self.input_file,".")[-1] == 'pdb' and type(self.models) == int:
            self.database = 0
            conf1 = PDBConfiguration(self.input_file,model=self.models)
            factory = PDBMoleculeFactory(conf1,None, lambda r,a: a.name in ['C','O','N'])
            #self.pdb = factory.retrieveAsymmetricUnit()
            self.pdb = factory.retrieveMolecules()
        elif string.split(self.input_file,".")[-1] == 'pdb' and type(self.models) == list:
            self.database = 0            
        elif string.split(self.input_file,".")[-1] == 'list':
            f=file(self.input_file,'r')
            self.database=map(lambda x: x[:6],f.readlines())
            f.close()
            self.db_name = string.split(self.input_file,".")[0]
            self.pdb = 0
        elif string.split(self.input_file,".")[-1] == 'nc':
            traj=Trajectory(None, self.input_file,mode='r')
            self.trajectory=traj
            self.universe=self.trajectory.universe           
        else:
            self.database = 0
            conf1 = PDBConfiguration(self.input_file)
            factory = PDBMoleculeFactory(conf1,None, lambda r,a: a.name in ['C','O','N'])
            self.pdb = factory.retrieveMolecules()

        if self.models and self.database:
            raise TypeError("You cannot define models for databse analysis")
        
def is_universe(data):
    if data.universe:
        return 1
    else:
        return None
        
            
