from MMTK import *
from MMTK.PDB import PDBConfiguration
from MMTK.PDBMoleculeFactory import PDBMoleculeFactory
from Scientific.Geometry import Quaternion, Vector, Transformation
from Scientific.IO.ArrayIO import *
from Scientific.Geometry.Objects3D import Line 
import string
import os
from numpy import zeros, sign, sqrt

from QuaternionFactory import *
from Utils import *


def screwFitAnalysis(pdb_str, nmr_model = None):
	try:
		traj_flag=pdb_str.is_universe	
		if len(pdb_str) == 1: pdb_str = pdb_str[0]
	except:
		traj_flag=None
# input model counting in "mmtk-style"
	results = {}
	if nmr_model is None: nmr_model = 0
# screwfit is calculated on every protein chain in the universe!
	not_avail = []
	for mol in pdb_str: 
		molname = string.replace(str(mol),' ','_')
		if molname[:5] == 'Chain':
			[sen, ang_dist, err], control = angularDistance(mol, traj_flag)
			# la variabile control per il momento non e' usata. 
			num_seq, num_pdb_s, num_pdb_h, h_radius, hand, straig = ChaselsAnalysis(mol, traj_flag)
			results[molname] = num_seq, sen, hand, err, ang_dist, h_radius, straig, nmr_model+1
                else:
			results[molname] = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
                        not_avail.append(mol)
# this option is needed to compute on databases...
#	if len(results) == 1:
#		return results.values()[0]
#	else:
	return results, not_avail

	

def angularDistance(molecule, traj_flag=None):

        proteina=molecule
	control = []
	
        """ Reference frame """
        Crefpos=Vector(0.00,0.00,0.00)
        Orefpos=Vector(0.00,0.00,0.123)
        Nrefpos=Vector(0.111,0.00,-0.0728)

        Catom=Atom('C',position=Crefpos)
        Oatom=Atom('O',position=Orefpos)
        Natom=Atom('N',position=Nrefpos)
        piano=Collection()
        piano.addObject(Catom)
        piano.addObject(Oatom)
        piano.addObject(Natom)		

        universe2=InfiniteUniverse()
        universe2.addObject(piano)
        refconf =copy(universe2.configuration())
        lung=len(proteina) - 1
        Dist=zeros((lung,))
        sen = zeros((lung,))
	
	molname = string.replace(str(molecule),' ','_')
	control.append("%s " %molname)
	control.append("   \n")
        control.append("   \n")

	if traj_flag is None:
		prim=str(proteina[0])
		lpr=string.split(prim,'_')
		j=int(lpr[-1])-1
	else:
		prim=str(proteina[0])
		j=int(prim[-1])
#	catena=string.split(prim)[1][0]
        list_emin = []

        for i in range(lung):

                j += 1
                sen[i]=j
#		nam=string.split(str(proteina[i]))[1][2:5]
                if traj_flag:
			Cpos = proteina[i].peptide.C.position()
			Opos = proteina[i].peptide.O.position()
			Npos = proteina[i+1].peptide.N.position()
		else:
			Cpos = proteina[i].C.position()
			Opos = proteina[i].O.position()
			Npos = proteina[i+1].N.position()

                Catom.setPosition(Cpos)
                Oatom.setPosition(Opos)
                Natom.setPosition(Npos)

                distV=universe2.distanceVector(Cpos,Crefpos)

                Catom.translateBy(distV)
                Oatom.translateBy(distV)
                Natom.translateBy(distV)

                Qm, v, e, emin, emax, rms = findQuaternionMatrix(piano,Catom,  \
                                                                 refconf)

                control.append("\t ---%d-- Peptidic Plane ---- Fit Error %f \n" %(i,emin))
                for q in v:
                        control.append("%f \n" %q)

                rmsD = piano.rmsDifference(refconf)
                AngDist= rmsD/sqrt(emax)
                control.append("\t AngDist : %f\n" %AngDist)
                control.append("            \n")
                control.append("            \n")

                Dist[i] = AngDist
                list_emin.append(emin)

                refconf = copy(universe2.configuration())

        return [sen, Dist, list_emin], control


def ChaselsAnalysis(molecule,traj_flag):
#Length to be considered is the length of the protein less the last
# residue (from which we get just the atom N). The latest peptide plane
# is not involved in calculations for screw motion. So the last residue
# to which apply the screw motion is |len(protein)-2|

        proteina=molecule

        """ Reference frame """
#        Crefpos=Vector(0.00,0.00,0.00)
#        Orefpos=Vector(0.00,0.00,0.123)
#        Nrefpos=Vector(0.111,0.00,-0.0728)

        if traj_flag:
		Crefpos = proteina[0].peptide.C.position()
		Orefpos = proteina[0].peptide.O.position()
		Nrefpos = proteina[1].peptide.N.position()
	else:
		Crefpos = proteina[0].C.position()
		Orefpos = proteina[0].O.position()
		Nrefpos = proteina[1].N.position()

        Catom=Atom('C',position=Crefpos)
        Oatom=Atom('O',position=Orefpos)
        Natom=Atom('N',position=Nrefpos)
        piano=Collection()
        piano.addObject(Catom)
        piano.addObject(Oatom)
        piano.addObject(Natom)

        universe2=InfiniteUniverse()
        universe2.addObject(piano)
        refconf =copy(universe2.configuration())
        lung=len(proteina) - 1
        Dist=zeros((lung,))
        sen=zeros((lung,))

	if traj_flag is None:
		prim=str(proteina[0])
		lpr=string.split(prim,'_')
		j=int(lpr[-1])-1
	else:
		prim=str(proteina[0])
		j=int(prim[-1])
		
	lung2= lung -1
        DatiAsse=[]
        u_orto=[]

        for i in range(lung2):

		if traj_flag:
			Cpos = proteina[i+1].peptide.C.position()
			Opos = proteina[i+1].peptide.O.position()
			Npos = proteina[i+2].peptide.N.position()
		else:
			Cpos = proteina[i+1].C.position()
			Opos = proteina[i+1].O.position()
			Npos = proteina[i+2].N.position()

                Catom.setPosition(Cpos)
                Oatom.setPosition(Opos)
                Natom.setPosition(Npos)

                newconf=copy(universe2.configuration())

                TL1,ROT,TL2, rms = findGenericTransformation(piano,Catom,refconf,newconf)
                Tr = TL1*ROT*TL2
                DatiAsse.append(Tr.screwMotion())
                refconf = copy(universe2.configuration())

        asseIP=[]
        for i in range(len(DatiAsse)):
                if DatiAsse[i][1].length() != 0.:
                        asseIP.append(Line(DatiAsse[i][0],DatiAsse[i][1]))
                else:
                        asseIP.append(Line(DatiAsse[i][0],DatiAsse[i+1][1]))


	num_seq = []
        num_pdb_s = []
        num_pdb_h = []
        h_radius = []
        hand = []
        straig = []

        for i in range(len(proteina)-1):

                j += 1
                if i <= len(proteina)-3:
			if traj_flag:
				Ca = proteina[i].peptide.C.position()
				Ca1 = proteina[i+1].peptide.C.position()
				Ca2 = proteina[i+2].peptide.C.position()
			else:
				Ca = proteina[i].C.position()
				Ca1 = proteina[i+1].C.position()
				Ca2 = proteina[i+2].C.position()
                        num_seq.append(i+1)
                        num_pdb_h.append(j)
                        h_radius.append(asseIP[i].distanceFrom(Ca))
                        hand.append(sign(DatiAsse[i][3]))
                else:
                        num_seq.append(i+1)
                        num_pdb_h.append(j)
                        h_radius.append(0.0)
                        hand.append(0.0)
                if i <= len(DatiAsse)-3:
                        num_pdb_s.append(j)
                        pscal=(asseIP[i].projectionOf(Ca)-asseIP[i+1].projectionOf(Ca1))*  \
                             (asseIP[i+1].projectionOf(Ca1)-asseIP[i+2].projectionOf(Ca2))
                        mod1=(asseIP[i].projectionOf(Ca)-asseIP[i+1].projectionOf(Ca1)).length()
                        mod2=(asseIP[i+1].projectionOf(Ca1)-asseIP[i+2].projectionOf(Ca2)).length()
                        straig.append(pscal/(mod1*mod2))
                else:
                        num_pdb_s.append(j)
                        straig.append(0.0)

        return num_seq, num_pdb_s, num_pdb_h, h_radius, hand, straig

