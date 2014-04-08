import SFit
from sys import argv
from copy import deepcopy
from getopt import *
from time import sleep,asctime,localtime,time
import string
import datetime
import os
import numpy as N



def header(PDBfile):
        inputfile = string.split(os.path.basename(PDBfile),".")[0]
        time_string = datetime.datetime.now().strftime("%y-%m-%d_%Hh%M")
	filename = string.split(inputfile,'.')[0]+"_"+time_string+".log"
        a= file(filename, "w")

	start_time = time()

	header = """\n
        =================== 
	ScrewFit version """+str(SFit.__version__)+""" 
	===================
	Started at: """+asctime(localtime(start_time))+"""
	\n
	Processing input file....

        Logfile : """+filename+""" 
	--- \n """

	print header
	a.write(header)
	a.close()
        return filename

def dump_logfile(logfile, not_avail):
	a= file(logfile, "a")
	for i in not_avail:
 	  a.write("ScrewFit analysis not available for molecule %s \n" %i)
	a.close()

def tail(logfile, PDBfile,results,nmr_model=None):
        
	a = file(logfile, "a")
        infile=string.split(PDBfile,'.')[0]
        end_time = time()
        
        if nmr_model:
                outputModel(infile,results,nmr_model)
        else:
                strs = output(infile,results)
              
        end = """
        Done {"""+asctime(localtime(end_time))+"""}.
        \n"""

	
	for item in strs:
            print item
	    a.write("%s \n" %item)
        print end 
	a.write(end)
	a.close()


def printHelp():
        print 	" \n\n"
	print "ScrewFit v.%s" % SFit.__version__ 
        print "\n"
	print "Usage:     screwfit filename"	
        print "Input file must be a python script or a single crystallographic PDBfile.\n\n"


def openFileNameFromPDBfile(flag, PDBfile, molname = None,nmr_model = None):

	if flag == 'control':
		prefix='AngDist_Control_'
		extension = '.log'
	elif flag == 'SMdata':
		prefix='SM_data_'
		extension = '.dat'
	elif flag == 'Assign':
		prefix='SecStr_'
		extension = '.fasta'
	else:
		prefix='ScrewFit'
		extension = '.dat'

	
	current = os.getcwd()
	basename = string.split(os.path.basename(PDBfile),".")[0]
	if nmr_model:
		filename = os.path.join(current, prefix+basename+"-model_"+str(nmr_model)+extension)
	else:
		if molname:
			filename = os.path.join(current, prefix+basename+"-"+molname+extension)
		else:
			filename = os.path.join(current, prefix+"_"+basename+extension)
			
	return open(filename,'w')

	
def output(input_file,results):
	
	#     num_seq, sen, hand, err, ang_dist, h_radius, straig
	stringhe = []
	for mol in results:
		params = results[mol]
		molname = string.replace(str(mol),' ','_')
		if molname[:5] == 'Chain':
			stringhe.append('             Screw Motion Parameters output file:     SM_data_%s-%s.dat\n' \
			      % (input_file,molname))
			f = openFileNameFromPDBfile('SMdata',input_file, molname)
			for i in range(len(params[0])):
				f.write("%s %d %d %s %s %s %s\n" \
				        %(params[0][i],params[1][i],params[2][i],\
				          params[3][i],params[4][i],params[5][i],params[6][i]))
			f.close()		
	return stringhe
		

def RulerSecStr(length):
	ten = "....:....|"
	resto = 0
	lrint = int((length)/10.)
	resto = resto + (length - lrint*10)
	ruler= lrint*ten
	if 0 < resto < 5:
		ruler=ruler+resto*"."
	elif resto == 5:
		ruler=ruler+"....:"
	elif 5 < resto < 10:
		ruler=ruler+"....:"+(resto-5)*"."
	else:
		pass
	
	return ruler

def ThreeToOne(data):
    transl={'ALA':'A', 'ARG':'R', \
            'ASN':'N', 'ASP':'D', \
            'CYS':'C', 'CYX':'C', \
            'GLU':'E', 'GLN':'Q', \
            'GLY':'G', 'HIS':'H', \
            'ILE':'I', 'LEU':'L', \
            'LYS':'K', 'MET':'M', \
            'PHE':'F', 'PRO':'P', \
            'SER':'S', 'THR':'T', \
            'TRP':'W', 'TYR':'Y', \
            'VAL':'V', 'HSD':'H'}

    sequences={}
    for structure in data.pdb:
	try:
	    structure.sequence
	    new_seq=[]
	    for item in structure.sequence:
		new_seq.append(transl[item[:3]])
	    new_string=string.join(new_seq,'')
	    new_string=new_string.replace(' ','')
	    sequences[str(structure)]=new_string
	except:
	    pass	    
    return sequences
    
def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def FASTAout(data, PDBfile,seq):	
	sequences=ThreeToOne(data)
	for item in seq.keys():
	  new_item = item.replace('_',' ')
          try:
              le = len(seq[item])
              ln = len(item)
              title = "> "+PDBfile+" | "+str(item)
              a=openFileNameFromPDBfile('Assign',PDBfile, item)
              a.write(title+"\n")
              length2=split_len(sequences[new_item],100)
              for i in range(len(length2)):
                  a.write(length2[i]+"\n")
              a.write("\n")
              a.write("\n")
              title = "> Secondary Structure | "+str(item)
              a.write(title+"\n")
              length1=split_len(seq[item],100)
              for i in range(len(length1)):
                  a.write(length1[i]+"\n")
              a.close()
          except:
              print "pass", item

def Sequence_out(data, PDBfile,seq):	
	sequences=ThreeToOne(data)
	for item in seq.keys():
	    new_item = item.replace('_',' ')
	    try:
              le = len(seq[item])
              ln = len(item)
              title = "> "+PDBfile+" | "+str(item)
              print title
              length1=split_len(seq[item],100)
              length2=split_len(sequences[new_item],100)
              length3=split_len(RulerSecStr(le),100)
              for i in range(len(length1)):
                  print length1[i]
                  print length2[i]
                  print length3[i]
                  print " "+"\n"
	    except:
              print "pass", item
