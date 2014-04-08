import numpy
import string

def DEF_print(helix):
    a=list(helix)
    b=map(lambda x: str(int(x)), a)
    st = ''
    for item in b:
        st=st+item
    print st


def DefineSecondaryStructure(results, model=None):
    
    results_sstructure = {}    

    if model:
        num_seq, sen, Madist, Mhrad, Mstraig, mod_ref = results
        indx1=mod_ref.index(model)
        angdist = Madist[indx1]
        radius  = Mhrad[indx1]
        straigt = Mstraig[indx1]
        return SSanalysis(sen,angdist,radius,straigt)
    else:
        results_structure = {}
        for mol in results:
	    # Don't perform analysis on not protein-like molecules
	    if type(results[mol][0]) == float:
		pass
	    else:
		angdist = numpy.array(results[mol][4])
		radius  = numpy.array(results[mol][5])
		straigt = numpy.array(results[mol][6])
		sen = numpy.array(results[mol][1])
		results_structure[mol]=SSanalysis(sen,angdist,radius,straigt)
	return results_structure
    

def SSanalysis(sen,angdist,radius,straigt):

    ang_lowers = [0.446, 0.579, 0.380, 0.721]
    ang_uppers = [0.628, 0.761, 0.562, 0.979]
    rad_lowers = [0.113, 0.091, 0.123, 0.001]
    rad_uppers = [0.223, 0.201, 0.233, 0.081]
    
    if type(sen) is float:
        pass
    else:
        set_vector = numpy.ones((len(sen),),dtype=float)

        helix_Ang_lower = ang_lowers[0]*set_vector
        helix_Ang_upper = ang_uppers[0]*set_vector
        helix_Rad_lower = rad_lowers[0]*set_vector
        helix_Rad_upper = rad_uppers[0]*set_vector

        helixtre_Ang_lower = ang_lowers[1]*set_vector
        helixtre_Ang_upper = ang_uppers[1]*set_vector
        helixtre_Rad_lower = rad_lowers[1]*set_vector
        helixtre_Rad_upper = rad_uppers[1]*set_vector

        helixpi_Ang_lower = ang_lowers[2]*set_vector
        helixpi_Ang_upper = ang_uppers[2]*set_vector
        helixpi_Rad_lower = rad_lowers[2]*set_vector
        helixpi_Rad_upper = rad_uppers[2]*set_vector

        strand_Ang_lower = ang_lowers[3]*set_vector
        strand_Ang_upper = ang_uppers[3]*set_vector
        strand_Rad_lower = rad_lowers[3]*set_vector
        strand_Rad_upper = rad_uppers[3]*set_vector

        strs = numpy.array((0.5 <= straigt)*(straigt <= 1.0), dtype=float)        

        helix_ang = numpy.array((helix_Ang_lower <= angdist)*(angdist <= helix_Ang_upper), dtype=float)
        helix_rad = numpy.array((helix_Rad_lower <= radius)*(radius <= helix_Rad_upper), dtype=float)

        helix = numpy.array( helix_ang*helix_rad, dtype= float) 
        
        helixI_ang = numpy.array((helixpi_Ang_lower <= angdist)*(angdist <= helixpi_Ang_upper), dtype=float)
        helixI_rad = numpy.array((helixpi_Rad_lower <= radius)*(radius <= helixpi_Rad_upper), dtype=float)

        helixI = numpy.array( helixI_ang*helixI_rad, dtype= float) 

        helixG_ang = numpy.array((helixtre_Ang_lower <= angdist)*(angdist <= helixtre_Ang_upper), dtype=float)
        helixG_rad = numpy.array((helixtre_Rad_lower <= radius)*(radius <= helixtre_Rad_upper), dtype=float)

        helixG = numpy.array( helixG_ang*helixG_rad, dtype= float) 
        
        strand_ang = numpy.array((strand_Ang_lower <= angdist)*(angdist <= strand_Ang_upper), dtype=float)
        strand_rad = numpy.array((strand_Rad_lower<=radius)*(radius<=strand_Rad_upper), dtype=float)

        strand = numpy.array( strand_ang * strand_rad, dtype= float) 
        
        SS_structure = []
        for j in range(len(set_vector)):
            SS_structure.append('-')

# ASSIGN ALPHA-HELICES AND ALPHA-TURNS
        for j in range(1,len(set_vector)-2):
            if helix[j-1] == 1. and  helix[j] == 1. and  helix[j+1] == 1. and  helix[j+2] == 1.:
                SS_structure[j-1] = 'H'
                SS_structure[j] = 'H'
                SS_structure[j+1] = 'H'
                SS_structure[j+2] = 'H'
            elif SS_structure[j-1] == 'H' and helix[j] == 1.:
                SS_structure[j] = 'H'
            elif helix[j-1] != 1. and helix[j] == 1. and helix[j+1] != 1.:
                SS_structure[j] = 'T'
            else:
                pass
# ASSIGN PI-HELICES
        for j in range(1,len(set_vector)-3):
            if SS_structure[j] == '-'  and SS_structure[j-1] != 'H':
                if helixI[j-1] == 1. and helixI[j] == 1. and helixI[j+1] == 1. and helixI[j+2] == 1. and helixI[j+3] == 1.:
                    SS_structure[j-1] = 'I'
                    SS_structure[j] = 'I'
                    SS_structure[j+1] = 'I'
                    SS_structure[j+2] = 'I'
                    SS_structure[j+3] = 'I'
                else:
                    pass

#ASSIGN 3-10 HELICES
        for j in range(1,len(set_vector)-1):
            if SS_structure[j] == '-' and SS_structure[j-1] != 'H':
                if helixG[j-1] == 1. and helixG[j] == 1. and helixG[j+1] == 1.:
                    SS_structure[j-1] = 'G'
                    SS_structure[j] = 'G'
                    SS_structure[j+1] = 'G'
                else:
                    pass

# ASSIGN BETA STRANDS
        for j in range(1,len(set_vector)-1):
            if SS_structure[j] == '-':
                if strand[j-1] == 1. and strand[j] == 1.:
                    SS_structure[j-1] = 'E'
                    SS_structure[j] = 'E'
                elif strand[j-1] != 1. and strand[j] == 1. and strand[j+1] != 1.: 
                    SS_structure[j] = 'B'
                else:
                    #                    SS_structure[j] = '-'
                    pass

    
    return string.join(SS_structure,'')+"-"

