########################################################################
#
#	University of Southampton, James Matthews, 130625
#
#	"extrapolate_sub.py"
#		
#	subroutines for extrapolate.py. individually described
#
#########################################################################

## import modules
from math import sqrt, fabs
import numpy as np



############################################################################
## CHOOSEMODE FUNCTION



def choose_mode(arg_array):

	Plot = False
	fmode=1
	MAX_ENERGY=10000

	if len(arg_array)>1:
		for i in range( len(arg_array) ):
			if arg_array[i] == '-e':
				MAX_ENERGY = float ( arg_array[ i + 1 ] )
			if arg_array[i] == '-f':
				fmode= int ( arg_array[ i + 1 ] )
			if arg_array[i] == 'plot':
				Plot = True
	else:
		fmode=1
		MAX_ENERGY=10000

	return fmode, MAX_ENERGY, Plot


############################################################################
## MINIMISATION FUNCTION

def gmin(f,a,c,tol=3.0e-8):
    """Golden section minimisation.

    Find minimum of f(x) known to have a single minimum
    between a and c, using golden section search. Returns
    minimum position, xmin, and minimum value, f(xmin).
    Suggest set tol to 1.0e-4 for single precision, 3.0e-8
    for double precision."""

    if (c<=a):
        raise ValueError, 'gmin(f,a,c) requires c>a'

    w=(3.0-sqrt(5.0))/2.0 # golden ratio

    b=a+w*(c-a) # set initial abscissas in interval

    bp=c-w*(c-a)

    while(fabs(c-a)>tol*(c+a)/2.0): # loop to refine
        if(f(b)<f(bp)):            # bracketing of minimum
            c=bp
            bp=b
            b=a+w*(c-a)
        else:
            a=b
            b=bp
            bp=c-w*(c-a)

    # finished looping, return best answer
    fb=f(b)
    fbp=f(bp)

    if (fb<fbp): return b, fb
    else: return bp, fbp


############################################################################

def read_topbase(filename):
	'''
	read in XS info from Topbase XS data in Python format

	:INPUT:
		filename 		string
						atomic data filename e.g. topbase_h1_phot.py
	:OUTPUT:
		top 			topbase class instance
						topbase class instance containing information 
						for this filename
	'''

	# read in summary records
	Z,ion,islp,l, E0, num_records = sum_records = np.loadtxt( filename, 
			dtype={'names': ('Z', 'ion', 'islp', 'l', 'E0', 'np'), 
			'formats': ('i4', 'i4', 'i4', 'i4', 'float', 'i4')}, 
                        comments='PhotTop ', delimiter=None, converters=None, 
                        skiprows=0, usecols=(1,2,3,4,5,6), unpack=True, ndmin=0)
	
	# then read the actual cross sections
	energy, XS = np.loadtxt(filename, dtype='float', 
                        comments='PhotTopS', delimiter=None, converters=None, 
                        skiprows=0, usecols=(1,2), unpack=True, ndmin=0)

	# create blank array to store topbase class instances
	top = np.ndarray( len(Z),dtype=np.object)

	nline = 0
	
	for i in range(len(top)):

		n_p = int(num_records[i])
		nmax = nline + n_p
		
		top[i] = topbase_class (Z[i], int(ion[i]), int(islp[i]), int(l[i]), E0[i], n_p, energy[nline:nmax], XS[nline:nmax])

		nline = nmax
		
	# top is an array of class instances like the topbase_ptr in PYTHONRT
	return top


def write_topbase(top, filename):

	'''write array of topbase class instances to file'''

	import numpy as np

	file_write = open( filename, 'w')

	for i in range(len(top)):

		## write the summary records
		file_write.write('PhotTopS  %1d  %1d %3d    %1d     %.6f  %2d\n' %
                                    ( top[i].Z, top[i].ion, top[i].islp, top[i].l, top[i].E0, top[i].np ))

		n_p = top[i].np

		## write the actual XSs
		for j in range( n_p ):
			
			file_write.write('PhotTop     %.6f %.3e\n' % ( top[i].energy[j], top[i].XS[j]) )


	file_write.close()
		
	return 0


class topbase_class:
	'''
	This is a class for topbase photoionization data

	Z 		atomic number
	ion 	ion stage
	islp 	islp number (2s+1, l, parity)
	E0 		threshold energy eV 
	l 		level
	np 		number of entries
	energy 	energies 
	XS 		cross sections
	'''	

	def __init__(self, nz, ne, islp_init, linit, E0_init, np_init, energies, cross_sections):
		self.Z = nz
		self.ion = ne
		self.islp = islp_init
		self.l = linit
		self.E0 = E0_init 
		self.np = np_init
		self.energy = energies
		self.XS = cross_sections



def check_odd_XS(XS_array, top_record):
	'''
	Routine which checks if the top_record class in question
	is in the flagged list

	:INPUT:
		XS_array 		topbase class instance
						lists of unusual XSections

		top_record 		topbase class instance
						record to compare

	:OUTPUT:
		totalmatch		Bool
						1/True if in list, 0/False if not
	'''

	import numpy as np

	# create boolean arrays for Z, islp and ion
	Zmatch = ( XS_array.Z == top_record.Z )
	ionmatch = ( XS_array.ion == top_record.ion )
	islpmatch = ( XS_array.islp == top_record.islp )
	levmatch = ( XS_array.l == top_record.l )

	# create total boolean array (note * means 'and')
	totalmatch = np.sum(ionmatch * islpmatch * Zmatch * levmatch)


	# error check
	if totalmatch != 1 and totalmatch!=0:
		print "ERROR: Does not equal 1 or zero, exiting"
		print totalmatch
		sys.exit()

	return totalmatch




def get_odd_XS():
	'''
	get_odd_XS is a real kluge, it just returns an array of XS identified by eye
	as being odd. These XS are then fudged with a nu**-3 extrapolation.

	:OUTPUT:

		XS_array 		topbase class instance
						array of XS to fudge in topbase_class format

	'''



	# Manually identified list of unusual / anomalous XSections
	XSlist = ['6 2 220 3', \
	'6 2 220 4', \
	'6 2 231 4', \
	'6 2 231 5', \
	'6 2 231 6', \
	'6 5 111 6', \
	'6 5 311 6', \
	'6 5 311 7', \
	'7 5 200 6', \
	'7 6 111 6', \
	'7 6 311 7', \
	'8 3 141 1', \
	'8 3 141 2', \
	'8 3 141 3', \
	'8 3 141 4', \
	'8 3 141 5', \
	'8 3 150 1', \
	'8 3 150 2', \
	'8 3 150 3', \
	'8 3 151 1', \
	'8 3 151 3', \
	'8 3 151 4', \
	'8 3 341 1', \
	'8 3 341 2', \
	'8 3 341 3', \
	'8 3 341 4', \
	'8 3 341 5', \
	'8 3 350 1', \
	'8 3 350 2', \
	'8 3 350 3', \
	'8 3 351 1', \
	'8 3 351 2', \
	'8 3 351 3', \
	'8 3 351 4', \
	'8 3 351 5', \
	'8 3 351 6', \
	'8 3 351 7', \
	'8 3 351 8', \
	'8 4 200 9', \
	'8 4 211 14', \
	'8 4 220 11', \
	'8 4 231 1', \
	'8 4 231 7', \
	'8 4 411 4', \
	'8 4 421 2', \
	'8 6 200 8']

	import numpy as np

	# number of odd XSections
	n_odd = len(XSlist)

	# create empty object array to store topbase class instances
	XS_array = np.ndarray( n_odd,dtype=np.object)

	Z, ion, islp, l = [],[],[], []

	# cycle over each XS and place in class instance
	for i in range(n_odd):


		# split string
		data = XSlist[i].split()

		
		# get information
		Z.append (int( data[0]))
		ion.append( int(data[1]))
		islp.append(int(data[2]))
		l.append(int(data[3]))

	
		# get topbase class instance
	XS_array = topbase_class (np.array(Z), np.array(ion), np.array(islp), np.array(l), 0, 0, [], [])


	# all done, return array
	return XS_array








