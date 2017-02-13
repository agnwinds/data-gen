#! /Library/Frameworks/Python.framework/Versions/2.7/Resources/Python.app/Contents/MacOS/Python
'''
University of Southampton, James Matthews, 130625

	"extrapolate_log.py"
		
This code is intended to extrapolate topbase to higher energies

The basic procedure for this code is as follows:
1. list the filenames you want to create
2. work out a fit to the topbase date
3. extrapolate this fit to a higher energy

usage:
	python extrapolate_log.py [-e max_energy] [-f function mode]
returns:
	topbase files in form topbase_h1_phot_extrap.py
'''
	
## import modules
import os, sys
import numpy as np
import pylab
import extrapolate_sub as sub


## we read the functional form and the maximum energy from the command line
## they are set to defaults if not read
func, NEW_MAX_ENERGY, Plot= sub.choose_mode(sys.argv)
legend = True


## open the standard73 atomic data file in order to read topbase filenames
atomic_filename = 'data/standard73'

## read standard73 into array
atomic_data_files = np.loadtxt(atomic_filename, dtype='string', comments='#', delimiter=None, converters=None, skiprows=0, usecols=None, unpack=False, ndmin=0)


## write topbase filenames to an array
print '''Your filenames are:
	Old       |   	    New
------------------------------------------------'''
topbase_data_files, new_data_files=[],[]
for i in range(len(atomic_data_files)):
	oldfile=atomic_data_files[i]

	## is it a topbase photoionization record?
	if 'topbase' in oldfile and 'phot' in oldfile:

		## if so, create new filename
		newfile=atomic_data_files[i][14:-3]+'_extrap.py'
		topbase_data_files.append(oldfile)
		new_data_files.append(newfile)
		print oldfile[14:],' | ',newfile

odd_XS = sub.get_odd_XS()


## we now have 2 arrays of old and new filenames
## let's do the fit

n_files=len(topbase_data_files)

## quick check
if n_files!=len(new_data_files): 
	print 'Error: differing array lengths'
	sys.exit(0)


## define our new maximum energy in eV
NEW_MAX_ENERGY=100000.0	#use 100KeV for the moment


for i_file in range(n_files):

	## do each filename in turn
	filename_read = topbase_data_files[i_file]

	filename_write = new_data_files[i_file]

	sys.stderr.write('Reading '+filename_read+'...')

	top = sub.read_topbase ( filename_read )
	topnew = sub.read_topbase ( filename_read )

	sys.stderr.write('done.\n')
	sys.stderr.write('Extrapolating data...')


	## Now do the fit

	nline=0	## counter for records

	
	## we now iterate 
	
	for i in range(len(top)):

		n_p = top[i].np

		

		Emin, Emax = top[i].energy[0], top[i].energy[-1]
		XSmax = top[i].XS[-1]

		#coefficient = XSmax / (Emax ** -3) 
		if NEW_MAX_ENERGY < Emax: 
			new_Emax = Emax 
			nmore=0

		else:
			new_Emax = NEW_MAX_ENERGY
			nmore=50	## the number of entries we extrapolate by

		topnew[i].np = n_p + nmore



		try: 
			dlogE = np.log10( (new_Emax - Emax) ) / nmore 
			dE = (new_Emax - Emax) / nmore
		except RuntimeWarning:
			print '0 divide in Logs, sometimes happens with Fe'
			dlogE, dE = 0,0



		ldx = np.log10 ( top[i].energy[-1] ) - np.log10 ( top[i].energy[-2] )
		ldy = np.log10 ( top[i].XS[-1] ) - np.log10 ( top[i].XS[-2] )


		# gradient in log space
		grad = ldy / ldx




		fudge = sub.check_odd_XS ( odd_XS, top[i] )


		if fudge and grad > -2.5:
			print "odd XS, fudging!, grad = ", grad

			ldx = np.log10 ( top[i].energy[-1] ) - np.log10 ( top[i].energy[-5] )
			ldy = np.log10 ( top[i].XS[-1] ) - np.log10 ( top[i].XS[-5] )


			# gradient in log space
			grad = ldy / ldx

			if grad > -2.5:
				grad = -3


		logXSmax = np.log10 ( XSmax ) 
		logEmax = np.log10 ( Emax ) 
		


		## now extrapolate using gradient in log log space at last 2 points in topbase entry
		for i_energy in range(1, nmore+1):

			E = Emax + 10.0 ** ( i_energy * dlogE) 

			logXSfit =  logXSmax + ( grad * ( np.log10 (E) - logEmax ) )
			XSfit = 10.0 ** logXSfit
			
			## append the values to temporary arrays for this topbase record
			topnew[i].XS = np.append ( topnew[i].XS, XSfit )
			topnew[i].energy = np.append ( topnew[i].energy, E )

		nfive = ( top[i].Z == "7" and top[i].ion == "5"  )
		ofour = ( top[i].Z == "8" and top[i].ion == "4"  )

		## plot up to check fits
		if Plot:

			'''First_time_round = i_file == 0 and i == 0

			if First_time_round == False and i_file == 3 and fudge_last:

				if islp_last != top[i].islp or ion_last != top[i].ion or Z_last != top[i].Z or l_last != top[i].l:

					print 'newsave Dropbox/XS_77/figures_%s/XS_%i_%i_%i_%i.png' % ( filename_last, Z_last, ion_last, islp_last, l_last)
					pylab.xlabel('Energy eV')
					pylab.ylabel('XS cm^-2')
					pylab.savefig('/Users/jmatthews/Dropbox/XS_77/figures_%s/XS_%i_%i_%i.png' %
				                  ( filename_last, Z_last, ion_last, islp_last))

					#if fudge_last: pylab.show()

					pylab.clf()
					legend = True'''

			if fudge:
				pylab.loglog(top[i].energy, top[i].XS, 'r', label = 'standard73')
				pylab.loglog(topnew[i].energy, topnew[i].XS, 'k--', label = 'standard77')
				print 'newsave Dropbox/XS_77/figures_%s/XS_%i_%i_%i_%i.png' % ( filename_last, Z_last, ion_last, islp_last, l_last)
				pylab.xlabel('Energy eV')
				pylab.ylabel('XS cm^-2')
				pylab.savefig('/Users/jmatthews/Dropbox/XS_77/figures_%s/XS_%i_%i_%i_%i.png' %
				                  ( filename_read[14:-3], topnew[i].Z, topnew[i].ion, topnew[i].islp, top[i].l))
				pylab.clf()

			#if fudge: pylab.show()
			
			if fudge: print 'summary', top[i].Z, top[i].ion, top[i].islp, top[i].l, top[i].E0

			if legend and fudge:
				pylab.legend()
				legend = False

		

		islp_last = top[i].islp
		ion_last = top[i].ion
		Z_last = top[i].Z
		filename_last = filename_read[14:-3]
		l_last = top[i].l

		fudge_last = fudge
	




	## Finally, write the files- no numpy modules for flexibility

	sys.stderr.write ('Writing to file...')

	sub.write_topbase ( top, filename_write)

	sys.stderr.write ('done.\n')


## now everything is done!
		

	
		
		






