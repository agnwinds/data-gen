#!/usr/bin/env python  

'''
                    Southampton University


Synopsis:  
	This code parses data from verner et al (1996) 
	and verner and yakovlev (1995) into python
	input files
	

Description:  

	This code expects two data files to be present in
	the directory. 
	Firstly, it uses the old VFKY data file present 
	in the atomic directory since time immorial.
	To aid input, a title line is added to the top
	
	title z nelec eth emax e0 sigma0 ya P yw y0 y1
	
	This data file should be titled photo_vfky.data
	
	This is supplemented by data downloaded from vizir
	accessed via the adsabs page for verner and yakovlev.
	It should be cut an pasted into a file called 
	photo_partial_vy.data and a title line appended
	
    recno z nelec n l eth e0 sigma0 ya P yw
	
	
	
	

Arguments:  

	none

Returns:

	Two tabulated data files of cross sections. one
	for outer shell ionization which should in principle
	be similar to the file produced by JM, but with 
	actual data past the inner shell threshold rather
	than extrapolations, and a new one which had inner 
	shell ionization cross sections.

Notes:
									   
History:
15aug	nsh	Coded


'''









import sys
import numpy as np
# from pyhdf import SD
from astropy import constants as consts
from astropy import units as u
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
# import astroconst as ac
# import atomic_data_sub as ads
import subprocess

data1=ascii.read('photo_vfky.data')
data2=ascii.read('photo_partial_vy.data')


emax=50000.0  #The maximum energy to be used
npoints=100   #The number of points in the interpolation

# set up a rather complex array which defines the subshells which are added together to make the outer shell cross section
# the numbers are general up to N=18 for all elements, however for heavier elements, it becomes much more complex.
# the vfky database ignores most heavier elements, only dealing with z=20 (Calcium) and z=26 (Iron)


outer_shell=[	[[1,0]],
				[[1,0]],
				[[2,0]],
				[[2,0]],
				[[2,1],[2,0]],
				[[2,1],[2,0]],
				[[2,1],[2,0]],
				[[2,1],[2,0]],
				[[2,1],[2,0]],
				[[2,1],[2,0]],
				[[3,0]],
				[[3,0]],
				[[3,1],[3,0]],
				[[3,1],[3,0]],
				[[3,1],[3,0]],
				[[3,1],[3,0]],
				[[3,1],[3,0]],
				[[3,1],[3,0]]
			]
			
#We make another couple of arrays which have the outer shells for Calcium and Iron.			
#If we want to have more complex elements, we will need to do some research to find out what
#the outer shells actually are. These two are from looking at the data in VFKY along with the text.
				
outer_shell_20=[[[4,0]],[[4,0]]]
outer_shell_26=[[[3,2]],[[3,2]],[[3,2]],[[3,2]],[[3,2]],[[3,2]],[[4,0],[3,2]],[[4,0],[3,2]]]


#This first loop generates outer shell cross sections. It uses the VFKY data up to the threshold for inner 
#sheel ionization, then supplements this data with the outer shell data from VY. It knows which of the 
#shell/subshell data to use to supplement from the arrays above.

cross_secs=[]   #an array of arrays for cross sections for each ion
energies=[]     #the energy relating to each cross section

for i in range(len(data1['z'])):
	e=np.logspace(np.log10(data1['eth'][i]+0.00001),np.log10(emax),npoints)  #Se up an energy array. Start just before threshold to make sure we get the threshold
	csec=np.zeros(len(e))    #Set all the cross sections to zero
	for j in range(len(e)):   #Loop over all energies
		if e[j]>=data1['eth'][i] and e[j]<=data1['emax'][i]:  #If we are between the threshold and max for VFKY, use that
			x=e[j]/data1['e0'][i]-data1['y0'][i]
			y=np.sqrt(x*x+data1['y1'][i]*data1['y1'][i])
			temp=(x-1.0)**2.0+data1['yw'][i]**2.0
			temp=temp*y**(0.5*data1['P'][i]-5.5)
			temp=temp*(1+np.sqrt(y/data1['ya'][i]))**(-1.0*data1['P'][i])
			csec[j]=csec[j]+(data1['sigma0'][i] * temp)
	for j in range(len(data2['z'])): #We now have the cross sections up to the limit of the VFKY data - search for relevant data in the VY data file to append
		if 	data2['z'][j]==data1['z'][i] and data2['nelec'][j]==data1['nelec'][i]: #We have the correct element
			if data1['nelec'][i]<19:    #We can check for the correct outer level using the generic table
				outer_shell_list=outer_shell[data1['nelec'][i]-1]
			elif data1['z'][i]==20: #We have calcuim, so we need to consult the special table for the outer two shells
				outer_shell_list=outer_shell_20[data1['nelec'][i]-19]
			elif data1['z'][i]==26: #We have iron.
				outer_shell_list=outer_shell_26[data1['nelec'][i]-19]
			else:
				print "Panic - we are being asked for outer shells for element ",data1['z'][i]
				exit()
			for ii in range(len(outer_shell_list)):  #We have found the correct list of outer shells, loop over them
				if data2['n'][j]==outer_shell_list[ii][0] and data2['l'][j]==outer_shell_list[ii][1]: #If the n and l for the current VY record, use it
					for k in range(len(e)):   #Loop over energy
						if e[k]>data1['emax'][i]:  #We only want more data if we are past the VFKY maximum energy
							y=e[k]/data2['e0'][j]
							temp=(y-1.0)**2.0+data2['yw'][j]**2.0
							q=5.5+data2['l'][j]-0.5*data2['P'][j]
							temp=temp*y**(-1.0*q)
							temp=temp*(1+np.sqrt(y/data2['ya'][j]))**(-1.0*data2['P'][j])
							csec[k]=csec[k]+data2['sigma0'][j] * temp    #Add the new data
	energies.append(e)   #Append the new array
	cross_secs.append(csec)  #Append the new array
	

#Output the arrays as a python friendly file
out1=open("vfky_outershell_tab.data",'w')
#Write out some preamble
out1.write("# This is outer shell photoionization data from Verner et al (1996) 1996ApJ...465..487V (VFKY)\n")
out1.write("# supplemented by outer shell data from Verner and Yakovlev (1995) 	1995A&AS..109..125V (VY)\n")
out1.write("# The data is from VFKY up to the inner shell edge (if any) then from VY past this point\n")
out1.write("# It is processed via the python script verner_2_py.py from simple text dumps of the online\n")
out1.write("# data from these two papers\n")
out1.write("#\n")
out1.write("# The data is read in using the topbase mechanics, so level and level_label are just placholders\n")
out1.write("Label z state islp(placeholder) level(placeholder) theshold_energy n_points\n")


for i in range(len(data1['z'])):
	istate=data1['z'][i]-data1['nelec'][i]+1
	out1.write('PhotVfkyS %d %d %d %d %5.3e %d\n'%(data1['z'][i],istate,-999,-999,data1['eth'][i],npoints))
	for j in range(npoints):
		out1.write("PhotVfky %7.4e %7.4e\n"%(energies[i][j],cross_secs[i][j]*1e-18))
out1.close()


#We are now going to produce inner shell cross sections
		
cross_secs=[]   #As before, an array of arrays for the cross sections
energies=[]     #An array of arrays for energies 
z=[]            #We need to save the atomic number
nelec=[]        #And the number of electrons
n=[]            #And the inner shell which is being ionized
l=[]		    #and the subshell
		
		
		

	
for i in range(len(data2['z'])):
	e=np.logspace(np.log10(data2['eth'][i]+0.00001),np.log10(emax),npoints)   #Generate the energy array, running from the threshold to the global maximum energy
	csec=np.zeros(len(e))   #And produce an ampty cross section array
	outer=0  #We are about to check to see if we have an inner shell
	if data2['nelec'][i]<19:    #We can check for the correct outer level using the generic table
		outer_shell_list=outer_shell[data2['nelec'][i]-1]
	elif data2['z'][i]==20: #We have calcuim, so we need to consult the special table for the outer two shells
		outer_shell_list=outer_shell_20[data2['nelec'][i]-19]
	elif data2['z'][i]==26: #We have iron.
		outer_shell_list=outer_shell_26[data2['nelec'][i]-19]
	else:
		outer=1 #This flag stops this shell being computed and output.
	if outer==0:  #This could still be an outer shell, we need to check
		for ii in range(len(outer_shell_list)):
			if data2['n'][i]==outer_shell_list[ii][0] and data2['l'][i]==outer_shell_list[ii][1]:
				outer=1   #It is an outer shell, flag this up				
	if outer==0:  #If we got through the whole check without changing the flag, this must be an inner shell of an ion that we understand
		for k in range(len(e)):  #Run over all energies
			y=e[k]/data2['e0'][i]
			temp=(y-1.0)**2.0+data2['yw'][i]**2.0
			q=5.5+data2['l'][i]-0.5*data2['P'][i]
			temp=temp*y**(-1.0*q)
			temp=temp*(1+np.sqrt(y/data2['ya'][i]))**(-1.0*data2['P'][i])
			csec[k]=csec[k]+data2['sigma0'][i] * temp    #This is the cross section, we will now append to all out arrays
		energies.append(e)
		cross_secs.append(csec)
		z.append(data2['z'][i])
		nelec.append(data2['nelec'][i])
		n.append(data2['n'][i])
		l.append(data2['l'][i])
		

#output the arrays.
				
out2=open("vy_innershell_tab.data",'w')

#write a little pramble
out2.write("# This is inner shell data from Verner and Yakovlev (1995) 	1995A&AS..109..125V (VY)\n")
out2.write("# It is processed via the python script verner_2_py.py from simple text dumps of the online\n")
out2.write("# data from these two papers\n")
out2.write("#\n")
out2.write("label z state n_shell l_subshell theshold_energy n_points\n")


for i in range(len(energies)):
	istate=z[i]-nelec[i]+1  #We use the astromonical state
	out2.write('InnerVYS %d %d %d %d %5.3e %d\n'%(z[i],istate,n[i],l[i],energies[i][0],npoints))
	for j in range(npoints):
		out2.write("InnerVY %7.4e %7.4e\n"%(energies[i][j],cross_secs[i][j]*1e-18))
	
out2.close()
	
	

