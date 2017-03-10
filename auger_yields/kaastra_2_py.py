#!/usr/bin/env python  -i

'''
                    Southampton University


Synopsis:  
	This code parses electron and photon yield data
	from Kaastra and Mewe into new data files which
	link it to inner shell PI rates from Verner and
	Vakovlev.
	

Description:  

	This code expects three data files to be present in
	the directory. 
	Firstly, it uses the inner shell photoionization 
	rate data file vy_innershell_tab.data which is in 
	turn parsed from VFKY and VY filed.
	
	
	This is supplemented by data downloaded from vizir
	accessed via the adsabs page for kaastra and mewe.
	The data for electron yield should be in a file called
	electron_yield.data and the photon yield in
	fluorescent_yield.data. 
	
	The code matches up the n/l code used in K+M to the 
	n and l level information in VFKY. 
	
	

Arguments:  

	none

Returns:

	Two tabulated data files of electron and photon 
	yields. Each contains a cross reference to the n
	and l numbers in the relevant VY innershall record.

Notes:
									   
History:
15aug	nsh	Coded


'''









import sys
import numpy as np
from pyhdf import SD
from astropy import constants as consts
from astropy import units as u
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import astroconst as ac
import atomic_data_sub as ads
import subprocess

data_yield=ascii.read('electron_yield.data')
data_fluor=ascii.read('fluorescent_yield.data')

prob=(data_yield['prej1'].data,data_yield['prej2'].data,data_yield['prej3'].data,data_yield['prej4'].data,
     data_yield['prej5'].data,data_yield['prej6'].data,data_yield['prej7'].data,data_yield['prej8'].data,
	 data_yield['prej9'].data,data_yield['prej10'].data)

prob=np.array(prob)




#The shell and subshells in kaastra and mewe are coded. The array below converts their 's' parameter
#into n and l to match the inner shell cross sections. 

levels=[[1,0],[2,0],[2,1],[2,1],[3,0],[3,1],[3,1]]



input=open('vy_innershell_tab.data','r')

z=[]
state=[]
n=[]
l=[]

for line in input.readlines():
	data=line.split()
	if data[0]=='InnerVYS':
		z.append(int(data[1]))
		state.append(int(data[2]))
		n.append(int(data[3]))
		l.append(int(data[4]))
	

#We now have a list of inner shell transitions that we require yields for. Firstly, we will find which
#records match

yield_records=[]
fluor_records=[]

for i in range(len(z)):
	yr_temp=[]
	yf_temp=[]
	for j in range(len(data_yield['z'])):
		n_temp=levels[data_yield['s'][j]-1][0]   #Convert the state into n and l to match inner shell data
		l_temp=levels[data_yield['s'][j]-1][1]
		if data_yield['z'][j]==z[i] and data_yield['state'][j]==state[i] and n_temp==n[i] and l_temp==l[i]:
			yr_temp.append(j)
	yield_records.append(yr_temp)  #Append the records that relate to this inner shell ejection
	for j in range(len(data_fluor['z'])):
		n_temp=levels[data_fluor['s'][j]-1][0]
		l_temp=levels[data_fluor['s'][j]-1][1]

		if data_fluor['z'][j]==z[i] and data_fluor['state'][j]==state[i] and n_temp==n[i] and l_temp==l[i]:
			yf_temp.append(j)
	fluor_records.append(yf_temp)
	
	
#And now for fluorescent yields
	
	


z_out=[]
state_out=[]
n_out=[]
l_out=[]
I_out=[]
EA_out=[]
prob_out=[]

print "outputting"

out=open("kaastra_electron_yield.data",'w')
out1=open("kaastra_fluorescent_yield.data",'w')


out.write("#This is electron yield data from Kaastra & Mewe (1993) -1993A&AS...97..443K\n")
out1.write("#This is fluorescent yield data from Kaastra & Mewe (1993) -1993A&AS...97..443K\n")
out.write("#It is processed from data downloaded from http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+AS/97/443\n")
out1.write("#It is processed from data downloaded from http://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+AS/97/443\n")
out.write("#Via the python script kaastra_2_py.py\n")
out1.write("#Via the python script kaastra_2_py.py\n")
out.write("#\n")
out1.write("#\n")
out.write("#Label z state n l IP mean_electron_energy Prob_of_1e Prob_of_2e Prob_of_3e ....\n")
out1.write("#Label z state n l photon_energy yield\n")



for i in range(len(z)):
	I_temp=0.0
	EA_temp=0.0
	prob_temp=np.zeros(10)
	if len(yield_records[i])>0:
		print "electrons",i,len(yield_records[i])
		for j in range(len(yield_records[i])):
			I_temp=I_temp+data_yield['I'][yield_records[i][j]]
			EA_temp=EA_temp+data_yield['EA'][yield_records[i][j]]
			prob_temp=prob_temp+prob[:,yield_records[i][j]]
		out.write("Kelecyield %d %d %d %d %6.4e %5.3e"%(z[i],state[i],n[i],l[i],
						I_temp/len(yield_records[i]),EA_temp/len(yield_records[i])))
		for j in range(10):
			out.write(" %5.3e"%(prob_temp[j]/len(yield_records[i])))
		out.write("\n")
			
		if len(fluor_records[i])>0:
			print "photons",i,len(fluor_records[i])
			for j in range(len(fluor_records[i])):
				out1.write("Kphotyield %d %d %d %d %5.3e %5.3e\n"%(z[i],state[i],n[i],l[i],
				data_fluor['E'][fluor_records[i][j]],data_fluor['omega'][fluor_records[i][j]]))
		
out.close()
out1.close()
			
			
			
			
			
