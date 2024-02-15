#!/usr/bin/env python 

'''
                    Southampton University

Synopsis:  

Program to extract data from a standard Chianti directory tree and put it into files that python can digest 


Command line usage (if any):

	None.  The program is run from within a python shell                     

Description:  

	The program expects to be run from a directory including a directory called Chianti.
	This directory contains subdirectories for all the elements.

Primary routines:

Notes:

	
									   
History:

19/7/12 NSH coding begun
240125 ksl Updated to Python3, and made to run on current version of Chianti

'''

import os
from types import *
from configparser import *
import pickle
from datetime import date
import numpy as np
# from FortranFormat import *
# import chianti.constants as const
# import chianti
import scipy


El = ['h','he','li','be','b','c','n','o','f','ne','na', \
    'mg','al','si','p','s','cl','ar','k','ca','sc','ti', \
    'v','cr','mn','fe','co','ni','cu','zn']





KB=1.380658e-16
H=6.6260755e-27
C=2.99792458e10
ME=9.1093897e-28
MP=1.6726231e-24
SB=5.6704373e-5
HEV=4.13620e-15
PI=3.1415

out=open("chianti_rr.dat","w")

for i in range(len(El)):
    for j in range(i+1):
        ion=El[i]+'_'+str(j+2)
        fname='Chianti/'+El[i]+'/'+ion+'/'+ion+'.rrparams'
        print(fname)
        if os.path.isfile(fname):
            input=open(fname,'r')
            #  need to read first line and see how many elements
            lines=input.readlines()
            input.close()
            rrtype=int(lines[0])
            ref=lines[3:-2]
            #
            if rrtype == 1:
                # a 4 parameter Badnell type
                params=lines[1].split()
                out.write("RR_BADNL ")
                for item in params:
                    out.write(str(item)+" ")
                out.write("\n")
            elif rrtype == 2:
                # a 6 parameter Badnell type
                params=lines[1].split()
                out.write("RR_BADNL ")
                for item in params:
                    out.write(str(item)+" ")
                out.write("\n")
            elif rrtype == 3:
                # a Shull type
                params=lines[1].split()
                out.write("RR_SHULL ")
                for item in params:
                    out.write(str(item)+" ")
                out.write("\n")
            else:
                RrParams=None
                print(' for ion %5s unknown RR type = %5i' %(ion, rrtype))
        else:
            RrParams=None
            print(' for ion %5s no RR file' %(ion))



            inp2=open(fname,'r')
            inp2.close()

out.close()
out=open("chianti_dr.dat","w")

for i in range(len(El)-1):
    print(i)
    for j in range(i+1):
        ion=El[i+1]+'_'+str(j+2)
        fname='Chianti/'+El[i+1]+'/'+ion+'/'+ion+'.drparams'
        print(fname)
        if os.path.isfile(fname):
            input=open(fname,'r')
            #  need to read first line and see how many elements
            lines=input.readlines()
            input.close()
            drtype=int(lines[0])
            ref=lines[4:-1]
            #
            if drtype == 1:
                # a Badnell type
                eparams=lines[1].split()
                cparams=lines[2].split()
                out.write("DR_BADNL E ")
                for item in eparams:
                    out.write(str(item)+" ")
                out.write("\n")
                out.write("DR_BADNL C ")
                for item in cparams:
                    out.write(str(item)+" ")
                out.write("\n")
            elif drtype == 2:
                # shull type
                params=lines[1].split()
                out.write("DR_SHULL ")
                for item in params:
                    out.write(str(item)+" ")
                out.write("\n")
            else:
                DrParams = None
                print(' for ion %5s unknown DR type = %5i' %(ion, drtype))
        else:
            DrParams=None
            print(' for ion %5s no DR file' %(ion))


out.close()

	


