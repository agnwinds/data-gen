2184 Fri Nov 16 17:25:19 EST 2001 muskie.stsci.edu
This directory contains awk scripts required to create atomic data for python from the 
data retrieved from Dina Verner's web site.  Some of the programs, notable py_read_verner
and py_link associated with this are in the directory py_atomic.  

Normally the procedure is to run the awkscript Awk.verner.  This creates from
abun.dat and ions.dat the file elem_ions_ver.py.  As written the file is
commented so that elements and ions with abundances > 5 are marked for use

Then one would run py_read_verner and py_link on atom1.dat to produce a line list.

170211 - ksl copied into data-gen.
