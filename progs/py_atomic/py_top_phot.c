
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis: parse_phot

  Description:	This program parses the results of topbase regarding
	photoionization crossections and reduces the total number of
	xsections in the list to a fixed value.  Boxcar smoothing of
	the original xsections are carried out so as to smooth over
	some of the top base resonances.  



  Arguments:  The program uses rdpar to get parameters

  rdstr ("Root.for.all.files", root);      Rootname for all of the input files.  The program assumes
					   that the input data file will have the extension .dat
						and therefore the output file will have the extension
						.py
  rdint ("No.of.output.points", &nout);    number of output x-sections for each input x section
  rdflo ("Fractional.smoothing", &delta);  The fractional size, e.g. 0.5 for the boxcar smooth

  Returns:

  Notes:

 	The program suppresses xsections that are below threshold, but does not at least
	at present adjust the thresholds themselves to experimental values.  This is not
	done here on the argument that if one is going to use these cross sections then one will
	also have a levels file for this data (in Python) and that it would be easier to fix
	the levels file and then do something in get_atomicdata


  History:
	01oct08 ksl	Began work.  Adapted from parse_levels
	01nov	ksl	Modified printout so more like python normally
			expects, e.g. ion number instead of number of
			electrons, and eV for energy.  Added boxcar 
			smoothing


 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include 	<strings.h>
#define		LINELENGTH 132

#define RYD     13.605698	//Rydberg in eV

#define NCROSS 10000		// The maximum number of cross sections in the input data files
int ntop_phot;

struct topbase_phot
{
  int nlev;			/* Internal index to the config structure for this x-section */
  int nion;			/* Internal index to the    ion structure for this x-section */
  int z, istate;
  int isp, ilv;
  int np;
  double ethresh;
  double freq[NCROSS], x[NCROSS];
  double f, sigma;		/*last freq, last x-section */
}
photin, photout;



main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *optr, *fopen ();
  char infile[80], outfile[80];
  char line[LINELENGTH];
  char firstword[LINELENGTH];
  int nz, ne, isp, ilv;
  float e, x;
  int np;
  int npp;
  char root[80];
  int nout;
  float delta;


// Initialize3 inputs

  strcpy (root, "topbase_phot_hhe");
  nout = 50;
  delta = 0.05;

// Get inputs
  rdstr ("Root.for.all.files", root);
  rdint ("No.of.output.points", &nout);
  rdflo ("Fractional.smoothing", &delta);

// Complete the file names

  strcpy (infile, root);
  strcat (infile, ".dat");
  strcpy (outfile, root);
  strcat (outfile, ".py");




/* Open the input file.  Exit if it is not opened successfully */
  if ((fptr = fopen (infile, "r")) == NULL)
    {
      printf ("Failed to open file %s\n", infile);
      exit (0);
    }

/* Open the output file */
  if ((optr = fopen (outfile, "w")) == NULL)
    {
      printf ("Failed to open file %s\n", outfile);
      exit (0);
    }



/* Read and print the input file */

  while (fgets (line, LINELENGTH, fptr) != NULL)
    {
      sscanf (line, "%s", firstword);
      if (firstword[0] == '=' || firstword[0] == 'I')
	{
	  printf ("#%s", line);
	}
      else
	if (sscanf
	    (line, "%*d %d %d %d %d %e %d", &nz, &ne, &isp, &ilv, &e,
	     &np) == 6)
	{
	  photin.z = nz;
	  photin.istate = nz - ne + 1;
	  photin.isp = isp;
	  photin.ilv = ilv;
	  photin.np = np;
	  photin.ethresh = -e * RYD;
	  printf ("PhotTopS %2d %2d %2d %4d %12.5e %3d\n", nz, nz - ne + 1,
		  isp, ilv, e * RYD, nout);
	  npp = 0;
	  if (np > NCROSS)
	    {
	      Error ("Too many input cross sections for this level %d>%d\n",
		     np, NCROSS);
	      exit (0);
	    }
	}
      else
	{
	  sscanf (line, "%e %e", &e, &x);
	  printf ("PhotTop %12.6e %8.3e\n", e, x);
	  photin.freq[npp] = e * RYD;
	  photin.x[npp] = x * 1e-18;
	  npp++;
	  if (npp == np)
	    {
	      resample_phot (nout, delta);
	      print_phot (optr);
	    }
	}
    }
}

/* resample_phot both resamples and smooths the function slightly.  This is
to account for resonances in the top_base database  */

#define MDELTA  11		// Number of points in boxcar smoothing

int
resample_phot (nout, delta)
     int nout;
     float delta;
{
  float emin, emax;
  int n, m;
  float frac, value, original_x ();
  float freq;
  float del;
  int mdelta;

  photout.ethresh = emin = photin.ethresh;
  emax = photin.freq[photin.np - 1];
  photout.z = photin.z;
  photout.istate = photin.istate;
  photout.isp = photin.isp;
  photout.ilv = photin.ilv;
  photout.np = nout;


  for (n = 0; n < nout; n++)
    {
      photout.freq[n] = emin + (emax - emin) / (nout - 1) * n;
    }

  for (n = 0; n < nout; n++)
    {
      value = 0;
      del = 2. * (1. - photout.ethresh / photout.freq[n]);
      if (del > delta)
	del = delta;
      for (m = -MDELTA / 2; m <= MDELTA / 2; m++)
	{
	  freq = photout.freq[n] * (1. + m * del / MDELTA);
	  value += original_x (freq);
	}
      photout.x[n] = value / MDELTA;
    }
}

/* Get the value of the x-section as it was read in by interpolation */

float
original_x (freq)
     float freq;
{
  float x, frac;
  int m;

  m = 0;
  while (photin.freq[m] < freq && m < photin.np - 1)
    m++;
  if (m == 0)
    {
      x = photin.x[0];
    }
  else
    {
      frac =
	(freq - photin.freq[m - 1]) / (photin.freq[m] - photin.freq[m - 1]);
      x = frac * photin.x[m] + (1. - frac) * photin.x[m - 1];
    }

  return x;

}

int
print_phot (optr)
     FILE *optr;
{
  int n;
  fprintf (optr, "PhotTopS %2d %2d %2d %4d %12.6f %3d\n",
	   photout.z, photout.istate, photout.isp, photout.ilv,
	   photout.ethresh, photout.np);
  for (n = 0; n < photout.np; n++)
    {
      fprintf (optr, "PhotTop %12.6f %8.3e\n", photout.freq[n], photout.x[n]);
    }
}
