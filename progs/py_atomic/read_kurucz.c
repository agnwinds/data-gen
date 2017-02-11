/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

This program reads the Kurucz linelist, and creates
files that can be read by get_atomic_data.  

  Description:	

read_kuruc  creates two files, levels.out, which contains the energy levels for all 
ions of interest, and, lines.out, which contains the lines associated 
with each ion.  

In python_32, no use is made of levels.out, but instead I generally edit
lines.out to produce a limited set of excited state lines to include in the
code.

  Arguments:

At present, there are no arguments.  The input file is hardwired

  Returns:

  Notes:


	The energy levels in the kurucz line lists are in cm**-1;
        The multiplicities are expressed in j so g = 2j +1

  History:
	00oct	ksl	Coded
	00nov25	ksl	Began trying to speed this up.  Originally, this
			took several hours to run because gfall.dat
			was read multiple times.  Now it is read
                        only once, and all the lines of "interest"
			are stored internally.  It takes about a minute.
	01apr13	ksl	Changed levels.out format to better match 
			existing get_atomic_data input requirements.
	01jul11	ksl	Modified the levels.out format so that it now
			matches the first part of the the atomic datafile.
			At present it produces a fully usable python input
			file for elements, ions, and levels, although it
			should be pointed out that this file is quite
			larege (with 7000 levels and 41,000 lines)


************************************************************************/

#include	<math.h>
#include	 <stdio.h>
#include        <strings.h>

#define HC	1.98587e-16
#define C	2.997925e10
#define HEV 	4.13620e-15	/* Planck's constant in eV */
#define BOLTZMANN 1.38062e-16
#define 	EPSILON 1e-5



struct kurucz
{
  double freq, gf;
  int elem, ion;
  double e1, e2;
  double j1, j2;
}
kz;

#define NIONS       500
#define NLEVELS   1000
#define NLINES    100000
#define LINELENGTH 132
#define KLINELENGTH 160

struct ions
{
  int z;
  char ele[20];
  int istate, g;
  int firstlevel;
  int nlevels;
  float xp;
}
ion[NIONS];
int nion;

struct configurations
{
  int q_num;			/* principal quantum number */
  double g;			/* multiplicity for level */
  double ex;			/*excitation energy of level */
}
xconfig[NLEVELS], config[NLEVELS];

struct lines
{
  int z, i, g1, g2;
  double w, f, e1, e2;
  int n1, n2;
}
xl[NLINES];
int nlines;




main (argc, argv)
     int argc;
     char *argv[];
{

  char kline[KLINELENGTH], jline[KLINELENGTH];
  int lineno;
  double eion;
  FILE *iptr, *fptr, *optr, *lptr, *kptr, *fopen ();
  int n, nlev, l1, l2;
  int ielem, istate;
  float t, x, z;
  char word1[10], word2[10];
  int iielem, iion, ig;
  float xp;
  int nn;
  int nline;
  char dummy[LINELENGTH];
  char root[LINELENGTH], input[LINELENGTH];
  char ionfile[LINELENGTH], linefile[LINELENGTH], levelfile[LINELENGTH];



  printf
    ("This program produces atomic data from the kurucz list of lines\n");

  /* Determine whether input data should be read from a file or from the terminal.  */
  if (argc == 2)
    {
      strcpy (dummy, argv[1]);
    }
  else
    {
      printf ("Input file (interactive=stdin):");
      fgets (dummy, LINELENGTH, stdin);
    }
  get_root (root, dummy);

  if (strncmp (root, "dummy", 5) == 0)
    {
      printf
	("Proceeding to create rdpar file in dummy.pf, but will not run prog\n");
    }
  else if (strncmp (root, "stdin", 5) == 0 || strncmp (root, "rdpar", 5) == 0
	   || root[0] == ' ' || strlen (root) == 0)
    {
      strcpy (root, "mod");
      printf
	("Proceeding in interactive mode\n Output files will have rootname mod\n");
    }
  else
    {
      strcpy (input, root);
      strcat (input, ".pf");

      opar (input);
      printf ("Reading data from file %s\n", input);

    }


/* Initialize variables */

  t = 100000;
  nion = 0;
  strcpy (ionfile, "atomic/atomic.data");
  strcpy (linefile, "lines.out");
  strcpy (levelfile, "levels.out");


/* Get inputs */

  rdstr ("atom_ion_infile", ionfile);
  rdstr ("atom_ion_level_outfile", levelfile);
  rdstr ("line_outfile", linefile);
  rdflo ("Temperature", &t);

// Got aall the inputs
  cpar ("read_kurucz.pf");

// Open the remainging input and output files
  iptr = fopen (ionfile, "r");
  lptr = fopen (levelfile, "w");
  optr = fopen (linefile, "w");

/* First read the file which contains the ions of interest */


  while (fgets (jline, KLINELENGTH, iptr) != NULL)
    {
      sscanf (jline, "%s %s %d %d %d %e ", word1, word2, &iielem, &iion,
	      &ig, &xp);
// Echo the elements line of the input file to the levelfile
      if (strncmp (word1, "Element", 7) == 0)
	fprintf (lptr, "%s", jline);

      if (strncmp (word1, "Ion", 3) == 0)
	{
//	  fputs (jline, lptr);
	  strcpy (ion[nion].ele, word2);
	  ion[nion].z = ielem = iielem;
	  ion[nion].istate = iion;
	  ion[nion].xp = xp;
	  ion[nion].g = ig;
	  nion++;
	}
    }

/* Create an intermediate file which contains all of the lines
   one wishes to use out of the total line list  */

  fptr = fopen ("linelist/gfall.dat", "r");
  while (fgets (kline, KLINELENGTH, fptr) != NULL)
    {
      lineno++;
      sscanf (kline, "%11lf%7lf%4d.%2d%12lf%5lf", &kz.freq, &kz.gf, &kz.elem,
	      &kz.ion, &kz.e1, &kz.j1);
      sscanf (&kline[52], "%12lf%5lf", &kz.e2, &kz.j2);

      for (nn = 0; nn < nion; nn++)
	{
	  ielem = ion[nn].z;
	  istate = ion[nn].istate - 1;
	  if (kz.elem == ielem && kz.ion == istate && ion[nn].xp > 13.4)
	    {
	      z = pow (10., (double) kz.gf);
	      x = z * exp ((double) (-HC * (double) kz.e2 / (BOLTZMANN * t)));
	      if (x > 1e-3 && kz.freq < 1000.)
		{
		  xl[nlines].z = kz.elem;
		  xl[nlines].i = kz.ion;
		  xl[nlines].f = pow (10., (double) kz.gf);
		  xl[nlines].g1 = 2 * kz.j1 + 1;
		  xl[nlines].f /=xl[nlines].g1; // Convert from gf to f as python expects
		  xl[nlines].g2 = 2 * kz.j2 + 1;
		  xl[nlines].w = 10. * kz.freq;	// Convert from nanometers to Angstroms

		  xl[nlines].e1 = HEV * C * kz.e1;
		  xl[nlines].e2 = HEV * C * kz.e2;
//if(xl[nlines].z==1 && xl[nlines].i==0) puts(kline);
		  nlines++;
		  if (nlines == NLINES)
		    {
		      printf
			("Too many lines selected at outset. Increase NLINES??\n");
		    }
		}
	    }
	}

    }
  fclose (fptr);

// Now work on only the lines we want.  Main loop is over the ions

  for (nn = 0; nn < nion; nn++)
    {
// Write out the ion line
      fprintf (lptr, "\nIon   %5s %8d %8d %8d %8.5g\n", ion[nn].ele, ion[nn].z,
	       ion[nn].istate, ion[nn].g, ion[nn].xp);

      ielem = ion[nn].z;
      istate = ion[nn].istate - 1;
/* Read and print the input file */
      lineno = 0;
      nlev = 0;
      for (nline = 0; nline < nlines; nline++)
	{
	  if (xl[nline].z == ielem && xl[nline].i == istate)
	    {

	      fprintf (optr,
		       "Line %2d %2d %11.6f %11.6f %3d %3d %11.6f %11.6f\n",
		       ion[nn].z, ion[nn].istate, xl[nline].w, xl[nline].f,
		       xl[nline].g1, xl[nline].g2, xl[nline].e1,
		       xl[nline].e2);
	      if (nlev == 0 && xl[nline].e1 >= 0.0)
		{
		  xconfig[0].ex = xl[nline].e1;
		  xconfig[0].g = xl[nline].g1;
		  nlev = 1;
		}

	      l1 = l2 = 0;
	      for (n = 0; n < nlev; n++)
		{
//                if (xl[nline].e1 == xconfig[n].ex && xl[nline].g1==xconfig[n].g)
		  if (xl[nline].e1 == xconfig[n].ex)
		    l1 = 1;
//                if (xl[nline].e2 == xconfig[n].ex && xl[nline].g2==xconfig[n].g)
		  if (xl[nline].e2 == xconfig[n].ex)
		    l2 = 1;
		}
	      if (l1 == 0 && xl[nline].e1 >= 0.0)
		{
		  xconfig[nlev].ex = xl[nline].e1;
		  xconfig[nlev].g = xl[nline].g1;
		  nlev++;
		}
	      if (l2 == 0 && xl[nline].e2 >= 0.0)
		{
		  xconfig[nlev].ex = xl[nline].e2;
		  xconfig[nlev].g = xl[nline].g2;
		  nlev++;
		}
	    }

	}
/* Now print out the results */

      if (nlev > 0)
	sortem (nlev);

      fprintf (optr, "# There are %d unique levels for %d %d\n", nlev,
	       ion[nn].z, ion[nn].istate);
      fprintf (lptr, "# There are %d unique levels for %d %d\n", nlev,
	       ion[nn].z, ion[nn].istate);


      for (n = 0; n < nlev; n++)
	{
//	  fprintf (lptr, "Level %3d  %3d %10.6f %3.0f %3d\n", ion[nn].z,
//		   ion[nn].istate, config[n].ex, config[n].g, n);
	  fprintf (lptr, "Level %3d %3d %3d  %3.0f %10.6f\n", ion[nn].z,ion[nn].istate,n, config[n].g,
		   config[n].ex);
	}
    }
}


int
sortem (nlev)
     int nlev;
{
  float x[NLEVELS + 1];
  int i, iorder[NLEVELS + 1];

  for (i = 0; i < nlev; i++)
    {
      x[i + 1] = xconfig[i].ex;
    }
  indexx (nlev, x, iorder);

  for (i = 0; i < nlev; i++)
    {
      config[i].ex = xconfig[iorder[i + 1] - 1].ex;
      config[i].g = xconfig[iorder[i + 1] - 1].g;
    }

}
