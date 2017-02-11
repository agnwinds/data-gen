/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

	This program reads the Kurucz linelist, and creates
	files that can be read by get_atomic_data.  

  Description:	

	read_kurucz creates two files, levels.out, which contains the energy levels for all 
	ions of interest, and, lines.out, which contains the lines associated 
	with each ion.  

  Arguments:


  Returns:

  Notes:


	The energy levels in the kurucz line lists are in cm**-1;
        The multiplicities are expressed in j so g = 2j +1

	The Kurucz list has lines in which the first level has a higher energy than the
	second level.  According to Ivan this is because Kurucz chose to order the
	levels by parity, even or odd in gfall.dat.

	0808 - According to the kurucz linelist format, there are allowances
	made for different isotopes.  I could not any indication that this
	field contained anythinng but zeros.

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
			large (with 7000 levels and 41,000 lines)
	01nov15	ksl	Wrote analogous program for Verner linelist, and
			mad a few fixes to this as a result.
	01nov23	ksl	Fixed so that under normal circumstances the
			levels file does not also contain elements and
			ions.  This can be changed using the variable integ
	0807	ksl	Checked this program on bluegill after not running
			it for many years, and made a few minor changes
			to check for files it is trying to read, instead of
			simply returning a bus error
	0808	ksl	Gave user more control over what is produced.  Added
			print out of radiative decay rate for uppoer level
			where available.  1.0 in the output for this indicates
			no known value.  It's not entirely clear how to 
			like this to the levels.


************************************************************************/

#include	<math.h>
#include	 <stdio.h>
#include        <strings.h>
#include        <stdlib.h>

#define HC	1.98587e-16
#define C	2.997925e10
#define HEV 	4.13620e-15	/* Planck's constant in eV */
#define BOLTZMANN 1.38062e-16
#define 	EPSILON 1e-5



struct kurucz
{
  double freq, gf;
//  int elem, ion,iso_no;
  int elem, ion;
  double e1, e2;
  double j1, j2;
  double gamma_rad;
}
kz;

#define NELEMENTS    50
#define NIONS       1000
#define NLEVELS   10000
#define NLINES    100000
#define LINELENGTH 132
#define KLINELENGTH 160

struct elements
{				/* Element contains physical parameters that apply to element as a whole and
				   provides and index to the atomic data */
  char name[20];		/* Element name */
  int z;			/* Atomic number */
  int firstion;			/*Index into struct ions  ion[firstion] is the lowest ionization state of this ion */
  int nions;			/*The number of ions actually read in from the data file for this element */
  double abun;			/*Abundance */
}
ele[NELEMENTS];
int nelements;			// The actual number of elements

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
  double gamma_rad;
//  int iso_no;
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
  char root[LINELENGTH], input[LINELENGTH],kuruczfile[LINELENGTH];
  char ionfile[LINELENGTH], linefile[LINELENGTH], levelfile[LINELENGTH];
  int integ;
  float wavmin,wavmax;  // wavelenght ranges for lines that are produced
  float ipmin,ipmax;   //  mim and maximum ioonization of ions to include
  float gf_boltz;         //


  integ = 1;			//0 implies writing out elements and ions into the level file, otherwise don't

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
      printf ("Proceeding in interactive mode\n ");
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
  nelements = nion = 0;

  strcpy(kuruczfile,"linelist/gfall.dat");
  strcpy (ionfile, "elem_ions_ver.py");
  strcpy (linefile, "lines_simple_kur.py");
  strcpy (levelfile, "levels_kur.py");

  wavmin=0.0;
  wavmax=10000;

  ipmin=13.4;  
  ipmax=1000.;
  gf_boltz=1.e-3;


/* Get inputs */

  rdstr ("atom_ion(infile)", ionfile);
  rdstr ("Kurucz.linelist",kuruczfile);
  rdstr ("levels(outfile)", levelfile);
  rdstr ("lines_simple(outfile)", linefile);
  rdflo ("Temperature", &t);
  rdflo ("wave.min(Ang)",&wavmin);
  rdflo ("wave.max(Ang)",&wavmax);
  rdflo ("ion.ipot.min(eV)",&ipmin);
  rdflo ("ion.ipot.max(ev)",&ipmax);
  rdflo ("gf_boltzmann_min",&gf_boltz);

wavmin/=10.; // convert to nm
wavmax/=10.;

// Got all the inputs
  cpar ("py_read_kurucz.pf");

// Open the remainging input and output files
  if((iptr = fopen (ionfile, "r"))==NULL){
	  Error("Could not open %s\n",ionfile);
	  exit(0);
  }

  lptr = fopen (levelfile, "w");
  optr = fopen (linefile, "w");

/* First read the file which contains the ions of interest */


  while (fgets (jline, KLINELENGTH, iptr) != NULL)
    {
      strcpy (word1, "");
      strcpy (word2, "");
      sscanf (jline, "%s %d %s", word1, &iielem, word2);
// Echo the elements line of the input file to the levelfile
// Note that the position of the element name is different for an ion record and element rec
      if (strncmp (word1, "Element", 7) == 0)
	{
	  strcpy (ele[nelements].name, word2);
	  ele[nelements].z = iielem;
	  nelements++;
	  if (nelements > NELEMENTS)
	    {
	      Error ("Too many elements to store %d\n", nelements);
	      exit (0);
	    }
	  if (integ == 0)
	    fprintf (lptr, "%s", jline);
	}
      else if (strncmp (word1, "Ion", 3) == 0)
	{
	  sscanf (jline, "%s %s %d %d %d %e ", word1, word2, &iielem, &iion,
		  &ig, &xp);
	  n = 0;
	  while (ele[n].z != iielem && n < nelements)
	    n++;
	  if (n < nelements)
	    {
	      strcpy (ion[nion].ele, word2);
	      ion[nion].z = ielem = iielem;
	      ion[nion].istate = iion;
	      ion[nion].xp = xp;
	      ion[nion].g = ig;
	      nion++;
	    }
	  else
	    {
	      printf ("Warning: There is no element for this ion:\n %s\n",
		      jline);
	    }
	}
    }

/* Create an intermediate file which contains all of the lines
   one wishes to use out of the total line list  */

  if ((fptr = fopen (kuruczfile, "r")) ==NULL){
	  Error("Could not open linelist %s \n",kuruczfile);
	  exit(0);
  }

  while (fgets (kline, KLINELENGTH, fptr) != NULL)
    {
      lineno++;
      sscanf (kline, "%11lf%7lf%4d.%2d%12lf%5lf", &kz.freq, &kz.gf, &kz.elem,
	      &kz.ion, &kz.e1, &kz.j1);
      sscanf (&kline[52], "%12lf%5lf", &kz.e2, &kz.j2);
      sscanf (&kline[80], "%6lf",&kz.gamma_rad);        //080810 - get the raditive constant
//      sscanf (&kline[106], "%3d",&kz.iso_no);        //080810 - iso_no

      for (nn = 0; nn < nion; nn++)
	{
	  ielem = ion[nn].z;
	  istate = ion[nn].istate - 1;


	  /* The next few lines "select" lines that are likely to be of interest.  
	   * The first constraint is that the ionization potential of the ion be
	   * greater than 13.4 eV
	   *
	   * The second is that the gf of the line times and occupation factor
	   * be greater than a certain value here, 1e-3
	   *
	   * */

	  if (kz.elem == ielem && kz.ion == istate && (ipmin <= ion[nn].xp) && (ion[nn].xp <= ipmax))
	    {
	      z = pow (10., (double) kz.gf);
	      x = z * exp ((double) (-HC * (double) kz.e2 / (BOLTZMANN * t)));
	      if (x > gf_boltz && wavmin< kz.freq && kz.freq < wavmax)
		{
		  xl[nlines].z = kz.elem;
		  xl[nlines].i = kz.ion;
		  xl[nlines].w = 10. * kz.freq;	// Convert from nanometers to Angstroms

		  if (kz.e2 > kz.e1)
		    {
		      xl[nlines].e1 = HEV * C * kz.e1;
		      xl[nlines].g1 = 2 * kz.j1 + 1;
		      xl[nlines].e2 = HEV * C * kz.e2;
		      xl[nlines].g2 = 2 * kz.j2 + 1;
		    }
		  else
		    {		// Switch the order
		      xl[nlines].e1 = HEV * C * kz.e2;
		      xl[nlines].g1 = 2 * kz.j2 + 1;
		      xl[nlines].e2 = HEV * C * kz.e1;
		      xl[nlines].g2 = 2 * kz.j1 + 1;
		    }
		  xl[nlines].f = pow (10., (double) kz.gf);
		  xl[nlines].f /= xl[nlines].g1;	// Convert from gf to f as python expects
		      xl[nlines].gamma_rad=pow(10.,(double) kz.gamma_rad);
//		      xl[nlines].iso_no=kz.iso_no;

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
      if (integ == 0)
	fprintf (lptr, "\nIon   %5s %8d %8d %8d %8.5g\n", ion[nn].ele,
		 ion[nn].z, ion[nn].istate, ion[nn].g, ion[nn].xp);
      else
	{
	  fprintf (lptr, "# Comment-- Ion %s %2d\n", ion[nn].ele,
		   ion[nn].istate);


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
			   "Line %2d %2d %11.6f %11.6f %3d %3d %11.6f %11.6f %8.3e\n",
			   ion[nn].z, ion[nn].istate, xl[nline].w,
			   xl[nline].f, xl[nline].g1, xl[nline].g2,
			   xl[nline].e1, xl[nline].e2, xl[nline].gamma_rad);
		  if (nlev == 0 && xl[nline].e1 >= 0.0)
		    {
		      xconfig[0].ex = xl[nline].e1;
		      xconfig[0].g = xl[nline].g1;
		      nlev = 1;
		    }

		  /* Create the configurations */

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
	}
/* Now print out the results */

      if (nlev > 0)
	sortem (nlev);

      fprintf (lptr,
	       "# There are %d unique levels for %s z %d istate %d\n",
	       nlev, ion[nn].ele, ion[nn].z, ion[nn].istate);


      for (n = 0; n < nlev; n++)
	{
	  fprintf (lptr, "LevKur %3d %3d %3d  %3.0f %10.6f\n", ion[nn].z,
		   ion[nn].istate, n, config[n].g, config[n].ex);
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
