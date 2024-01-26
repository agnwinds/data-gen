
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  py_read_verner  name.pf

  Description:	

This program reads the file atom1.dat containing permitted resonance lines 
compiled by Verner, Verner and Ferland Atomic Data Nucl Data Tables in press
and reformats it so that it can be used with get_atomic.data.  A file
containing all the levels is also produced. 

The results are nominal written to lines_ver.py, although this is can be 
modified when running the program.  

Originally, the elements and ions selected were controlled by the structure 
"element", and this can still be done by responding "none" to the 
"Atom_ion_input_filename(none=oldstyle)" question.   


  Arguments:		

The program uses rdpar:

"Atom_ion_input_file(none=oldstyle)	ionfile  
	file which controls what elements and ions are included on output.  
	none implies elements and ions are are controlled internally 
	within the program using the element structure.
"Line_output_file" 			linefile);                  
	file containing lines selected in python format
"Atom_ion_level_output_file" 		levelfile);        
	file containing lines selected in python format

  Returns:

  Notes:

Unlike the Kurucz linelist, the Verner list does not include the electron 
number.

The output line list produced by py_read_verner does not contain the links
to the level numbers. To produce these links (which are/will be neccessary)
for sophisticated line treatment the program py_link must be run successfully.


  History:
	1/97	ksl 	Source for this modification Written and debugged. 
	01nov	ksl	I have lost the version of this program that
			was used to create lines.data.  Began to modify
			the only version of the program I could find to
			first reproduce lines.data as it exists. 
				(a) fixed element struct to produce same lines for
					same ions as previously
				(b) eliminated lines with unknown wavelength (which
					Verner had set to zero.
	01nov	ksl	Modify output to include uper and lower level energies of
			line.  This will enable get_atomicdata.c to be simplified,
			and a more accurate treatment of the lines in python.
	01nov	ksl	Modify to also generate a level file in the same format
			as produced by py_kurucz
	01nov	ksl	Added a currently hardwired variable integ to control
			whether element and ion lines were printed out along
			with the levels.
	12oct	ksl	Fixed a problem with parsing lines with elements which
			had two characters and ions that had high numbers
	

 ************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "templates.h"
#include "log.h"

#define LINELENGTH 132
#define KLINELENGTH 160
#define WORDLENGTH	25
#define WORDS		11

#define H   				6.6262e-27
#define HC  				1.98587e-16
#define HEV				4.13620e-15	/* Planck's constant in eV */
#define C   				2.997925e10

#define NIONS       1000
#define NLEVELS   1000
#define NLINES    100000



/* This is used to locate which character a parameter is in the fixed format file */
int first[] = { 0, 5, 39, 47, 57, 73, 89, 92, 95, 104, 113 };

/* This is the structure that tells one which elements and ions to select from the 
 * table */

#define NELEM 14
struct element
{
  char name[5];
  int z;
  int ion_min, ion_max;
}
e[] =
{
"H", 1, 1, 2,
    "He", 2, 1, 2,
    "C", 6, 2, 6,
    "N", 7, 3, 7,
    "O", 8, 3, 7,
    "Ne", 10, 1, 10,
    "Na", 11, 2, 11,
    "Mg", 12, 2, 8,
    "Al", 13, 2, 8,
    "Si", 14, 2, 14,
    "S", 16, 2, 16, 
    "Ar", 18, 1, 18, 
    "Ca", 20, 1, 20, 
    "Fe", 26, 2, 26};

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


int
main (argc, argv)
     int argc;
     char *argv[];
{

  char dummy[LINELENGTH];
  char root[LINELENGTH], input[LINELENGTH];
  char ionfile[LINELENGTH], linefile[LINELENGTH], levelfile[LINELENGTH];
  int integ;



  integ=1;  //0 implies writing out elements and ions into the level file

  printf
    ("This program produces atomic data from the Verner list of lines\n");

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

  nion = 0;
  strcpy (ionfile, "elem_ions_ver.py");
  strcpy (linefile, "lines_simple_ver.py");
  strcpy (levelfile, "levels_ver.py");


/* Get inputs */

  rdstr ("Atom_ion_input_file(none=oldstyle)", ionfile);
  rdstr ("Line_output_file", linefile);
  rdstr ("Atom_ion_level_output_file", levelfile);

// Got all the inputs
  cpar ("py_read_verner.pf");

  if (strncmp (ionfile, "none", 4) == 0)
    {
      populate_ions ();
      remove (levelfile);
    }

  else
    read_ions (ionfile, levelfile,0);

  getlines ();
  putlines (linefile);
  putlevels (levelfile,integ);

}

int
read_ions (ionfile, levelfile,integ)
     char *ionfile, *levelfile;
     int integ;
{
  FILE *fopen (), *iptr, *lptr;
  char jline[KLINELENGTH];
  char word1[10], word2[10];
  int iielem, iion, ig;
  float xp;


  lptr = fopen (levelfile, "w");
  iptr = fopen (ionfile, "r");

  while (fgets (jline, KLINELENGTH, iptr) != NULL)
    {
      strcpy (word1, "");

      sscanf (jline, "%s %s %d %d %d %e ", word1, word2, &iielem, &iion,
	      &ig, &xp);
// Echo the elements line of the input file to the levelfile
      if (strncmp (word1, "Element", 7) == 0 && integ==0)
	{
	  fprintf (lptr, "%s", jline);
	}
      else if (strncmp (word1, "Ion", 3 ) == 0)
	{
	  strcpy (ion[nion].ele, word2);
	  ion[nion].z = iielem;
	  ion[nion].istate = iion;
	  ion[nion].xp = xp;
	  ion[nion].g = ig;
//xxx          if(integ==0) fputs (jline, lptr);
//xxx         else fprintf(lptr,"# Comment-- Ion %s %2d\n",word2,iion);
	  nion++;
	  if (nion > NIONS)
	    {
	      Error
		("read_ions: nion %d > NIONS (Need to increase NIONS in program) %d\n",
		 nion, NIONS);
	      exit (0);
	    }
	}
    }
  fclose (lptr);
  fclose (iptr);
}

int
populate_ions ()
{
  int nelem, m;

  nion = 0;
  for (nelem = 0; nelem < NELEM; nelem++)
    {
      m = e[nelem].ion_min;
      while (m <= e[nelem].ion_max)
	{
	  ion[nion].z = e[nelem].z;
	  strcpy (ion[nion].ele, e[nelem].name);
	  ion[nion].istate = m;
	  nion++;
	  m++;


	}
      printf ("nelem,nion %d %d\n", nelem, nion);

    }
  printf ("nion %d\n", nion);
}

int
getlines ()
{
  FILE *fptr, *fopen ();
  char line[LINELENGTH];
  int lengths[WORDS];
  char words[WORDS][WORDLENGTH];
  char element[WORDS], old_element[WORDS];
  int istate, old_istate;
  int nelem;
  int nnion;

  int nrecords;
  int i, j, k, l, m, n, mmax;
  int z, ielec;
  float f, lambda;
  int g_i, g_k;
  float energy_i, energy_k;
  int icount;


  for (n = 0; n < WORDS - 1; n++)
    {
      lengths[n] = first[n + 1] - first[n] - 1;
    }
  lengths[WORDS - 1] = 2;

  fptr = fopen ("atom1.dat", "r");

  nrecords = nlines = 0;
  while (fgets (line, LINELENGTH, fptr) != NULL)
    {
      // printf ("%s\n", line);
      nrecords++;
      for (n = 0; n < WORDS; n++)
	{
	  strncpy (&words[n][0], &line[first[n]], lengths[n]);
	}
      strcpy (element, "");
      /* The next line fails when the ion sstate is greater then 10 and element has two characters so
       * made modfications to capture this state*/

      printf("Test x%sx\n",words[0]);
      if (icount=sscanf (&words[0][0], "%s %d", element, &istate) == 1)
	{
		sscanf(&words[0][0],"%2s%d",element, &istate);
		printf("Found %s %d\n",element,istate);
	}
      else if (icount!=2) {
	      printf ("%d Skipping record: %s \n",icount,line);
	      continue;

     }
      else {
		printf("Found %s %d\n",element,istate);
      }

      nnion = 0;
      while (nnion < nion &&
	     (strncmp (element, ion[nnion].ele, 2) != 0
	      || istate != ion[nnion].istate))
	nnion++;

      if (nnion < nion)		// Then we have located an ion that we want
	{
	  sscanf (&words[10][0], "%d", &mmax);
	  printf ("Found %s %d %d \n", element, istate, mmax);
	  for (m = 0; m < mmax; m++)
	    {
	      if (fgets (line, LINELENGTH, fptr) == NULL)
		exit (0);
	      for (n = 0; n < WORDS; n++)
		{
		  strncpy (&words[n][0], &line[first[n]], lengths[n]);
		}

	      sscanf (&words[3][0], "%f", &lambda);
	      sscanf (&words[4][0], "%f", &energy_i);
	      sscanf (&words[5][0], "%f", &energy_k);
	      energy_i *= (HEV * C);
	      energy_k *= (HEV * C);
	      sscanf (&words[6][0], "%d", &g_i);
	      sscanf (&words[7][0], "%d", &g_k);
	      sscanf (&words[9][0], "%f", &f);
	      if (lambda > 0)
		{
		  printf
		    ("Line %2d %2d %11.6f %11.6f %3d %3d %11.6f %11.6f\n",
		     ion[nnion].z, istate, lambda, f, g_i, g_k, energy_i,
		     energy_k);
		  xl[nlines].z = ion[nnion].z;
		  xl[nlines].i = istate;
		  xl[nlines].g1 = g_i;
		  xl[nlines].g2 = g_k;
		  xl[nlines].w = lambda;
		  xl[nlines].f = f;
		  xl[nlines].e1 = energy_i;
		  xl[nlines].e2 = energy_k;
		  nlines++;
		}
	    }


	}

    }

  printf ("Grand total nlines %d\n", nlines);
  return EXIT_SUCCESS;
}



int
putlines (linefile)
     char *linefile;
{
  FILE *fopen (), *optr;
  int n;

  optr = fopen (linefile, "w");
  for (n = 0; n < nlines; n++)
    {
      fprintf (optr,
	       "Line %2d %2d %11.6f %11.6f %3d %3d %11.6f %11.6f\n",
	       xl[n].z, xl[n].i, xl[n].w, xl[n].f,
	       xl[n].g1, xl[n].g2, xl[n].e1, xl[n].e2);

    }
}

int
putlevels (levelfile,integ)
     char *levelfile;
     int integ;
{
  FILE *fopen (), *lptr;
  int nn;
  int ielem, istate, lineno, nlev, nline, l1, l2, n;

  lptr = fopen (levelfile, "a");	// We may already have created this


  /* The outermoost loop is over ions */

  for (nn = 0; nn < nion; nn++)
    {
// Write out the ion line
      if(integ==0) fprintf (lptr, "\nIon   %5s %8d %8d %8d %8.5g\n", ion[nn].ele,
	       ion[nn].z, ion[nn].istate, ion[nn].g, ion[nn].xp);
         else fprintf(lptr,"# Comment-- Ion %s %2d\n",ion[nn].ele,ion[nn].istate);

      ielem = ion[nn].z;
//      istate = ion[nn].istate - 1;  ?? This was the line in py_kurucz.  Why ??
      istate = ion[nn].istate;

/* Read and print the input file */
      lineno = 0;
      nlev = 0;
      for (nline = 0; nline < nlines; nline++)
	{
		/* The next lines selects out those lines associated with a particular
		 * ion and ion state
		 */

	  if (xl[nline].z == ielem && xl[nline].i == istate)
	    {

		    /* By definition the next line has to be populated*/
	      if (nlev == 0 && xl[nline].e1 >= 0.0)
		{
		  xconfig[0].ex = xl[nline].e1;
		  xconfig[0].g = xl[nline].g1;
		  nlev = 1;
		}

	      l1 = l2 = 0;
	      for (n = 0; n < nlev; n++)
		{
		  if (xl[nline].e1 == xconfig[n].ex)
		    l1 = 1;
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


      for (n=1;n<nlev;n++){
	      // printf("Check %8.6f\n", xconfig[n].ex);
      }
	  /* Line above is the end of the loop over a particular ion and state
	   * At this point xconfig is complete but it needs to be sorted */
/* Now print out the results */

      /* Resort the levels into ascending order and number them */
      if (nlev > 0)
	sortem (nlev);

      fprintf (lptr, "# There are %d unique levels for %d %d\n", nlev,
	       ion[nn].z, ion[nn].istate);
      for (n = 0; n < nlev; n++)
	{
//        fprintf (lptr, "Level %3d  %3d %10.6f %3.0f %3d\n", ion[nn].z,
//                 ion[nn].istate, config[n].ex, config[n].g, n);
	  fprintf (lptr, "Level %3d %3d %3d  %3.0f %10.6f\n", ion[nn].z,
		   ion[nn].istate, n, config[n].g, config[n].ex);
	}
    }
}


/* this routines reads the unsorted levels, and sorts them 
 * into ascending order putting the results in the structure 
 * xcomfig.  It returns the number of levels
 */

int
sortem (nlev)
     int nlev;
{
  float x[NLEVELS + 1];
  // int i, iorder[NLEVELS + 1];
  int i;
  unsigned long nnlev,iorder[NLEVELS + 1];

  // printf ("Check: Entered sortem %d \n",nlev);
  for (i = 0; i < nlev; i++)
    {
      x[i + 1] = xconfig[i].ex;

      // printf("Check: sortem %f\n",xconfig[i].ex);
    }

  nnlev=nlev;
  indexx (nnlev, x, iorder);

  for (i = 0; i < nlev; i++)
    {
      config[i].ex = xconfig[iorder[i + 1] - 1].ex;
      config[i].g = xconfig[iorder[i + 1] - 1].g;
      // printf("Check2: sortem %d  %f\n",iorder[i+1],config[i].ex);
     
    }

}
