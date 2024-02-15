
/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:

  Description:	This program parses the results of topbase regarding 
energy levels.  


  Arguments:  parse_levels 

  Returns:

  Notes:
	To get the appropriate input file, be sure to check
	all output options, except for quantum defect.

	I have supressed the first column because it is just the
	line number in the query

	Topbase does not use ionization state to identify and ion; 
	instead it uses the number of electrons.  Python uses the
	astronomical ionization state, e.g. C I is neutral.

  History:
	01sep25	ksl	Began work
	01nov	ksl	Make modifications so that insofar as possible
			the output from this program is in python units,
			e.g electron volts for energy levels.  Also 
			correct to experimental energies.
	01dec	ksl	Allowed for lines without radiative lifetimes.
			This was some kind of problem I was having
			with Topbase retrievals and if that gets fixed
			by them, it should be put back.  For now, the
			program will read files with and without the
			radiative lifetime.
	

 ************************************************************************/
#include	<math.h>
#include	<stdio.h>
#include	<strings.h>
#include        <stdlib.h>
#include        "templates.h"
#include        "log.h"

#define RYD     13.605698	//Rydberg in eV

#define LINELENGTH 132
#define NLEVELS 10000

struct configurations
{
  int z;			/*The element associated with this configuration */
  int istate;			/*The ion associated with the configuration */
  int nion;			/*Internal index to ion structure */
  int nden;			/*Internal index to levden array in wind structure. -1 if no such level */
  int isp;			/*Topbase description of angular momentum (2s+1)*100+l*10+p */
  int ilv;			/*Topbase level number */
  int g;			/* multiplicity for level */
  double q_num;			/* principal quantum number.  In Topbase this has non-integer values */
  double e, ex;			/*excitation energy of level */
  double eqn, rl;		/* effective quantum number and radiative lifetime */
  char name[20];
}
config[NLEVELS], *zz, *yy[NLEVELS], *yhold[NLEVELS];
int nlev;

main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *optr, *fopen ();
  char infile[80];
  char line[LINELENGTH];
  char firstword[LINELENGTH];
  char root[80], outfile[80];
  int nz, ne, islp, ilv;
  char name[16], xconf[16];
  float e, ex, g, eqn, rl;
  float extop, ebot;
  float qmax;
  int n;
  int i, imin, imax;
  int j, jmin, jmax;
  int nwords;



  strcpy (name, "");
  strcpy (root, "topbase_levels_hhe");
  extop = 4;
  ebot = 0.1;
  qmax = 4;
  nlev = 0;
  imax = -1;
  imin = 1000;

  rdstr ("Root.for.all.files", root);
  rdflo ("Max.energy.above.ground(Ryd)", &extop);
  rdflo ("Min.energy.below.continuum(Ryd)", &ebot);
  rdflo ("Max.qnum", &qmax);
  if (ebot > 0)
    ebot *= -1;


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

  optr = fopen (outfile, "w");

/* Read and print the input file */
  fprintf (optr,
	   "# Maximum excitation above ground (in Ryd)  for inclusion  %f\n",
	   extop);
  fprintf (optr,
	   "# Miniumum excitation below continuum (in Ryd) for inclusion %f\n",
	   ebot);
  fprintf (optr,
	   "# Topbase levels: Order changed to move config to end of line\n");
  fprintf (optr,
	   "#LevTop  z ion  iSLP iLV iCONF   E(eV) Te(eV) gi RL(s) eqn RL(s)\n");

  while (fgets (line, LINELENGTH, fptr) != NULL)
    {
      sscanf (line, "%s", firstword);
      if (firstword[0] == '=' || firstword[0] == 'i')
	{
	  fprintf (optr, "#%s", line);
	}
      else
	{
	  strcpy (name, "");
	  nwords =
	    sscanf (line, "%*d %d %d %d %d %15c %e %e %e %e %e", &nz, &ne,
		    &islp, &ilv, name, &e, &ex, &g, &eqn, &rl);
	  if (nwords < 9)
	    {
	      printf ("Error: only %d words in line: %s\n", nwords, line);
	      exit (0);
	    }
	  strcpy (xconf, "");
	  name[15] = '\0';
	  strncpy (xconf, name, 15);
	  if (ex <= extop && e < ebot && eqn <= qmax)
	    {
	      config[nlev].z = i = nz;
	      config[nlev].istate = nz - ne + 1;
	      config[nlev].isp = islp;
	      config[nlev].ilv = ilv;
	      config[nlev].e = e * RYD;
	      config[nlev].ex = ex * RYD;
	      config[nlev].g = g;
	      config[nlev].eqn = eqn;
//Allow for the possibility of no rl  01dec12 ksl
	      if (nwords == 10)
		config[nlev].rl = rl * 1.e-9;
	      else
		config[nlev].rl = -1;
	      strcpy (config[nlev].name, name);
	      if (i < imin)
		imin = i;
	      if (i > imax)
		imax = i;
	      nlev++;
	    }
	}
    }

  sortem (nlev);

  for (n = 0; n < nlev; n++)
    {
      {
	zz = yy[n];
	fprintf (optr,
		 "LevTop %2d %2d %4d %2d %10f %10f %2d %7.4f %.2e %s \n",
		 zz->z, zz->istate, zz->isp, zz->ilv, zz->e, zz->ex,
		 zz->g, zz->eqn, zz->rl, zz->name);
      }

    }
}


/* Sort the levels into normal python order

Note!! indexx the Numerical Recipes routine does not preserve the order in the event
that two values (of x) are identical, and thus one cannot search first on enrgy level,
then on istate and then on z, unless you add small increments to istate and z to 
keep the current order of the array.  Seems like a major flaw in indexx...and
I have adopted a Kluge which should work as long as there are not too many levels

*/


int
sortem (nlev)
     int nlev;
{
  float x[NLEVELS + 1];
  int i, iorder[NLEVELS + 1];
  int jorder[NLEVELS + 1];

  for (i = 0; i < nlev; i++)
    {
      yhold[i] = yy[i] = &config[i];	// Start by setting up pointers to config
    }
  for (i = 0; i < nlev; i++)
    {
      x[i + 1] = yy[i]->ex;

    }
  indexx (nlev, x, iorder);

// Now rearrange the indices
  for (i = 0; i < nlev; i++)
    {
      yy[i] = yhold[iorder[i + 1] - 1];	// Copy the new order of the pointers into yy
    }
  for (i = 0; i < nlev; i++)
    {
      yhold[i] = yy[i];
    }


//OK now both yy and yhold have excitation order


  for (i = 0; i < nlev; i++)
    {
      x[i + 1] = yy[i]->istate + 1e-4 * i;	// See note
    }
  indexx (nlev, x, iorder);

// Now rearrange the indices
  for (i = 0; i < nlev; i++)
    {
      yy[i] = yhold[iorder[i + 1] - 1];	// Copy the new order of the pointers into yy
    }
  for (i = 0; i < nlev; i++)
    {
      yhold[i] = yy[i];
    }
//OK now both yy and yhold have istate & excitation order


  for (i = 0; i < nlev; i++)
    {
      x[i + 1] = yy[i]->z + 1e-4 * i;	// See note
    }
  indexx (nlev, x, iorder);
// Now rearrange the indices
  for (i = 0; i < nlev; i++)
    {
      yy[i] = yhold[iorder[i + 1] - 1];	// Copy the new order of the pointers into yy
    }
  for (i = 0; i < nlev; i++)
    {
      yhold[i] = yy[i];
    }

//OK now both yy and yhold have istate & excitation and element order

}
