


/**************************************************************************
                    Space Telescope Science Institute


  Synopsis:  

  Description:	

  Arguments:  (Input via .pf file)		


  Returns:

  Notes:

  This is roughly what the Makefile should look like. Uncomment the mv 
  statement to store binary in ~/bin

CC = gcc
CFLAGS = -c -g -I$$HOME/include -Wall
LDFLAGS= -L$$HOME/lib -lm -lkpar
BIN = $$HOME/bin

simple: simple.o
	gcc simple.o  $(LDFLAGS) -o simple
#	mv $@ $(BIN)

 

  History:
2004	ksl	Coded as better.c

 ************************************************************************/

#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


#include "log.h"

#define LINELENGTH 132

int
main (argc, argv)
     int argc;
     char *argv[];
{

  FILE *fptr, *fopen ();
  char parfile[LINELENGTH];
  char infile[LINELENGTH];
  char line[LINELENGTH];
  char firstword[LINELENGTH];

/* Parse the command line */
  if (argc == 1)
    strcpy (parfile, "simple.pf");
  else if (argc == 2)
    {
      strcpy (parfile, argv[1]);
      if ((strstr (parfile, ".pf")) == NULL)
	strcat (parfile, ".pf");
    }
  else
    {
      printf ("Usage: simple [file.pf]\n");
      exit (0);
    }

/* Initialize variables that will be read by rdpar */

  strcpy (infile, "file.ls");

/* Get data from a parameter file */


  opar (parfile);


  rdstr ("filename", infile);

  cpar (parfile);

/* End of input section */


  /* Open the input file.  Exit if it is not opened successfully */
  if ((fptr = fopen (infile, "r")) == NULL)
    {
      printf ("Failed to open file %s\n", infile);
      exit (0);
    }

/* Read and print the input file */
  while (fgets (line, LINELENGTH, fptr) != NULL)
    {
      sscanf (line, "%s", firstword);
      Error (line);
    }


/* OK now end program and go on to the next important task. */
  error_summary ("End of program");
  exit (0);

}
