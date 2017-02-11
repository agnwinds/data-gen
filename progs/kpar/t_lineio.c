#include <stdio.h>
main ()
{
  FILE *ptr, *fopen (), *dptr;
  char line[132];

  ptr = fopen ("test.dat", "r");

  while (get_line (ptr, line) != EOF)
    printf ("%s", line);

  printf ("Now test a failure to open a file\n");


  dptr = fopen ("test.bog", "r");
  get_line (dptr, line);
  printf ("Should have gotten an error message\n");
}
