#include <stdio.h>
#include "validation_functions.c"

int main()
{

	//printf("%d\n", dow(2012,12,29));

	printf("Julian Date: %f\n", greg_to_jul(2012,12,30,20,28,0));

	printf("Time Mod: %f\n", time_mod(greg_to_jul(2012,12,30,20,28,0)));

	return 0;
}
