//#include <iostream.h>
//#include <fstream.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

main () {
  FILE *fin/*, *fout*/;
  int i, j, numBins = 400;
  double RDF[numBins][4], Lhalf = 1.327, n = 80.22, Integral = 0.0, Integral2 = 0.0, deltaR;
  double pi = 3.14159;
  fin = fopen ("RDF.out", "r");
 // fout = fopen ("DielConst.out", "w");

for(i=0; i<numBins; i++)
//fscanf(fin, "%lf   %lf    %lf    %lf    %lf   %lf    %lf    %lf     %lf   %lf    %lf", &Induction[i][0], &Induction[i][1], &Induction[i][2], &Induction[i][3], &Induction[i][4], &Induction[i][5], &Induction[i][6], &Induction[i][7], &Induction[i][8], &Induction[i][9], &Induction[i][10]); 
fscanf(fin, "%lf   %lf    %lf    %lf", &RDF[i][0], &RDF[i][1], &RDF[i][2], &RDF[i][3]); 

//for(i=0;i<numBins; i++)
//printf("%f   %f    %f    %f    %f\n", Induction[i][0], Induction[i][1], Induction[i][2], Induction[i][3]);

deltaR = Lhalf/numBins;

for (i=1; i<=118; i++)
	{
	Integral = Integral + RDF[i][0]*RDF[i][0]*RDF[i][3];
	}

for (i=1; i<=75; i++)
	{
	Integral2 = Integral2 + RDF[i][0]*RDF[i][0]*RDF[i][2];
	}

	Integral = 4*pi*(n/3)*Integral*deltaR;
	Integral2 = 4*pi*(n)*Integral2*deltaR;

//printf ("deltaR = Lhalf/numBins = %f \n", deltaR);
//printf ("O-O Integration range - 1st or 2nd minima = %d \n", i);
printf ("Coordination Number O-O = %f \n", Integral);
printf ("Coordination Number O-H = %f \n", Integral2);

fclose ( fin );
//fclose ( fout );

}
