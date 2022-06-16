//#include <iostream.h>
//#include <fstream.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

main () {
  FILE *fin/*, *fout*/;
  int i, j, numsteps = 600, numMolec = 500;
  double DipoleSum[3]={0,0,0}, Induction[numsteps][5], dipole=0, Dipole2Sum=0, L = 2.6541;
  double T = 650, pi = 3.14159, DipoleAver2, V = L*L*L;
  fin = fopen ("induction.out", "r");
 // fout = fopen ("DielConst.out", "w");
  
for(i=0; i<numsteps; i++)
fscanf(fin, "%lf   %lf    %lf    %lf    %lf", &Induction[i][0], &Induction[i][1], &Induction[i][2], &Induction[i][3], &Induction[i][4]); 
 

//for(i=0;i<numsteps; i++)
//printf("%f   %f    %f    %f    %f\n", Induction[i][0], Induction[i][1], Induction[i][2], Induction[i][3], Induction[i][4]); 

for (i=0; i<numsteps; i++)
	{
	dipole += Induction[i][1];
	DipoleSum[0] += Induction[i][2];
	DipoleSum[1] += Induction[i][3];
	DipoleSum[2] += Induction[i][4];
	Dipole2Sum += Induction[i][2]*Induction[i][2] + Induction[i][3]*Induction[i][3] + Induction[i][4]*Induction[i][4];
	DipoleAver2 = (DipoleSum[0]/(numsteps))*(DipoleSum[0]/(numsteps)) + \
		    (DipoleSum[1]/(numsteps))*(DipoleSum[1]/(numsteps)) + \
		    (DipoleSum[2]/(numsteps))*(DipoleSum[2]/(numsteps));
        }

	dipole = dipole/numsteps;
	double Dipole2Aver = Dipole2Sum/numsteps;

double Epsilon = 880031.1376*(Dipole2Aver - DipoleAver2)/(T*V);
//	DielConst = 1 + (4*3.14159*Rho*m_i_permanent2*g_k)/(3*T*k_b);
//	DielConst = 1 + (4*pi*Rho*numMolec*m_i*m_i)/(3*T*k_b);
//printf ("Square of the total system dipole <M><M> = %f \n", DipoleAver2);
//printf ("Average of the square total system dipole <M^2>= %f \n", Dipole2Aver);
printf ("Dielectric constant = %f \n", Epsilon);
printf ("Dipole moment m = 2.19366D + m_induced = %f \n", dipole);
//fprintf ( fout, "Dielectric Constant = %f \n", DielConst);

fclose ( fin );
//fclose ( fout );

}
