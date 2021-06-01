#include <iostream>

#include "tnt_array2d.h"
using namespace TNT;

#include <iostream>     // std::cout

#include <cstddef>      // std::size_t

#include <cmath>        // std::atan2

#include <valarray>     // std::valarray, std::atan2#include <fstream>

#include <sstream>

#include <cstdlib>

#include <ctime>

#include <fstream>

#include <string>

#include <stdio.h>

#include "math.h"

using namespace std;

const double M_PI = 3.14159365368;

double L;
double stDev;
double b;
double angleUpperBound;
double dx;
double dy;
double dz;
double dL;
double maxDY;
double maxDZ;
double wu;
double matrixLength;
double matrixHeight;
double matrixWidth;
double thetaXY;
double thetaXZ;
double omegaXY;
double omegaXZ;
double phiXY;
double phiXZ;


int aa, n, j, ii, jj, i, inc, nsim, nCNT, nCNTElements;

double unif() {
	return rand() / double(RAND_MAX);
}

double unif(double a, double b) {
	return (b - a)*unif() + a;
}

long unif(long a) {
	if (a < 0) a = -a;
	if (a == 0) return 0;
	return long(unif()*a) + 1;
}

double randn() 
{
	static bool has_saved = false;
	static double saved;
	
	if (has_saved) 
	{
		has_saved = false;
		return saved;	
	}
	double x, y, r;
	do 
	{
		x = unif(-1., 1.);
		y = unif(-1., 1.);
		r = x * x + y * y;
		
	} while (r >= 1.);
		
	double mu = sqrt(-2.0 * log(r) / r);
		
	saved = mu * y;
	has_saved = true;
	return mu * x;
}

template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

double SDF1(double mr, double cl, double st)
{
	return 1. / 4. * st*st * cl*cl*cl * mr*mr * exp(-cl * mr);
}
// double SDF2(double p, double q, double r, double s)
// {
//	 if (r < 0.)
//	 {
//		 return 0;
//	 }
//	 double Sp = 0.;
//	 double dz = q / 1000.;
//	 double z=dz/2.;
//	 for (int ii =1; ii <= 1000; ii++)
//	 {
//		double nf = sqrt(2.) * s * cos(p * z + r);
//		Sp += nf*cos(p*z)*dz;
//		z += dz;
//	 }
//	 Sp = .5 / M_PI / q * Sp*Sp;
//	 return Sp;
//}
 double SDF3(double kx, double ky, double bx, double by, double s)
{ 
	 return .25 / M_PI * s*s * bx*by *exp(-.25*((bx*kx*bx*kx) + (by*ky*by*ky)));
}
int main()
{
	const char* input_file_name = "geometry_generation_parameters.dat";
	ifstream my_in(input_file_name);
	if (my_in.fail()) {
		cerr << "unable to open the file " << input_file_name
			<< " for input." << endl;
		return 1;
	}
	
	string str;
	//cout << "DEFINE NUMBER OF SIMULATIONS" << endl;
	my_in >> nsim;

	//cout << "DEFINE NUMBER OF CNT ELEMENTS" << endl;
	my_in >> nCNTElements;

	//cout << "DEFINE CNT LENGTH" << endl;
	my_in >> L;
	dL = L / double(nCNTElements);

	//cout << "DEFINE TOTAL CNT NUMBER" << endl;
	my_in >> nCNT;

	//cout << "DEFINE STOCHASTIC FIELD STANDARD DEVIATION" << endl;
	my_in >> stDev;
	stDev = abs(stDev);

	//cout << "DEFINE CORRELATION LENGTH PARAMETER" << endl;
	my_in >> b;

	//cout << "DEFINE ANGLE UPPER BOUND ABSOLUTE VALUE" << endl;
	my_in >> angleUpperBound;
	angleUpperBound = abs(angleUpperBound);

	//cout << "DEFINE MATRIX DIMENSIONS" << endl;
	my_in >> matrixLength;

	my_in >> matrixHeight; // 

	my_in >> matrixWidth; // 

	Array2D<double> localCNTcoordinates(nCNTElements + 1, 3);
	Array2D<double> globalCNTcoordinates(nCNT * (nCNTElements + 1), 3);
	const int num = 200;
	double ku = 1.22718464 * 10/(0.8*b);
	double dk = ku / float(num);
	double *phi1 = new double[num];
	double* phi2 = new double[num];

	for (int isim = 1; isim <= nsim; isim++)
	{

		string s = to_string(isim);
		string output_file_name_1, output_file_name_2;
		output_file_name_1 = "nodes_.txt";
		output_file_name_1.insert(6, s);
		ofstream my_out_1(output_file_name_1);
		output_file_name_2 = "connectivity.txt";
		ofstream my_out_2(output_file_name_2);

		for (int i = 1; i <= nCNT; i++)
		{

			for (j = 1; j <= nCNTElements + 1; j++)
			{

				for (jj = 1; jj <= 3; jj++)
				{
					globalCNTcoordinates[(i - 1) * (nCNTElements + 1) + j - 1][jj - 1] = 0;
				}
			}

		}
		for (int iCNT = 1; iCNT <= nCNT; iCNT++)
		{
			
			for (j = 1; j <= nCNTElements + 1; j++)
			{
				for (jj = 1; jj <= 3; jj++)
				{
					localCNTcoordinates[j - 1][jj - 1] = 0;
				}
			}

			for (int i = 1; i <= num; i++)
			{
				phi1[i - 1] = unif(0., (2. * M_PI));
			}
			for (int i = 1; i <= num; i++)
			{
				phi2[i - 1] = unif(0., (2. * M_PI));
			}

			int iNode;
			int iNode0 = (iCNT - 1) * (nCNTElements + 1) + 1;		
			globalCNTcoordinates[iNode0 - 1][0] = unif(0., matrixLength);
			globalCNTcoordinates[iNode0 - 1][1] = unif(0., matrixHeight);
			globalCNTcoordinates[iNode0 - 1][2] = unif(0., matrixWidth);
			double x0 = globalCNTcoordinates[iNode0 - 1][0];
			double y0 = globalCNTcoordinates[iNode0 - 1][1];
			double z0 = globalCNTcoordinates[iNode0 - 1][2];


			for (int iCNTElement = 1; iCNTElement <= nCNTElements; iCNTElement++)
			{
				iNode = (iCNT - 1)*(nCNTElements + 1) + iCNTElement;
				if (iCNTElement == 1)
				{
					thetaXY = unif(-M_PI,  M_PI);
					thetaXZ = unif(-M_PI, M_PI);
					dx = cos(thetaXY) * cos(thetaXZ) * dL;
					dy = sin(thetaXY) * cos(thetaXZ) * dL;
					dz = sin(thetaXZ) * dL;
				}
				else
				{
					thetaXY = 0;
					thetaXZ = 0;
					double k = dk / 2.;
					double x = (double(iCNTElement)-1) * dL;
					for (int i = 1; i <= num; i++)
					{
						thetaXY += sqrt(2.) * sqrt(2. * SDF1(k, b, stDev) * dk) * double(cos(k * x + phi1[i - 1]));
						thetaXZ += sqrt(2.) * sqrt(2. * SDF1(k, b, stDev) * dk) * double(cos(k * x + phi2[i - 1]));
						k += dk;
					}

					if (abs(thetaXY) > angleUpperBound) thetaXY = sgn(thetaXY) * angleUpperBound;
					if (abs(thetaXZ) > angleUpperBound) thetaXZ = sgn(thetaXZ) * angleUpperBound;
					phiXY = atan2(localCNTcoordinates[iCNTElement- 1][1] - localCNTcoordinates[iCNTElement - 2][1],
						localCNTcoordinates[iCNTElement - 1][0] - localCNTcoordinates[iCNTElement - 2][0]);
					phiXZ = atan2(localCNTcoordinates[iCNTElement - 1][2] - localCNTcoordinates[iCNTElement - 2][2]
						,sqrt(pow(localCNTcoordinates[iCNTElement - 1][0] - localCNTcoordinates[iCNTElement - 2][0],2)
							+pow(localCNTcoordinates[iCNTElement - 1][1] - localCNTcoordinates[iCNTElement - 2][1],2)));
					dx = cos(thetaXY + phiXY) * cos(thetaXZ + phiXZ) * dL;
					dy = sin(thetaXY + phiXY) * cos(thetaXZ + phiXZ) * dL;
					dz = sin(thetaXZ + phiXZ) * dL;
				}
				double dist = sqrt(dx * dx + dy * dy + dz * dz);
				localCNTcoordinates[iCNTElement][0] = localCNTcoordinates[iCNTElement - 1][0] + dx;
				localCNTcoordinates[iCNTElement][1] = localCNTcoordinates[iCNTElement - 1][1] + dy;
				localCNTcoordinates[iCNTElement][2] = localCNTcoordinates[iCNTElement - 1][2] + dz;
				globalCNTcoordinates[iNode][0] = globalCNTcoordinates[iNode0 - 1][0] + localCNTcoordinates[iCNTElement][0];
				globalCNTcoordinates[iNode][1] = globalCNTcoordinates[iNode0 - 1][1] + localCNTcoordinates[iCNTElement][1];
				globalCNTcoordinates[iNode][2] = globalCNTcoordinates[iNode0 - 1][2] + localCNTcoordinates[iCNTElement][2];

				if ((globalCNTcoordinates[iNode][0] > matrixLength) || (globalCNTcoordinates[iNode][0] < 0.)
					|| (globalCNTcoordinates[iNode][1] > matrixHeight) || (globalCNTcoordinates[iNode][1] < 0.)
					|| (globalCNTcoordinates[iNode][2] > matrixWidth) || (globalCNTcoordinates[iNode][2] < 0.))
				{
					iCNTElement = 0; 
					globalCNTcoordinates[iNode0 - 1][0] = unif(0., matrixLength);
					globalCNTcoordinates[iNode0 - 1][1] = unif(0., matrixHeight);
					globalCNTcoordinates[iNode0 - 1][2] = unif(0., matrixWidth);
					double x0 = globalCNTcoordinates[iNode0 - 1][0];
					double y0 = globalCNTcoordinates[iNode0 - 1][1];
					double z0 = globalCNTcoordinates[iNode0 - 1][2];
					continue;
				}

			}
			my_out_1 << iNode0 << ", " << globalCNTcoordinates[iNode0 - 1][0] << ", " << globalCNTcoordinates[iNode0 - 1][1]
				<< ", " << globalCNTcoordinates[iNode0 - 1][2] << endl;
			for (int iCNTElement = 1; iCNTElement <= nCNTElements; iCNTElement++)
			{
				iNode = (iCNT - 1) * (nCNTElements + 1) + iCNTElement;
				my_out_1 << iNode + 1 << ", " << globalCNTcoordinates[iNode][0] << ", " << globalCNTcoordinates[iNode][1] << ", " <<
					globalCNTcoordinates[iNode][2] << endl;
				my_out_2 << (iCNT - 1) * nCNTElements + iCNTElement << ", " << iNode << ", " << iNode + 1 << endl;
			}
		}
	}
	delete[] phi1, phi2;
	return 0;

}