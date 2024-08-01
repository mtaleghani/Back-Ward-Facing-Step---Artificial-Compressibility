#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>
#include <math.h>
using namespace std;

#define  imax 141 //30 and 42
#define  jmax 51

// BLOCK TRIDIAG SUBROUTINE FOR X-SWEEP
const int IM = imax - 2; // number of points to solve using block tridiag FOR X-SWEEP

void solve1x(long double P[][3][IM + 1], long double Q[][3][IM + 1], long double R[][3][IM + 1], int N)
{
	long double det = 0, Temp[3][3];
	det = (P[0][0][N] * P[1][1][N] * P[2][2][N] + P[0][1][N] * P[1][2][N]
		* P[2][0][N] + P[0][2][N] * P[1][0][N] * P[2][1][N]) -
		(P[0][2][N] * P[1][1][N] * P[2][0][N] + P[0][1][N] * P[1][0][N]
			* P[2][2][N] + P[0][0][N] * P[1][2][N] * P[2][1][N]);
	det = 1 / det;
	Temp[0][0] = det * (P[1][1][N] * P[2][2][N] - P[1][2][N] * P[2][1][N]);
	Temp[1][0] = -det * (P[1][0][N] * P[2][2][N] - P[1][2][N] * P[2][0][N]);
	Temp[2][0] = det * (P[1][0][N] * P[2][1][N] - P[1][1][N] * P[2][0][N]);
	Temp[0][1] = -det * (P[0][1][N] * P[2][2][N] - P[0][2][N] * P[2][1][N]);
	Temp[1][1] = det * (P[0][0][N] * P[2][2][N] - P[0][2][N] * P[2][0][N]);
	Temp[2][1] = -det * (P[0][0][N] * P[2][1][N] - P[0][1][N] * P[2][0][N]);
	Temp[0][2] = det * (P[0][1][N] * P[1][2][N] - P[0][2][N] * P[1][1][N]);
	Temp[1][2] = -det * (P[0][0][N] * P[1][2][N] - P[0][2][N] * P[1][0][N]);
	Temp[2][2] = det * (P[0][0][N] * P[1][1][N] - P[0][1][N] * P[1][0][N]);
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
		{
			Q[i][j][N] = 0;
			for (int k = 0; k <= 2; k++)
				Q[i][j][N] += Temp[i][k] * R[k][j][N];
		}
}

void solve2x(long double P[][3][IM + 1], long double Q[][IM + 1], long double R[][IM + 1], int N)
{
	long double det = 0, Temp[3][3];
	det = (P[0][0][N] * P[1][1][N] * P[2][2][N] + P[0][1][N] * P[1][2][N]
		* P[2][0][N] + P[0][2][N] * P[1][0][N] * P[2][1][N]) -
		(P[0][2][N] * P[1][1][N] * P[2][0][N] + P[0][1][N] * P[1][0][N]
			* P[2][2][N] + P[0][0][N] * P[1][2][N] * P[2][1][N]);
	det = 1 / det;
	Temp[0][0] = det * (P[1][1][N] * P[2][2][N] - P[1][2][N] * P[2][1][N]);
	Temp[1][0] = -det * (P[1][0][N] * P[2][2][N] - P[1][2][N] * P[2][0][N]);
	Temp[2][0] = det * (P[1][0][N] * P[2][1][N] - P[1][1][N] * P[2][0][N]);
	Temp[0][1] = -det * (P[0][1][N] * P[2][2][N] - P[0][2][N] * P[2][1][N]);
	Temp[1][1] = det * (P[0][0][N] * P[2][2][N] - P[0][2][N] * P[2][0][N]);
	Temp[2][1] = -det * (P[0][0][N] * P[2][1][N] - P[0][1][N] * P[2][0][N]);
	Temp[0][2] = det * (P[0][1][N] * P[1][2][N] - P[0][2][N] * P[1][1][N]);
	Temp[1][2] = -det * (P[0][0][N] * P[1][2][N] - P[0][2][N] * P[1][0][N]);
	Temp[2][2] = det * (P[0][0][N] * P[1][1][N] - P[0][1][N] * P[1][0][N]);
	for (int i = 0; i <= 2; i++)
	{
		Q[i][N] = 0;
		for (int j = 0; j <= 2; j++)
			Q[i][N] += Temp[i][j] * R[j][N];
	}
}

void BTRIDx(long double aa[][3][IM + 1], long double bb[][3][IM + 1], long double cc[][3][IM + 1], long double dd[][IM + 1], long double x[][IM + 1])
{
	long double m[3][3][IM + 1] = { 0 }, n[3][3][IM + 1] = { 0 }, y[3][IM + 1] = { 0 };
	int N = IM;
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
			m[i][j][1] = bb[i][j][1];
	solve1x(m, n, cc, 1);
	solve2x(m, y, dd, 1);
	for (int ND = 2; ND <= N - 1; ND++)
	{
		for (int i = 0; i <= 2; i++)
			for (int j = 0; j <= 2; j++)
			{
				for (int k = 0; k <= 2; k++)
					m[i][j][ND] += -aa[i][k][ND] * n[k][j][ND - 1];
				m[i][j][ND] += bb[i][j][ND];
			}
		solve1x(m, n, cc, ND);
	}

	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
		{
			for (int k = 0; k <= 2; k++)
				m[i][j][N] += -aa[i][k][N] * n[k][j][N - 1];
			m[i][j][N] += bb[i][j][N];
		}
	for (int ND = 2; ND <= N; ND++)
	{
		for (int i = 0; i <= 2; i++)
			for (int j = 0; j <= 2; j++)
				dd[i][ND] += -aa[i][j][ND] * y[j][ND - 1];
		solve2x(m, y, dd, ND);
	}
	for (int i = 0; i <= 2; i++)
		x[i][N] = y[i][N];
	for (int ND = N - 1; ND >= 1; ND--)
		for (int i = 0; i <= 2; i++)
		{
			for (int j = 0; j <= 2; j++)
				y[i][ND] += -n[i][j][ND] * x[j][ND + 1];
			x[i][ND] = y[i][ND];
		}
}


// BLOCK TRIDIAG SUBROUTINE FOR Y-SWEEP
const int IN = jmax - 2; // number of points to solve using block tridiag FOR Y-SWEEP

void solve1y(long double P[][3][IN + 1], long double Q[][3][IN + 1], long double R[][3][IN + 1], int N)
{
	long double det = 0, Temp[3][3];
	det = (P[0][0][N] * P[1][1][N] * P[2][2][N] + P[0][1][N] * P[1][2][N]
		* P[2][0][N] + P[0][2][N] * P[1][0][N] * P[2][1][N]) -
		(P[0][2][N] * P[1][1][N] * P[2][0][N] + P[0][1][N] * P[1][0][N]
			* P[2][2][N] + P[0][0][N] * P[1][2][N] * P[2][1][N]);
	det = 1 / det;
	Temp[0][0] = det * (P[1][1][N] * P[2][2][N] - P[1][2][N] * P[2][1][N]);
	Temp[1][0] = -det * (P[1][0][N] * P[2][2][N] - P[1][2][N] * P[2][0][N]);
	Temp[2][0] = det * (P[1][0][N] * P[2][1][N] - P[1][1][N] * P[2][0][N]);
	Temp[0][1] = -det * (P[0][1][N] * P[2][2][N] - P[0][2][N] * P[2][1][N]);
	Temp[1][1] = det * (P[0][0][N] * P[2][2][N] - P[0][2][N] * P[2][0][N]);
	Temp[2][1] = -det * (P[0][0][N] * P[2][1][N] - P[0][1][N] * P[2][0][N]);
	Temp[0][2] = det * (P[0][1][N] * P[1][2][N] - P[0][2][N] * P[1][1][N]);
	Temp[1][2] = -det * (P[0][0][N] * P[1][2][N] - P[0][2][N] * P[1][0][N]);
	Temp[2][2] = det * (P[0][0][N] * P[1][1][N] - P[0][1][N] * P[1][0][N]);
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
		{
			Q[i][j][N] = 0;
			for (int k = 0; k <= 2; k++)
				Q[i][j][N] += Temp[i][k] * R[k][j][N];
		}
}

void solve2y(long double P[][3][IN + 1], long double Q[][IN + 1], long double R[][IN + 1], int N)
{
	long double det = 0, Temp[3][3];
	det = (P[0][0][N] * P[1][1][N] * P[2][2][N] + P[0][1][N] * P[1][2][N]
		* P[2][0][N] + P[0][2][N] * P[1][0][N] * P[2][1][N]) -
		(P[0][2][N] * P[1][1][N] * P[2][0][N] + P[0][1][N] * P[1][0][N]
			* P[2][2][N] + P[0][0][N] * P[1][2][N] * P[2][1][N]);
	det = 1 / det;
	Temp[0][0] = det * (P[1][1][N] * P[2][2][N] - P[1][2][N] * P[2][1][N]);
	Temp[1][0] = -det * (P[1][0][N] * P[2][2][N] - P[1][2][N] * P[2][0][N]);
	Temp[2][0] = det * (P[1][0][N] * P[2][1][N] - P[1][1][N] * P[2][0][N]);
	Temp[0][1] = -det * (P[0][1][N] * P[2][2][N] - P[0][2][N] * P[2][1][N]);
	Temp[1][1] = det * (P[0][0][N] * P[2][2][N] - P[0][2][N] * P[2][0][N]);
	Temp[2][1] = -det * (P[0][0][N] * P[2][1][N] - P[0][1][N] * P[2][0][N]);
	Temp[0][2] = det * (P[0][1][N] * P[1][2][N] - P[0][2][N] * P[1][1][N]);
	Temp[1][2] = -det * (P[0][0][N] * P[1][2][N] - P[0][2][N] * P[1][0][N]);
	Temp[2][2] = det * (P[0][0][N] * P[1][1][N] - P[0][1][N] * P[1][0][N]);
	for (int i = 0; i <= 2; i++)
	{
		Q[i][N] = 0;
		for (int j = 0; j <= 2; j++)
			Q[i][N] += Temp[i][j] * R[j][N];
	}
}

void BTRIDy(long double aa[][3][IN + 1], long double bb[][3][IN + 1], long double cc[][3][IN + 1], long double dd[][IN + 1], long double x[][IN + 1])
{
	long double m[3][3][IN + 1] = { 0 }, n[3][3][IN + 1] = { 0 }, y[3][IN + 1] = { 0 };
	int N = IN;
	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
			m[i][j][1] = bb[i][j][1];
	solve1y(m, n, cc, 1);
	solve2y(m, y, dd, 1);
	for (int ND = 2; ND <= N - 1; ND++)
	{
		for (int i = 0; i <= 2; i++)
			for (int j = 0; j <= 2; j++)
			{
				for (int k = 0; k <= 2; k++)
					m[i][j][ND] += -aa[i][k][ND] * n[k][j][ND - 1];
				m[i][j][ND] += bb[i][j][ND];
			}
		solve1y(m, n, cc, ND);
	}

	for (int i = 0; i <= 2; i++)
		for (int j = 0; j <= 2; j++)
		{
			for (int k = 0; k <= 2; k++)
				m[i][j][N] += -aa[i][k][N] * n[k][j][N - 1];
			m[i][j][N] += bb[i][j][N];
		}
	for (int ND = 2; ND <= N; ND++)
	{
		for (int i = 0; i <= 2; i++)
			for (int j = 0; j <= 2; j++)
				dd[i][ND] += -aa[i][j][ND] * y[j][ND - 1];
		solve2y(m, y, dd, ND);
	}
	for (int i = 0; i <= 2; i++)
		x[i][N] = y[i][N];
	for (int ND = N - 1; ND >= 1; ND--)
		for (int i = 0; i <= 2; i++)
		{
			for (int j = 0; j <= 2; j++)
				y[i][ND] += -n[i][j][ND] * x[j][ND + 1];
			x[i][ND] = y[i][ND];
		}
}


double MatMult(long double A[imax][jmax][3][3], long double Q[imax][jmax][3], int i, int j, int n, int row)
{
	long double result;

	result = A[i][j][row][0] * Q[i][j][0] + A[i][j][row][1] * Q[i][j][1] + A[i][j][row][2] * Q[i][j][2];

	return result;
}

void Mesh(double x[imax], double y[jmax], double dx, double dy)
{
	for (int i = 0; i < imax; i++)
	{
		x[i] = i * dx;
	}

	for (int j = 0; j < jmax; j++)
	{
		y[j] = j * dy;
	}
}

void InitAndBC(double Q[imax][jmax][3], double x[imax], double y[jmax], double ubar, double H, double h)
{
	int i, j;

	// left Boundary
	i = 0;

	for (j = 0; j < jmax; j++)
	{
		Q[i][j][1 - 1] = 1.0; //p
		Q[i][j][3 - 1] = 0.0; //v
	}

	// u inlet
	for (j = 0; j < (jmax - 1) / 2; j++)
	{
		Q[i][j][2 - 1] = 0.0; //u
	}
	for (j = (jmax - 1) / 2; j < jmax; j++)
	{
		Q[i][j][2 - 1] = 6.0 * ubar * (H - y[j]) * (y[j] - h) / pow(H - h, 2); //u
	}

	// interior
	for (i = 1; i < imax - 1; i++)
	{
		for (j = 0; j < jmax; j++)
		{
			Q[i][j][1 - 1] = 1.0; // p
			Q[i][j][2 - 1] = 0.0;// Q[i - 1][j][2 - 1]; // u inlet as initial guess
			Q[i][j][3 - 1] = 0.0; // v
		}
	}

	// right Boundary
	i = imax - 1;
	for (j = 0; j < jmax; j++)
	{
		Q[i][j][1 - 1] = 1.0; //POUT
		Q[i][j][2 - 1] = 0.0;// Q[i - 1][j][2 - 1];// u inlet as initial guess
		Q[i][j][3 - 1] = 0.0; //v
	}
}

// WARNING WARNING WARNING WARNING
void TecplotOutput(const char* Name, const char* var, double Q[imax][3], double x[imax], int row)
{
	ofstream MyFile(Name);

	/* Tecplot Output*/

	/* Header */
	MyFile << "Variables = x , " << var << "\n";

	/* Header */

	/* Data */
	MyFile << "Zone T = \"Test\" , I = " << imax << "\n";

	for (int i = 0; i < imax; i++)
	{
		MyFile << x[i];
		MyFile << "		";
		MyFile << Q[i][row];
		MyFile << "\n";
	}

	/* Data */

	MyFile.close();
}

void TecplotOutputAll(const char* Name, double Q[imax][jmax][3], double x[imax], double y[jmax])
{
	FILE* stream;
	fopen_s(&stream, Name, "w");

	/* Tecplot Output*/

	/* Header */
	fprintf_s(stream, "Variables = x, y, p, u, v\n");

	/* Header */

	/* Data */
	fprintf_s(stream, "Zone T = \"Test\" , I = %d , J = %d\n", imax, jmax);

	for (int j = 0; j < jmax; j++)
	{
		for (int i = 0; i < imax; i++)
		{
			fprintf_s(stream, "%-17.10f", x[i]);
			fprintf_s(stream, "		");
			fprintf_s(stream, "%-17.10f", y[j]);
			fprintf_s(stream, "		");
			fprintf_s(stream, "%-17.10f", Q[i][j][0]);
			fprintf_s(stream, "		");
			fprintf_s(stream, "%-17.10f", Q[i][j][1]);
			fprintf_s(stream, "		");
			fprintf_s(stream, "%-17.10f", Q[i][j][2]);
			fprintf_s(stream, "			");
			fprintf_s(stream, "\n");
		}
	}

	/* Data */

	fclose(stream);
}

double maxvec(double arr[], int n)
{
	int i;

	// Initialize maximum element
	long double max = arr[0];

	// Traverse array elements
	// from second and compare
	// every element with current max
	for (i = 1; i < n; i++)
		if (arr[i] > max)
			max = arr[i];

	return max;
}

double maxarr(double arr[imax][jmax], int n, int m)
{
	int i, j;

	// Initialize maximum element
	long double max = arr[0][0];

	// Traverse array elements
	// from second and compare
	// every element with current max
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			if (arr[i][j] > max)
			{
				max = arr[i][j];
			}
		}
	}
	return max;
}

double maxmat(double arr[imax][jmax][3], int n, int m, int row)
{
	int i, j;

	// Initialize maximum element
	long double max = arr[0][0][row];

	// Traverse array elements
	// from second and compare
	// every element with current max
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < m; j++)
		{
			if (arr[i][j][row] > max)
			{
				max = arr[i][j][row];
			}
		}
	}
	return max;
}

int main()
{
	double xmax = 30.0;
	double ymax = 1.0;

	// problem parameters
	double H = ymax, h = H / 2.0, ubar = 1.0;

	double cfl = 1.0, beta = 8.0, Re = 1., err = 1e-7, epsi = 0.0, epse = 0.0, L_inft0 = 0.0, L_inft1 = 0.0, L_inft2 = 0.0, L_inft3 = 0.0, divV = 0.0, umax = 0.0, vmax = 0.0;
	double dx = (xmax) / (imax - 1.0), dy = (ymax) / (jmax - 1.0), dt = 0.1;
	double x[imax] = { 0.0 }, y[jmax] = { 0.0 }, diff0[imax][jmax] = { 0.0 }, diff1[imax][jmax] = { 0.0 }, diff2[imax][jmax] = { 0.0 }, DivV[imax][jmax][3] = { 0.0 }, F[imax][jmax][3] = { 0.0 }, G[imax][jmax][3] = { 0.0 };

	double Dissip[imax][jmax][3] = { 0.0 }, dQ[imax][jmax][3] = { 0.0 }, dQstar[imax][jmax][3] = { 0.0 }, Q[imax][jmax][3] = { 0.0 };

	// Coeff. Matrixes
	double Dbar[imax][jmax][3] = { 0.0 }, D2bar[imax][jmax][3] = { 0.0 };
	double A[imax][jmax][3][3] = { 0.0 }, B[imax][jmax][3][3] = { 0.0 }, C[imax][jmax][3][3] = { 0.0 }, N[imax][jmax][3][3] = { 0.0 },
	   	Abar[imax][jmax][3][3] = { 0.0 }, Bbar[imax][jmax][3][3] = { 0.0 }, Cbar[imax][jmax][3][3] = { 0.0 };

	int i, j;
	int iter;

	iter = 0;

	Mesh(x, y, dx, dy);
	InitAndBC(Q, x, y, ubar, H, h);
	TecplotOutputAll("BC.plt", Q, x, y);
	
	// Laplacian Terms in x and y directions.
	double Q2x[imax][jmax][3] = { 0.0 };
	double Q2y[imax][jmax][3] = { 0.0 };

	// 4th order Terms in x and y directions.
	double Q4x[imax][jmax][3] = { 0.0 };
	double Q4y[imax][jmax][3] = { 0.0 };

	long double AX[3][3][IM + 1], BX[3][3][IM + 1], CX[3][3][IM + 1], DX[3][IM + 1], XX[3][IM + 1];
	long double AY[3][3][IN + 1], BY[3][3][IN + 1], CY[3][3][IN + 1], DY[3][IN + 1], XY[3][IN + 1];



	FILE* convergence;
	fopen_s(&convergence,"convergence.plt", "w");
	/* Tecplot Output*/
	/* Header */

	fprintf_s(convergence, "Variables = Iter, \"L<sub>2</sub> Norm\", \"L<sub>2</sub> Prev\"\n");
	fprintf_s(convergence, "Zone T = \"Convergence History of div(V)\"");
	fprintf_s(convergence, "\n");


	while (1)
	{

		// update Fluxes and Matrixes
		for (i = 0; i < imax; i++)
		{
			for (j = 0; j < jmax; j++)
			{
				// Fluxes
				F[i][j][1 - 1] = beta * Q[i][j][2 - 1];
				F[i][j][2 - 1] = pow(Q[i][j][2 - 1], 2) + Q[i][j][1 - 1];
				F[i][j][3 - 1] = Q[i][j][2 - 1] * Q[i][j][3 - 1];

				G[i][j][1 - 1] = beta * Q[i][j][3 - 1];
				G[i][j][2 - 1] = Q[i][j][2 - 1] * Q[i][j][3 - 1];
				G[i][j][3 - 1] = pow(Q[i][j][3 - 1], 2) + Q[i][j][1 - 1];
				//

				// Jacobian Matrixes
				A[i][j][0][0] = 0.0;		A[i][j][0][1] = beta;						A[i][j][0][2] = 0.0;
				A[i][j][1][0] = 1.0;		A[i][j][1][1] = 2.0 * Q[i][j][2 - 1]; 		A[i][j][1][2] = 0.0;
				A[i][j][2][0] = 0.0;		A[i][j][2][1] = Q[i][j][3 - 1]; 			A[i][j][2][2] = Q[i][j][2 - 1];

				B[i][j][0][0] = 0.0;		B[i][j][0][1] = 0.0;						B[i][j][0][2] = beta;
				B[i][j][1][0] = 0.0;		B[i][j][1][1] = Q[i][j][3 - 1]; 			B[i][j][1][2] = Q[i][j][2 - 1];
				B[i][j][2][0] = 1.0;		B[i][j][2][1] = 0.0; 						B[i][j][2][2] = 2.0 * Q[i][j][3 - 1];
				//
			}
		}

		/*umax = abs(maxmat(Q, imax, jmax, 1));
		vmax = abs(maxmat(Q, imax, jmax, 2));
		
		double dty = cfl * dy / (vmax + sqrt(vmax * vmax + beta));
		double dtx = cfl * dx / (umax + sqrt(umax * umax + beta));

		if (dty >= dtx)
		{
			dt = dtx;
		}
		else
		{
			dt = dty;
		}*/

		//dt = cfl * dy / (1.8 + sqrt(1.8 * 1.8 + beta));


	  // Dissipation terms ( zero for i,j = {0,1} , i = {imax-1,imax} , j = {jmax-1,max} )

		for (i = 2; i < imax - 2; i++)          // for (i = 0; i < imax; i++)                              
		{
			for (j = 2; j < jmax - 2; j++)		// for (j = 0; j < jmax; j++)

			{

				// commented since already satisfied by Dissip[imax][jmax][3]

				//if ((j < 2 & i < 2) || (j < 2 & i > imax - 3) || (j > jmax - 3 & i > imax - 3) || (j > jmax - 3 & i < 2))
				//{
				//	Dissip[i][j][0] = 0.0; //epse;
				//	Dissip[i][j][1] = 0.0; //epse;
				//	Dissip[i][j][2] = 0.0; //epse;
				//}

				//else {
					Q4x[i][j][0] = Q[i - 2][j][0] - 4.0 * Q[i - 1][j][0] + 6.0 * Q[i][j][0] - 4.0 * Q[i + 1][j][0] + Q[i + 2][j][0];
					Q4x[i][j][1] = Q[i - 2][j][1] - 4.0 * Q[i - 1][j][1] + 6.0 * Q[i][j][1] - 4.0 * Q[i + 1][j][1] + Q[i + 2][j][1];
					Q4x[i][j][2] = Q[i - 2][j][2] - 4.0 * Q[i - 1][j][2] + 6.0 * Q[i][j][2] - 4.0 * Q[i + 1][j][2] + Q[i + 2][j][2];

					Q4y[i][j][0] = Q[i][j - 2][0] - 4.0 * Q[i][j - 1][0] + 6.0 * Q[i][j][0] - 4.0 * Q[i][j + 1][0] + Q[i][j + 2][0];
					Q4y[i][j][1] = Q[i][j - 2][1] - 4.0 * Q[i][j - 1][1] + 6.0 * Q[i][j][1] - 4.0 * Q[i][j + 1][1] + Q[i][j + 2][1];
					Q4y[i][j][2] = Q[i][j - 2][2] - 4.0 * Q[i][j - 1][2] + 6.0 * Q[i][j][2] - 4.0 * Q[i][j + 1][2] + Q[i][j + 2][2];

					Dissip[i][j][0] = (Q4x[i][j][0] + Q4y[i][j][0]) * epse;
					Dissip[i][j][1] = (Q4x[i][j][1] + Q4y[i][j][1]) * epse;
					Dissip[i][j][2] = (Q4x[i][j][2] + Q4y[i][j][2]) * epse;
				//}
			}
		}

		// Filling the COEFF MATRIXES
		for (i = 1; i < imax - 1; i++)
		{
			for (j = 1; j < jmax - 1; j++)
			{
				Abar[i][j][0][0] = -A[i - 1][j][0][0] * dt / (2.0 * dx) + epsi;
				Abar[i][j][0][1] = -A[i - 1][j][0][1] * dt / (2.0 * dx) + epsi;
				Abar[i][j][0][2] = -A[i - 1][j][0][2] * dt / (2.0 * dx) + epsi;
				Abar[i][j][1][0] = -A[i - 1][j][1][0] * dt / (2.0 * dx) + epsi;
				Abar[i][j][1][1] = -A[i - 1][j][1][1] * dt / (2.0 * dx) + epsi - (1.0 / Re) * dt / pow(dx, 2);
				Abar[i][j][1][2] = -A[i - 1][j][1][2] * dt / (2.0 * dx) + epsi;
				Abar[i][j][2][0] = -A[i - 1][j][2][0] * dt / (2.0 * dx) + epsi;
				Abar[i][j][2][1] = -A[i - 1][j][2][1] * dt / (2.0 * dx) + epsi;
				Abar[i][j][2][2] = -A[i - 1][j][2][2] * dt / (2.0 * dx) + epsi - (1.0 / Re) * dt / pow(dx, 2);

				Bbar[i][j][0][0] = 1.0 - 2.0 * epsi;
				Bbar[i][j][0][1] = -2.0 * epsi;
				Bbar[i][j][0][2] = -2.0 * epsi;
				Bbar[i][j][1][0] = -2.0 * epsi;
				Bbar[i][j][1][1] = 1.0 + (2.0 / Re) * dt / pow(dx, 2) - 2.0 * epsi;
				Bbar[i][j][1][2] = -2.0 * epsi;
				Bbar[i][j][2][0] = -2.0 * epsi;
				Bbar[i][j][2][1] = -2.0 * epsi;
				Bbar[i][j][2][2] = 1.0 + (2.0 / Re) * dt / pow(dx, 2) - 2.0 * epsi;

				Cbar[i][j][0][0] = A[i + 1][j][0][0] * dt / (2.0 * dx) + epsi;
				Cbar[i][j][0][1] = A[i + 1][j][0][1] * dt / (2.0 * dx) + epsi;
				Cbar[i][j][0][2] = A[i + 1][j][0][2] * dt / (2.0 * dx) + epsi;
				Cbar[i][j][1][0] = A[i + 1][j][1][0] * dt / (2.0 * dx) + epsi;
				Cbar[i][j][1][1] = A[i + 1][j][1][1] * dt / (2.0 * dx) + epsi - (1.0 / Re) * dt / pow(dx, 2);
				Cbar[i][j][1][2] = A[i + 1][j][1][2] * dt / (2.0 * dx) + epsi;
				Cbar[i][j][2][0] = A[i + 1][j][2][0] * dt / (2.0 * dx) + epsi;
				Cbar[i][j][2][1] = A[i + 1][j][2][1] * dt / (2.0 * dx) + epsi;
				Cbar[i][j][2][2] = A[i + 1][j][2][2] * dt / (2.0 * dx) + epsi - (1.0 / Re) * dt / pow(dx, 2);



				//Q2x[i][j][0] = (Q[i + 1][j][0] - 2.0 * Q[i][j][0] + Q[i - 1][j][0]);
				Q2x[i][j][1] = (Q[i + 1][j][1] - 2.0 * Q[i][j][1] + Q[i - 1][j][1]);
				Q2x[i][j][2] = (Q[i + 1][j][2] - 2.0 * Q[i][j][2] + Q[i - 1][j][2]);

				//Q2y[i][j][0] = (Q[i][j + 1][0] - 2.0 * Q[i][j][0] + Q[i][j - 1][0]);
				Q2y[i][j][1] = (Q[i][j + 1][1] - 2.0 * Q[i][j][1] + Q[i][j - 1][1]);
				Q2y[i][j][2] = (Q[i][j + 1][2] - 2.0 * Q[i][j][2] + Q[i][j - 1][2]);

				Dbar[i][j][0] = (-(F[i + 1][j][0] - F[i - 1][j][0]) / (2.0 * dx) - (G[i][j + 1][0] - G[i][j - 1][0]) / (2.0 * dy)) * dt - Dissip[i][j][0];
				Dbar[i][j][1] = (-(F[i + 1][j][1] - F[i - 1][j][1]) / (2.0 * dx) - (G[i][j + 1][1] - G[i][j - 1][1]) / (2.0 * dy) + (Q2x[i][j][1] / (dx * dx) + Q2y[i][j][1] / (dy * dy)) * 1.0 / Re) * dt - Dissip[i][j][1];
				Dbar[i][j][2] = (-(F[i + 1][j][2] - F[i - 1][j][2]) / (2.0 * dx) - (G[i][j + 1][2] - G[i][j - 1][2]) / (2.0 * dy) + (Q2x[i][j][2] / (dx * dx) + Q2y[i][j][2] / (dy * dy)) * 1.0 / Re) * dt - Dissip[i][j][2];
			}
		}

		// x-sweep starts
		for (j = 1; j < jmax - 1; j++)
		{
			//Filling the BLOCK MATRIXES.

			// Modified Bbar and Cbar for LEFT INLET Boundary conditions.
			i = 1;

			BX[0][0][i] = Bbar[i][j][0][0] + Abar[i][j][0][0] * (4.0 / 3.0);		BX[0][1][i] = Bbar[i][j][0][1];			BX[0][2][i] = Bbar[i][j][0][2];
																					 										 
			BX[1][0][i] = Bbar[i][j][1][0] + Abar[i][j][1][0] * (4.0 / 3.0);		BX[1][1][i] = Bbar[i][j][1][1];			BX[1][2][i] = Bbar[i][j][1][2];
																					 										 
			BX[2][0][i] = Bbar[i][j][2][0] + Abar[i][j][2][0] * (4.0 / 3.0);		BX[2][1][i] = Bbar[i][j][2][1];			BX[2][2][i] = Bbar[i][j][2][2];
																					 										
																					 										
			CX[0][0][i] = Cbar[i][j][0][0] - Abar[i][j][0][0] / 3.0;				CX[0][1][i] = Cbar[i][j][0][1];			CX[0][2][i] = Cbar[i][j][0][2];
			 																		 										 
			CX[1][0][i] = Cbar[i][j][1][0] - Abar[i][j][1][0] / 3.0;				CX[1][1][i] = Cbar[i][j][1][1];			CX[1][2][i] = Cbar[i][j][1][2];
			 																		 										 
			CX[2][0][i] = Cbar[i][j][2][0] - Abar[i][j][2][0] / 3.0;				CX[2][1][i] = Cbar[i][j][2][1];			CX[2][2][i] = Cbar[i][j][2][2];


			DX[0][i] = Dbar[i][j][0];
			DX[1][i] = Dbar[i][j][1];
			DX[2][i] = Dbar[i][j][2];

			for (i = 2; i < imax - 2; i++)
			{
				AX[0][0][i] = Abar[i][j][0][0];
				AX[0][1][i] = Abar[i][j][0][1];
				AX[0][2][i] = Abar[i][j][0][2];
				AX[1][0][i] = Abar[i][j][1][0];
				AX[1][1][i] = Abar[i][j][1][1];
				AX[1][2][i] = Abar[i][j][1][2];
				AX[2][0][i] = Abar[i][j][2][0];
				AX[2][1][i] = Abar[i][j][2][1];
				AX[2][2][i] = Abar[i][j][2][2];

				BX[0][0][i] = Bbar[i][j][0][0];
				BX[0][1][i] = Bbar[i][j][0][1];
				BX[0][2][i] = Bbar[i][j][0][2];
				BX[1][0][i] = Bbar[i][j][1][0];
				BX[1][1][i] = Bbar[i][j][1][1];
				BX[1][2][i] = Bbar[i][j][1][2];
				BX[2][0][i] = Bbar[i][j][2][0];
				BX[2][1][i] = Bbar[i][j][2][1];
				BX[2][2][i] = Bbar[i][j][2][2];

				CX[0][0][i] = Cbar[i][j][0][0];
				CX[0][1][i] = Cbar[i][j][0][1];
				CX[0][2][i] = Cbar[i][j][0][2];
				CX[1][0][i] = Cbar[i][j][1][0];
				CX[1][1][i] = Cbar[i][j][1][1];
				CX[1][2][i] = Cbar[i][j][1][2];
				CX[2][0][i] = Cbar[i][j][2][0];
				CX[2][1][i] = Cbar[i][j][2][1];
				CX[2][2][i] = Cbar[i][j][2][2];

				DX[0][i] = Dbar[i][j][0];
				DX[1][i] = Dbar[i][j][1];
				DX[2][i] = Dbar[i][j][2];
			}

			// Modified Abar and Bbar for RIGHT OUTLET Boundary conditions.
			i = (imax - 1) - 1;

			AX[0][0][i] = Abar[i][j][0][0];			AX[0][1][i] = Abar[i][j][0][1] - Cbar[i][j][0][1] / 3.0;			AX[0][2][i] = Abar[i][j][0][2] - Cbar[i][j][0][2] / 3.0;
																														
			AX[1][0][i] = Abar[i][j][1][0];			AX[1][1][i] = Abar[i][j][1][1] - Cbar[i][j][1][1] / 3.0;			AX[1][2][i] = Abar[i][j][1][2] - Cbar[i][j][1][2] / 3.0;
																														
			AX[2][0][i] = Abar[i][j][2][0];			AX[2][1][i] = Abar[i][j][2][1] - Cbar[i][j][2][1] / 3.0;			AX[2][2][i] = Abar[i][j][2][2] - Cbar[i][j][2][2] / 3.0;
							

			BX[0][0][i] = Bbar[i][j][0][0];			BX[0][1][i] = Bbar[i][j][0][1] + 4.0 * Cbar[i][j][0][1] / 3.0;		BX[0][2][i] = Bbar[i][j][0][2] + 4.0 * Cbar[i][j][0][2] / 3.0;
																														
			BX[1][0][i] = Bbar[i][j][1][0];			BX[1][1][i] = Bbar[i][j][1][1] + 4.0 * Cbar[i][j][1][1] / 3.0;		BX[1][2][i] = Bbar[i][j][1][2] + 4.0 * Cbar[i][j][1][2] / 3.0;
																														
			BX[2][0][i] = Bbar[i][j][2][0];			BX[2][1][i] = Bbar[i][j][2][1] + 4.0 * Cbar[i][j][2][1] / 3.0;		BX[2][2][i] = Bbar[i][j][2][2] + 4.0 * Cbar[i][j][2][2] / 3.0;


			DX[0][i] = Dbar[i][j][0];
			DX[1][i] = Dbar[i][j][1];
			DX[2][i] = Dbar[i][j][2];

			// Solve Process Begins for x-sweep
			BTRIDx(AX, BX, CX, DX, XX);

			for (i = 1; i <= imax - 2; i++)
			{
				dQstar[i][j][0] = XX[0][i];
				dQstar[i][j][1] = XX[1][i];
				dQstar[i][j][2] = XX[2][i];

				//cout << XX[0][i] << "\n";
				//cout << XX[1][i] << "\n";
				//cout << XX[2][i] << "\n";
			}

			//TecplotOutputAll("dQstar.dat", dQstar, x, y);
			// Solve Ends

			// Step 1 BCs Completion

			//Left Wall and Inlet
			dQstar[0][j][1 - 1] = (4.0 * dQstar[1][j][1 - 1] - dQstar[2][j][1 - 1]) / 3.0;

			//Right Outlet
			dQstar[imax - 1][j][2 - 1] = (4.0 * dQstar[imax - 2][j][2 - 1] - dQstar[imax - 3][j][2 - 1]) / 3.0;
			dQstar[imax - 1][j][3 - 1] = (4.0 * dQstar[imax - 2][j][3 - 1] - dQstar[imax - 3][j][3 - 1]) / 3.0;
		}

		// end of x-sweep

		// Filling the COEFF MATRIXES
		for (i = 1; i < imax - 1; i++)
		{
			for (j = 1; j < jmax - 1; j++)
			{
				Abar[i][j][0][0] = -B[i][j - 1][0][0] * dt / (2.0 * dy) + epsi;
				Abar[i][j][0][1] = -B[i][j - 1][0][1] * dt / (2.0 * dy) + epsi;
				Abar[i][j][0][2] = -B[i][j - 1][0][2] * dt / (2.0 * dy) + epsi;
				Abar[i][j][1][0] = -B[i][j - 1][1][0] * dt / (2.0 * dy) + epsi;
				Abar[i][j][1][1] = -B[i][j - 1][1][1] * dt / (2.0 * dy) + epsi - (1.0 / Re) * dt / pow(dy, 2);
				Abar[i][j][1][2] = -B[i][j - 1][1][2] * dt / (2.0 * dy) + epsi;
				Abar[i][j][2][0] = -B[i][j - 1][2][0] * dt / (2.0 * dy) + epsi;
				Abar[i][j][2][1] = -B[i][j - 1][2][1] * dt / (2.0 * dy) + epsi;
				Abar[i][j][2][2] = -B[i][j - 1][2][2] * dt / (2.0 * dy) + epsi - (1.0 / Re) * dt / pow(dy, 2);

				Bbar[i][j][0][0] = 1.0 - 2.0 * epsi;
				Bbar[i][j][0][1] = -2.0 * epsi;
				Bbar[i][j][0][2] = -2.0 * epsi;
				Bbar[i][j][1][0] = -2.0 * epsi;
				Bbar[i][j][1][1] = 1.0 + (2.0 / Re) * dt / pow(dy, 2) - 2.0 * epsi;
				Bbar[i][j][1][2] = -2.0 * epsi;
				Bbar[i][j][2][0] = -2.0 * epsi;
				Bbar[i][j][2][1] = -2.0 * epsi;
				Bbar[i][j][2][2] = 1.0 + (2.0 / Re) * dt / pow(dy, 2) - 2.0 * epsi;

				Cbar[i][j][0][0] = B[i][j + 1][0][0] * dt / (2.0 * dy) + epsi;
				Cbar[i][j][0][1] = B[i][j + 1][0][1] * dt / (2.0 * dy) + epsi;
				Cbar[i][j][0][2] = B[i][j + 1][0][2] * dt / (2.0 * dy) + epsi;
				Cbar[i][j][1][0] = B[i][j + 1][1][0] * dt / (2.0 * dy) + epsi;
				Cbar[i][j][1][1] = B[i][j + 1][1][1] * dt / (2.0 * dy) + epsi - (1.0 / Re) * dt / pow(dy, 2);
				Cbar[i][j][1][2] = B[i][j + 1][1][2] * dt / (2.0 * dy) + epsi;
				Cbar[i][j][2][0] = B[i][j + 1][2][0] * dt / (2.0 * dy) + epsi;
				Cbar[i][j][2][1] = B[i][j + 1][2][1] * dt / (2.0 * dy) + epsi;
				Cbar[i][j][2][2] = B[i][j + 1][2][2] * dt / (2.0 * dy) + epsi - (1.0 / Re) * dt / pow(dy, 2);

				Dbar[i][j][0] = dQstar[i][j][0];
				Dbar[i][j][1] = dQstar[i][j][1];
				Dbar[i][j][2] = dQstar[i][j][2];
			}
		}

		//TecplotOutputAll("dQ.dat", dQstar, x, y);

		// y-sweep starts
		for (i = 1; i < imax - 1; i++)
		{
			//Filling the BLOCK MATRIXES.

			// Modified Bbar and Cbar for LOWER WALL Boundary conditions.
			j = 1;

			BY[0][0][j] = Bbar[i][j][0][0] + 4.0 / 3.0 * Abar[i][j][0][0];		BY[0][1][j] = Bbar[i][j][0][1];		BY[0][2][j] = Bbar[i][j][0][2];
																				 									 
			BY[1][0][j] = Bbar[i][j][1][0] + 4.0 / 3.0 * Abar[i][j][1][0];		BY[1][1][j] = Bbar[i][j][1][1];		BY[1][2][j] = Bbar[i][j][1][2];
																				 									 
			BY[2][0][j] = Bbar[i][j][2][0] + 4.0 / 3.0 * Abar[i][j][2][0];		BY[2][1][j] = Bbar[i][j][2][1];		BY[2][2][j] = Bbar[i][j][2][2];
																				 									 
																				 									 
			CY[0][0][j] = Cbar[i][j][0][0] - Abar[i][j][0][0] / 3.0;			CY[0][1][j] = Cbar[i][j][0][1];		CY[0][2][j] = Cbar[i][j][0][2];
			 																	 									 
			CY[1][0][j] = Cbar[i][j][1][0] - Abar[i][j][1][0] / 3.0;			CY[1][1][j] = Cbar[i][j][1][1];		CY[1][2][j] = Cbar[i][j][1][2];
			 																	 									 
			CY[2][0][j] = Cbar[i][j][2][0] - Abar[i][j][2][0] / 3.0;			CY[2][1][j] = Cbar[i][j][2][1];		CY[2][2][j] = Cbar[i][j][2][2];


			DY[0][j] = Dbar[i][j][0];
			DY[1][j] = Dbar[i][j][1];
			DY[2][j] = Dbar[i][j][2];

			for (j = 2; j < jmax - 2; j++)
			{
				AY[0][0][j] = Abar[i][j][0][0];
				AY[0][1][j] = Abar[i][j][0][1];
				AY[0][2][j] = Abar[i][j][0][2];
				AY[1][0][j] = Abar[i][j][1][0];
				AY[1][1][j] = Abar[i][j][1][1];
				AY[1][2][j] = Abar[i][j][1][2];
				AY[2][0][j] = Abar[i][j][2][0];
				AY[2][1][j] = Abar[i][j][2][1];
				AY[2][2][j] = Abar[i][j][2][2];

				BY[0][0][j] = Bbar[i][j][0][0];
				BY[0][1][j] = Bbar[i][j][0][1];
				BY[0][2][j] = Bbar[i][j][0][2];
				BY[1][0][j] = Bbar[i][j][1][0];
				BY[1][1][j] = Bbar[i][j][1][1];
				BY[1][2][j] = Bbar[i][j][1][2];
				BY[2][0][j] = Bbar[i][j][2][0];
				BY[2][1][j] = Bbar[i][j][2][1];
				BY[2][2][j] = Bbar[i][j][2][2];

				CY[0][0][j] = Cbar[i][j][0][0];
				CY[0][1][j] = Cbar[i][j][0][1];
				CY[0][2][j] = Cbar[i][j][0][2];
				CY[1][0][j] = Cbar[i][j][1][0];
				CY[1][1][j] = Cbar[i][j][1][1];
				CY[1][2][j] = Cbar[i][j][1][2];
				CY[2][0][j] = Cbar[i][j][2][0];
				CY[2][1][j] = Cbar[i][j][2][1];
				CY[2][2][j] = Cbar[i][j][2][2];

				DY[0][j] = Dbar[i][j][0];
				DY[1][j] = Dbar[i][j][1];
				DY[2][j] = Dbar[i][j][2];
			}

			// Modified Abar and Bbar for UPPER WALL Boundary conditions.
			j = (jmax - 1) - 1;

			AY[0][0][j] = Abar[i][j][0][0] + 4.0 / 3.0 * Cbar[i][j][0][0];		AY[0][1][j] = Abar[i][j][0][1];		AY[0][2][j] = Abar[i][j][0][2];
			 																	 									 
			AY[1][0][j] = Abar[i][j][1][0] + 4.0 / 3.0 * Cbar[i][j][1][0];		AY[1][1][j] = Abar[i][j][1][1];		AY[1][2][j] = Abar[i][j][1][2];
			 																	 									 
			AY[2][0][j] = Abar[i][j][2][0] + 4.0 / 3.0 * Cbar[i][j][2][0];		AY[2][1][j] = Abar[i][j][2][1];		AY[2][2][j] = Abar[i][j][2][2];


			BY[0][0][j] = Bbar[i][j][0][0] - Cbar[i][j][0][0] / 3.0;			BY[0][1][j] = Bbar[i][j][0][1];		BY[0][2][j] = Bbar[i][j][0][2];
															  														 
			BY[1][0][j] = Bbar[i][j][1][0] - Cbar[i][j][1][0] / 3.0;			BY[1][1][j] = Bbar[i][j][1][1];		BY[1][2][j] = Bbar[i][j][1][2];
															  														 
			BY[2][0][j] = Bbar[i][j][2][0] - Cbar[i][j][2][0] / 3.0;			BY[2][1][j] = Bbar[i][j][2][1];		BY[2][2][j] = Bbar[i][j][2][2];


			DY[0][j] = Dbar[i][j][0];
			DY[1][j] = Dbar[i][j][1];
			DY[2][j] = Dbar[i][j][2];

			// Solve Process Begins FOR Y-SWEEP
			BTRIDy(AY, BY, CY, DY, XY);

			for (j = 1; j <= jmax - 2; j++)
			{
				dQ[i][j][0] = XY[0][j];
				dQ[i][j][1] = XY[1][j];
				dQ[i][j][2] = XY[2][j];
			}
			// Solve Ends FOR Y-SWEEP

			// Step 2 BCs Completion

			//LOWER WALL
			dQ[i][0][1 - 1] = (4.0 * dQ[i][1][1 - 1] - dQ[i][2][1 - 1]) / 3.0;

			//UPPER WALL
			dQ[i][jmax - 1][1 - 1] = (4.0 * dQ[i][jmax - 2][1 - 1] - dQ[i][jmax - 3][1 - 1]) / 3.0;
		}
		// end of Y-SWEEP
		///////////////////////////////////////////  HHHHHHHHEEEEEEEEEEEEERRRRRRRRRRRRREEEEEEEEEEEE !!!!!!!!!!!!!!!!!!!
		for (j = 1; j < jmax - 1; j++)
		{
			dQ[0	   ][j][1 - 1] = (4.0 * dQ[1	   ][j][1 - 1] - dQ[2	    ][j][1 - 1]) / 3.0;

			dQ[imax - 1][j][2 - 1] = (4.0 * dQ[imax - 2][j][2 - 1] - dQ[imax - 3][j][2 - 1]) / 3.0;

			dQ[imax - 1][j][3 - 1] = (4.0 * dQ[imax - 2][j][3 - 1] - dQ[imax - 3][j][3 - 1]) / 3.0;
		}

		// STEP 3
		for (i = 0; i < imax; i++)
		{
			for (j = 0; j < jmax; j++)
			{
				// UPDATES Q TO NEXT TIME LEVEL.
				Q[i][j][0] = dQ[i][j][0] + Q[i][j][0];
				Q[i][j][1] = dQ[i][j][1] + Q[i][j][1];
				Q[i][j][2] = dQ[i][j][2] + Q[i][j][2];

				diff0[i][j] = abs(dQ[i][j][0]);
				diff1[i][j] = abs(dQ[i][j][1]);
				diff2[i][j] = abs(dQ[i][j][2]);

				DivV[i][j][0] = DivV[i][j][1];

				/*dQ[i][j][0] = 0.0;
				dQ[i][j][1] = 0.0;
				dQ[i][j][2] = 0.0;

				dQstar[i][j][0] = 0.0;
				dQstar[i][j][1] = 0.0;
				dQstar[i][j][2] = 0.0;*/
			}
		}



		for (i = 1; i < imax - 1; i++)
		{
			//LOWER WALL
			Q[i][0][1 - 1] = (4.0 * Q[i][1][1 - 1] - Q[i][2][1 - 1]) / 3.0;

			//UPPER WALL
			Q[i][jmax - 1][1 - 1] = (4.0 * Q[i][jmax - 2][1 - 1] - Q[i][jmax - 3][1 - 1]) / 3.0;
		}


		for (j = 1; j < jmax - 1; j++)
		{
			Q[0][j][1 - 1] = (4.0 * Q[1][j][1 - 1] - Q[2][j][1 - 1]) / 3.0; // inlet pressure

			Q[imax - 1][j][2 - 1] = (4.0 * Q[imax - 2][j][2 - 1] - Q[imax - 3][j][2 - 1]) / 3.0; // outlet u

			Q[imax - 1][j][3 - 1] = (4.0 * Q[imax - 2][j][3 - 1] - Q[imax - 3][j][3 - 1]) / 3.0; // outlet v
		}

		// END of STEP 3

		iter++;

		for (i = 1; i < imax - 1; i++)
		{
			for (j = 1; j < jmax - 1; j++)
			{
				DivV[i][j][1] = (Q[i + 1][j][2 - 1] - Q[i - 1][j][2 - 1]) / (2.0 * dx) + (Q[i][j + 1][3 - 1] - Q[i][j - 1][3 - 1]) / (2.0 * dy);
			}
		}

		L_inft0 = maxarr(diff0, imax, jmax);
		L_inft1 = maxarr(diff1, imax, jmax);
		L_inft2 = maxarr(diff2, imax, jmax);
		divV = maxmat(DivV, imax, jmax, 1);
		L_inft3 = abs(divV - maxmat(DivV, imax, jmax, 0));


		double sum = 0;

		for (i = 0; i < imax; i++)
		{
			for (j = 0; j < jmax; j++)
			{
				sum = sum + pow((DivV[i][j][1]), 2);
			}
		}
		double L2div = pow(sum / (imax * jmax), 0.5);

		sum = 0;

		for (i = 0; i < imax; i++)
		{
			for (j = 0; j < jmax; j++)
			{
				sum = sum + pow((DivV[i][j][1]- DivV[i][j][0]), 2);
			}
		}
		double L2prev = pow(sum / (imax * jmax), 0.5);


		fprintf_s(convergence, "%-17.10d", iter);
		fprintf_s(convergence, "		");
		fprintf_s(convergence, "%-17.10f", L2div);
		fprintf_s(convergence, "		");
		fprintf_s(convergence, "%-17.10f", L2prev);
		fprintf_s(convergence, "		");
		fprintf_s(convergence, "\n");


		if (iter % 1000 == 0)
		{
			cout << "Iter is : " << iter << "\n\n";
			cout << "dt : " << dt << "\n\n";
			cout << "L_inft0 : " << L_inft0 << "\n\n";
			cout << "L_inft1 : " << L_inft1 << "\n\n";
			cout << "L_inft2 : " << L_inft2 << "\n\n";
			cout << "L_inft3 : " << L_inft3 << "\n\n";
			cout << "div(V) : " << divV << "\n\n\n\n";

			//TecplotOutputAll("all.dat", Q, x, y);
		}

		if (L_inft3 < err  /* && L_inft0 < err*/ && L_inft1 < err && L_inft2 < err)
		{
			cout << "Iter is : " << iter << "\n\n";

			double tawup[imax], tawlow[imax];

			for (int i = 0; i < imax; i++)
			{
				int j = 0;
				tawlow[i] = (-3. * Q[i][j][2 - 1] + 4. * Q[i][j + 1][2 - 1] - Q[i][j + 2][2 - 1]) / (2. * dy);
				j = jmax - 1;
				tawup[i] = -(3. * Q[i][j][2 - 1] - 4. * Q[i][j - 1][2 - 1] + Q[i][j - 2][2 - 1])/ (2. * dy);
			}


			FILE* shear;
			fopen_s(&shear, "shear.plt", "w");
			/* Tecplot Output*/
			/* Header */
			fprintf_s(shear, "Variables = x, du/dyL, du/dyU\n");
			fprintf_s(shear, "Zone T = \"Wall Shear Stress\"");
			fprintf_s(shear, "\n");
			for (int i = 0; i < imax; i++)
			{
				fprintf_s(shear, "%-17.10f", x[i]);
				fprintf_s(shear, "		");
				fprintf_s(shear, "%-17.10f", tawlow[i]);
				fprintf_s(shear, "		");
				fprintf_s(shear, "%-17.10f", tawup[i]);
				fprintf_s(shear, "		");
				fprintf_s(shear, "\n");
			}
			fclose(shear);
			fclose(convergence);


			break;
		}
	}
	// end of while
	TecplotOutputAll("all.plt", Q, x, y);

	//*/

	return 0;
}
