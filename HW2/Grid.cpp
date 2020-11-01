#include "pch.h"
#include "Grid.h"
#include <iostream>

#define PI 3.14159265

Grid::Grid()
{
}

Grid::Grid(int m)
{
	m_iMeshM = m;
	m_dH = 1. / (m_iMeshM + 1);
	m_iNtotal = (m_iMeshM + 2)*(m_iMeshM + 2);

	m_iT = 0;

	m_pUij = new double *[m_iMeshM + 2];   //m+2*m+2    0~m+1
	for (int i = 0; i <= m_iMeshM + 1; i++) m_pUij[i] = new double[m_iMeshM + 2];

	m_pU0ij = new double *[m_iMeshM + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) m_pU0ij[i] = new double[m_iMeshM + 2];

	m_pRij = new double *[m_iMeshM + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) m_pRij[i] = new double[m_iMeshM + 2];


	m_pU00ij = new double *[m_iMeshM + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) m_pU00ij[i] = new double[m_iMeshM + 2];



	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pUij[i][j] = 0.;
				m_pU0ij[i][j] = 0.;
				m_pRij[i][j] = 0.;

				m_pU00ij[i][j] = 0.;
			}
		}

}

Grid::~Grid()
{
	if (m_pUij != NULL)
	{
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			delete[] m_pUij[i];
			m_pUij[i] = NULL;
		}
		delete[] m_pUij;
		m_pUij = NULL;
	}

	if (m_pU0ij != NULL)
	{
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			delete[] m_pU0ij[i];
			m_pU0ij[i] = NULL;
		}
		delete[] m_pU0ij;
		m_pU0ij = NULL;
	}

	if (m_pRij != NULL)
	{
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			delete[] m_pRij[i];
			m_pRij[i] = NULL;
		}
		delete[] m_pRij;
		m_pRij = NULL;
	}

	if (m_pU00ij != NULL)
	{
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			delete[] m_pU00ij[i];
			m_pU00ij[i] = NULL;
		}
		delete[] m_pU00ij;
		m_pU00ij = NULL;
	}
}

void Grid::TwoDWaveStep2()
{

	int iCount;
	double Error, MaxError;
	double x, y;
	double La;


	// BC  u(x,0) = 0
	double dudy, f;
	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pUij[i][0] = 0.;
	}

	// BC  u(x,1) = 0
	for (int i = 0; i <= m_iMeshM + 1; i++)
	{
		m_pUij[i][m_iMeshM + 1] = 0.;
	}

	// BC  u(0,y) = f(y)
	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		y = j * m_dH;
		m_pUij[0][j] = 0.;

	}

	// BC  u(1,y) = 0
	for (int j = 0; j <= m_iMeshM + 1; j++)
	{
		m_pUij[m_iMeshM + 1][j] = 0.;
	}

	// Inner grid
	for (int i = 1; i <= m_iMeshM; i++)
		for (int j = 1; j <= m_iMeshM; j++)
		{
			x = i * m_dH;
			y = j * m_dH;

			La = 1. / m_dH / m_dH * (m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);


			m_pUij[i][j] = 2.* m_pU0ij[i][j]- m_pU00ij[i][j] + m_dt * m_dt *La;;

		}


}


void Grid::SolveTwoDWave()
{
	m_dt = 1.e-4;  //m_dt must smaller than dx
	
	int iMax = (1. / m_dt) - 1;
	//iMax = 1;

	TwoDWaveStep1();
	UpateTime();
	OutputResultAsVTK();

	for (int i = 1; i <= iMax; i++)
	{
		TwoDWaveStep2();
		UpateTime();
		//if ((i + 1) % 1000 == 0) OutputResultAsVTK();
	}

	OutputResultAsVTK();
}


void Grid::TwoDWaveStep1( )
{

	int iCount;
	double Error, MaxError;
	double x,y;
	double La;

	
		// BC  u(x,0) = 0
		double dudy, f;
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			m_pUij[i][0] = 0.;
		}

		// BC  u(x,1) = 0
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			m_pUij[i][m_iMeshM + 1] = 0.;
		}

		// BC  u(0,y) = f(y)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			y = j * m_dH;
			m_pUij[0][j] = 0.;

		}

		// BC  u(1,y) = 0
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			m_pUij[m_iMeshM + 1][j] = 0.;
		}

		// Inner grid
		for (int i = 1; i <= m_iMeshM; i++)
			for (int j = 1; j <= m_iMeshM; j++)
			{
				x = i * m_dH;
				y = j * m_dH;

				La = 1./ m_dH/ m_dH*(m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] - 4.*m_pU0ij[i][j] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);


				m_pUij[i][j] = m_pU0ij[i][j] + FunctionX(x)*FunctionX(y)*m_dt + m_dt * m_dt / 2.*La;;

			}

	

}

void Grid::UpateTime()
{
	m_iT = m_iT + 1;

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				m_pU00ij[i][j] = m_pU0ij[i][j];
				m_pU0ij[i][j] = m_pUij[i][j];
			}
		}
}


void Grid::SolveByJacobi(int iLoopMax, double dRelax)
{
	// Jacobi Loop
	//int iLoopMax = 1e6;

	m_iJacobiLoopMax = iLoopMax;
	double TOL = 1e-8;

	int iCount;
	double Error, MaxError;
	double y;

	for (iCount = 1; iCount <= m_iJacobiLoopMax; iCount++)
	{

		// BC  dudy(x,0) = 0
		double dudy, f;
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			dudy = 0.;
			f = 0.;
			m_pUij[i][0] = m_pU0ij[i][1] - m_dH * dudy + m_dH * m_dH / 2.*f;
		}

		// BC  dudy(x,1) = 0
		for (int i = 0; i <= m_iMeshM + 1; i++)
		{
			dudy = 0.;
			f = 0.;
			m_pUij[i][m_iMeshM + 1] = m_pU0ij[i][m_iMeshM] - m_dH * dudy + m_dH * m_dH / 2.*f;
		}

		// BC  u(0,y) = f(y)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			y = j * m_dH;
			m_pUij[0][j] = cos(2.*PI*y);

			/*f = cos(2.*PI*y);
			if (f >= 0)
			{
				m_pUij[0][j] = 1;
			}
			else
			{
				m_pUij[0][j] = -1;
			}*/
			
		}

		// BC  u(1,y) = 0
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			m_pUij[m_iMeshM + 1][j] = 0.;
		}

		// Inner grid
		for (int i = 1; i <= m_iMeshM; i++)
			for (int j = 1; j <= m_iMeshM; j++)
			{
				{
					m_pUij[i][j] = 0.25*(m_pU0ij[i + 1][j] + m_pU0ij[i][j + 1] + m_pU0ij[i - 1][j] + m_pU0ij[i][j - 1]);

					m_pUij[i][j] = m_pU0ij[i][j] + dRelax * (m_pUij[i][j] - m_pU0ij[i][j]);
				}
			}

		MaxError = 1.e-30;
		for (int i = 0; i <= m_iMeshM + 1; i++)
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				{
					Error = fabs(m_pU0ij[i][j] - m_pUij[i][j]);
					if (Error > MaxError) MaxError = Error;

					m_pU0ij[i][j] = m_pUij[i][j];
				}
			}

		//printf("iCount=%d, MaxError=%lg \n", iCount, MaxError);

		if (MaxError < TOL) break;
	}

	m_dMaxError = MaxError;
	m_iIterations = iCount;

	if(iLoopMax >100) printf("    J MeshM=%d, Iterations=%d, MaxError=%lg \n", m_iMeshM, m_iIterations, m_dMaxError);
}

void Grid::OutputResultAsVTK()
{
	//output result
	double x, y;

	char StrFileName[100] = { 0 };

	//sprintf(StrFileName, "file_%d.vtk", m_iMeshM);

	sprintf(StrFileName, "file_Grid%d_Time%g.vtk", m_iMeshM, m_dt*m_iT);

	FILE *fp;
	//fp = fopen("file.vtk", "w");
	fp = fopen(StrFileName, "w");
	fprintf(fp, "# vtk DataFile Version 3.0\n");
	fprintf(fp, "Result Field\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_GRID\n");
	fprintf(fp, "DIMENSIONS %d %d %d\n", m_iMeshM + 2, m_iMeshM + 2, 1);
	fprintf(fp, "POINTS %d double\n", m_iNtotal);


	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				x = i * m_dH;
				y = j * m_dH;

				fprintf(fp, "%f %f %f\n", x, y, 0.);
			}
		}

	fprintf(fp, "\nPOINT_DATA %d\n", m_iNtotal);
	fprintf(fp, "SCALARS U double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");

	for (int i = 0; i <= m_iMeshM + 1; i++)
		for (int j = 0; j <= m_iMeshM + 1; j++)
		{
			{
				fprintf(fp, "%f\n", m_pUij[i][j]);

			}
		}

	fclose(fp);//closing file.
}

void Grid::SolveByMultiGrid(int iLoopMax)
{
	// 2h grid
	int m_2h = (m_iMeshM - 1) / 2;

	double ** Rij_2h, ** R0ij_2h;
	Rij_2h = new double *[m_2h + 2];
	for (int i = 0; i <= m_2h + 1; i++) Rij_2h[i] = new double[m_2h + 2];

	R0ij_2h = new double *[m_2h + 2];
	for (int i = 0; i <= m_2h + 1; i++) R0ij_2h[i] = new double[m_2h + 2];

	for (int i = 0; i <= m_2h + 1; i++)
		for (int j = 0; j <= m_2h + 1; j++)
		{
			{
				Rij_2h[i][j] = 0.;
				R0ij_2h[i][j] = 0.;
			}
		}

	double TOL = 1e-8;

	int iCount;
	double Error, MaxError;
	
	for ( iCount = 1; iCount <= iLoopMax; iCount++)
	{

			for (int i = 0; i <= m_iMeshM + 1; i++)
				for (int j = 0; j <= m_iMeshM + 1; j++)
				{
					{

						m_pU00ij[i][j] = m_pU0ij[i][j];  //m_pU00ij used for error comparison
					}
				}


			//
			SolveByJacobi(3, 0.6666);

			// Compeute the residual of h grid


				for (int i = 1; i <= m_iMeshM; i++)     // only do inner
					for (int j = 1; j <= m_iMeshM; j++)
					{
						{
							m_pRij[i][j] = 0. - (m_pUij[i + 1][j] + m_pUij[i][j + 1] - 4.* (m_pUij[i][j]) + m_pUij[i - 1][j] + m_pUij[i][j - 1]);

							//printf("m_pRij[%d][%d]=%lg \n", i, j, m_pRij[i][j]);
						}
					}

				
				//for (int i = 0; i <= m_iMeshM + 1; i++) // BC  dudy(x,0) = 0
				//{
				//	m_pRij[i][0] = 0. - (-m_pUij[i][0] + m_pUij[i][1]);
				//}

				//
				//for (int i = 0; i <= m_iMeshM + 1; i++) // BC  dudy(x,1) = 0
				//{
				//	m_pRij[i][m_iMeshM + 1] = 0. - (-m_pUij[i][m_iMeshM + 1] + m_pUij[i][m_iMeshM]);
				//}

				// BC Dirichlet , residul =0


			// Coarsen the residual 
			for (int i = 0; i <= m_2h+1; i++)    
				for (int j = 0; j <= m_2h+1; j++)
				{
					{
						Rij_2h[i][j] = m_pRij[i * 2][j * 2];
						R0ij_2h[i][j] = Rij_2h[i][j];

						//printf("Rij_2h[%d][%d]=%lg \n", i, j, Rij_2h[i][j]);
					}
				}

			//Solve residual vector on 2h grid
			 for (int k = 1; k <= 3; k++)
			{

				for (int i = 1; i <= m_2h; i++)     // only do inner
					for (int j = 1; j <= m_2h; j++)
					{
						{
							Rij_2h[i][j] = 0.25*(R0ij_2h[i][j]  + R0ij_2h[i + 1][j] + R0ij_2h[i][j + 1] + R0ij_2h[i - 1][j] + R0ij_2h[i][j - 1]);
						}
					}
				
					//for (int i = 0; i <= m_2h + 1; i++) // BC  dudy(x,0) = 0
					//{
					//	Rij_2h[i][0] = (R0ij_2h[i][0] + R0ij_2h[i][1]);
					//}


					//for (int i = 0; i <= m_2h + 1; i++) // BC  dudy(x,1) = 0
					//{
					//	Rij_2h[i][m_2h + 1] = (R0ij_2h[i][m_2h + 1] + R0ij_2h[i][m_2h]);
					//}

				for (int i = 0; i <= m_2h+1; i++)    
					for (int j = 0; j <= m_2h+1; j++)
					{
						{

							R0ij_2h[i][j] = Rij_2h[i][j];

							//printf(" after Solve , Rij_2h[%d][%d]=%lg \n", i, j, Rij_2h[i][j]);
						}
					}
			}

			int iX, iY;
			// Interpolate back to h grid
			for (int i = 1; i <= m_iMeshM; i++)     // only do inner
				for (int j = 1; j <= m_iMeshM; j++)
				{

					if ((i % 2 == 0) && (j % 2 == 0))
					{
						iX = int(i / 2);
						iY = int(j / 2);
						m_pRij[i][j] = Rij_2h[iX][iY];
					}
					else if ((i % 2 == 0) && (j % 2 != 0))
					{
						iX = int(i / 2);
						iY = int(j / 2);
						m_pRij[i][j] = 0.5*(Rij_2h[iX][iY] + Rij_2h[iX][iY + 1]);
					}
					else if ((i % 2 != 0) && (j % 2 == 0))
					{
						iX = int(i / 2);
						iY = int(j / 2);
						m_pRij[i][j] = 0.5*(Rij_2h[iX][iY] + Rij_2h[iX + 1][iY]);
					}
					else
					{
						iX = int(i / 2);
						iY = int(j / 2);
						m_pRij[i][j] = 0.25*(Rij_2h[iX][iY] + Rij_2h[iX + 1][iY] + Rij_2h[iX][iY + 1] + Rij_2h[iX + 1][iY + 1]);
					}

					m_pUij[i][j] -= m_pRij[i][j];

					//printf(" Interpolate , m_pU0ij[%d][%d]=%lg \n", i, j, m_pU0ij[i][j]);
					//printf(" Interpolate , m_pUij[%d][%d]=%lg \n", i, j, m_pUij[i][j]);


					m_pU0ij[i][j] = m_pUij[i][j];
				}

			// iterations on h grid to smooth out error induced by interpolation
			SolveByJacobi(3, 0.6666);



			MaxError = 1.e-30;
			for (int i = 0; i <= m_iMeshM + 1; i++)
				for (int j = 0; j <= m_iMeshM + 1; j++)
				{
					{
						Error = fabs(m_pU00ij[i][j] - m_pUij[i][j]);
						if (Error > MaxError) MaxError = Error;
				
					}
				}

			printf("iCount=%d, MaxError=%lg \n", iCount, MaxError);

			if (MaxError < TOL) break;
	}

	m_dMaxError = MaxError;
	m_iIterations = iCount;

	//printf("MeshM=%d, Iterations=%d, MaxError=%lg \n", m_iMeshM, m_iIterations, m_dMaxError);
}

/*
void Grid::SolveByMultiGrid(int iLoopMax)
{
	// 2h grid
	int m_2h = (m_iMeshM - 1) / 2;

	double ** Rij_2h, ** R0ij_2h;
	Rij_2h = new double *[m_2h + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) Rij_2h[i] = new double[m_2h + 2];

	R0ij_2h = new double *[m_2h + 2];
	for (int i = 0; i <= m_iMeshM + 1; i++) R0ij_2h[i] = new double[m_2h + 2];

	for (int i = 0; i <= m_2h + 1; i++)
		for (int j = 0; j <= m_2h + 1; j++)
		{
			{
				Rij_2h[i][j] = 0.;
				R0ij_2h[i][j] = 0.;
			}
		}

	double TOL = 1e-8;

	int iCount;
	double Error, MaxError;

	for (iCount = 1; iCount <= iLoopMax; iCount++)
	{

		for (int i = 0; i <= m_iMeshM + 1; i++)
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				{

					m_pU00ij[i][j] = m_pU0ij[i][j];  //m_pU00ij used for error comparison
				}
			}


		//
		SolveByJacobi(3, 0.6666);

		// Compeute the residual of h grid


		for (int i = 1; i <= m_iMeshM; i++)     // only do inner
			for (int j = 1; j <= m_iMeshM; j++)
			{
				{
					m_pRij[i][j] = 0. - (m_pUij[i + 1][j] + m_pUij[i][j + 1] - 4.* (m_pUij[i][j]) + m_pUij[i - 1][j] + m_pUij[i][j - 1]);

					//printf("m_pRij[%d][%d]=%lg \n", i, j, m_pRij[i][j]);
				}
			}



		// Coarsen the residual 
		for (int i = 1; i <= m_2h; i++)     // only do inner
			for (int j = 1; j <= m_2h; j++)
			{
				{
					Rij_2h[i][j] = m_pRij[i * 2][j * 2];
					R0ij_2h[i][j] = Rij_2h[i][j];

					//printf("Rij_2h[%d][%d]=%lg \n", i, j, Rij_2h[i][j]);
				}
			}

		//Solve residual vector on 2h grid
		for (int k = 1; k <= 3; k++)
		{

			for (int i = 1; i <= m_2h; i++)     // only do inner
				for (int j = 1; j <= m_2h; j++)
				{
					{
						Rij_2h[i][j] = 0.25*(R0ij_2h[i][j] + R0ij_2h[i + 1][j] + R0ij_2h[i][j + 1] + R0ij_2h[i - 1][j] + R0ij_2h[i][j - 1]);
					}
				}

			for (int i = 1; i <= m_2h; i++)     // only do inner
				for (int j = 1; j <= m_2h; j++)
				{
					{

						R0ij_2h[i][j] = Rij_2h[i][j];

						//printf(" after Solve , Rij_2h[%d][%d]=%lg \n", i, j, Rij_2h[i][j]);
					}
				}
		}

		int iX, iY;
		// Interpolate back to h grid
		for (int i = 1; i <= m_iMeshM; i++)     // only do inner
			for (int j = 1; j <= m_iMeshM; j++)
			{

				if ((i % 2 == 0) && (j % 2 == 0))
				{
					iX = int(i / 2);
					iY = int(j / 2);
					m_pRij[i][j] = Rij_2h[iX][iY];
				}
				else if ((i % 2 == 0) && (j % 2 != 0))
				{
					iX = int(i / 2);
					iY = int(j / 2);
					m_pRij[i][j] = 0.5*(Rij_2h[iX][iY] + Rij_2h[iX][iY + 1]);
				}
				else if ((i % 2 != 0) && (j % 2 == 0))
				{
					iX = int(i / 2);
					iY = int(j / 2);
					m_pRij[i][j] = 0.5*(Rij_2h[iX][iY] + Rij_2h[iX + 1][iY]);
				}
				else
				{
					iX = int(i / 2);
					iY = int(j / 2);
					m_pRij[i][j] = 0.25*(Rij_2h[iX][iY] + Rij_2h[iX + 1][iY] + Rij_2h[iX][iY + 1] + Rij_2h[iX + 1][iY + 1]);
				}

				m_pUij[i][j] -= m_pRij[i][j];

				//printf(" Interpolate , m_pU0ij[%d][%d]=%lg \n", i, j, m_pU0ij[i][j]);
				//printf(" Interpolate , m_pUij[%d][%d]=%lg \n", i, j, m_pUij[i][j]);


				m_pU0ij[i][j] = m_pUij[i][j];
			}

		// iterations on h grid to smooth out error induced by interpolation
//		SolveByJacobi(3, 0.6666);



		MaxError = 1.e-30;
		for (int i = 0; i <= m_iMeshM + 1; i++)
			for (int j = 0; j <= m_iMeshM + 1; j++)
			{
				{
					Error = fabs(m_pU00ij[i][j] - m_pUij[i][j]);
					if (Error > MaxError) MaxError = Error;

				}
			}

		printf("iCount=%d, MaxError=%lg \n", iCount, MaxError);

		if (MaxError < TOL) break;
	}

	m_dMaxError = MaxError;
	m_iIterations = iCount;

	//printf("MeshM=%d, Iterations=%d, MaxError=%lg \n", m_iMeshM, m_iIterations, m_dMaxError);
}

*/