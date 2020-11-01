#pragma once
class Grid
{
public:
	Grid();
	Grid(int m);
	~Grid();

	void SolveByJacobi(int iLoopMax, double dRelax); // HW1
	void SolveByMultiGrid(int iLoopMax);             // HW2
    void OutputResultAsVTK();                        // HW1
    
	void SolveTwoDWave();
	void TwoDWaveStep1();
	void TwoDWaveStep2();
	void UpateTime();


	int m_iMeshM;
	int m_iNtotal;
	int m_iT;
	int m_iIterations;
	int m_iJacobiLoopMax;
	double m_dH;
	double m_dMaxError;
	double m_dt;

	double **m_pUij, **m_pU0ij , **m_pRij;
	double **m_pU00ij;
};

