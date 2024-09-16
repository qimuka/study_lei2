/**
 * @Author: Your name
 * @Date:   2024-09-16 11:02:49
 * @Last Modified by:   Your name
 * @Last Modified time: 2024-09-16 14:58:28
 */
#ifndef DG_H  
#define DG_H 
#pragma once

#include"../head/head.h"
#include"baseSolver.h"

class DG:public BaseSolver {
public:
	DG(int Nxx, int Nyy, int Nzz):BaseSolver(Nxx, Nyy, Nzz){};
	~DG() override;

	void get_kspace(double** ProjMatrix, double** pkspace, double* pk2space,
		double* ktspace, double* ktmpspace, double* kt_gradientspace) override;
	void get_initial_value(fftw_complex* Imag) override;
    void get_ProjMatrix(double **ProjMatrx);
    double getTau() const override;
    double getgamma() const override;
    double getc() const override;
    double getdt() const override;
    int getNx() const override;
    int getNy() const override;
    int getNz() const override;
    int getdim() const override;
    int getDim() const override;
    
protected:
    double c = 1.0;
	double dt = 0.1;
	double tau = -0.4;
	double gamma = 0.4;

};

DG::~DG()
{
	printf("调用了DG的析构函数\n");
}

inline void DG::get_ProjMatrix(double **ProjMatrix)
{
    for (int j1 = 0; j1 < dim; j1++)
		for (int j2 = 0; j2 < Dim; j2++)
			ProjMatrix[j1][j2] = 0.0;
	for (int j1 = 0; j1 < dim; j1++)
		ProjMatrix[j1][j1] = 1.0 / sqrt_6;
}

inline void DG::get_kspace(double** ProjMatrix, double** pkspace, double* pk2space,
	double* ktspace, double* ktmpspace, double* kt_gradientspace)
{
	int* k = (int*)malloc(sizeof(int) * Dim);
	double* pk1 = (double*)malloc(sizeof(double) * Dim);
	int index;
	for (int j1 = 0; j1 < Nx; j1++) {
		k[0] = j1 < Nx / 2 ? j1 : j1 - Nx;
		for (int j2 = 0; j2 < Ny; j2++) {
			k[1] = j2 < Ny / 2 ? j2 : j2 - Ny;
			for (int j3 = 0; j3 < Nz / 2 + 1; j3++) {
				k[2] = j3;
				index = j1 * Ny * (Nz / 2 + 1) + j2 * (Nz / 2 + 1) + j3;
				pk2space[index] = 0.0;
				for (int i1 = 0; i1 < dim; i1++) {
					pk1[i1] = 0.0;
					for (int i2 = 0; i2 < Dim; i2++)
						pk1[i1] += ProjMatrix[i1][i2] * double(k[i2]);
					pkspace[index][i1] = pk1[i1];
					pk2space[index] = pk2space[index] + pk1[i1] * pk1[i1];
				}
				ktspace[index] = 1.0 - pk2space[index];
				ktmpspace[index] = 1.0 + dt * c * ktspace[index] * ktspace[index];
				kt_gradientspace[index] = c * ktspace[index] * ktspace[index];
			}
		}
	}
	free(k);
	free(pk1);
}

inline void DG::get_initial_value(fftw_complex* Imag)
{
	fftw_complex* temp = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Ny * Nz);
	int index;
	for (int j1 = 0; j1 < Nx; j1++)
	{
		for (int j2 = 0; j2 < Ny; j2++)
		{
			for (int j3 = 0; j3 < Nz; j3++)
			{
				index = j1 * Ny * Nz + j2 * Nz + j3;
				**(temp + index) = 0.0;
				*(*(temp + index) + 1) = 0.0;
			}
		}
	}

	//211
	index = 2 * Ny * Nz + Nz + 1;
	**(temp + index) = -0.3;
	//21-1
	index = 2 * Ny * Nz + Nz + (Nz - 1);
	**(temp + index) = 0.3;
	//-121
	index = (Nx - 1) * Ny * Nz + 2 * Nz + 1;
	**(temp + index) = 0.3;
	//-12-1
	index = (Nx - 1) * Ny * Nz + 2 * Nz + (Nz - 1);
	**(temp + index) = 0.3;
	//-2-11
	index = (Nx - 2) * Ny * Nz + (Ny - 1) * Nz + 1;
	**(temp + index) = 0.3;
	//-2-1-1
	index = (Nx - 2) * Ny * Nz + (Ny - 1) * Nz + (Nz - 1);
	**(temp + index) = -0.3;
	//1-21
	index = Ny * Nz + (Ny - 2) * Nz + 1;
	**(temp + index) = 0.3;
	//1-2-1
	index = Ny * Nz + (Ny - 2) * Nz + (Nz - 1);
	**(temp + index) = 0.3;

	//121
	index = Ny * Nz + 2 * Nz + 1;
	**(temp + index) = -0.3;
	//12-1
	index = Ny * Nz + 2 * Nz + (Nz - 1);
	**(temp + index) = -0.3;
	//-211
	index = (Nx - 2) * Ny * Nz + Nz + 1;
	**(temp + index) = 0.3;
	//-21-1
	index = (Nx - 2) * Ny * Nz + Nz + (Nz - 1);
	**(temp + index) = -0.3;
	//-1-21
	index = (Nx - 1) * Ny * Nz + (Ny - 2) * Nz + 1;
	**(temp + index) = -0.3;
	//-1-2-1
	index = (Nx - 1) * Ny * Nz + (Ny - 2) * Nz + (Nz - 1);
	**(temp + index) = -0.3;
	//2-11
	index = 2 * Ny * Nz + (Ny - 1) * Nz + 1;
	**(temp + index) = -0.3;
	//2-1-1
	index = 2 * Ny * Nz + (Ny - 1) * Nz + (Nz - 1);
	**(temp + index) = 0.3;

	//112
	index = Ny * Nz + Nz + 2;
	**(temp + index) = -0.3;
	//11-2
	index = Ny * Nz + Nz + (Nz - 2);
	**(temp + index) = 0.3;
	//-112
	index = (Nx - 1) * Ny * Nz + Nz + 2;
	**(temp + index) = -0.3;
	//-11-2
	index = (Nx - 1) * Ny * Nz + Nz + (Nz - 2);
	**(temp + index) = 0.3;
	//-1-12
	index = (Nx - 1) * Ny * Nz + (Ny - 1) * Nz + 2;
	**(temp + index) = 0.3;
	//-1-1-2
	index = (Nx - 1) * Ny * Nz + (Ny - 1) * Nz + (Nz - 2);
	**(temp + index) = -0.3;
	//1-12
	index = Ny * Nz + (Ny - 1) * Nz + 2;
	**(temp + index) = 0.3;
	//1-1-2
	index = Ny * Nz + (Ny - 1) * Nz + (Nz - 2);
	**(temp + index) = -0.3;

	for (int j1 = 0; j1 < Nx; j1++) {
		for (int j2 = 0; j2 < Ny; j2++) {
			for (int j3 = 0; j3 < Nz / 2 + 1; j3++) {
				int id = j1 * Ny * Nz + j2 * Nz + j3;
				int id1 = j1 * Ny * (Nz / 2 + 1) + j2 * (Nz / 2 + 1) + j3;
				Imag[id1][0] = temp[id][0];
				Imag[id1][1] = temp[id][1];
			}
		}
	}
	fftw_free(temp);
}

double DG::getTau() const {  
    return tau;  
} 
double DG::getgamma() const {  
    return gamma;  
} 
double DG::getc() const  {  
    return c;  
} 
double DG::getdt() const {  
    return dt;  
} 
int DG::getNx() const {  
    return Nx;  
} 
int DG::getNy() const {  
    return Ny;  
} 
int DG::getNz() const {  
    return Nz;  
} 
int DG::getdim() const {  
    return dim;  
} 
int DG::getDim() const {  
    return Dim;  
} 

#endif // DG_H
