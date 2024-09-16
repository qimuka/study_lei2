/**
 * @Author: Your name
 * @Date:   2024-09-16 11:02:49
 * @Last Modified by:   Your name
 * @Last Modified time: 2024-09-16 14:58:28
 */
#ifndef BCC_H  
#define BCC_H  
 


#pragma once

#include"../head/head.h"
#include"baseSolver.h"

class BCC:public BaseSolver {
public:
	BCC(int Nxx, int Nyy, int Nzz):BaseSolver(Nxx, Nyy, Nzz){};
	~BCC() override;

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
	double tau = -0.05;
	double gamma = 1.0;
};

BCC::~BCC()
{
	printf("调用了BCC的析构函数\n");
}

inline void BCC::get_ProjMatrix(double **ProjMatrix)
{
    for (int j1 = 0; j1 < dim; j1++)
		for (int j2 = 0; j2 < Dim; j2++)
			ProjMatrix[j1][j2] = 0.0;
	for (int j1 = 0; j1 < dim; j1++)
		ProjMatrix[j1][j1] = 1.0 / sqrt_2;
}

inline void BCC::get_kspace(double** ProjMatrix, double** pkspace, double* pk2space,
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

inline void BCC::get_initial_value(fftw_complex* Imag)
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

	//110
	index = Ny * Nz + Nz;
	temp[index][0] = 0.3;
	//101
	index = Ny * Nz + 1;
	temp[index][0] = 0.3;
	//011
	index = Nz + 1;
	temp[index][0] = 0.3;
	//-110
	index = (Nx - 1) * Ny * Nz + Nz;
	temp[index][0] = 0.3;
	//-101
	index = (Nx - 1) * Ny * Nz + 1;
	temp[index][0] = 0.3;
	//0-11
	index = (Ny - 1) * Nz + 1;
	temp[index][0] = 0.3;
	//1-10
	index = Ny * Nz + (Ny - 1) * Nz;
	temp[index][0] = 0.3;
	//10-1
	index = Ny * Nz + (Nz - 1);
	temp[index][0] = 0.3;
	//01-1
	index = Nz + (Nz - 1);
	temp[index][0] = 0.3;
	//-1-10
	index = (Nx - 1) * Ny * Nz + (Ny - 1) * Nz;
	temp[index][0] = 0.3;
	//-10-1
	index = (Nx - 1) * Ny * Nz + (Nz - 1);
	temp[index][0] = 0.3;
	//0-1-1
	index = (Ny - 1) * Nz + (Nz - 1);
	temp[index][0] = 0.3;

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

double BCC::getTau() const {  
	return tau;  
} 
double BCC::getgamma() const {  
    return gamma;  
} 
double BCC::getc() const  {  
    return c;  
} 
double BCC::getdt() const {  
    return dt;  
} 
int BCC::getNx() const {  
    return Nx;  
} 
int BCC::getNy() const {  
    return Ny;  
} 
int BCC::getNz() const {  
    return Nz;  
} 
int BCC::getdim() const {  
    return dim;  
} 
int BCC::getDim() const {  
    return Dim;  
} 

// 头文件内容  
 
#endif // DG_H