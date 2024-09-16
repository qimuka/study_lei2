/**
 * @Author: Your name
 * @Date:   2024-09-16 11:02:49
 * @Last Modified by:   Your name
 * @Last Modified time: 2024-09-16 14:58:28
 */
#ifndef BASE_H  
#define BASE_H  
 


#pragma once

#include"../head/head.h"
#ifndef Pi
#define Pi 3.1415926535897932384626	
#endif // !Pi

class BaseSolver {
public:
	BaseSolver(int Nxx, int Nyy, int Nzz);
	virtual ~BaseSolver(){}

	virtual void get_kspace(double** ProjMatrix, double** pkspace, double* pk2space,
		double* ktspace, double* ktmpspace, double* kt_gradientspace) = 0;
	virtual void get_initial_value(fftw_complex* Imag) = 0;
	virtual void get_ProjMatrix(double **ProjMatrx) = 0;
	virtual double getTau() const = 0;
	virtual double getgamma() const = 0;
	virtual double getc() const = 0;
	virtual double getdt() const = 0;
	virtual int getNx() const = 0;
	virtual int getNy() const = 0;
	virtual int getNz() const = 0;
	virtual int getdim() const = 0;
	virtual int getDim() const = 0;
protected:
	int Nx, Ny, Nz;
	int dim = 3;
	int Dim = 3;
	double sqrt_2 = 1.41421356237309504880;
 	double sqrt_3 = 1.73205080756887729352;
	double sqrt_5 = 2.23606797749978969640;
	double sqrt_6 = 2.44948974278317809819;
};

BaseSolver::BaseSolver(int Nxx, int Nyy, int Nzz)
{
	Nx = Nxx;
	Ny = Nyy;
	Nz = Nzz;
}

// 头文件内容  
 
#endif // DG_H