/**
 * @Author: Your name
 * @Date:   2024-09-16 11:03:46
 * @Last Modified by:   Your name
 * @Last Modified time: 2024-09-16 15:00:42
 */
#pragma once
#include"../head/head.h"
#include"../base_fun/basefun.h"
#include"../crystal_lei/baseSolver.h"
class LBSolver {
public:
	LBSolver(BaseSolver *base1): base(base1){}
	~LBSolver()
	{
		fn_mat_free<double>(ProjMatrix);
		fn_mat_free<double>(pk);
		fn_vec_free<double>(pk2);
		fn_vec_free<double>(kt);
		fn_vec_free<double>(ktmp);
		fn_vec_free<double>(kt_gradient);
		fn_vec_free<double>(Phi);
		fn_vec_free<double>(Phi2);
		fn_vec_free<double>(Phi3);
		fn_vec_free<double>(Phi4);
		fn_vec_free<double>(LB_Nonlinear_terms);
		fn_vec_free<double>(LB_e2);
		fn_vec_free<double>(LB_e);
		fn_vec_free<double>(Gradient);
		fn_vec_free<double>(Energy);
		fn_vec_free_cplx<fftw_complex>(Phi_hat);
		fn_vec_free_cplx<fftw_complex>(Phi_hat_t);
		fn_vec_free_cplx<fftw_complex>(LB_Nonlinear_terms_hat);
		fn_vec_free_cplx<fftw_complex>(LB_gradient);
		fn_vec_free_cplx<fftw_complex>(LB_e2_hat);
		fn_vec_free_cplx<fftw_complex>(LB_e_hat);
	}
	void fftw_plan_initial();
	void set_para();
	void solver();

private:
	BaseSolver *base;
	double energy, gradient = 1.0, energy0 = 0.0, energy_err = 1.0;
	int iter = 0;
	double tol = 1e-14;
	int Iter_max = 3000;
	fftw_plan plan_r2c_forward, plan_c2r_backward;
	int Nx = base->getNx();
	int Ny = base->getNy();
	int Nz = base->getNz();
	int dim = base->getdim();
	int Dim = base->getDim();
	double gamma = base -> getgamma();
	double tau = base -> getTau();
	double dt = base -> getdt();
	double c = base -> getc();
	int RealSize = Nx * Ny * Nz;
	int ImagSize = Nx * Ny * (Nz / 2 + 1);
	fftw_complex* Phi_hat = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* Phi_hat_t = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* LB_Nonlinear_terms_hat = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* LB_gradient = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* LB_e2_hat = fn_vec_init_cplx<fftw_complex>(ImagSize);
	fftw_complex* LB_e_hat = fn_vec_init_cplx<fftw_complex>(ImagSize);
	double* Phi = fn_vec_init<double>(RealSize);
	double* Phi2 = fn_vec_init<double>(RealSize);
	double* Phi3 = fn_vec_init<double>(RealSize);
	double* Phi4 = fn_vec_init<double>(RealSize);
	double* LB_Nonlinear_terms = fn_vec_init<double>(RealSize);
	double* LB_e2 = fn_vec_init<double>(RealSize);
	double* LB_e = fn_vec_init<double>(RealSize);
	double* Gradient = fn_vec_init<double>(Iter_max);
	double* Energy = fn_vec_init<double>(Iter_max);
	double** ProjMatrix = fn_mat_init<double>(dim, Dim);
	double** pk = fn_mat_init<double>(ImagSize, dim);
	double* pk2 = fn_vec_init<double>(ImagSize);
	double* kt = fn_vec_init<double>(ImagSize);
	double* ktmp = fn_vec_init<double>(ImagSize);
	double* kt_gradient = fn_vec_init<double>(ImagSize);
};

void LBSolver::fftw_plan_initial()
{
	char filename[50]; // 分配足够的空间来存储文件名和可能的路径    
    sprintf(filename, "fftw_wisdom_%d.dat", Nx); 
    clock_t t1 = clock();
    if (access(filename, F_OK) != -1) {  
        // 文件存在  
        if (fftw_import_wisdom_from_filename(filename) == 0) {
            printf("Failed to import wisdom\n");
        }
        plan_r2c_forward = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, Phi, Phi_hat, FFTW_ESTIMATE);
		plan_c2r_backward = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, Phi_hat, Phi, FFTW_ESTIMATE);
    } else {
        // 文件不存在  
        plan_r2c_forward = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, Phi, Phi_hat, FFTW_PATIENT);
		plan_c2r_backward = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, Phi_hat, Phi, FFTW_PATIENT);

        fftw_export_wisdom_to_filename(filename); 
    }
    clock_t time = clock() - t1;
	std::cout << "creat plan time: " << time << std::endl;
}

inline void LBSolver::set_para()
{
	fftw_plan_initial();
	base->get_ProjMatrix(ProjMatrix);
	base->get_kspace(ProjMatrix, pk, pk2, kt, ktmp, kt_gradient);
	base->get_initial_value(Phi_hat);
	fftw_execute_dft_c2r(plan_c2r_backward, Phi_hat, Phi);
	base->get_initial_value(Phi_hat);
	
	multiply_xy(Phi, Phi, Phi2, Nx * Ny * Nz);
	multiply_xy(Phi, Phi2, Phi3, Nx * Ny * Nz);
	multiply_xy(Phi, Phi3, Phi4, Nx * Ny * Nz);
	for (int j1 = 0; j1 < Nx * Ny * Nz; j1++) {
		LB_Nonlinear_terms[j1] = 0.5 * gamma * Phi2[j1] - Phi3[j1] / 6.0 - tau * Phi[j1];
	}
	
	fftw_execute_dft_r2c(plan_r2c_forward, LB_Nonlinear_terms, LB_Nonlinear_terms_hat);
	average_fourier(LB_Nonlinear_terms_hat, Nx * Ny * (Nz / 2 + 1), Nx * Ny * Nz);
}

inline void LBSolver::solver()
{
	clock_t t1 = clock();
	while (energy_err > tol && iter < Iter_max)
	{
		iter = iter + 1;
		
		for (int j1 = 0; j1 < Nx * Ny * (Nz / 2 + 1); j1++) {
			Phi_hat_t[j1][0] = (Phi_hat[j1][0] + dt * LB_Nonlinear_terms_hat[j1][0]) / ktmp[j1];
			Phi_hat_t[j1][1] = (Phi_hat[j1][1] + dt * LB_Nonlinear_terms_hat[j1][1]) / ktmp[j1];
		}
		Phi_hat_t[0][0] = 0.0;
		Phi_hat_t[0][1] = 0.0;
		
		for (int j1 = 0; j1 < Nx * Ny * (Nz / 2 + 1); j1++) {
			LB_gradient[j1][0] = fabs(Phi_hat_t[j1][0] - Phi_hat[j1][0]) / dt;
			LB_gradient[j1][1] = fabs(Phi_hat_t[j1][1] - Phi_hat[j1][1]) / dt;
		}
		for (int j1 = 0; j1 < Nx * Ny * (Nz / 2 + 1); j1++) {
			Phi_hat[j1][0] = Phi_hat_t[j1][0];
			Phi_hat[j1][1] = Phi_hat_t[j1][1];
		}

		gradient = Max(LB_gradient, Nx * Ny * (Nz / 2 + 1));
		Gradient[iter - 1] = gradient;
		
		for (int j1 = 0; j1 < Nx * Ny * (Nz / 2 + 1); j1++) {
			LB_e2_hat[j1][0] = kt[j1] * Phi_hat_t[j1][0];
			LB_e2_hat[j1][1] = kt[j1] * Phi_hat_t[j1][1];
		}
		fftw_execute_dft_c2r(plan_c2r_backward, Phi_hat_t, Phi);
		multiply_xy(Phi, Phi, Phi2, Nx * Ny * Nz);
		multiply_xy(Phi2, Phi, Phi3, Nx * Ny * Nz);
		multiply_xy(Phi3, Phi, Phi4, Nx * Ny * Nz);
		for (int j1 = 0; j1 < Nx * Ny * Nz; j1++) {
			LB_Nonlinear_terms[j1] = 0.5 * gamma * Phi2[j1] - Phi3[j1] / 6.0 - tau * Phi[j1];
		}
		
		fftw_execute_dft_r2c(plan_r2c_forward, LB_Nonlinear_terms, LB_Nonlinear_terms_hat);
		average_fourier(LB_Nonlinear_terms_hat, Nx * Ny * (Nz / 2 + 1), Nx * Ny * Nz);

		fftw_execute_dft_c2r(plan_c2r_backward, LB_e2_hat, LB_e2);
		multiply_xy(LB_e2, LB_e2, LB_e, Nx * Ny * Nz);
		for (int j1 = 0; j1 < Nx * Ny * Nz; j1++) {
			LB_e[j1] = c * LB_e[j1] / 2.0 + tau * Phi2[j1] / 2.0 - gamma * Phi3[j1] / 6.0 + Phi4[j1] / 24.0;
		}
		fftw_execute_dft_r2c(plan_r2c_forward, LB_e, LB_e_hat);
		average_fourier(LB_e_hat, Nx * Ny * (Nz / 2 + 1), Nx * Ny * Nz);
		
		energy = **LB_e_hat;
		*(Energy + iter - 1) = energy;
		energy_err = fabs(energy - energy0);
		energy0 = energy;
		std::cout << "iter = " << iter << " ";
		std::cout << "gradient = " << std::scientific << gradient << " ";
		std::cout << std::setprecision(15) << "energy = " << energy << " ";
		std::cout << std::setprecision(6) << "energy_err = " << std::scientific << energy_err << std::endl;
	}
	clock_t time = clock() - t1;
	std::cout << "run time: " << time << std::endl;
}