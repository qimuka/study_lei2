/**
 * @Author: Your name
 * @Date:   2024-09-16 11:02:48
 * @Last Modified by:   Your name
 * @Last Modified time: 2024-09-16 15:01:00
 */
#include"./head/head.h"
#include"./LB_solver/LBSolver.h"
#include"./crystal_lei/BCC.h"
#include"./crystal_lei/DG.h"
int main()
{
	int Nx = 64;
	int Ny = 64;
	int Nz = 64;
	int flag;
	std::cout << "please input crystallographic groups number: ";
	std::cin >> flag;

	BaseSolver* solverPtr = nullptr; 
	if(flag == 229)
	{
		solverPtr = new BCC(Nx,Ny,Nz);
	}
	else if(flag == 230) 
	{
		solverPtr = new DG(Nx,Ny,Nz);
	}
	else
	{
		std::cout<<"just solve BCC and DG"<<std::endl;
		return 1;
	}
	
	LBSolver LB(solverPtr);
	LB.set_para();
	LB.solver();

	delete solverPtr;
	return 0;
}