/*********************************************************************************
pi_calc.cpp calculates value of pi and compares with actual
value (to 25digits) of pi to give error. Integrates function f(x)=4/(1+x^2).
July 6, 2001 K. Spry CSCI3145
**********************************************************************************/
#include <math.h> 
#include <iostream>
#include <string>
#include <limits>
#include <cmath>
#include <chrono>
#include <ctime>
#include "mpi.h"

int main(int argc, char *argv[])
{
	if(argc < 2){
		std::cout<< "Ussage:\n\tPlease provide the number of intervals" <<std::endl;
		return -1;
	}
	double actual_pi = 3.141592653589793238462643;
	//for comparison later
	int n, rank, num_proc, i, p;
	double temp_pi, calc_pi, int_size, part_sum, x;
	double max =std::numeric_limits<double>::max();
	char response = 'y', resp1 = 'y';

	MPI_Init(&argc, &argv); 
		
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	if (rank == 0) { /*I am root node*/
		//std::cout <<"__________________________________" <<std::endl;
		//std::cout <<"\nEnter the number of intervals: (0 will exit)" << std::endl;
		//std::cin >> n;
		 p = std::stoi(argv[1]);
		 n = std::pow(2,p);

		std::cout << "Processor" <<"," << "intervals" << "," << "Calculated pi" << "," << "Error 25 places" <<"," << "Cpu Time"  << std::endl;
	}
	while (n < max ){
		auto start = std::chrono::system_clock::now();
		auto start_cpu = std::clock();
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); //broadcast n
		if( n < 0)
		{
						
			break;
		}
		int_size = 1.0 / (double) n; //calcs interval size
		part_sum = 0.0;

		for ( i = rank + 1; i <= n; i += num_proc )
		{
			//calcs partial sums
			x = int_size * ((double)i - 0.5);
			part_sum += (4.0 / (1.0 + x*x));
		}
		temp_pi = int_size * part_sum;
		//collects all partial sums computes pi
		MPI_Reduce(&temp_pi,&calc_pi, 1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		if (rank == 0) {		
		
			auto end = std::chrono::system_clock::now();
			auto end_cpu = std::clock();	
			std::chrono::duration<double> diff = end-start;
			std::cout << num_proc <<"," << n << "," << calc_pi << "," << fabs(calc_pi - actual_pi) <<"," << diff.count()  << std::endl;
				n = std::pow(2,++p);	
		}
		
	}//while

	MPI_Finalize(); //terminate MPI
	return 0;
} //end main
