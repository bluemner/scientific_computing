/*********************************************************************************
pi_calc.cpp calculates value of pi and compares with actual
value (to 25digits) of pi to give error. Integrates function f(x)=4/(1+x^2).
July 6, 2001 K. Spry CSCI3145
**********************************************************************************/
#include <math.h> //include files
#include <iostream>
#include "mpi.h"

void print_it(); //function prototypes


int main(int argc, char *argv[])
{
	double actual_pi = 3.141592653589793238462643;
	//for comparison later
	int n, rank, num_proc, i;
	double temp_pi, calc_pi, int_size, part_sum, x;
	char response = 'y', resp1 = 'y';

	MPI_Init(&argc, &argv); 
		
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	print_it(); /* I am root node, print out welcome */
	while (response == 'y') {
	if (resp1 == 'y') {
		if (rank == 0) { /*I am root node*/
			std::cout <<"__________________________________" <<std::endl;
			std::cout <<"\nEnter the number of intervals: (0 will exit)" << std::endl;
			std::cin >> n;
		}
	}
	else
		n = 0;
	
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); //broadcast n
	
	if (n==0)
		break; //check for quit condition
	else {
		int_size = 1.0 / (double) n; //calcs interval size
		part_sum = 0.0;
		for (i = rank + 1; i <= n; i += num_proc)
		{ //calcs partial sums
		x = int_size * ((double)i - 0.5);
		part_sum += (4.0 / (1.0 + x*x));
		}
		temp_pi = int_size * part_sum;
		//collects all partial sums computes pi
		MPI_Reduce(&temp_pi,&calc_pi, 1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (rank == 0) { /*I am server*/
			std::cout << "pi is approximately " << calc_pi
			<< ". Error is " << fabs(calc_pi - actual_pi)
			<< std::endl
			<<"_______________________________________"
			<< std::endl;
			}
	} //end else
		if (rank == 0) { /*I am root node*/
		std::cout << "\nCompute with new intervals? (y/n)" << std::endl;
		std::cin >> resp1;
		}
	}//end while
		MPI_Finalize(); //terminate MPI
		return 0;
} //end main
//functions
void print_it()
{
	std::cout << "\n*********************************" << std::endl;
	std::cout << "Welcome to the pi calculator!" << std::endl;
	std::cout << "Programmer: K. Spry" << std::endl;
	std::cout << "You set the number of divisions \nfor estimating the integral:\n\tf(x)=4/(1+x^2)" << std::endl;
	std::cout << "*********************************" << std::endl;;
} //end print_it