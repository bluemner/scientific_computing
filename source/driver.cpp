/** MIT License

	Copyright (c) 2017 Brandon Bluemner

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.
 */
#include <iostream>
#include "../headers/matrix.hpp"
#include "../headers/cubic_spline.hpp"
#include "../headers/lagrangian.hpp"
#include "../headers/golden_section_search.hpp"
#include <chrono>
#include <cmath>

/**
 * Test Driver Class
 */
int main(int argc, char * argv[]){

	#pragma region Region_Matrix
	{
		double d22[2][2] ={
						{4,3.0},
						{6,3},					
						};
		betacore::Matrix<double> d23(d22);
		betacore::Matrix<double> L(2);
		betacore::Matrix<double> U(2);
		d23.lu(L,U);
		L.print();
		U.print();
		double d[3][2] ={
						{1.0,2.0},
						{2.1,3.4},
						{4.0,2.0}
						};
		betacore::Matrix<double> d1(d);				
		betacore::Matrix<double> d2(d);
		d2.print();
		d2+= d2;
		d2.print();
		d2-=d1;
		d2.print();
		d2=d2+d2;
		d2.print();
		d2=d2-d1;
		d2.print();
		d1.transpose();
		d2*=d1;
		d2.print();
			double _d[2][1] ={
						{2.0},
						{-1},
					
						};
		betacore::Matrix<double> d3(_d);
		betacore::Matrix<double> d4(_d);
		d4.transpose();
		d3*=d4;
		//d3/=d3;
		d3.print();
		std::cout<<"done with matrix" <<std::endl;
		// betacore::Matrix<double> matrix(3);
		// std::cout<<"___________________"<<std::endl;
		// std::cout<<"Matrix:["<<matrix.row_count()<<","<<matrix.col_count()<<"]"<<std::endl;
		// std::cout<<"___________________"<<std::endl;
		// matrix.print();
		// matrix.transpose();
		// std::cout<<"___________________"<<std::endl;
		// std::cout<<"Matrix:["<<matrix.row_count()<<","<<matrix.col_count()<<"]"<<std::endl;
		// std::cout<<"___________________"<<std::endl;
		// matrix.print();
		// std::cout<<"___________________"<<std::endl;
		// double d[3][2] ={
		// 				{1.0,2.0},
		// 				{2.1,3.4},
		// 				{4.0,2.0}
		// 				};
		// betacore::Matrix<double> d2(d);
		// std::cout<<"___________________"<<std::endl;
		// std::cout<<"Matrix:["<<d2.row_count()<<","<<d2.col_count()<<"]"<<std::endl;
		// std::cout<<"___________________"<<std::endl;
		// d2.print();
		// d2.transpose();
		// std::cout<<"___________________"<<std::endl;
		// std::cout<<"Matrix:["<<d2.row_count()<<","<<d2.col_count()<<"]"<<std::endl;
		// std::cout<<"___________________"<<std::endl;
		// d2.print();
		// std::cout<<"___________________"<<std::endl;
		
		// d2.transpose();
		// matrix+=d2;

		// std::cout<<"______________________________________"<<std::endl;
		// std::cout<<"Matrix:["<<matrix.row_count()<<","<<matrix.col_count()<<"]"<<std::endl;
		// std::cout<<"______________________________________"<<std::endl;
		// matrix.print();
		// std::cout<<"______________________________________"<<std::endl;
		// char x;
		// std::cin>>x;		
	}
	#pragma endregion Region_Matrix

	#pragma region Linear_System
	{
		std::cout << "______________________________________"<<std::endl;
		std::cout <<"linear_system"<<std::endl;
		std::cout << "______________________________________"<<std::endl;
		double A[3][3] ={
						{2.0,4.0,-2.0},
						{0.0,1.0,1.0}, 
						{0.0,0.0,4.0} 
					};
		double x[3]={1.0,1.0,1.0};
		double b[3]={2.0,4.0,8.0};

		try{
			betacore::Linear_System<double> ls(A,x,b);
			ls.print();
			ls.backward_substitution();
			ls.print();
		}catch(const std::exception& e){
			std::cout<<"Problem\n"<< e.what()<<std::endl;
		}


		std::cout << "______________________________________"<<std::endl;
		std::cout <<"linear_system"<<std::endl;
		std::cout << "______________________________________"<<std::endl;
		double A2[3][3] ={
						{ 2.0,0.0,0.0},
						{ 4.0,1.0,0.0}, 
						{-2.0,1.0,4.0} 
						// { 2.0, 4.0,-2.0},
						// { 4.0, 9.0,-3.0}, 
						// {-2.0,-3.0, 7.0} 
					};
		double x2[3]={1.0,1.0, 1.0};
		double b2[3]={4.0,10.0,-6.0}; //{2.0,8.0,10.0};

		try{
			betacore::Linear_System<double> ls2(A2,x2,b2);
			ls2.print();
			ls2.forward_substitution();
			ls2.print();
		}catch(const std::exception& e){
			std::cout<<"Problem\n"<< e.what()<<std::endl;
		}
	}
	#pragma endregion Linear_System

	#pragma region Cubic_Spline
	{
		auto start = std::chrono::system_clock::now();
		
		std::cout << "______________________________________"<<std::endl;
		std::cout <<"Cubic_Spline"<<std::endl;
		std::cout << "______________________________________"<<std::endl;

		//(1,4)(2,2)(3,4)(4,2)
		double data[6][2] ={
						{1.0,4.0},
						{2.0,2.0},
						{3.0,4.0},
						{4.0,2.0},
						{6.0,0.0},
						{8.0,5.0}
					};
		std::cout<<"Input" <<std::endl;
	    for(int i=0; i< 6 ; i++){  	
  			
			std::cout<<"("<< data[i][0] <<","<< data[i][1]  << ")" <<std::endl;
		}
		std::cout<<"Output Natural Cubic Spline" <<std::endl;
		try{
			betacore::Cubic_Spline<double>cs(data);			
			cs.print();
		}catch(const std::exception& e){
			std::cout<<"Problem\n"<< e.what()<<std::endl;
		}

		std::cout << "______________________________________"<<std::endl;
		std::cout <<"Cubic_Spline 2"<<std::endl;
		std::cout << "______________________________________"<<std::endl;
		//(1,4)(2,2)(3,4)(4,2)
		double data2[3][2] ={
						{-1.0,1.0},
						{0.0,0.0},
						{1.0,1.0}
					};
		std::cout<<"Input" <<std::endl;
	    for(int i=0; i< 3 ; i++){  	
  			
			std::cout<<"("<< data2[i][0] <<","<< data2[i][1]  << ")" <<std::endl;
		}
		std::cout<<"Output Natural Cubic Spline" <<std::endl;
		try{
			betacore::Cubic_Spline<double>cs2(data2);
			cs2.print();
		}catch(const std::exception& e){
			std::cout<<"Problem\n"<< e.what()<<std::endl;
		}

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = end-start;
		std::cout<< "Running time:\t" << diff.count() <<"s" << std::endl;
	}
	#pragma endregion Cubic_Spline

	#pragma region lagrangian
	{
		
	}
	#pragma endregion lagrangian

	double a = 0.0, b =1.0, tol = 0.001;
	betacore::GoldenSectionSearch<double> d;
	d.search(a, b, tol);


	a = -1.0, b =2.0, tol = 0.001;
	d.search(a, b, tol);

	// #ifdef _WIN32
	// 	std::cout<<"Press any key to exit"<<std::endl;
	// 	char x;
	// 	std::cin>>x;
	// #endif
	return 0;
}
//*/
