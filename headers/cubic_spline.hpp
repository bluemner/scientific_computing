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
#ifndef _BETACORE_CUBIC_SPLINE_H_
#define _BETACORE_CUBIC_SPLINE_H_

#ifndef _BETACORE_LINAR_SYSTEM_H_
#include "./linar_system.hpp"
#endif
namespace betacore{
struct Cubic_Spline_Exception : public std::exception {
   const char * what () const throw () {
      return "Linar System Exception";
   }
};
	template <typename T>
	class Cubic_Spline{
		private:
			Linar_System<T> * linar_system;
			T left_bound_value;
			T right_bound_value;
			int left_bound_order;
			int right_bound_order;

			template<size_t rows, size_t cols>
			void make_linar_system(T (&data)[rows][cols], int mode){
				const size_t n = rows-1;

				T matrix[n][n];
				T right_hand_side[n];
				T b_val[n];
				//Only neeed to populate the upper matrix
				if(mode ==0){			
					for(int i=0; i<n;i++){
						for(int j=0; j<n; j++){
							if(i==j){
								if(i==0 || i== n-1){
									matrix[i][j]=4;
								}else{
									matrix[i][j]=2;
								}
							}else{
								if(j-1 ==i || j+1 ==i){
									matrix[i][j]=1;
								}
								else{
									matrix[i][j]=0;
								}
							}
						}
					}
				}else if(mode ==1){
					for(int i=1; i<n+1; i++) {
						matrix[i][i-1] = (T)(1.0/3.0) * (data[i][1]-data[i-1][1]);
						matrix[i][i]   = (T)(2.0/3.0) * (data[i+1][1]-data[i-1][1]);
						matrix[i][i+1] = (T)(1.0/3.0) * (data[i+1][1]-data[i][1]);
						right_hand_side[i]=(data[i+1][2]-data[i][2])/(data[i+1][1]-data[i][1]) - (data[i][2]-data[i-1][2])/(data[i][1]-data[i-1][1]);
					}
					if(left_bound_order == FIRST_DIRIVATE){
							matrix[n-1][n-1]=(T)2.0*(data[n-1][1]-data[n-2][1]);
							matrix[n-1][n-2]=(T)1.0*(data[n-1][1]-data[n-2][1]);
							right_hand_side[n-1]=3.0*(right_bound_value-(data[n-1][2]-data[n-2][2])/(data[n-1][1]-data[n-2][1]));
					}else if(left_bound_order == SECOND_DERIVATIVE){
						matrix[0][0] = 2.0;
						matrix[0][1] = 0.0;
						right_hand_side[0]=left_bound_value;
					}
				}else{
					//TODO;
				}
				T  b[n];
				T  x[n];
				for(int i=0; i<n; i++){
					b[i]= 3 * (data[i+1][2]-data[i][2]);
					x[i]= 1;// just to have a value;
				}
				T L[n][n];
				T U[n][n]; 
				Matrix<T>::lu(matrix,L,U);

				std::cout<<"Cublic spline matrix A:"<<std::endl;
				Matrix<T>(matrix).print();
				std::cout<<"LU:"<<std::endl;
				std::cout<<"L:"<<std::endl;
				Matrix<T>(L).print();
				std::cout<<"U:"<<std::endl;
				Matrix<T>(U).print();
				std::cout<<"Loading U into linar_system Pre backward_substitution"<<std::endl;
				this->linar_system=new Linar_System<T>( U,x,b);
				this->linar_system->print();
			}
		public:
			
			static const int FIRST_DIRIVATE=1;
			static const int SECOND_DERIVATIVE=2;
			
			Cubic_Spline(){	}
			~Cubic_Spline(){
				delete this->linar_system;
				this->linar_system= nullptr;
			}
			template<size_t rows, size_t cols>
			Cubic_Spline(T (&data)[rows][cols]){
				this->make_linar_system(data, 1);
				this->linar_system->backward_substitution();

				this->left_bound_value= (T) 0;
				this->left_bound_order= SECOND_DERIVATIVE;

				this->right_bound_value= (T) 0;
				this->right_bound_order= SECOND_DERIVATIVE;
				;
			}
			void print(){
				std::cout<<"Loading U into linar_system Post backward_substitution"<<std::endl;
				this->linar_system->print();
			}
	};
}
#endif