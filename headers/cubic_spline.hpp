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

			template<size_t rows, size_t cols>
			void make_linar_system(T (&data)[rows][cols], int open){
				const size_t n = rows-1;
				T matrix[n][n];
				//Only neeed to populate the upper matrix
				if(open ==0){			
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
				}else{
					//TODO;
				}
				T  b[n];
				T  x[n];
				for(int i=0; i<n; i++){
					b[i]=  3 * (data[i+1]-data[i]);
					x[i]= 1;// just to have a value;
				}
				T L[n][n];
				T U[n][n];
				Matrix<T>::lu(matrix,L,U);
				this->linar_system=new Linar_System<T>( U,x,b);
			}
		public:
			
			Cubic_Spline(){	}
			~Cubic_Spline(){
				delete this->linar_system;
				this->linar_system= nullptr;
			}
			template<size_t rows, size_t cols>
			Cubic_Spline(T (&data)[rows][cols]){
				this->make_linar_system(data, 0);
			}
			void print(){
				this->linar_system->backward_substitution();
				this->linar_system->print();
			}
	};
}
#endif