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
#ifndef _BETACORE_HOUSEHOLDER_H_
#define _BETACORE_HOUSEHOLDER_H_

#ifndef _BETACORE_MATRIX_H_
#include "./matrix.hpp"
#endif

#include <math.h>

namespace betacore{
	struct Householder_Exception : public std::exception {
		const char * what () const throw () {
			return "Householder  Exception";
		}
	};
	template<typename T>
	class Householder{
		private:
			Matrix<T> * A;
			Matrix<T> * R;
			unsigned int _rows;
			unsigned int _cols;
		public:
			template<size_t rows, size_t cols>
			Householder(T (&A)[rows][cols]){
				this->A = new Matrix<T>(A);
				this->R = new Matrix<T>(A);
				this->_rows = rows;
				this->_cols = cols;
			}
			~Householder(){
				delete A;
				A = nullptr;
				delete R;
				R = nullptr;
			}
			void run(){
				T magnitude;
				T alpha;
				Matrix<T> v(this->_rows,1);
				Matrix<T> vt;
				Matrix<T> u(this->_rows,1);
				Matrix<T> Identity((const unsigned int) this->_cols);
				Matrix<T> P(this->_rows, this->_rows);
				Matrix<T> Q(this->_rows, this->_rows);
				for ( unsigned int i = 0 ; i < this->_cols ; i++){
					u.zero();
					v.zero();
					
					//GET ALPHA
					magnitude = (T) 0.0;
					for(unsigned int j=0; j < this->_rows; j++){
						u[j] = R->at(i,j);
						magnitude += (T) u[j]*u[j];
					}
					
					alpha = u[i] < 0 ? magnitude : -magnitude;
					magnitude= (T) 0.0;
					for(unsigned int j=i; j< this->_rows; j++){
						v[j]= j==i? u[j]+ alpha : u[j];
						magnitude += v[j] * v[j];
					}
					magnitude = sqrt(magnitude);
					
					if(magnitude == 0){
						continue;
					}

					for(size_t j=i; j < this->_rows; j++){
						v[j] /= magnitude;
					}
					
					vt = v;
					vt.transpose();
					P = Identity - (v*vt) * 2.0;
					P.print();
					(*R) = P * (*R);
					R->print();
					Q = Q * P;
					Q.print();
				}
			}
			void print(){
			
			}
	};
}
#endif