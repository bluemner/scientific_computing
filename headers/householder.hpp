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

#include <cmath>

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
				Matrix<T> Identity((const unsigned int) this->_rows);
				Matrix<T> P(this->_rows, this->_rows);
				Matrix<T> Q((const unsigned int) this->_rows);
				for ( unsigned int i = 0 ; i < this->_cols ; i++){
					u.erase();
					v.erase();
					
					//GET ALPHA
					magnitude = (T) 0.0;
					for(unsigned int j=i; j < this->_rows; j++){
						auto temp =this-> R->at(j,i);
						u.set(j,0,temp );
						magnitude += (T) u.at(j,0)*u.at(j,0);
					}
					u.print();
					alpha = u[i] < 0 ? (T) sqrt(magnitude) : (T) -sqrt(magnitude);
					magnitude= (T) 0.0;
					
					for(unsigned int j=i; j< this->_rows; j++){
						T temp  = (j==i)? u.at(j,0)+ alpha : u.at(j,0);
						v.set(j,0,temp);
						magnitude += v.at(j,0) * v.at(j,0);
					}
					magnitude = sqrt(magnitude);
					
					if(magnitude == 0){
						continue;
					}

					for(size_t j=i; j < this->_rows; j++){
						T temp = (T)  v.at(j,0) / magnitude;
						v.set(j,0, temp);
					}
					std::cout<<"____________________________________________________________________________"<<"\n";
					std::cout<<"\tIteration:\t"<< i<<"\n";
					std::cout<<"____________________________________________________________________________"<<"\n";

					vt = v;
					vt.print();
					vt.transpose();
					vt.print();
					vt = (v*vt);
					vt.print();
					vt = vt *2.0;
					vt.print();
					
					P = Identity - vt;
					P.print();
					R = P * R;
					R->print();
					Q = Q * P;
					Q.print();
				}
			}
			void print(){
				A->print();
			}
	};
}
#endif