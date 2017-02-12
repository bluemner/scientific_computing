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
#ifndef _BETACORE_MATRIX_H_
#define _BETACORE_MATRIX_H_
#include <iostream>
#include <exception>
namespace betacore
{
struct MatrixException : public std::exception {
   const char * what () const throw () {
      return "Matrix Exception";
   }
};
/**
 * @desc instance of a NxM matrix
 */
template <typename T>
class Matrix
{
  private:
	T **M;
	unsigned int _rows;
	unsigned int _cols;
	template<size_t rows, size_t cols>
	void copy(T (&M)[rows][cols] ){
		this->M = new T*[rows];
		this->_rows = rows;
		this->_cols = cols;
		for(int i = 0; i < rows; ++i)
			this->M[i] = new T[cols];
		for(size_t i = 0; i < rows; ++i)
		{
			for(size_t j = 0; j < cols; ++j){
				this->M[i][j] = M[i][j];
			}
		}
	}
	void copy(T** M,const int rows, const int cols){
		this->_rows=rows;
		this->_cols=cols;
		this->M =new T*[rows];
		for(int i=0; i<rows;++i)
			this->M[i]=new T[cols];
		for(size_t i =0; i< rows; ++i)
		{
			for(size_t j=0; j<cols; ++j){
				this->M[i][j] = M[i][j];
			}
		}
	}
	T** identity(unsigned int size){
		T** t  =new T*[size];
		for(size_t i=0; i<size;++i)
			t[i]=new T[size];
		for(size_t i=0; i< size; ++i){
			for(size_t j=0; j < size; ++j){
				t[i][j]= (i==j)?1:0;
			}
		}
		return t;
	}
  public:
  	/*
	 * Create null matrix
	 */
	Matrix(){
		this->M = NULL;
		this->rows = 0;
		this->cols = 1 ;
	}
	/**
	 * @desc Creates an instace of MxN Matrix of Type <T>
	 * @param {T} matrix with type <T> and size = rows x cols
	 */
	template<size_t rows, size_t cols>
	Matrix(T (&matrix) [rows][cols]){
		copy(matrix);
	}

	/*
	 * @desc Create identiy matrix
	 * @param {int} size Create identity matrix of M x M st. M = Size
	 */
	Matrix(const unsigned int size){
		auto id = identity(size);
		copy(id,size,size);
		delete id;
		id = nullptr;
	}
	/**
	 * Deconstructor for matrix
	 */
	~Matrix(){
		for(int i=0; i< this->_rows; i++){
			delete this->M[i];
			this->M[i]=nullptr;
		}
		delete this->M;
		this->M  = nullptr;
		this->_rows =0;
		this->_cols =0;
	}
	/**
	 * @desc prints out matrix to terminal
	 */
	void print(){
		// for(unsigned int i=0; i< this->_rows; ++i){
		// 	for(unsigned int j =0; j < this->_cols; ++j){
		// 		std::cout<<this->M[i][j]<<"\t";
		// 	}
		// 	std::cout<<std::endl;
		// }

		size_t i = 0;
		size_t j = 0;
		size_t precision = 5;
		size_t width = 10 + precision;
		char fill =' ';
		std::cout<< "┌";
		for(i=0; i< this->_cols; i++){
			std::cout<<" ";
			for(j=0; j< width; ++j){
				std::cout<<" ";
			}			
			
		}			
		std::cout<<" ┐"<<std::endl;

		for(i=0; i< this->_rows; i++){
			std::cout<<"│";
			
			for(j=0; j< this->_cols; j++){
				std::cout<<" "<<std::fixed<<std::setprecision(precision)<<std::setfill(fill)<<std::setw(width);
				std::cout<<this->M[i][j];
			}
			std::cout<<" │";	
							
			std::cout<<std::endl;
		}
		std::cout<< "└";
		for(i=0; i< this->_cols; i++){	
			std::cout<<" ";				
			for(j=0; j< width; ++j){
				std::cout<<" ";
			}	
		}
		std::cout<<" ┘"<<std::endl;
		
		std::cout<<std::endl;
	
	}
	/**
	 *
	 */
	void transpose(){
		T** transpose = new T*[this->_cols];
		for (unsigned int j =0; j< this->_cols; ++j)
			transpose[j] = new T[this->_rows];
		for( unsigned int i=0; i< this->_rows; ++i){
			for( unsigned int j=0; j < this->_cols; ++j){
				transpose[j][i] = this->M[i][j];
			}
		}
		for (unsigned int j =0; j< this->_cols; ++j)
			 delete this->M[j];
		delete this->M;
		this->M = transpose;
		auto temp = this->_rows;
		this->_rows = this->_cols;
		this->_cols = temp;
	}

	/*
	 * Matrix row Count
	 */
	unsigned int row_count(){
		return this->_rows;
	}
	/*
	 * Matrix Column Count
	 */
	unsigned int col_count(){
		return this->_cols;
	}
	T at(size_t &row, size_t &col){
		if( (unsigned int) row > this->_rows -1 ){
			throw MatrixException();
		}
		if( (unsigned int) col >  this->_cols -1){
			throw MatrixException();
		}
		return this->M[row][col];
	}

	template <size_t size>
	static void lu(T (&A)[size][size], T (&L)[size][size], T (&U)[size][size]){
		size_t i=0;
		size_t j=0;
		size_t k=0;
		
		for (i = 0; i < size; i++)
		{
			for (j = 0; j < size; j++)
			{
				if (j < i)
					L[j][i] = 0;
				else
				{
					L[j][i] = A[j][i];
					for (k = 0; k < i; k++)
					{
						L[j][i] = L[j][i] - L[j][k] * U[k][i];
					}
				}
			}
			for (j = 0; j < size; j++)
			{
				if (j < i)
					U[i][j] = 0;
				else if (j == i)
					U[i][j] = 1;
				else
				{
					U[i][j] = A[i][j] / L[i][i];
					for (k = 0; k < i; k++)
					{
						U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
					}
				}
			}
		}
	}
	/**
	 *
	 */
	Matrix& operator+=(const Matrix& rhs){
		if(rhs._rows != this->_rows){
			throw  MatrixException();
		}
		T ** union_matrix = new T*[this->_rows];
		size_t count = this->_cols+rhs._cols;
		for(unsigned int i=0; i < this->_rows; ++i){
			union_matrix[i]= new T[count];
		}
		for(unsigned int i=0; i < this->_rows; ++i){
			for(unsigned int j=0; j<count; ++j){
				union_matrix[i][j]= (j<this->_rows)? this->M[i][j] : rhs.M[i][j-this->_rows];
			}
		}
		delete this->M;
		this->M = union_matrix;
		this->_cols=count;
		return * this;
	}
	

	friend Matrix operator+(Matrix lhs, const Matrix& rhs){
		lhs += rhs; // reuse compound assignment
		return lhs; // return the result by value (uses move constructor)
	}
};
}
#endif