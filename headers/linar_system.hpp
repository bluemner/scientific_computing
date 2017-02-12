#ifndef _BETACORE_LINAR_SYSTEM_H_
#define _BETACORE_LINAR_SYSTEM_H_

#ifndef _BETACORE_MATRIX_H_
#include "./matrix.hpp"
#endif
#include <exception>
#include <iostream>
#include <iomanip>
namespace betacore{
	struct Linar_System_Exception : public std::exception {
   const char * what () const throw () {
      return "Linar System Exception";
   }
};
	template <typename T>
	class Linar_System{
		private:
			Matrix<T> * A;
			size_t size;
			T* x;
			T* b;
			template<size_t rows, size_t cols>
			void copy(T (&A)[rows][cols], T (&x)[rows], T (&b)[rows]){
				this->size= rows;
				this->A = new Matrix<T>(A);
				this->x = new T[rows];
				this->b = new T[rows];
				size_t i=0;
				//vectors have the same size;
				for(; i< rows; i++){
					this->x[i] = x[i];
					this->b[i] = b[i];
				}
			}
		public:
			Linar_System(){}
			template<size_t rows, size_t cols>
			Linar_System(T (&A)[rows][cols], T (&x)[rows], T (&b)[rows]){
				if(rows != cols){
					throw Linar_System_Exception();	
				}
				copy(A,x,b);
			}
			~Linar_System(){
				delete this->A;
				this->A = nullptr;
				delete this->x;
				this->x = nullptr;
				delete this->b;
				this->b = nullptr;
			}
			void forward_substitution(){
				size_t i=0;
				size_t j=0;
			

				size_t n = this->size;
				for ( ; j < n ; j++ ){
					if (this->A->at(i,j) == 0){
						break;
					}
					std::cout<<"i:"<<i<<"\tj:"<<j<<std::endl;
					this->x[j] = this->b[j] / this->A->at(i,j);//[i][j];
					for (i=j+1; i < j-1 ; i++  ){
						this->b[i]= this->b[i] - this->A->at(i,j) * this->x[j];
					}
				}

			}
						
			void backward_substitution( ){

				int n = (int) this->size;
				int i=0;
				int j=this->size-1;
				for ( ; j >=0 ; j-- ){					
					size_t uj = (size_t) j;
					T Ajj=this->A->at(uj,uj);				
					if (Ajj == 0.0){					
						break; 
					}
					T xj = this->b[j] / Ajj;
					this->x[j] =xj;				
					for (i=0; i<=(int)j-1 ; i++  ){
						size_t ui = (size_t) i;
						this->b[i]= this->b[i] -  this->A->at(ui,uj) * this->x[j];
					}
					
				}
			}
			void print(){
			
				size_t i = 0;
				size_t j = 0;
				size_t precision = 5;
				size_t width = 10 + precision;
				char fill =' ';
				std::cout<< "┌";
				for(i=0; i< this->size; i++){
					std::cout<<" ";
					for(j=0; j< width; ++j){
						std::cout<<" ";
					}			
					
				}
				std::cout<<" ┐   ┌ ";
				for(j=0; j< width; ++j){
						std::cout<<" ";
					}
				std::cout<<" ┐   ┌ ";
				for(j=0; j< width; ++j){
						std::cout<<" ";
					}		
				std::cout<<" ┐"<<std::endl;

				for(i=0; i< this->size; i++){
					
						std::cout<<"│";
					
					for(j=0; j< this->size; j++){
						std::cout<<" "<<std::fixed<<std::setprecision(precision)<<std::setfill(fill)<<std::setw(width);
						std::cout<<this->A->at(i,j);
					}					
					
					std::cout<<" │";
					if( i == (this->size/2)   ){
						std::cout<<" X "<<"│"<<" "<<std::setfill(fill) << std::setw(width)<<this->x[i]<<" "<<"│";
						std::cout<<" = "<<"│"<<" "<< std::setw(width)<< this->b[i]<<" "<<"│";
					}else{
						std::cout<<"   "<<"│"<<" "<< std::setw(width)<< this->x[i]<<" "<<"│";
						std::cout<<"   "<<"│"<<" "<< std::setw(width)<< this->b[i]<<" "<<"│";
					}
							
					std::cout<<std::endl;
				}


				std::cout<< "└";
				for(i=0; i< this->size; i++){	
					std::cout<<" ";				
					for(j=0; j< width; ++j){
						std::cout<<" ";
					}	
				}
				std::cout<<" ┘   └ ";
				for(j=0; j< width; ++j){
						std::cout<<" ";
					}	
				std::cout<<" ┘   └ ";
				for(j=0; j< width; ++j){
						std::cout<<" ";
					}	
				std::cout<<" ┘"<<std::endl;
				for(i=0; i< this->size; i++){
					std::cout << "("<<x[i]<<","<<b[i]<<")";
					
				}
				std::cout<<std::endl;
			}
			
	};
}



#endif