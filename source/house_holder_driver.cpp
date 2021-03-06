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
#include "../headers/householder.hpp"

typedef double T;

int main(int argc, char * argv[]){
	T data[5][4] = {
						{ 2.0,  1.0, 3.0, 5.0},
						{-1.0,  0.0, 7.0, 1.0},
						{ 0.0, -1.0, -1.0, 3.0},
						{-3.0,  7.0, 4.0, 3.0},
						{ 1.0,  6.0, 4.0, -3.0}
					};
	betacore::Householder<T> HH(data);
	HH.run();

	// T ex2[5][3] = {
	// 					{1,-1.0,1.0},
	// 					{1,-0.5,0.25},
	// 					{1, 0.0,0.0},
	// 					{1, 0.5,0.24},
	// 					{1, 1.0,1.0}
	// 				};
	// betacore::Householder<T> HH2(ex2);
	// HH2.run();
	return 0;
}