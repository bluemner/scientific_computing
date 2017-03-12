#ifndef _BETACORE_GOLDEN_SECTION_SEARCH_H_
#define _BETACORE_GOLDEN_SECTION_SEARCH_H_
#include <iostream>
#include <cmath>
namespace betacore{
	template<typename T>
	class GoldenSectionSearch{
		private:

		public:
		static void search(T a, T b, T tolerance){
			T t = sqrt(5) - 1 /2;
			T x1 = a +  (1 - t) ( b - 1 );
			T x2 = a + t * (b - a);
			T f1 = x1;
			T f2 = x2;
			while((b-a) > tolerance) {
				if( f1 > f2){
					a = x1;
					x1 = x2;
					f1 = f2;
					x2 =  a + t (b-a);
					f2 = x2;
				}else{
					b = x2;
					x2 = x1;
					f2 = f1;
					x1 = a + (1 - t )(b-a);
					f1 = x1;
				}
				std::cout<< "a:" <<  a << "\tb:" << b << "\tx1:" << x1 << "\tx2" << x2 << std::endl;
			}
		}
	}

}
#endif