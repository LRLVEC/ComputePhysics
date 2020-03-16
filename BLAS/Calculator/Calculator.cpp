#include <_BLAS.h>
#include <_Time.h>

int main()
{
	using namespace BLAS;
	mat L1{ {1},{2,1},{3,0,1} };
	mat L2{ {1},{0,1},{0,2,1} };
	L1(L2).print();
}