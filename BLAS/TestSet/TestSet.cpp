#include <_BLAS.h>

int main()
{
	BLAS::vec a({1,2,3,4,5});
	::printf("%f\n", a[2]);
	a += a;
	a.print();
	a.printInfo();
	(a / a).print();
	printf("%f\n", a.normInf());

	BLAS::mat b({ {1,2,3}, {4,5,6}, {7,8,9} });
	BLAS::mat c;
	b.print();
	b.printInfo();
	c = b;
	c.print();
}