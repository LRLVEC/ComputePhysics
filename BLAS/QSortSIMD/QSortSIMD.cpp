#include <_BLAS.h>
#include <_Time.h>

using namespace BLAS;

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(-1.0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);
	Timer timer;

	vec a(1024 * 1024, false);
	randomVec(a, mt, rd);
	timer.begin();
	//a.qsort();
	a.qsortAVX();
	timer.end();
	timer.print();
	//a.print();
}