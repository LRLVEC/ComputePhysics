#include <_BLAS.h>
#include <_Time.h>

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(-1.0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);
	Timer timer;

	using namespace BLAS;

	mat a(9, 9, false);
	randomMatSymmetric(a, mt, rd, 0.01);
	a.printToTableTxt("E:\\files\\C++\\ComputePhysics\\BLAS\\EigenvalueTest\\a.txt");
	a.print();
	mat t(a.symmetricMatDiagonalization());
	t.print();
}