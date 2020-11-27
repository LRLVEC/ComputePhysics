#include <_BLAS.h>
#include <_Time.h>

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(-1.0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);
	Timer timer;

	using namespace BLAS;

	mat a(8, 8, false);
	randomMatSymmetric(a, mt, rd, 0.5);
	a.printToTableTxt("E:\\files\\C++\\ComputePhysics\\BLAS\\EigenvalueTest\\a.txt");
	a.print();
	timer.begin();
	mat tH(a.tridiagonalizationHouseholder());
	mat tL(a.tridiagonalizationLanczos());
	timer.end();
	timer.print();

	tH.print();
	tL.print();

	timer.begin();
	vec lambdasH(8, false);
	vec lambdasL(8, false);
	tH.implicitSymmetricQR(1e-20, lambdasH);
	tL.implicitSymmetricQR(1e-20, lambdasL);
	timer.end();
	timer.print();
	//t.print();
	lambdasH.print();
	lambdasL.print();
	lambdasH.printToTableTxt("E:\\files\\C++\\ComputePhysics\\BLAS\\EigenvalueTest\\b.txt", false);
}