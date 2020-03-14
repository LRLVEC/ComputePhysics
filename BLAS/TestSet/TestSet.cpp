#include <_BLAS.h>
#include <_Time.h>
#include <random>

template<class T>void randomVec(BLAS::vec& a, std::mt19937& mt, T& rd)
{
	for (unsigned int c0(0); c0 < a.dim; ++c0)
		a.data[c0] = rd(mt);
}
template<class T>void randomMat(BLAS::mat& a, std::mt19937& mt, T& rd)
{
	for (unsigned int c0(0); c0 < a.width * a.height; ++c0)
		a.data[c0] = rd(mt);
}

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(0, 2);
	std::uniform_int_distribution<unsigned int> rduint(1, 10);
	Timer timer;

	using namespace BLAS;

	vec ta({ 1,2 });
	mat tb({ {1,2},{3,4} });
	mat tc({ {1,2},{3,4} });
	tb(ta).print();
	ta(tb).print();

	vec vecA(128 * 128);
	vec vecB(128 * 128);
	vec vecC(1024);
	randomVec(vecA, mt, rd);
	randomVec(vecB, mt, rd);
	randomVec(vecC, mt, rd);

	timer.begin();
	for (unsigned int c0(0); c0 < 100; ++c0)
		vecA += vecB;
	timer.end();
	timer.print("vec add:");

	timer.begin();
	for (unsigned int c0(0); c0 < 100; ++c0)
		vecA -= vecB;
	timer.end();
	timer.print("vec minus:");

	timer.begin();
	for (unsigned int c0(0); c0 < 100; ++c0)
		vecA *= vecB;
	timer.end();
	timer.print("vec mult:");

	timer.begin();
	for (unsigned int c0(0); c0 < 100; ++c0)
		vecA /= vecB;
	timer.end();
	timer.print("vec div:");

	timer.begin();
	double dots;
	dots = (vecA, vecB);
	timer.end();
	::printf("vec dot:%.4f\t", dots);
	timer.print();

	timer.begin();
	double norm(vecA.norm1());
	timer.end();
	::printf("vec norm1:%.4f\t", norm);
	timer.print();

	timer.begin();
	norm = vecA.norm2();
	timer.end();
	::printf("vec norm2:%.4f\t", norm);
	timer.print();

	timer.begin();
	norm = vecA.normInf();
	timer.end();
	::printf("vec normInf:%.4f\t", norm);
	timer.print();

	mat matA(1024, 1024, false);
	mat matB(1024, 1024, false);
	mat matC(3, 2, false);
	mat matD(64, 3, false);
	randomMat(matA, mt, rd);
	randomMat(matB, mt, rd);

	timer.begin();
	matA += matB;
	timer.end();
	timer.print("mat add:");

	timer.begin();
	matA -= matB;
	timer.end();
	timer.print("mat minus:");

	timer.begin();
	matA *= matB;
	timer.end();
	timer.print("mat mult:");

	timer.begin();
	matA /= matB;
	timer.end();
	timer.print("mat div:");

	randomMat(matA, mt, rd);
	randomMat(matB, mt, rd);

	timer.begin();
	vec vecD(matA(vecC));
	timer.end();
	timer.print("mat mult vec:");

	/*timer.begin();
	mat matE(matA(matB));
	timer.end();
	timer.print("mat mult mat:");*/

	randomMat(matC, mt, rd);
	randomMat(matD, mt, rd);

	timer.begin();
	mat matE(matC(matD));
	timer.end();
	timer.print("mat mult mat:");

	matC.printToTxt("./matC.txt");
	matD.printToTxt("./matD.txt");
	matE.printToTxt("./matE.txt");

	//::printf("%f\n", a[2]);
	//a += a;
	//a.print();
	//a.printInfo();
	//(a / a).print();
	//printf("%f\n", a.normInf());

}