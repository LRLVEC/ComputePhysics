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
	for (unsigned int c0(0); c0 < a.height; ++c0)
		for (unsigned int c1(0); c1 < a.width; ++c1)
			a(c0, c1) = rd(mt);
}
template<class T>void randomMatGood(BLAS::mat& a, std::mt19937& mt, T& rd)
{
	for (unsigned int c0(0); c0 < a.height; ++c0)
	{
		unsigned int c1(0);
		for (; c1 < a.width && c1 < c0; ++c1)
			a(c0, c1) = 0.1 * rd(mt);
		a(c0, c1++) = 1 + 0.01 * rd(mt);
		for (; c1 < a.width; ++c1)
			a(c0, c1) = 0.1 * rd(mt);
	}
}
template<class T>void randomMatL(BLAS::mat& a, std::mt19937& mt, T& rd)
{
	for (unsigned int c0(0); c0 < a.height; ++c0)
	{
		unsigned int c1(0);
		for (; c1 < a.width && c1 < c0; ++c1)
			a(c0, c1) = 0.1 * rd(mt);
		a(c0, c1++) = 1 + 0.2 * rd(mt);
		for (; c1 < a.width; ++c1)
			a(c0, c1) = 0;
	}
}
template<class T>void randomMatU(BLAS::mat& a, std::mt19937& mt, T& rd)
{
	for (unsigned int c0(0); c0 < a.height; ++c0)
	{
		unsigned int c1(0);
		for (; c1 < a.width && c1 < c0; ++c1)
			a(c0, c1) = 0;
		a(c0, c1++) = 1 + 0.2 * rd(mt);
		for (; c1 < a.width; ++c1)
			a(c0, c1) = 0.1 * rd(mt);
	}
}
template<class T>void randomMatBandL(BLAS::mat& a, std::mt19937& mt, T& rd, unsigned band)
{
	for (unsigned int c0(0); c0 < a.height; ++c0)
	{
		unsigned int start(int(c0 - band) > 0 ? c0 - band : 0);
		unsigned int c1(0);
		for (; c1 < start; ++c1)
			a(c0, c1) = 0;
		for (; c1 <= c0; ++c1)
			a(c0, c1) = rd(mt);
		for (; c1 < a.width; ++c1)
			a(c0, c1) = 0;
	}
}

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(-1.0, 1.0);
	std::uniform_int_distribution<unsigned int> rduint(1, 10);
	Timer timer;

	using namespace BLAS;

	//vec ta({ 1,2 });
	//mat tb({ {1,2},{3,4} });
	//mat tc({ {1,2},{3,4} });
	//tb(ta).print();
	//ta(tb).print();

	vec vecA(128 * 128);
	vec vecB(128 * 128);

	vec vecC(512, false);
	vec vecD(512, false);
	vec vecE(512, false);
	mat matA(512, 512, false);
	//mat matB(1024, 1024, false);
	//mat matC(7, 7, false);
	//mat matD(64, 3, false);
	randomVec(vecA, mt, rd);
	randomVec(vecB, mt, rd);
	//randomVec(vecC, mt, rd);
	randomVec(vecE, mt, rd);
	randomMatGood(matA, mt, rd);
	//randomMat(matB, mt, rd);
	//randomMat(matC, mt, rd);

	//vec
	{
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
		for (unsigned int c0(0); c0 < 100; ++c0)
			vecA.fmadd(1.0, vecB);
		timer.end();
		timer.print("vec fmadd:");

		randomVec(vecA, mt, rd);
		randomVec(vecB, mt, rd);

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

		timer.begin();
		norm = vecA.normP(3.6);
		timer.end();
		::printf("vec normP:%.4f\t", norm);
		timer.print();
	}

	//mat

	matA(vecE, vecC);
	//matA.printToTableTxt("./matA.txt");
	//vecC.printToTableTxt("./vecC.txt", false);
	//vecE.printToTableTxt("./vecE.txt", false);

	timer.begin();
	matA.solveJacobiIter(vecC, vecD, 1e-5);
	timer.end();

	vec delta(vecE - vecD);
	::printf("solveJacobiIter delta norm:%e \t", delta.norm2());
	timer.print();
	//vecD.printToTableTxt("./vecD.txt", false);
	//delta.printToTableTxt("./delta.txt", false);

	//timer.begin();
	//for (unsigned int c0(0); c0 < 100; ++c0)
	//	matA(vecC);
	//timer.end();
	//timer.print("mat mult vec:");

	//timer.begin();
	//for (unsigned int c0(0); c0 < 100; ++c0)
	//	vecC(matA);
	//timer.end();
	//timer.print("vec mult mat:");

	//timer.begin();
	//mat matE(matA(matB));
	//timer.end();
	//timer.print("mat mult mat:");

	/*timer.begin();
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
	timer.print("mat div:");*/

	//timer.begin();
	//mat matE(matA(matB));
	//timer.end();
	//timer.print("mat mult mat:");

	//timer.begin();
	//vec vecE(matA(vecC));
	//timer.end();
	//timer.print("mat mult vec:");

	//timer.begin();
	//vec vecE(vecC(matA));
	//timer.end();
	//timer.print("vec mult mat:");

	//timer.begin();
	//vec vecE(matA(vecD));
	//timer.end();
	//timer.print("mat mult vec:");


	//matA.printToTableTxt("./matA.txt");
	//matB.printToTableTxt("./matB.txt");
	//matE.printToTableTxt("./matE.txt");
	//matD.printToTxt("./matD.txt");
	//matE.printToTxt("./matE.txt");
	//vecA.printToTableTxt("./vecA.txt");
	//vecB.printToTableTxt("./vecB.txt");

	//timer.begin();
	//vecA += vecB;
	//timer.end();
	//timer.print("vec add vec:");

	//vecD.printToTableTxt("./vecC.txt", false);
	//vecE.printToTableTxt("./vecD.txt", false);
	//vecE.printToTableTxt("./vecE.txt");
}
/*matA := Import["D:\\files\\C++\\ComputePhysics\\BLAS\\TestSet\\matA.txt", "Table"];
vecC := Import["D:\\files\\C++\\ComputePhysics\\BLAS\\TestSet\\vecC.txt", "Table"];
vecD := Import["D:\\files\\C++\\ComputePhysics\\BLAS\\TestSet\\vecD.txt", "Table"]*/