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
			a(c0, c1) = 0.01 * rd(mt);
		a(c0, c1++) = 1 + 0.1 * rd(mt);
		for (; c1 < a.width; ++c1)
			a(c0, c1) = 0.01 * rd(mt);
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
template<class T>void randomMatSymmetric(BLAS::mat& a, std::mt19937& mt, T& rd)
{
	//make sure square matrix!
	for (unsigned int c0(0); c0 < a.height; ++c0)
	{
		unsigned int c1(0);
		for (; c1 < a.width && c1 < c0; ++c1)
			a(c1, c0) = a(c0, c1) = 0.3 * rd(mt);
		a(c0, c1++) = 1 + 0.2 * rd(mt);
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
template<class T>void randomMatBandL(BLAS::mat& a, std::mt19937& mt, T& rd)
{
	a.clear();
	for (unsigned int c0(0); c0 < a.height; ++c0)
	{
		unsigned int c1(a.LBandBeginOffset(c0));
		unsigned int ending(c0 <= a.halfBandWidth ? c0 + 1 : a.halfBandWidth + 1);
		ending += c1;
		for (; c1 < ending - 1; ++c1)
			a.data[c0 * a.width4d + c1] = 0.1 * rd(mt);
		a.data[c0 * a.width4d + c1] = 1 + 0.1 * rd(mt);
	}
}
template<class T>void randomMatBandU(BLAS::mat& a, std::mt19937& mt, T& rd)
{
	a.clear();
	for (unsigned int c0(0); c0 < a.height; ++c0)
	{
		unsigned int c1(c0 % 4);
		unsigned int ending(a.height - c0 <= a.halfBandWidth ? a.height - c0 : a.halfBandWidth + 1);
		ending += c1;
		a.data[c0 * a.width4d + c1++] = 1 + 0.1 * rd(mt);
		for (; c1 < ending; ++c1)
			a.data[c0 * a.width4d + c1] = 0.1 * rd(mt);
	}
}
BLAS::mat& transBandToNormalMat(BLAS::mat const& a, BLAS::mat& b)
{
	if (b.width != a.height || b.height != a.height)
		b.reconstruct(a.height, a.height, false);
	b.clear();
	if (a.matType == BLAS::MatType::LBandMat)
	{
		for (unsigned int c0(0); c0 < a.height; ++c0)
		{
			unsigned int bgn(a.LBandBeginOffset(c0));
			unsigned int ost(c0 <= a.halfBandWidth ? c0 : a.halfBandWidth);
			memcpy(b.data + c0 * b.width4d + c0 - ost,
				a.data + c0 * a.width4d + bgn, (ost + 1) * sizeof(double));
		}
	}
	else if (a.matType == BLAS::MatType::UBandMat)
	{
		for (unsigned int c0(0); c0 < a.height; ++c0)
		{
			unsigned int ost(a.height - c0 <= a.halfBandWidth ? a.height - c0 : a.halfBandWidth + 1);
			memcpy(b.data + c0 * b.width4d + c0,
				a.data + c0 * a.width4d + c0 % 4, ost * sizeof(double));
		}
	}
	return b;
}
BLAS::mat& transLBandToSymmetricMat(BLAS::mat const& a, BLAS::mat& b)
{
	if (b.width != a.height || b.height != a.height)
		b.reconstruct(a.height, a.height, false);
	b.clear();
	for (unsigned int c0(0); c0 < a.height; ++c0)
	{
		unsigned int bgn(a.LBandBeginOffset(c0));
		unsigned int ost(c0 <= a.halfBandWidth ? c0 : a.halfBandWidth);
		memcpy(b.data + c0 * b.width4d + c0 - ost,
			a.data + c0 * a.width4d + bgn, (ost + 1) * sizeof(double));
	}
	for (unsigned int c0(0); c0 < a.height - 1; ++c0)
		for (unsigned int c1(c0 + 1); c1 < a.height; ++c1)
			b.data[c0 * b.width4d + c1] = b.data[c1 * b.width4d + c0];
	return b;
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

	unsigned int bd(128);
	unsigned int h(129 * 129 - 1);
	mat matLB(bd, h, MatType::LBandMat, false);
	//mat matUB(bd, h, MatType::UBandMat, false);
	mat matSym(h, h, false);
	mat matCho(h, h, false);
	vec vv(h, false);
	randomMatBandL(matLB, mt, rd);
	//randomMatBandU(matUB, mt, rd);
	transLBandToSymmetricMat(matLB, matSym);
	//transBandToNormalMat(matUB, matU);
	randomVec(vv, mt, rd);
	//matUB.print();
	//matU.print();
	//vec s, r;
	//matLB(vv, s);
	//matL(vv, r);
	//vec det(s - r);
	//::printf("differ LBand: %f\n", det.norm2());
	//matUB(vv, s);
	//matU(vv, r);
	//det = s;
	//det -= r;
	//::printf("differ UBand: %f\n", det.norm2());

	//timer.begin();
	//for (unsigned int c0(0); c0 < 100; ++c0)
	//	matLB(vv);
	//timer.end();
	//timer.print("Lband mat mult vec:");

	//timer.begin();
	//for (unsigned int c0(0); c0 < 100; ++c0)
	//	matL(vv);
	//timer.end();
	//timer.print("mat mult vec:");

	vec rr(h, false);
	vec ans(h, false);
	vec ansB(h, false);
	vec tp(vv);
	matSym(vv, rr);

	timer.begin();
	matLB.solveCholeskyBand(rr, ansB);
	timer.end();
	timer.print("solveCholeskyBand:");

	transLBandToSymmetricMat(matLB, matCho);
	tp -= ansB;
	printf("ans - solveCholeskyBand: %e\n", tp.norm2());
	//tp.printToTableTxt("./det.txt", false);
	//ansB.printToTableTxt("./ans.txt", false);
	//matSym.printToTableTxt("./mat.txt");
	//matCho.printToTableTxt("./matCho.txt");

	//timer.begin();
	//matSym.solveL(rr, ans);
	//timer.end();
	//timer.print("solveL:");
	//tp -= ans;
	//printf("ans - solveL: %e\n", tp.norm2());

	//timer.begin();
	//matLB.solveL(rr, ansB);
	//timer.end();
	//timer.print("solveL band:");
	//tp = vv;
	//tp -= ansB;
	//printf("ans - solveL band: %e\n", tp.norm2());
	//tp = ans;
	//tp -= ansB;
	//printf("solveL - solveL band: %e\n", tp.norm2());

	//matU(vv, rr);

	//timer.begin();
	//matU.solveU(rr, ans);
	//timer.end();
	//timer.print("solveU:");
	//tp = vv; tp -= ans;
	//printf("ans - solveU: %e\n", tp.norm2());

	//timer.begin();
	//matUB.solveU(rr, ansB);
	//timer.end();
	//timer.print("solveU band:");
	//tp = vv; tp -= ansB;
	//printf("ans - solveU band: %e\n", tp.norm2());
	//tp = ans;
	//tp -= ansB;
	//printf("solveU - solveU band: %e\n", tp.norm2());
	//tp.printToTableTxt("./det.txt", false);

	vec vecA(128 * 128);
	vec vecB(128 * 128);

	vec vecC(1025, false);
	vec vecD(1025, false);
	vec vecE(1025, false);
	mat matA(1025, 1025, false);
	//mat matB(1024, 1024, false);
	//mat matC(7, 7, false);
	//mat matD(64, 3, false);
	randomVec(vecA, mt, rd);
	randomVec(vecB, mt, rd);
	//randomVec(vecC, mt, rd);
	randomVec(vecE, mt, rd);
	randomMatSymmetric(matA, mt, rd);
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

	//matA.solveCholeskyBand(vecC, vecD);
	timer.begin();
	matA.solveCholesky(vecC, vecD);
	timer.end();
	//matA.printToTableTxt("./matACho.txt");

	//::printf("%d\n", int(1u - 4u) / 4);

	vec delta(vecE - vecD);
	::printf("solveCholesky delta norm:%e \t", delta.norm2());
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