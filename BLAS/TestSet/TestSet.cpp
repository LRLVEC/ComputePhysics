#include <_BLAS.h>
#include <_Time.h>
#include <random>

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(-1.0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);
	Timer timer;

	using namespace BLAS;

	//vec ta({ 1,2 });
	//mat tb({ {1,2},{3,4} });
	//mat tc({ {1,2},{3,4} });
	//tb(ta).print();
	//ta(tb).print();

	//vec jj({ 1,2,3,4 });
	//mat gg({ {0,1},{0,0,1},{0,0,0,1},{1} });
	//mat matSp(MatType::SparseMat, 4);
	//matSp.data[0] = 1;
	//matSp.data[1] = 1;
	//matSp.data[2] = 1;
	//matSp.data[3] = 1;
	//matSp.rowIndice[0] = 0;
	//matSp.rowIndice[1] = 1;
	//matSp.rowIndice[2] = 2;
	//matSp.rowIndice[3] = 3;
	//matSp.colIndice[0] = 1;
	//matSp.colIndice[1] = 2;
	//matSp.colIndice[2] = 3;
	//matSp.colIndice[3] = 0;
	//gg(jj).print();
	//matSp(jj).print();

	vecCplx vecCplxA(3);
	vecCplx vecCplxB(3);
	//vecCplx vecCplxC(1024 * 1024);
	matCplx matCplxA(MatType::SparseMat, 4, 4);

	matCplxA.re.data[0] = 1;
	matCplxA.re.data[1] = 2;
	matCplxA.re.data[2] = 3;
	matCplxA.re.data[3] = 4;
	matCplxA.im.data[0] = 8;
	matCplxA.im.data[1] = 7;
	matCplxA.im.data[2] = 6;
	matCplxA.im.data[3] = 5;

	matCplxA.re.rowIndice[0] = 0;
	matCplxA.re.rowIndice[1] = 1;
	matCplxA.re.rowIndice[2] = 2;
	matCplxA.re.rowIndice[3] = 3;
	matCplxA.im.rowIndice[0] = 0;
	matCplxA.im.rowIndice[1] = 1;
	matCplxA.im.rowIndice[2] = 2;
	matCplxA.im.rowIndice[3] = 3;

	matCplxA.re.colIndice[0] = 0;
	matCplxA.re.colIndice[1] = 1;
	matCplxA.re.colIndice[2] = 2;
	matCplxA.re.colIndice[3] = 3;
	matCplxA.im.colIndice[0] = 1;
	matCplxA.im.colIndice[1] = 2;
	matCplxA.im.colIndice[2] = 3;
	matCplxA.im.colIndice[3] = 0;

	matCplxA.printSparse();

	randomVecCplx(vecCplxA, mt, rduint);
	//randomVecCplx(vecCplxB, mt, rduint);

	vecCplxA.print();
	//vecCplxB.print();
	timer.begin();
	//vecCplxA /= vecCplxB;
	//cplx t(vecCplxA.normSquare());
	matCplxA(vecCplxA, vecCplxB);
	timer.end();
	timer.print();
	vecCplxB.print();

	//vecCplxA.print();


	//unsigned long long bd(64);
	//unsigned long long h(64 * 64 - 1);
	//mat matB(bd, h, MatType::BandMat, false);
	//mat matLB(bd, h, MatType::LBandMat, false);
	//mat matUB(bd, h, MatType::UBandMat, false);
	//mat matSym(h, h, false);
	//mat matCho(h, h, false);
	//vec vv(h, false);
	//randomMatBandL(matLB, mt, rd, 1.0 / bd);
	//randomMatBandU(matUB, mt, rd);
	//transLBandToSymmetricMat(matLB, matSym);
	//transLBandToSymmetricBandMat(matLB, matB);
	//transNormalMatToBand(matSym, matB);
	//transBandToNormalMat(matUB, matU);
	//randomVec(vv, mt, rd);

	//vec rr(h, false);
	//vec ans(h, false);
	//vec ansB(h, false);
	//vec tp(h, false);
	//matB(vv, rr);
	//::printf("started\n");

	//timer.begin();
	//matLB.solveCholeskyBand(rr, ansB);
	//timer.end();
	//tp = vv; tp -= ansB;
	//printf("Cholesky (band):\t\t%e", tp.norm2());
	//timer.print();

	//timer.begin();
	//matB.solveConjugateGradient(rr, ansB, 1.0e-10);
	//timer.end();
	//tp = vv; tp -= ansB;
	//printf("Conjugate Gradient (band):\t%e", tp.norm2());
	//timer.print();

	//timer.begin();
	//matB.solveSteepestDescent(rr, ansB, 5.0e-9);
	//timer.end();
	//tp = vv; tp -= ansB;
	//printf("Steepest Descent (band):\t%e", tp.norm2());
	//timer.print();

	//timer.begin();
	//matSym.solveSteepestDescent(rr, ansB, 1e-14);
	//timer.end();
	//tp = vv; tp -= ansB;
	//printf("Steepest Descent:\t\t%e", tp.norm2());
	//timer.print();

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

	//vec vecA(128 * 128);
	//vec vecB(128 * 128);

	//vec vecC(h, false);
	//vec vecD(h, false);
	//vec vecE(h, false);
	//mat matA(h, h, false);
	//mat matB(1024, 1024, false);
	//mat matC(7, 7, false);
	//mat matD(64, 3, false);
	//randomVec(vecA, mt, rd);
	//randomVec(vecB, mt, rd);
	//randomVec(vecC, mt, rd);
	//randomVec(vecE, mt, rd);
	//randomMatSymmetric(matA, mt, rd);
	//randomMat(matB, mt, rd);
	//randomMat(matC, mt, rd);

	//vec
	/*{
		timer.begin();
		for (unsigned long long c0(0); c0 < 100; ++c0)
			vecA += vecB;
		timer.end();
		timer.print("vec add:");

		timer.begin();
		for (unsigned long long c0(0); c0 < 100; ++c0)
			vecA -= vecB;
		timer.end();
		timer.print("vec minus:");

		timer.begin();
		for (unsigned long long c0(0); c0 < 100; ++c0)
			vecA *= vecB;
		timer.end();
		timer.print("vec mult:");

		timer.begin();
		for (unsigned long long c0(0); c0 < 100; ++c0)
			vecA /= vecB;
		timer.end();
		timer.print("vec div:");

		timer.begin();
		for (unsigned long long c0(0); c0 < 100; ++c0)
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
	}*/

	//mat

	//matA(vecE, vecC);
	//matA.printToTableTxt("./matA.txt");
	//vecC.printToTableTxt("./vecC.txt", false);
	//vecE.printToTableTxt("./vecE.txt", false);

	//timer.begin();
	//matA.solveSteepestDescent(vecC, vecD, 1e-20);
	//timer.end();
	////matA.printToTableTxt("./matACho.txt");

	//vec delta(vecE - vecD);
	//::printf("solveSteepestDescent delta norm:%e \t", delta.norm2());
	//timer.print();
	//vecD.printToTableTxt("./vecD.txt", false);
	//vecE.printToTableTxt("./vecE.txt", false);
	//delta.printToTableTxt("./delta.txt", false);

	//timer.begin();
	//for (unsigned long long c0(0); c0 < 100; ++c0)
	//	matA(vecC);
	//timer.end();
	//timer.print("mat mult vec:");

	//timer.begin();
	//for (unsigned long long c0(0); c0 < 100; ++c0)
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