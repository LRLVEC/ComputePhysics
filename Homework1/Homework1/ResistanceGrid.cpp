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
	unsigned long long bd(64);
	unsigned long long h(64 * 64 - 1);
	mat matB(bd, h, MatType::BandMat, false);
	mat matLB(bd, h, MatType::LBandMat, false);
	vec vv(h, false);
	randomMatBandL(matLB, mt, rd, 1.0 / bd);
	transLBandToSymmetricBandMat(matLB, matB);
	randomVec(vv, mt, rd);

	vec rr(h, false);
	vec ans(h, false);
	vec ansB(h, false);
	vec tp(h, false);
	matB(vv, rr);
	::printf("started\n");

	//timer.begin();
	//matLB.solveCholeskyBand(rr, ansB);
	//timer.end();
	//tp = vv; tp -= ansB;
	//printf("Cholesky (band):\t\t%e", tp.norm2());
	//timer.print();

	timer.begin();
	matB.solveConjugateGradient(rr, ansB, 1.0e-10);
	timer.end();
	tp = vv; tp -= ansB;
	printf("Conjugate Gradient (band):\t%e", tp.norm2());
	timer.print();

	timer.begin();
	matB.solveSteepestDescent(rr, ansB, 5.0e-9);
	timer.end();
	tp = vv; tp -= ansB;
	printf("Steepest Descent (band):\t%e", tp.norm2());
	timer.print();
}