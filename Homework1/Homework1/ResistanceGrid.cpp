#include <_BLAS.h>
#include <_Time.h>
#include <random>

using namespace BLAS;

template<unsigned long long _dim>struct SquareGrid
{
	static constexpr unsigned long long dim = _dim + 1;

	mat matLBand;
	mat matBand;
	mat matSparse;
	vec u;
	vec i;
	unsigned long long clipA;
	unsigned long long clipB;
	unsigned long long clipID;
	double rCholesky;
	double rSteepestDescent;
	double rConjugateGradient;
	double rConjugateGradientSparse;
	Timer timer;
	SquareGrid(unsigned long long _clipA, unsigned long long _clipB)
		:
		matLBand(dim, dim* dim - 1, MatType::LBandMat, true),
		matBand(dim, dim* dim - 1, MatType::BandMat, true),
		matSparse(MatType::SparseMat, 5 * (dim * dim - 1)),
		u(dim* dim - 1, false),
		i(dim* dim - 1, true),
		clipA(_clipA),
		clipB(_clipB),
		clipID(_clipA* dim + _clipB),
		rCholesky(0),
		rSteepestDescent(0),
		rConjugateGradient(0),
		rConjugateGradientSparse(0)
	{
		setGrid();
	}
	double sumG(unsigned long long a, unsigned long long b)const
	{
		double s(4);
		if (a == 0)s--; if (a == _dim)s--;
		if (b == 0)s--; if (b == _dim)s--;
		return s;
	}
	void setGrid()
	{
		unsigned long long cnt(0);
		for (unsigned long long c0(0); c0 < dim; ++c0)
			for (unsigned long long c1(0); c1 < dim; ++c1)
			{
				unsigned long long n(c0 * dim + c1);
				if (n < clipID)
				{
					if (c0)
					{
						matLBand.LBandEleRef(n, n - dim) = matBand.BandEleRef(n, n - dim) = -1;
						matSparse.addSparse(n, n - dim, -1, cnt);
					}
					if (c1)
					{
						matLBand.LBandEleRef(n, n - 1) = matBand.BandEleRef(n, n - 1) = -1;
						matSparse.addSparse(n, n - 1, -1, cnt);
					}
					double ss(sumG(c0, c1));
					matLBand.LBandEleRef(n, n) = matBand.BandEleRef(n, n) = ss;
					matSparse.addSparse(n, n, ss, cnt);
					if (c1 != _dim && n + 1 != clipID)
					{
						matBand.BandEleRef(n, n + 1) = -1;
						matSparse.addSparse(n, n + 1, -1, cnt);
					}
					if (c0 != _dim && n + dim != clipID)
					{
						unsigned long long ns;
						if (n + dim < clipID)ns = n;
						else ns = n - 1;
						matBand.BandEleRef(n, ns + dim) = -1;
						matSparse.addSparse(n, ns + dim, -1, cnt);
					}
				}
				else if (n > clipID)
				{
					unsigned long long nd(n - 1);
					unsigned long long ns;
					if (c0 && n - dim != clipID)
					{
						if (n - dim < clipID)ns = n;
						else ns = nd;
						matLBand.LBandEleRef(nd, ns - dim) = matBand.BandEleRef(nd, ns - dim) = -1;
						matSparse.addSparse(nd, ns - dim, -1, cnt);
					}
					if (c1 && nd != clipID)
					{
						matLBand.LBandEleRef(nd, nd - 1) = matBand.BandEleRef(nd, nd - 1) = -1;
						matSparse.addSparse(nd, nd - 1, -1, cnt);
					}
					double ss(sumG(c0, c1));
					matLBand.LBandEleRef(nd, nd) = matBand.BandEleRef(nd, nd) = ss;
					matSparse.addSparse(nd, nd, ss, cnt);
					if (c1 != _dim)
					{
						matBand.BandEleRef(nd, nd + 1) = -1;
						matSparse.addSparse(nd, nd + 1, -1, cnt);
					}
					if (c0 != _dim)
					{
						matBand.BandEleRef(nd, nd + dim) = -1;
						matSparse.addSparse(nd, nd + dim, -1, cnt);
					}
				}
			}
		matSparse.elementNum = cnt;
		i.data[0] = 1;
	}
	double solveCholesky()
	{
		u = 0;
		timer.begin();
		matLBand.solveCholeskyBand(i, u);
		timer.end();
		rCholesky = u.data[0];
		::printf("SquareGrid<%llu> Cholesky:\t\t%.15e", _dim, rCholesky);
		timer.print();
		return rCholesky;
	}
	double solveSteepestDescent(double _esp)
	{
		//does not converge!
		u = 0;
		timer.begin();
		matBand.solveSteepestDescent(i, u, _esp);
		timer.end();
		rSteepestDescent = u.data[0];
		::printf("SquareGrid<%llu> SteepestDescent:\t\t%.15e", _dim, rSteepestDescent);
		timer.print();
		return rSteepestDescent;
	}
	double solveConjugateGradient(double _esp)
	{
		u = 0;
		timer.begin();
		matBand.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradient = u.data[0];
		::printf("SquareGrid<%llu> ConjugateGradient:\t%.15e", _dim, rConjugateGradient);
		timer.print();
		return rConjugateGradient;
	}
	double solveConjugateGradientSparse(double _esp)
	{
		u = 0;
		timer.begin();
		matSparse.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradientSparse = u.data[0];
		::printf("SquareGrid<%llu> ConjugateGradientSparse:\t%.15e", _dim, rConjugateGradientSparse);
		timer.print();
		return rConjugateGradientSparse;
	}
};

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(-1.0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);
	Timer timer;

	SquareGrid<16>ahh(0, 16);
	ahh.solveCholesky();
	ahh.solveConjugateGradient(1e-10);
	ahh.solveConjugateGradientSparse(1e-10);
	SquareGrid<64>bhh(0, 64);
	bhh.solveCholesky();
	bhh.solveConjugateGradient(1e-10);
	bhh.solveConjugateGradientSparse(1e-10);
	//SquareGrid<256>chh(0, 256);
	//chh.solveCholesky();
	//chh.solveConjugateGradient(1e-3);
	//chh.solveConjugateGradientSparse(1e-10);
	//SquareGrid<512>dhh(512, 512);
	//dhh.solveCholesky();
	//dhh.solveConjugateGradient(1e-3);
	//dhh.solveConjugateGradientSparse(1e-10);
	//SquareGrid<1024>ehh(1024, 1024);
	//ehh.solveCholesky();
	//ehh.solveConjugateGradient(1e-3);
	//ehh.solveConjugateGradientSparse(1e-10);
}
