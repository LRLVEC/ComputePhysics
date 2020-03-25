#include <_BLAS.h>
#include <_Time.h>
#include <random>

using namespace BLAS;

template<unsigned long long _dim>struct SquareGrid
{
	static constexpr unsigned long long dim = _dim + 1;
	static constexpr unsigned long long matDim = dim * dim - 1;

	mat matLBand;
	//mat matBand;
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
		matLBand(dim, matDim, MatType::LBandMat, true),
		//matBand(dim, matDim, MatType::BandMat, true),
		matSparse(MatType::SparseMat, 5 * matDim),
		u(matDim, false),
		i(matDim, true),
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
						matLBand.LBandEleRef(n, n - dim) = -1;
						//matBand.BandEleRef(n, n - dim) = -1;
						matSparse.addSparse(n, n - dim, -1, cnt);
					}
					if (c1)
					{
						matLBand.LBandEleRef(n, n - 1) = -1;
						//matBand.BandEleRef(n, n - 1) = -1;
						matSparse.addSparse(n, n - 1, -1, cnt);
					}
					double ss(sumG(c0, c1));
					matLBand.LBandEleRef(n, n) = ss;
					//matBand.BandEleRef(n, n) = ss;
					matSparse.addSparse(n, n, ss, cnt);
					if (c1 != _dim && n + 1 != clipID)
					{
						//matBand.BandEleRef(n, n + 1) = -1;
						matSparse.addSparse(n, n + 1, -1, cnt);
					}
					if (c0 != _dim && n + dim != clipID)
					{
						unsigned long long ns;
						if (n + dim < clipID)ns = n;
						else ns = n - 1;
						//matBand.BandEleRef(n, ns + dim) = -1;
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
						matLBand.LBandEleRef(nd, ns - dim) = -1;
						//matBand.BandEleRef(nd, ns - dim) = -1;
						matSparse.addSparse(nd, ns - dim, -1, cnt);
					}
					if (c1 && nd != clipID)
					{
						matLBand.LBandEleRef(nd, nd - 1) = -1;
						//matBand.BandEleRef(nd, nd - 1) = -1;
						matSparse.addSparse(nd, nd - 1, -1, cnt);
					}
					double ss(sumG(c0, c1));
					matLBand.LBandEleRef(nd, nd) = ss;
					//matBand.BandEleRef(nd, nd) = ss;
					matSparse.addSparse(nd, nd, ss, cnt);
					if (c1 != _dim)
					{
						//matBand.BandEleRef(nd, nd + 1) = -1;
						matSparse.addSparse(nd, nd + 1, -1, cnt);
					}
					if (c0 != _dim)
					{
						//matBand.BandEleRef(nd, nd + dim) = -1;
						matSparse.addSparse(nd, nd + dim, -1, cnt);
					}
				}
			}
		matSparse.elementNum = cnt;
		i.data[0] = 1;
	}
	double solveCholesky()
	{
		//u = 0;
		timer.begin();
		matLBand.solveCholeskyBand(i, u);
		timer.end();
		rCholesky = u.data[0];
		::printf("%.15e\tSquareGrid<%llu> Cholesky\t\t", rCholesky, _dim);
		timer.print();
		return rCholesky;
	}
	/*double solveSteepestDescent(double _esp)
	{
		u = 0;
		timer.begin();
		matBand.solveSteepestDescent(i, u, _esp);
		timer.end();
		rSteepestDescent = u.data[0];
		::printf("%.15e\tSquareGrid<%llu> SteepestDescent\t", rSteepestDescent, _dim);
		timer.print();
		return rSteepestDescent;
	}*/
	/*double solveConjugateGradient(double _esp)
	{
		u = 0;
		timer.begin();
		matBand.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradient = u.data[0];
		::printf("%.15e\tSquareGrid<%llu> ConjugateGradient\t", rConjugateGradient, _dim);
		timer.print();
		return rConjugateGradient;
	}*/
	double solveConjugateGradientSparse(double _esp)
	{
		u = 0;
		timer.begin();
		matSparse.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradientSparse = u.data[0];
		::printf("%.15e\tSquareGrid<%llu> ConjugateGradientSparse\t", rConjugateGradientSparse, _dim);
		timer.print();
		return rConjugateGradientSparse;
	}
};

template<unsigned long long _dim>struct TriangleGrid
{
	static constexpr unsigned long long dim = _dim + 1;
	static constexpr unsigned long long matDim = (dim * (dim + 1)) / 2 - 1;

	mat matLBand;
	//mat matBand;
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
	TriangleGrid(unsigned long long _clipA, unsigned long long _clipB)
		:
		matLBand(dim + 1, matDim, MatType::LBandMat, true),
		//matBand(dim + 1, matDim, MatType::BandMat, true),
		matSparse(MatType::SparseMat, 7 * matDim),
		u(matDim, false),
		i(matDim, true),
		clipA(_clipA),
		clipB(_clipB),
		clipID(id(_clipA, _clipB)),
		rCholesky(0),
		rSteepestDescent(0),
		rConjugateGradient(0),
		rConjugateGradientSparse(0)
	{
		setGrid();
	}
	unsigned long long id(unsigned long long a, unsigned long long b)const
	{
		return ((a + 1) * a) / 2 + b;
	}
	double sumG(unsigned long long a, unsigned long long b)const
	{
		double s(6);
		if (a == _dim)s -= 2;
		if (b == 0)s -= 2;
		if (b == a)s -= 2;
		return s;
	}
	void setGrid()
	{
		unsigned long long cnt(0);
		for (unsigned long long c0(0); c0 < dim; ++c0)
			for (unsigned long long c1(0); c1 <= c0; ++c1)
			{
				unsigned long long n(id(c0, c1));
				if (n < clipID)
				{
					if (c0)
					{
						if (c1)
						{
							matLBand.LBandEleRef(n, n - c0 - 1) = -1;
							//matBand.BandEleRef(n, n - c0 - 1) = -1;
							matSparse.addSparse(n, n - c0 - 1, -1, cnt);
						}
						if (c0 != c1)
						{
							matLBand.LBandEleRef(n, n - c0) = -1;
							//matBand.BandEleRef(n, n - c0) = -1;
							matSparse.addSparse(n, n - c0, -1, cnt);
						}
					}
					if (c1)
					{
						matLBand.LBandEleRef(n, n - 1) = -1;
						//matBand.BandEleRef(n, n - 1) = -1;
						matSparse.addSparse(n, n - 1, -1, cnt);
					}
					double ss(sumG(c0, c1));
					matLBand.LBandEleRef(n, n) = ss;
					//matBand.BandEleRef(n, n) = ss;
					matSparse.addSparse(n, n, ss, cnt);
					if (c1 != c0 && n + 1 != clipID)
					{
						//matBand.BandEleRef(n, n + 1) = -1;
						matSparse.addSparse(n, n + 1, -1, cnt);
					}
					if (c0 != _dim)
						for (unsigned long long ahh(1); ahh <= 2; ++ahh)
							if (n + c0 + ahh != clipID)
							{
								unsigned long long ns(n + c0 + ahh);
								if (ns > clipID)ns -= 1;
								//matBand.BandEleRef(n, ns) = -1;
								matSparse.addSparse(n, ns, -1, cnt);
							}
				}
				else if (n > clipID)
				{
					unsigned long long nd(n - 1);
					unsigned long long ns(n - c0 - 1);
					if (c1 && ns != clipID)
					{
						if (ns > clipID)ns -= 1;
						matLBand.LBandEleRef(nd, ns) = -1;
						//matBand.BandEleRef(nd, ns) = -1;
						matSparse.addSparse(nd, ns, -1, cnt);
					}
					ns = n - c0;
					if (c0 != c1 && ns != clipID)
					{
						if (ns > clipID)ns -= 1;
						matLBand.LBandEleRef(nd, ns) = -1;
						//matBand.BandEleRef(nd, ns) = -1;
						matSparse.addSparse(nd, ns, -1, cnt);
					}
					if (c1 && nd != clipID)
					{
						matLBand.LBandEleRef(nd, nd - 1) = -1;
						//matBand.BandEleRef(nd, nd - 1) = -1;
						matSparse.addSparse(nd, nd - 1, -1, cnt);
					}
					double ss(sumG(c0, c1));
					matLBand.LBandEleRef(nd, nd) = ss;
					//matBand.BandEleRef(nd, nd) = ss;
					matSparse.addSparse(nd, nd, ss, cnt);
					if (c1 != c0)
					{
						//matBand.BandEleRef(nd, nd + 1) = -1;
						matSparse.addSparse(nd, nd + 1, -1, cnt);
					}
					ns = nd + c0 + 1;
					if (c0 != _dim)
						for (unsigned long long ahh(ns); ahh <= ns + 1; ++ahh)
							if (ahh != clipID)
							{
								//matBand.BandEleRef(nd, ahh) = -1;
								matSparse.addSparse(nd, ahh, -1, cnt);
							}
				}
			}
		matSparse.elementNum = cnt;
		i.data[0] = 1;
	}
	double solveCholesky()
	{
		//u = 0;
		timer.begin();
		matLBand.solveCholeskyBand(i, u);
		timer.end();
		rCholesky = u.data[0];
		::printf("%.15e\tTriangleGrid<%llu> Cholesky\t\t", rCholesky, _dim);
		timer.print();
		return rCholesky;
	}
	/*double solveSteepestDescent(double _esp)
	{
		u = 0;
		timer.begin();
		matBand.solveSteepestDescent(i, u, _esp);
		timer.end();
		rSteepestDescent = u.data[0];
		::printf("%.15e\tTriangleGrid<%llu> SteepestDescent\t\t", rSteepestDescent, _dim);
		timer.print();
		return rSteepestDescent;
	}*/
	/*double solveConjugateGradient(double _esp)
	{
		u = 0;
		timer.begin();
		matBand.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradient = u.data[0];
		::printf("%.15e\tTriangleGrid<%llu> ConjugateGradient\t", rConjugateGradient, _dim);
		timer.print();
		return rConjugateGradient;
	}*/
	double solveConjugateGradientSparse(double _esp)
	{
		u = 0;
		timer.begin();
		matSparse.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradientSparse = u.data[0];
		::printf("%.15e\tTriangleGrid<%llu> ConjugateGradientSparse\t", rConjugateGradientSparse, _dim);
		timer.print();
		return rConjugateGradientSparse;
	}
};

template<unsigned long long _dim>struct HexagonGrid
{
	static constexpr unsigned long long dim = _dim + 1;
	static constexpr unsigned long long matDim = (dim * (dim + 1)) / 2 - 1;

	mat matLBand;
	//mat matBand;
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

	HexagonGrid(unsigned long long _clipA, unsigned long long _clipB)
		:
		matLBand(dim + 1, matDim, MatType::LBandMat, true),
		//matBand(dim + 1, matDim, MatType::BandMat, true),
		matSparse(MatType::SparseMat, 7 * matDim),
		u(matDim, false),
		i(matDim, true),
		clipA(_clipA),
		clipB(_clipB),
		clipID(id(_clipA, _clipB)),
		rCholesky(0),
		rSteepestDescent(0),
		rConjugateGradient(0),
		rConjugateGradientSparse(0)
	{
		setGrid();
	}
	unsigned long long id(unsigned long long a, unsigned long long b)const
	{
		return ((a + 1) * a) / 2 + b;
	}
	double sumG(unsigned long long a, unsigned long long b)const
	{
		unsigned long long flag(0);
		if (a == _dim)flag++;
		if (b == 0)flag++;
		if (b == a)flag++;
		switch (flag)
		{
		case 0: return 2;
		case 1: return 5.0 / 3;
		case 2: return 1;
		}
	}
	void setGrid()
	{
		unsigned long long cnt(0);
		for (unsigned long long c0(0); c0 < dim; ++c0)
			for (unsigned long long c1(0); c1 <= c0; ++c1)
			{
				unsigned long long n(id(c0, c1));
				if (n < clipID)
				{
					if (c0)
					{
						if (c1)
						{
							double s(c0 != c1 ? -1.0 / 3 : -1.0 / 2);
							matLBand.LBandEleRef(n, n - c0 - 1) = s;
							//matBand.BandEleRef(n, n - c0 - 1) = s;
							matSparse.addSparse(n, n - c0 - 1, s, cnt);
						}
						if (c0 != c1)
						{
							double s(c1 ? -1.0 / 3 : -1.0 / 2);
							matLBand.LBandEleRef(n, n - c0) = s;
							//matBand.BandEleRef(n, n - c0) = s;
							matSparse.addSparse(n, n - c0, s, cnt);
						}
					}
					if (c1)
					{
						double s(c0 != _dim ? -1.0 / 3 : -1.0 / 2);
						matLBand.LBandEleRef(n, n - 1) = s;
						//matBand.BandEleRef(n, n - 1) = s;
						matSparse.addSparse(n, n - 1, s, cnt);
					}
					double ss(sumG(c0, c1));
					matLBand.LBandEleRef(n, n) = ss;
					//matBand.BandEleRef(n, n) = ss;
					matSparse.addSparse(n, n, ss, cnt);
					if (c1 != c0 && n + 1 != clipID)
					{
						double s(c0 != _dim ? -1.0 / 3 : -1.0 / 2);
						//matBand.BandEleRef(n, n + 1) = s;
						matSparse.addSparse(n, n + 1, s, cnt);
					}
					if (c0 != _dim)
						for (unsigned long long ahh(1); ahh <= 2; ++ahh)
							if (n + c0 + ahh != clipID)
							{
								double s((ahh == 1 ? c1 : (c0 != c1)) ? -1.0 / 3 : -1.0 / 2);
								unsigned long long ns(n + c0 + ahh);
								if (ns > clipID)ns -= 1;
								//matBand.BandEleRef(n, ns) = s;
								matSparse.addSparse(n, ns, s, cnt);
							}
				}
				else if (n > clipID)
				{
					unsigned long long nd(n - 1);
					unsigned long long ns(n - c0 - 1);
					if (c1 && ns != clipID)
					{
						double s(c0 != c1 ? -1.0 / 3 : -1.0 / 2);
						if (ns > clipID)ns -= 1;
						matLBand.LBandEleRef(nd, ns) = s;
						//matBand.BandEleRef(nd, ns) = s;
						matSparse.addSparse(nd, ns, s, cnt);
					}
					ns = n - c0;
					if (c0 != c1 && ns != clipID)
					{
						double s(c1 ? -1.0 / 3 : -1.0 / 2);
						if (ns > clipID)ns -= 1;
						matLBand.LBandEleRef(nd, ns) = s;
						//matBand.BandEleRef(nd, ns) = s;
						matSparse.addSparse(nd, ns, s, cnt);
					}
					if (c1 && nd != clipID)
					{
						double s(c0 != _dim ? -1.0 / 3 : -1.0 / 2);
						matLBand.LBandEleRef(nd, nd - 1) = s;
						//matBand.BandEleRef(nd, nd - 1) = s;
						matSparse.addSparse(nd, nd - 1, s, cnt);
					}
					double ss(sumG(c0, c1));
					matLBand.LBandEleRef(nd, nd) = ss;
					//matBand.BandEleRef(nd, nd) = ss;
					matSparse.addSparse(nd, nd, ss, cnt);
					if (c1 != c0)
					{
						double s(c0 != _dim ? -1.0 / 3 : -1.0 / 2);
						//matBand.BandEleRef(nd, nd + 1) = s;
						matSparse.addSparse(nd, nd + 1, s, cnt);
					}
					ns = nd + c0 + 1;
					if (c0 != _dim)
						for (unsigned long long ahh(ns); ahh <= ns + 1; ++ahh)
							if (ahh != clipID)
							{
								double s((ahh == ns ? c1 : (c0 != c1)) ? -1.0 / 3 : -1.0 / 2);
								//matBand.BandEleRef(nd, ahh) = s;
								matSparse.addSparse(nd, ahh, s, cnt);
							}
				}
			}
		matSparse.elementNum = cnt;
		i.data[0] = 1;
	}
	double solveCholesky()
	{
		//u = 0;
		timer.begin();
		matLBand.solveCholeskyBand(i, u);
		timer.end();
		rCholesky = u.data[0];
		::printf("%.15e\tHexagonGrid<%llu> Cholesky\t\t", rCholesky, _dim);
		timer.print();
		return rCholesky;
	}
	/*double solveSteepestDescent(double _esp)
	{
		u = 0;
		timer.begin();
		matBand.solveSteepestDescent(i, u, _esp);
		timer.end();
		rSteepestDescent = u.data[0];
		::printf("%.15e\tHexagonGrid<%llu> SteepestDescent\t\t", rSteepestDescent, _dim);
		timer.print();
		return rSteepestDescent;
	}*/
	/*double solveConjugateGradient(double _esp)
	{
		u = 0;
		timer.begin();
		matBand.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradient = u.data[0];
		::printf("%.15e\tHexagonGrid<%llu> ConjugateGradient\t", rConjugateGradient, _dim);
		timer.print();
		return rConjugateGradient;
	}*/
	double solveConjugateGradientSparse(double _esp)
	{
		u = 0;
		timer.begin();
		matSparse.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradientSparse = u.data[0];
		::printf("%.15e\tHexagonGrid<%llu> ConjugateGradientSparse\t", rConjugateGradientSparse, _dim);
		timer.print();
		return rConjugateGradientSparse;
	}
};

template<unsigned long long _dim>struct TriangleGridCplx
{
	static constexpr unsigned long long dim = _dim + 1;
	static constexpr unsigned long long matDim = dim * (dim + 1) - 2;

	mat matLBand;
	mat matSparse;
	vec u;
	vec i;
	unsigned long long clipA;
	unsigned long long clipB;
	unsigned long long clipID;
	double omega;
	cplx rCholesky;
	cplx rConjugateGradientSparse;
	Timer timer;

	TriangleGridCplx(unsigned long long _clipA, unsigned long long _clipB)
		:
		matLBand((dim + 1) * 2, matDim, MatType::LBandMat, false),
		matSparse(MatType::SparseMat, 14 * matDim),
		u(matDim, false),
		i(matDim, true),
		clipA(_clipA),
		clipB(_clipB),
		clipID(id(_clipA, _clipB)),
		rCholesky({ 0 }),
		rConjugateGradientSparse({ 0 })
	{
		i.data[0] = 1;
	}
	unsigned long long id(unsigned long long a, unsigned long long b)const
	{
		return ((a + 1) * a) / 2 + b;
	}
	cplx sumG(unsigned long long a, unsigned long long b)
	{
		cplx s{ 2,2 * (omega - 1 / omega) };
		if (a == _dim) { s.re -= 1; s.im -= omega; }
		if (b == 0) { s.im -= omega - 1 / omega; }
		if (b == a) { s.re -= 1; s.im += 1 / omega; }
		return s;
	}
	void setGrid(double _omega)
	{
		matLBand.clear();
		omega = _omega;
		double divOmega(1 / omega);
		unsigned long long cnt(0);
		double tp[14];
		unsigned long long colIndices[14];
		for (unsigned long long c0(0); c0 < dim; ++c0)
		{
			for (unsigned long long c1(0); c1 <= c0; ++c1)
			{
				unsigned long long num(0);
				unsigned long long rowID(0);
				unsigned long long n(id(c0, c1));
				if (n < clipID)
				{
					unsigned long long row(2 * n);
					rowID = 2 * n + 1;
					if (c0)
					{
						if (c1)
						{
							matLBand.LBandEleRef(row, 2 * (n - c0 - 1)) = omega;
							matSparse.addSparse(row, 2 * (n - c0 - 1), omega, cnt);
							tp[num] = matLBand.LBandEleRef(rowID, 2 * (n - c0 - 1) + 1) = -omega;
							colIndices[num++] = 2 * (n - c0 - 1) + 1;
						}
						if (c0 != c1)
						{
							tp[num] = matLBand.LBandEleRef(rowID, 2 * (n - c0)) = matLBand.LBandEleRef(row, 2 * (n - c0) + 1) = -1;
							matSparse.addSparse(row, 2 * (n - c0) + 1, -1, cnt);
							colIndices[num++] = 2 * (n - c0);
						}
					}
					if (c1)
					{
						matLBand.LBandEleRef(row, 2 * (n - 1)) = -divOmega;
						matSparse.addSparse(row, 2 * (n - 1), -divOmega, cnt);
						tp[num] = matLBand.LBandEleRef(rowID, 2 * (n - 1) + 1) = divOmega;
						colIndices[num++] = 2 * (n - 1) + 1;
					}
					cplx ss(sumG(c0, c1));
					matLBand.LBandEleRef(row, 2 * n) = -ss.im;
					//matLBand.LBandEleRef(2 * n, 2 * n + 1) = ss.re;
					matSparse.addSparse(row, 2 * n, -ss.im, cnt);
					matSparse.addSparse(row, 2 * n + 1, ss.re, cnt);
					tp[num] = matLBand.LBandEleRef(rowID, 2 * n) = ss.re;
					colIndices[num++] = 2 * n;
					tp[num] = matLBand.LBandEleRef(rowID, 2 * n + 1) = ss.im;
					colIndices[num++] = 2 * n + 1;
					if (c1 != c0 && n + 1 != clipID)
					{
						matSparse.addSparse(row, 2 * (n + 1), -divOmega, cnt);
						tp[num] = divOmega;
						colIndices[num++] = 2 * (n + 1) + 1;
					}
					if (c0 != _dim)
					{
						if (n + c0 + 1 != clipID)
						{
							unsigned long long ns(n + c0 + 1);
							if (ns > clipID)ns -= 1;
							tp[num] = -1;
							matSparse.addSparse(row, 2 * ns + 1, -1, cnt);
							colIndices[num++] = 2 * ns;
						}
						if (n + c0 + 2 != clipID)
						{
							unsigned long long ns(n + c0 + 2);
							if (ns > clipID)ns -= 1;
							matSparse.addSparse(row, 2 * ns, omega, cnt);
							tp[num] = -omega;
							colIndices[num++] = 2 * ns + 1;
						}
					}
				}
				else if (n > clipID)
				{
					unsigned long long nd(n - 1);
					unsigned long long ns(n - c0 - 1);
					unsigned long long row(2 * nd);
					rowID = 2 * nd + 1;
					if (c1 && ns != clipID)
					{
						if (ns > clipID)ns -= 1;
						matLBand.LBandEleRef(row, 2 * ns) = omega;
						matSparse.addSparse(row, 2 * ns, omega, cnt);
						tp[num] = matLBand.LBandEleRef(rowID, 2 * ns + 1) = -omega;
						colIndices[num++] = 2 * ns + 1;
					}
					ns = n - c0;
					if (c0 != c1 && ns != clipID)
					{
						if (ns > clipID)ns -= 1;
						tp[num] = matLBand.LBandEleRef(rowID, 2 * ns) = matLBand.LBandEleRef(row, 2 * ns + 1) = -1;
						matSparse.addSparse(row, 2 * ns + 1, -1, cnt);
						colIndices[num++] = 2 * ns;
					}
					if (c1 && nd != clipID)
					{
						matLBand.LBandEleRef(row, 2 * (nd - 1)) = -divOmega;
						matSparse.addSparse(row, 2 * (nd - 1), -divOmega, cnt);
						tp[num] = matLBand.LBandEleRef(rowID, 2 * (nd - 1) + 1) = divOmega;
						colIndices[num++] = 2 * (nd - 1) + 1;
					}
					cplx ss(sumG(c0, c1));
					matLBand.LBandEleRef(row, 2 * nd) = -ss.im;
					//matLBand.LBandEleRef(2 * nd, 2 * nd + 1) = ss.re;
					matSparse.addSparse(row, 2 * nd, -ss.im, cnt);
					matSparse.addSparse(row, 2 * nd + 1, ss.re, cnt);
					tp[num] = matLBand.LBandEleRef(rowID, 2 * nd) = ss.re;
					colIndices[num++] = 2 * nd;
					tp[num] = matLBand.LBandEleRef(rowID, 2 * nd + 1) = ss.im;
					colIndices[num++] = 2 * nd + 1;
					if (c1 != c0)
					{
						matSparse.addSparse(row, 2 * (nd + 1), -divOmega, cnt);
						tp[num] = divOmega;
						colIndices[num++] = 2 * (nd + 1) + 1;
					}
					if (c0 != _dim)
					{
						ns = nd + c0 + 1;
						if (ns != clipID)
						{
							matSparse.addSparse(row, 2 * ns + 1, -1, cnt);
							tp[num] = -1;
							colIndices[num++] = 2 * ns;
						}
						++ns;
						if (ns != clipID)
						{
							matSparse.addSparse(row, 2 * ns, omega, cnt);
							tp[num] = -omega;
							colIndices[num++] = 2 * ns + 1;
						}
					}
				}
				if (rowID)
				{
					memcpy64d(matSparse.data + cnt, tp, num);
					memcpy64d(matSparse.colIndice + cnt, colIndices, num);
					for (unsigned long long cp(0); cp < num; ++cp)
						matSparse.rowIndice[cnt + cp] = rowID;
					cnt += num;
				}
			}
		}
		matSparse.elementNum = cnt;
	}
	cplx solveCholesky()
	{
		timer.begin();
		matLBand.solveCholeskyBand(i, u);
		timer.end();
		rCholesky.im = u.data[0];
		rCholesky.re = u.data[1];
		::printf("(%.15e, %.15e), %.15e\tTriangleGridCplx<%llu, omega=%.3e> Cholesky\t\t",
			rCholesky.re, rCholesky.im, rCholesky.norm(), _dim, omega);
		timer.print();
		return rCholesky;
	}
	cplx solveConjugateGradientSparse(double _esp)
	{
		timer.begin();
		matSparse.solveConjugateGradient(i, u, _esp);
		timer.end();
		rConjugateGradientSparse.im = u.data[0];
		rConjugateGradientSparse.re = u.data[1];
		::printf("(%.15e, %.15e), %.15e\tTriangleGridCplx<%llu, omega=%.3e> ConjugateGradientSparse\t",
			rConjugateGradientSparse.re, rConjugateGradientSparse.im, rConjugateGradientSparse.norm(), _dim, omega);
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

	//SquareGrid<4>test(4, 4);
	//mat ts(5 * 5 - 1, 5 * 5 - 1, false);
	//mat tr(5 * 5 - 1, 5 * 5 - 1, false);
	//transLBandToSymmetricMat(test.matLBand, ts);
	//ts(ts, tr);
	//ts.printToTableTxt("./ts.txt");
	//tr.printToTableTxt("./tr.txt");

	timer.begin();

	constexpr double eps(1e-18);

	SquareGrid<1>sq1_1(1, 1);
	sq1_1.solveCholesky();
	sq1_1.solveConjugateGradientSparse(eps);
	SquareGrid<4>sq4_1(4, 4);
	sq4_1.solveCholesky();
	sq4_1.solveConjugateGradientSparse(eps);
	SquareGrid<16>sq16_1(16, 16);
	sq16_1.solveCholesky();
	sq16_1.solveConjugateGradientSparse(eps);
	SquareGrid<64>sq64_1(64, 64);
	sq64_1.solveCholesky();
	sq64_1.solveConjugateGradientSparse(eps);

	::printf("\n");

	SquareGrid<1>sq1_2(0, 1);
	sq1_2.solveCholesky();
	sq1_2.solveConjugateGradientSparse(eps);
	SquareGrid<4>sq4_2(0, 4);
	sq4_2.solveCholesky();
	sq4_2.solveConjugateGradientSparse(eps);
	SquareGrid<16>sq16_2(0, 16);
	sq16_2.solveCholesky();
	sq16_2.solveConjugateGradientSparse(eps);
	SquareGrid<64>sq64_2(0, 64);
	sq64_2.solveCholesky();
	sq64_2.solveConjugateGradientSparse(eps);

	::printf("\n");

	TriangleGrid<1>tr1(1, 0);
	tr1.solveCholesky();
	tr1.solveConjugateGradientSparse(eps);
	TriangleGrid<4>tr4(4, 0);
	tr4.solveCholesky();
	tr4.solveConjugateGradientSparse(eps);
	TriangleGrid<16>tr16(16, 0);
	tr16.solveCholesky();
	tr16.solveConjugateGradientSparse(eps);
	TriangleGrid<64>tr64(64, 0);
	tr64.solveCholesky();
	tr64.solveConjugateGradientSparse(eps);

	::printf("\n");

	HexagonGrid<4>he4(3, 0);
	he4.solveCholesky();
	he4.solveConjugateGradientSparse(1e-18);
	HexagonGrid<16>he16(15, 0);
	he16.solveCholesky();
	he16.solveConjugateGradientSparse(1e-18);
	HexagonGrid<64>he64(63, 0);
	he64.solveCholesky();
	he64.solveConjugateGradientSparse(1e-18);

	timer.end();
	timer.print("Total time:");

	TriangleGridCplx<16>trCplx1(16, 0);
	trCplx1.setGrid(2);
	trCplx1.solveCholesky();
	trCplx1.solveConjugateGradientSparse(1e-18);
	trCplx1.setGrid(5);
	trCplx1.solveCholesky();
	trCplx1.solveConjugateGradientSparse(1e-18);
	trCplx1.setGrid(10);
	trCplx1.solveCholesky();
	trCplx1.solveConjugateGradientSparse(1e-18);
}
