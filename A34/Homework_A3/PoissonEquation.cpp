#include <_BLAS.h>
#include <_Time.h>

using namespace BLAS;

template<unsigned long long _dim>struct Grid
{
	static constexpr unsigned long long dim = _dim;
	static constexpr unsigned long long blockDim = dim - 1;
	static constexpr double h = 1.0 / dim;
	mat m;
	vec f;
	vec h2f;
	vec u;
	mat uij;
	double L2Error;
	double(*func)(double, double);
	double(*answer)(double, double);

	Grid(double(*_func)(double, double), double(*_answer)(double, double))
		:
		m(MatType::SparseMat, dim* (5 * dim - 14) + 9),
		f(blockDim* blockDim, false),
		h2f(blockDim* blockDim, false),
		u(blockDim* blockDim, false),
		uij(dim + 1, dim + 1, false),
		func(_func),
		answer(_answer)
	{
		setGrid();
		solveGrid();
		setAnswerGrid();
		L2Error = totalL2Error();
	}
	void setGrid()
	{
		unsigned long long n(0);
		unsigned long long cnt(0);
		for (unsigned long long c0(0); c0 < blockDim; ++c0)
			for (unsigned long long c1(0); c1 < blockDim; ++c1)
			{
				f[n] = func(double(1 + c0) / dim, double(1 + c1) / dim);
				if (c0)m.addSparse(n, n - blockDim, -1, cnt);
				if (c1)m.addSparse(n, n - 1, -1, cnt);
				m.addSparse(n, n, 4, cnt);
				if (c1 != dim - 2)m.addSparse(n, n + 1, -1, cnt);
				if (c0 != dim - 2)m.addSparse(n, n + blockDim, -1, cnt);
				n++;
			}
		m.elementNum = cnt;
		h2f = f;
		h2f *= (h * h);
	}
	void solveGrid()
	{
		m.solveConjugateGradient(h2f, u, 1e-15);
	}
	void setAnswerGrid()
	{
		memset64d(uij.data, 0, uij.width4d);
		for (unsigned long long c0(1); c0 < dim; ++c0)
		{
			uij(c0, 0) = 0;
			memcpy64d(uij.data + uij.width4d * c0 + 1,
				u.data + (c0 - 1) * blockDim, blockDim);
			uij(c0, dim) = 0;
		}
		memset64d(uij.data + uij.width4d * dim, 0, uij.width4d);
	}
	double check(unsigned long long c0, unsigned long long c1)
	{
		return abs(u[c0 * blockDim + c1] - answer(double(1 + c0) / dim, double(1 + c1) / dim));
	}
	double maxDelta()
	{
		double er(0);
		for (unsigned long long c0(0); c0 < blockDim; c0++)
			for (unsigned long long c1(0); c1 < blockDim; c1++)
			{
				double d(check(c0, c1));
				if (er < d)er = d;
			}
		return er;
	}
	double gaussPointValue(unsigned long long c0, unsigned long long c1, long long idx)
	{
		double rsq3(1 / sqrt(3.0));
		double x((2 * (idx & 1) - 1) * rsq3);
		double y((idx & 2 - 1) * rsq3);
		double r(uij(c0, c1) * (1 - x) * (1 - y));
		r += uij(c0 + 1, c1) * (1 + x) * (1 - y);
		r += uij(c0, c1 + 1) * (1 - x) * (1 + y);
		r += uij(c0 + 1, c1 + 1) * (1 + x) * (1 + y);
		return r / 4;
	}
	double blockL2Error(unsigned long long c0, unsigned long long c1)
	{
		double s(0);
		for (long long idx(0); idx < 4; ++idx)
		{
			double rsq3(1 / sqrt(3.0));
			double x((2 * (idx & 1) - 1) * rsq3);
			double y((idx & 2 - 1) * rsq3);
			x = (c0 + (1 + x) / 2) / dim;
			y = (c1 + (1 + y) / 2) / dim;
			double dt(gaussPointValue(c0, c1, idx) - answer(x, y));
			s += dt * dt;
		}
		return s / (dim * dim);
	}
	double totalL2Error()
	{
		double er(0);
		for (unsigned long long c0(0); c0 < dim; ++c0)
			for (unsigned long long c1(0); c1 < dim; ++c1)
				er += blockL2Error(c0, c1);
		return sqrt(er);
	}
};

double f(double x, double y)
{
	static constexpr double Pi2 = 2.0 * Pi * Pi;
	return  Pi2 * sin(Pi * x) * sin(Pi * y);
}
double answer(double x, double y)
{
	return  sin(Pi * x) * sin(Pi * y);
}

int main()
{
	Timer timer;

	::printf("Grid<16>:\t");
	timer.begin();
	Grid<16>grid16(f, answer);
	timer.end();
	timer.print("Time used: ");
	::printf("MaxError:\t%.2e\tL2Error:\t%.2e\n", grid16.maxDelta(), grid16.L2Error);

	::printf("Grid<32>:\t");
	timer.begin();
	Grid<32>grid32(f, answer);
	timer.end();
	timer.print("Time used: ");
	::printf("MaxError:\t%.2e\tL2Error:\t%.2e\n", grid32.maxDelta(), grid32.L2Error);

	::printf("Grid<64>:\t");
	timer.begin();
	Grid<64>grid64(f, answer);
	timer.end();
	timer.print("Time used: ");
	::printf("MaxError:\t%.2e\tL2Error:\t%.2e\n", grid64.maxDelta(), grid64.L2Error);

	::printf("Grid<128>:\t");
	timer.begin();
	Grid<128>grid128(f, answer);
	timer.end();
	timer.print("Time used: ");
	::printf("MaxError:\t%.2e\tL2Error:\t%.2e\n", grid128.maxDelta(), grid128.L2Error);

	::printf("Grid<256>:\t");
	timer.begin();
	Grid<256>grid256(f, answer);
	timer.end();
	timer.print("Time used: ");
	::printf("MaxError:\t%.2e\tL2Error:\t%.2e\n", grid256.maxDelta(), grid256.L2Error);

	::printf("Grid<512>:\t");
	timer.begin();
	Grid<512>grid512(f, answer);
	timer.end();
	timer.print("Time used: ");
	::printf("MaxError:\t%.2e\tL2Error:\t%.2e\n", grid512.maxDelta(), grid512.L2Error);

	::printf("Grid<1024>:\t");
	timer.begin();
	Grid<1024>grid1024(f, answer);
	timer.end();
	timer.print("Time used: ");
	::printf("MaxError:\t%.2e\tL2Error:\t%.2e\n", grid1024.maxDelta(), grid1024.L2Error);
}