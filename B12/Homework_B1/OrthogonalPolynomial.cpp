#include <_BLAS.h>
#include <_Time.h>

using namespace BLAS;

struct Legendre
{
	double value(unsigned long long n, double x)const
	{
		if (n == 0)return sqrt(1.0 / 2);
		if (n == 1)return x * sqrt(2.0 / 3.0);
		double x0(1), x1(x);
		double x2;
		for (unsigned long long c0(2); c0 <= n; ++c0)
		{
			x2 = ((2 * c0 - 1) * x * x1 - (c0 - 1) * x0) / c0;
			x0 = x1;
			x1 = x2;
		}
		return x2 * sqrt((2 * n + 1) / 2.0);
	}
};

struct Laguerre
{
	double value(unsigned long long n, double x)const
	{
		if (n == 0)return 1;
		if (n == 1)return x;
		double x0(1), x1(x);
		double x2;
		for (unsigned long long c0(2); c0 <= n; ++c0)
		{
			x2 = ((2 * c0 - 1 - x) * x1 - (c0 - 1) * x0) / c0;
			x0 = x1;
			x1 = x2;
		}
		return x2;
	}
};

struct Hermite
{
	unsigned long long n;
	mat m;
	vec zeroPoints;
	vec eigenVector;
	vec Aj;

	Hermite(unsigned long long _n)
		:
		n(_n),
		m(1, _n, MatType::BandMat, true),
		zeroPoints(_n, false),
		eigenVector(_n, false),
		Aj(_n, false)
	{
		for (unsigned long long c0(0); c0 < n - 1; ++c0)
			m.BandEleRef(c0, c0 + 1) = m.BandEleRef(c0 + 1, c0) = sqrt((c0 + 1) / 2.0);
	}
	double value(unsigned long long _n, double x)const
	{
		if (_n == 0)return pow(Pi, -0.25);
		if (_n == 1)return x * sqrt(2) * pow(Pi, -0.25);
		double x0(pow(Pi, -0.25)), x1(x * sqrt(2) * pow(Pi, -0.25));
		double x2;
		for (unsigned long long c0(2); c0 <= _n; ++c0)
		{
			x2 = sqrt(2.0 / c0) * x * x1 - sqrt(double(c0 - 1) / c0) * x0;
			x0 = x1;
			x1 = x2;
		}
		return x2;
	}
	vec& getZeros()
	{
		return m.implicitSymmetricQR(1e-40, zeroPoints);
	}
	vec& getAj()
	{
		for (unsigned long long c0(0); c0 < n; ++c0)
		{
			for (unsigned long long c1(0); c1 < n; ++c1)
				eigenVector[c1] = value(c1, zeroPoints[c0]);
			Aj[c0] = 1 / eigenVector.norm2Square();
		}
		return Aj;
	}
};

struct Chebyshev
{
	unsigned long long n;
	mat m;
	vec zeroPoints;
	vec eigenVector;
	vec Aj;

	Chebyshev(unsigned long long _n)
		:
		n(_n),
		m(1, _n, MatType::BandMat, true),
		zeroPoints(_n, false),
		eigenVector(_n, false),
		Aj(_n, false)
	{
		m.BandEleRef(0, 1) = m.BandEleRef(1, 0) = sqrt(1.0 / 2);
		for (unsigned long long c0(1); c0 < n - 1; ++c0)
			m.BandEleRef(c0, c0 + 1) = m.BandEleRef(c0 + 1, c0) = 1.0 / 2.0;
	}
	double value(unsigned long long _n, double x)const
	{
		if (_n == 0)return pow(Pi, -0.5);
		if (_n == 1)return x * sqrt(2 / Pi);
		double x0(1), x1(x);
		double x2;
		for (unsigned long long c0(2); c0 <= _n; ++c0)
		{
			x2 = 2 * x * x1 - x0;
			x0 = x1;
			x1 = x2;
		}
		return x2 * sqrt(2 / Pi);
	}
	vec& getZeros()
	{
		return m.implicitSymmetricQR(1e-40, zeroPoints);
	}
	vec& getAj()
	{
		for (unsigned long long c0(0); c0 < n; ++c0)
		{
			for (unsigned long long c1(0); c1 < n; ++c1)
				eigenVector[c1] = value(c1, zeroPoints[c0]);
			Aj[c0] = 1 / eigenVector.norm2Square();
		}
		return Aj;
	}
	bool checkZeros()
	{
		zeroPoints.qsort();
		double dt(Pi / (2 * zeroPoints.dim));
		for (unsigned long long c0(0); c0 < zeroPoints.dim; ++c0)
			if (abs(zeroPoints[c0] - cos(dt * (2 * (zeroPoints.dim - c0) - 1))) > 1e-14)
				return false;
		return true;
	}
};

void q2_1()
{
	Legendre legendre;
	Laguerre laguerre;
	for (unsigned long long c0(2); c0 <= 1024; c0 *= 8)
		::printf("%.10e\n", legendre.value(c0, 0.5));
	for (unsigned long long c0(2); c0 <= 1024; c0 *= 8)
		::printf("%.10e\n", laguerre.value(c0, 0.5));
}

void q2_2()
{
	Timer timer;
	//q2_1();
	for (unsigned long long c0(2); c0 <= 1024; c0 *= 8)
	{
		Hermite hermite(c0);
		timer.begin();
		hermite.getZeros();
		timer.end();
		printf("Hermite %llu: Time:", c0);
		timer.print();

		Chebyshev chebyshev(c0);
		timer.begin();
		chebyshev.getZeros();
		timer.end();
		printf("Chebyshev %llu: ", c0);
		::printf("Check: %d Time:", chebyshev.checkZeros());
		timer.print();
	}
}

void q2_3_1()
{
	Timer timer;
	for (unsigned long long c0(2); c0 <= 128; c0 *= 8)
	{
		Hermite hermite(c0);
		hermite.getZeros().qsort();
		timer.begin();
		hermite.getAj();
		timer.end();
		printf("Hermite %llu: Time:", c0);
		timer.print();
	}
	for (unsigned long long c0(2); c0 <= 1024; c0 *= 8)
	{
		Chebyshev chebyshev(c0);
		chebyshev.getZeros().qsort();
		timer.begin();
		chebyshev.getAj();
		timer.end();
		printf("Chebyshev %llu: Time:", c0);
		timer.print();
	}
}

void q2_3_2_1(unsigned long long n, unsigned long long j)
{
	Hermite hermite(n);
	hermite.getZeros().qsort();
	hermite.getAj();
	unsigned long long div(1000);
	vec f(div + 1, false);
	double a(1.5 * hermite.zeroPoints[0]), b(1.5 * hermite.zeroPoints[n - 1]);
	double dx((b - a) / div);
	double xj(hermite.zeroPoints[j]);
	vec q1(n, false), q2(n, false);
	for (unsigned long long c1(0); c1 < n; ++c1)
		q1[c1] = hermite.value(c1, xj);
	for (unsigned long long c1(0); c1 <= div; ++c1)
	{
		double x(a + dx * c1);
		for (unsigned long long c2(0); c2 < n; ++c2)
			q2[c2] = hermite.value(c2, x);
		f[c1] = exp((xj * xj - x * x) / 2) * (q1, q2);
	}
	::printf("%.12e %.12e\n", a, b);
	f.print();
}

void q2_3_2_2(unsigned long long n, unsigned long long j)
{
	Chebyshev chebyshev(n);
	chebyshev.getZeros().qsort();
	chebyshev.getAj();
	unsigned long long div(1000);
	vec f(div + 1, false);
	double a(chebyshev.zeroPoints[0]), b(chebyshev.zeroPoints[n - 1]);
	double dx((b - a) / div);
	double xj(chebyshev.zeroPoints[j]);
	vec q1(n, false), q2(n, false);
	for (unsigned long long c1(0); c1 < n; ++c1)
		q1[c1] = chebyshev.value(c1, xj);
	for (unsigned long long c1(0); c1 <= div; ++c1)
	{
		double x(a + dx * c1);
		for (unsigned long long c2(0); c2 < n; ++c2)
			q2[c2] = chebyshev.value(c2, x);
		f[c1] = sqrt(1 - xj * xj) * (q1, q2) / sqrt(1 - x * x);
	}
	::printf("%.12e %.12e\n", a, b);
	f.print();
}

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);
	Timer timer;
	//q2_1();
	//q2_2();
	//q2_3_1();
	//q2_3_2_1(128, 64);
	//q2_3_2_2(128, 64);
}