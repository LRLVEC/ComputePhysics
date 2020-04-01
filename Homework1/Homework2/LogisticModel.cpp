#include <_BLAS.h>
#include <random>

using namespace BLAS;
inline double logistic(double r, double x)
{
	return x * r * (1 - x);
}
inline double logisticNew(double r, double x)
{
	return r * sin(x);
}

double sections[13][3] =
{
	{0.0000000000,2.9999859514,2.9999859514},
	{2.9999859515,3.4494842807,0.4494983292},
	{3.4494842808,3.5440881112,0.0946038304},
	{3.5440881113,3.5644063456,0.0203182343},
	{3.5644063457,3.5687590436,0.0043526979},
	{3.5687590437,3.5696914566,0.0009324129},
	{3.5696914567,3.5698911971,0.0001997404},
	{3.5698911972,3.5699339931,0.0000427959},
	{3.5699339932,3.5699431659,0.0000091727},
	{3.5699431660,3.5699451332,0.0000019672},
	{3.5699451333,3.5699455557,0.0000004224},
	{3.5699455558,3.5699456467,0.0000000909},
	{3.5699456468,3.5699456663,0.0000000195}
};
double sectionsNew[13][3] =
{
	{0.0000000000,2.2616322271,2.2616322271},
	{2.2616322272,2.6176971311,0.3560649039},
	{2.6176971312,2.6973637854,0.0796666542},
	{2.6973637855,2.7145855603,0.0172217748},
	{2.7145855604,2.7182850522,0.0036994918},
	{2.7182850523,2.7190794259,0.0007943736},
	{2.7190794260,2.7192502852,0.0001708592},
	{2.7192502853,2.7192871702,0.0000368849},
	{2.7192871703,2.7192951884,0.0000080181},
	{2.7192951885,2.7192969531,0.0000017646},
	{2.7192969532,2.7192973498,0.0000003966},
	{2.7192973499,2.7192974421,0.0000000922}
};

void q1()
{
	FILE* temp(::fopen("q1.txt", "w+"));
	double x[9];
	unsigned long long num(9);
	unsigned long long iters(20);
	for (unsigned long long c0(0); c0 < num; ++c0)
		x[c0] = 0.1 * c0 + 0.1;
	for (unsigned long long c0(0); c0 < iters; ++c0)
	{
		for (unsigned long long c1(0); c1 < num; ++c1)
		{
			::fprintf(temp, "%.10e\t", x[c1]);
			x[c1] = logistic(0.5, x[c1]);
		}
		::fprintf(temp, "\n");
	}
	::fprintf(temp, "\n");
	for (unsigned long long c0(0); c0 < num; ++c0)
		x[c0] = 0.1 * c0 + 0.1;
	for (unsigned long long c0(0); c0 < iters; ++c0)
	{
		for (unsigned long long c1(0); c1 < num; ++c1)
		{
			::fprintf(temp, "%.10e\t", x[c1]);
			x[c1] = logistic(1.5, x[c1]);
		}
		::fprintf(temp, "\n");
	}
	::fclose(temp);
}

void q3()
{
	FILE* temp(::fopen("q3.txt", "w+"));
	double x[9];
	unsigned long long num(9);
	unsigned long long iters(20);
	for (unsigned long long c0(0); c0 < num; ++c0)
		x[c0] = 0.1 * c0 + 0.1;
	for (unsigned long long c0(0); c0 < iters; ++c0)
	{
		for (unsigned long long c1(0); c1 < num; ++c1)
		{
			::fprintf(temp, "%.10e\t", x[c1]);
			x[c1] = logistic(3.1, x[c1]);
		}
		::fprintf(temp, "\n");
	}
	::fclose(temp);
}

void q4()
{
	FILE* temp(::fopen("q4.txt", "w+"));
	double x, r;
	constexpr unsigned long long num(100);
	unsigned long long iters(2000);
	unsigned long long s(0);
	for (unsigned long long cc(0); cc < 2; ++cc)
	{
		for (unsigned long long c0(0); c0 < 10; ++c0)
		{
			x = 0.5;
			r = sections[cc][0] + sections[cc][2] * c0 / 1000;
			for (unsigned long long c1(0); c1 < iters; ++c1)
				x = logistic(r, x);
			::fprintf(temp, "%.10e", r);
			for (unsigned long long c1(0); c1 < s; ++c1)
				::fprintf(temp, "\t");
			for (unsigned long long c1(0); c1 < (1llu << cc); ++c1)
				::fprintf(temp, "\t%.10e", x = logistic(r, x));
			::fprintf(temp, "\n");
		}
		for (unsigned long long c0(1); c0 < num; ++c0)
		{
			x = 0.5;
			r = sections[cc][0] + sections[cc][2] * c0 / 100;
			for (unsigned long long c1(0); c1 < iters; ++c1)
				x = logistic(r, x);
			::fprintf(temp, "%.10e", r);
			for (unsigned long long c1(0); c1 < s; ++c1)
				::fprintf(temp, "\t");
			for (unsigned long long c1(0); c1 < (1llu << cc); ++c1)
				::fprintf(temp, "\t%.10e", x = logistic(r, x));
			::fprintf(temp, "\n");
		}
		s += (1llu << cc);
	}
	::fclose(temp);
}

void q5()
{
	FILE* temp(::fopen("q5_1.txt", "w+"));
	double x, r;
	constexpr unsigned long long num(100);
	unsigned long long preIters(1000);
	unsigned long long iters(65536);
	//4T
	unsigned long long s(0);
	for (unsigned long long cc(0); cc < 5; ++cc)
	{
		for (unsigned long long c0(0); c0 < 10; ++c0)
		{
			x = 0.5;
			r = sections[cc][0] + sections[cc][2] * c0 / 1000;
			for (unsigned long long c1(0); c1 < iters; ++c1)
				x = logistic(r, x);
			::fprintf(temp, "%.10e", r);
			for (unsigned long long c1(0); c1 < s; ++c1)
				::fprintf(temp, "\t");
			for (unsigned long long c1(0); c1 < (1llu << cc); ++c1)
				::fprintf(temp, "\t%.10e", x = logistic(r, x));
			::fprintf(temp, "\n");
		}
		for (unsigned long long c0(1); c0 < num; ++c0)
		{
			x = 0.5;
			r = sections[cc][0] + sections[cc][2] * c0 / 100;
			for (unsigned long long c1(0); c1 < iters; ++c1)
				x = logistic(r, x);
			::fprintf(temp, "%.10e", r);
			for (unsigned long long c1(0); c1 < s; ++c1)
				::fprintf(temp, "\t");
			for (unsigned long long c1(0); c1 < (1llu << cc); ++c1)
				::fprintf(temp, "\t%.10e", x = logistic(r, x));
			::fprintf(temp, "\n");
		}
		s += (1llu << cc);
	}
	::fclose(temp);
	temp = ::fopen("q5_2.txt", "w+");
	for (unsigned long long cc(0); cc < 12; ++cc)
	{
		x = 0.5;
		for (unsigned long long c0(0); c0 < 50; ++c0)
		{
			x = 0.5;
			r = sections[cc][0] + sections[cc][2] * c0 / 50;
			for (unsigned long long c1(0); c1 < iters; ++c1)
				x = logistic(r, x);
			double x1 = x + 1e-8;
			for (unsigned long long c1(0); c1 < (1llu << cc); ++c1)
				x1 = logistic(r, x1);
			if (abs(x1 - x) * 1e8 < 1)
				::fprintf(temp, "%.10e\t%.10e", r, (x1 - x) * 1e8);
			::fprintf(temp, "\n");
		}
	}
	::fclose(temp);
}

void q9()
{
	FILE* temp(::fopen("q9_1.txt", "w+"));
	/*double x, r;
	constexpr unsigned long long num(100);
	unsigned long long iters(65536);
	unsigned long long s(0);
	for (unsigned long long c0(0); c0 < 100; ++c0)
	{
		x = Pi / 2;
		r = 2.7192974419 + 1e-10 * c0;
		for (unsigned long long c1(0); c1 < iters; ++c1)
			x = logisticNew(r, x);
		double x1(logisticNew(r, x));
		for (unsigned long long c1(0); c1 < 2047; ++c1)
			x1 = logisticNew(r, x1);
		::printf("%.10e\t%.12e\n", r, (x1 - x));
	}*/
	double xt[9];
	unsigned long long num(9);
	unsigned long long iters(20);
	for (unsigned long long c0(0); c0 < num; ++c0)
		xt[c0] = (0.1 * c0 + 0.1) * Pi;
	for (unsigned long long c0(0); c0 < iters; ++c0)
	{
		for (unsigned long long c1(0); c1 < num; ++c1)
		{
			::fprintf(temp, "%.10e\t", xt[c1]);
			xt[c1] = logisticNew(0.5, xt[c1]);
		}
		::fprintf(temp, "\n");
	}
	::fprintf(temp, "\n");
	for (unsigned long long c0(0); c0 < num; ++c0)
		xt[c0] = (0.1 * c0 + 0.1) * Pi;
	for (unsigned long long c0(0); c0 < iters; ++c0)
	{
		for (unsigned long long c1(0); c1 < num; ++c1)
		{
			::fprintf(temp, "%.10e\t", xt[c1]);
			xt[c1] = logisticNew(1.5, xt[c1]);
		}
		::fprintf(temp, "\n");
	}
	::fprintf(temp, "\n");
	for (unsigned long long c0(0); c0 < num; ++c0)
		xt[c0] = (0.1 * c0 + 0.1) * Pi;
	for (unsigned long long c0(0); c0 < iters; ++c0)
	{
		for (unsigned long long c1(0); c1 < num; ++c1)
		{
			::fprintf(temp, "%.10e\t", xt[c1]);
			xt[c1] = logisticNew(sectionsNew[1][0] + 0.1, xt[c1]);
		}
		::fprintf(temp, "\n");
	}
	::fclose(temp);
	double x, r;
	iters = 65536, num = 100;
	temp = ::fopen("q9_4.txt", "w+");
	unsigned long long s(0);
	for (unsigned long long cc(0); cc < 5; ++cc)
	{
		for (unsigned long long c0(0); c0 < 10; ++c0)
		{
			x = 0.5*Pi;
			r = sectionsNew[cc][0] + sectionsNew[cc][2] * c0 / 1000;
			for (unsigned long long c1(0); c1 < iters; ++c1)
				x = logisticNew(r, x);
			::fprintf(temp, "%.10e", r);
			for (unsigned long long c1(0); c1 < s; ++c1)
				::fprintf(temp, "\t");
			for (unsigned long long c1(0); c1 < (1llu << cc); ++c1)
				::fprintf(temp, "\t%.10e", x = logisticNew(r, x));
			::fprintf(temp, "\n");
		}
		for (unsigned long long c0(1); c0 < num; ++c0)
		{
			x = 0.5*Pi;
			r = sectionsNew[cc][0] + sectionsNew[cc][2] * c0 / 100;
			for (unsigned long long c1(0); c1 < iters; ++c1)
				x = logisticNew(r, x);
			::fprintf(temp, "%.10e", r);
			for (unsigned long long c1(0); c1 < s; ++c1)
				::fprintf(temp, "\t");
			for (unsigned long long c1(0); c1 < (1llu << cc); ++c1)
				::fprintf(temp, "\t%.10e", x = logisticNew(r, x));
			::fprintf(temp, "\n");
		}
		s += (1llu << cc);
	}
	temp = ::fopen("q9_5.txt", "w+");
	for (unsigned long long cc(0); cc < 12; ++cc)
	{
		x = 0.5 * Pi;
		for (unsigned long long c0(0); c0 < 50; ++c0)
		{
			x = 0.5;
			r = sectionsNew[cc][0] + sectionsNew[cc][2] * c0 / 50;
			for (unsigned long long c1(0); c1 < iters; ++c1)
				x = logisticNew(r, x);
			double x1 = x + 1e-8;
			for (unsigned long long c1(0); c1 < (1llu << cc); ++c1)
				x1 = logisticNew(r, x1);
			if (abs(x1 - x) * 1e8 < 1)
				::fprintf(temp, "%.10e\t%.10e", r, (x1 - x) * 1e8);
			::fprintf(temp, "\n");
		}
	}
	::fclose(temp);
}

int main()
{
	//q1();
	//q3();
	//q4();
	//q5();
	q9();
}
