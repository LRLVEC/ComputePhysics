#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <immintrin.h>
//#include <ammintrin.h>
#include <_Time.h>

void testALU(double* a, double* b, double* c, unsigned int l)
{
	for (unsigned int c0(0); c0 < 10; ++c0)
	{
		for (unsigned int c1(0); c1 < l; ++c1)
			c[c1] = a[c1] * b[c1];
		for (unsigned int c1(0); c1 < l; ++c1)
			a[c1] = b[c1] * c[c1];
		for (unsigned int c1(0); c1 < l; ++c1)
			b[c1] = a[c1] * c[c1];
	}
}
void testAVX2(__m256d* a, __m256d* b, __m256d* c, unsigned int l)
{
	for (unsigned int c0(0); c0 < 10; ++c0)
	{
		for (unsigned int c1(0); c1 < l; ++c1)
			c[c1] = _mm256_mul_pd(a[c1], b[c1]);
		for (unsigned int c1(0); c1 < l; ++c1)
			a[c1] = _mm256_mul_pd(b[c1], c[c1]);
		for (unsigned int c1(0); c1 < l; ++c1)
			b[c1] = _mm256_mul_pd(a[c1], c[c1]);
	}
}

void clear(double* a, double* b, unsigned int l)
{
	for (unsigned int c0(0); c0 < l; ++c0)
	{
		a[c0] = double(rand()) / rand();
		b[c0] = double(rand()) / rand();
	}
}

int main()
{
	srand(time(nullptr));
	unsigned int l(20000000);
	double* a((double*)::malloc(l * sizeof(double)));
	double* b((double*)::malloc(l * sizeof(double)));
	double* c((double*)::malloc(l * sizeof(double)));
	unsigned int lv(l / 4);
	__m256d* av((__m256d*)a);
	__m256d* bv((__m256d*)b);
	__m256d* cv((__m256d*)c);

	Timer timer;

	for (unsigned int c0(0); c0 < 10; ++c0)
	{
		clear(a, b, l);
		timer.begin();
		testALU(a, b, c, l);
		timer.end();
		timer.print("ALU:");

	}
	for (unsigned int c0(0); c0 < 10; ++c0)
	{
		clear(a, b, l);
		timer.begin();
		testAVX2(av, bv, cv, lv);
		timer.end();
		timer.print("AVX2:");
	}
	::free(a);
	::free(b);
	::free(c);
}