#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <immintrin.h>
#include <random>
#include <chrono>

void clear(double* a, unsigned int length)
{
	for (unsigned int c0(0); c0 < length; ++c0)
		a[c0] = 0;
}
void randomMat(double* a, unsigned int length, std::mt19937& mt, std::uniform_real_distribution<double>& rd)
{
	for (unsigned int c0(0); c0 < length; ++c0)
		a[c0] = rd(mt);
}

void matMult(double* a, double* b, double* c, unsigned int widthA,
	unsigned int heightA, unsigned int widthB, unsigned int heightB)
{
	__m256d* aData((__m256d*)a);
	__m256d* cData((__m256d*)c);
	unsigned int widthA4(widthA / 4);
	unsigned int minDim(widthA > heightB ? heightB : widthA);
	constexpr unsigned int warp = 16;
	for (unsigned int c0(0); c0 < heightA; c0 += 2)
		for (unsigned int c1(0); c1 < widthA4; c1 += warp)
		{
			__m256d ans0[warp] = { 0 };
			__m256d ans1[warp] = { 0 };
			for (unsigned int c2(0); c2 < minDim; ++c2)
			{
				double s = a[c0 * widthA + c2];
				__m256d tp0 = { s,s,s,s };
				s = a[(c0 + 1) * widthA + c2];
				__m256d tp1 = { s,s,s,s };
#pragma unroll(4)
				for (unsigned int c3(0); c3 < warp; ++c3)
				{
					__m256d b = aData[widthA4 * c2 + c1 + c3];
					ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
					ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
				}
			}
#pragma unroll(4)
			for (unsigned int c3(0); c3 < warp; ++c3)
			{
				cData[c0 * widthA4 + c1 + c3] = ans0[c3];
				cData[(c0 + 1) * widthA4 + c1 + c3] = ans1[c3];
			}
		}
}

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(0, 2);
	unsigned int l(1024 * 1024);
	double* a((double*)_mm_malloc(l * sizeof(double), 32));
	double* b((double*)_mm_malloc(l * sizeof(double), 32));
	double* c((double*)_mm_malloc(l * sizeof(double), 32));
	__m256d* av((__m256d*)a);
	__m256d* bv((__m256d*)b);
	__m256d* cv((__m256d*)c);

	::printf("%p\n", a);
	randomMat(a, l, mt, rd);
	randomMat(b, l, mt, rd);
	clear(c, l);

	auto t1=std::chrono::system_clock::now();
	matMult(a, b, c, 1024, 1024, 1024, 1024);
	auto t2=std::chrono::system_clock::now();
	::printf("%d us\n", std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count());

	::free(a);
	::free(b);
	::free(c);
}