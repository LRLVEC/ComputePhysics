#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <immintrin.h>
#include <random>
//#include <_Time.h>
#include <chrono>


void clear(double* a, unsigned int length)
{
	for (unsigned int c0(0); c0 < length; ++c0)
		a[c0] = 0;
}
void randomIt(double* a, unsigned int length, std::mt19937& mt, std::uniform_real_distribution<double>& rd)
{
	for (unsigned int c0(0); c0 < length; ++c0)
		a[c0] = rd(mt);
}

void matMultmat(double* a, double* b, double* c, unsigned int widthA,
	unsigned int heightA, unsigned int widthB, unsigned int heightB)
{
	__m256d* aData((__m256d*)b);
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

void matMultvec(double* a, double* b, double* c, unsigned int widthA,
	unsigned int heightA, unsigned int lengthB)
{
	unsigned int minDim(widthA > lengthB ? lengthB : widthA);
	if (minDim && heightA)
	{
		constexpr unsigned int warp = 8;
		unsigned int minWidth4((minDim - 1) / 4 + 1);
		__m256d* aData((__m256d*)a);
		__m256d* bData((__m256d*)b);
		__m256d* rData((__m256d*)c);
		for (unsigned int c0(0); c0 < heightA; c0 += 4)
		{
			__m256d ans[4] = { 0 };
			__m256d tp[warp];
			for (unsigned int c1(0); c1 < minWidth4; c1 += warp)
			{
				__m256d* s(aData + minWidth4 * c0 + c1);
#pragma unroll(4)
				for (unsigned int c2(0); c2 < warp; ++c2)
					tp[c2] = bData[c1 + c2];
				for (unsigned int c2(0); c2 < 4; ++c2, s += minWidth4)
				{
#pragma unroll(4)
					for (unsigned int c3(0); c3 < warp; ++c3)
					{
						__m256d t = s[c3];
						ans[c2] = _mm256_fmadd_pd(t, tp[c3], ans[c2]);
					}
				}
			}
			__m256d s;
			for (unsigned int c1(0); c1 < 4; ++c1)
			{
				s.m256d_f64[c1] = ans[c1].m256d_f64[0];
				s.m256d_f64[c1] += ans[c1].m256d_f64[1];
				s.m256d_f64[c1] += ans[c1].m256d_f64[2];
				s.m256d_f64[c1] += ans[c1].m256d_f64[3];
			}
			rData[c0 >> 2] = s;
		}
	}
}

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(0, 2);
	unsigned int l(1024 * 1024);
	double* vecA((double*)_mm_malloc(1024 * sizeof(double), 32));
	double* vecB((double*)_mm_malloc(1024 * sizeof(double), 32));
	double* matA((double*)_mm_malloc(l * sizeof(double), 32));
	double* matB((double*)_mm_malloc(l * sizeof(double), 32));
	double* matC((double*)_mm_malloc(l * sizeof(double), 32));

	Timer timer;


	randomIt(vecA, 1024, mt, rd);
	//randomIt(vecB, 1024, mt, rd);
	randomIt(matA, l, mt, rd);
	randomIt(matB, l, mt, rd);
	//clear(matC, l);

	auto t1 = std::chrono::system_clock::now();
	for (unsigned int c0(0); c0 < 100; ++c0)
		matMultvec(matA, vecA, vecB, 1024, 1024, 1024);
	auto t2 = std::chrono::system_clock::now();
	::printf("%d us\n", std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());

	//timer.begin();
	//matMultmat(matA, matB, matC, 1024, 1024, 1024, 1024);
	//timer.end();
	//timer.print("AVX2:");

	//timer.begin();
	//timer.end();
	//timer.print("mat mult vec (AVX2): ");

	::_mm_free(vecA);
	::_mm_free(vecB);
	::_mm_free(matA);
	::_mm_free(matB);
	::_mm_free(matC);
}