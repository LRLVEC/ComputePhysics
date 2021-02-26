#include <_BLAS.h>
#include <_Time.h>
#include <random>
#ifdef _WIN32
#include <Windows.h>
#include <process.h>
#else
#include <pthread.h>
#include <sys/sysinfo.h>
#endif
BLAS::mat& matMultMT_(BLAS::mat const& ts, BLAS::mat const& a, BLAS::mat& b)
{
	using namespace BLAS;
	unsigned long long threadNum;
#ifdef _WIN32
	SYSTEM_INFO systemInfo;
	GetSystemInfo(&systemInfo);
	threadNum = systemInfo.dwNumberOfProcessors;
#else
	threadNum = get_nprocs_conf();
#endif
	unsigned long long minDim(ts.width > a.height ? a.height : ts.width);
	if (minDim)
	{
		bool overflow(ceiling4(a.width, ts.height) > b.width4d * b.height);
		if (overflow && b.type != Type::Native)return b;
		mat const* source(&ts);
		mat r;
		if (&b == &ts)
		{
			source = &r;
			r = ts;
		}
		if (overflow)b.reconstruct(a.width, source->height, false);
		__m256d* aData((__m256d*)a.data);
		__m256d* bData((__m256d*)b.data);
		constexpr unsigned long long warp = 16;
		unsigned long long aWidth256d(a.width4d / 4);
		unsigned long long aWidthWarpFloor((aWidth256d / warp) * warp);
		unsigned long long warpLeft(aWidth256d - aWidthWarpFloor);
		unsigned long long height2Floor(source->height & (-2));

		//ptr: A beginning
		//     A ending
		//     B
		//     C beginning
		//     width4d
		//     aWidthWarpFloor
		//     minDim
		//     aWidth256d
		//     warpLeft
#ifdef _WIN32
		void (*lambda)(void*) = [](void* ptr)
#else
		void* (*lambda)(void*) = [](void* ptr)->void*
#endif
		{
			void** ptrs((void**)ptr);
			double* source = ((double*)ptrs[0]);
			double* sourceEnding((double*)ptrs[1]);
			__m256d* aData((__m256d*)ptrs[2]);
			__m256d* bData((__m256d*)ptrs[3]);
			unsigned long long width4d((unsigned long long)ptrs[4]);
			unsigned long long aWidthWarpFloor((unsigned long long)ptrs[5]);
			unsigned long long minDim((unsigned long long)ptrs[6]);
			unsigned long long aWidth256d((unsigned long long)ptrs[7]);
			unsigned long long warpLeft((unsigned long long)ptrs[8]);

			constexpr unsigned long long warp = 16;
			for (; source < sourceEnding; source += width4d * 2, bData += aWidth256d * 2)
			{
				unsigned long long c1(0);
				for (; c1 < aWidthWarpFloor; c1 += warp)
				{
					__m256d ans0[warp] = { 0 };
					__m256d ans1[warp] = { 0 };
					for (unsigned long long c2(0); c2 < minDim; ++c2)
					{
						__m256d tp0 = _mm256_set1_pd(source[c2]);
						__m256d tp1 = _mm256_set1_pd(source[width4d + c2]);
#pragma unroll(4)
						for (unsigned long long c3(0); c3 < warp; ++c3)
						{
							__m256d b = aData[aWidth256d * c2 + c1 + c3];
							ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
							ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
						}
					}
#pragma unroll(4)
					for (unsigned long long c3(0); c3 < warp; ++c3)
					{
						bData[c1 + c3] = ans0[c3];
						bData[aWidth256d + c1 + c3] = ans1[c3];
					}
				}
				if (c1 < aWidth256d)
				{
					__m256d ans0[warp] = { 0 };
					__m256d ans1[warp] = { 0 };
					for (unsigned long long c2(0); c2 < minDim; ++c2)
					{
						__m256d tp0 = _mm256_set1_pd(source[c2]);
						__m256d tp1 = _mm256_set1_pd(source[width4d + c2]);
						for (unsigned long long c3(0); c3 < warpLeft; ++c3)
						{
							__m256d b = aData[aWidth256d * c2 + c1 + c3];
							ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
							ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
						}
					}
					for (unsigned long long c3(0); c3 < warpLeft; ++c3)
					{
						bData[c1 + c3] = ans0[c3];
						bData[aWidth256d + c1 + c3] = ans1[c3];
					}
				}
			}
#ifndef _WIN32
			return 0;
#endif
		};

		unsigned long long threadHeight;
		unsigned long long rowBeginning;
		unsigned long long leftHeight;
		if (source->height / threadNum >= 2)
		{
			threadHeight = (source->height / threadNum) & (-2);
		}
		else
		{
			threadHeight = 2;
			threadNum = source->height / 2;
		}
		rowBeginning = threadHeight * threadNum;
		leftHeight = source->height - rowBeginning;

#ifdef _WIN32
		HANDLE* threads(nullptr);
#else
		pthread_t* threads(nullptr);
#endif
		void** paras(nullptr);
		if (threadNum)
		{
#ifdef _WIN32
			threads = (HANDLE*)::malloc(threadNum * sizeof(HANDLE));
#else
			threads = (pthread_t*)::malloc(threadNum * sizeof(pthread_t));
#endif
			paras = (void**)::malloc(threadNum * 9 * sizeof(void*));
			for (unsigned long long c0(0); c0 < threadNum; ++c0)
			{
				paras[c0 * 9] = (void*)(source->data + threadHeight * source->width4d * c0);			//A beginning
				paras[c0 * 9 + 1] = (void*)(source->data + threadHeight * source->width4d * (c0 + 1));	//A ending
				paras[c0 * 9 + 2] = (void*)(a.data);													//B
				paras[c0 * 9 + 3] = (void*)(b.data + threadHeight * a.width4d * c0);					//C beginning
				paras[c0 * 9 + 4] = (void*)(source->width4d);											//width4d
				paras[c0 * 9 + 5] = (void*)(aWidthWarpFloor);											//aWidthWarpFloor
				paras[c0 * 9 + 6] = (void*)(minDim);													//minDim
				paras[c0 * 9 + 7] = (void*)(aWidth256d);												//aWidth256d
				paras[c0 * 9 + 8] = (void*)(warpLeft);													//warpLeft
#ifdef _WIN32
				threads[c0] = (HANDLE)_beginthread(lambda, 0, paras + c0 * 9);
#else
				pthread_create(threads + c0, 0, lambda, paras + c0 * 9);
				pthread_detach(threads[c0]);
#endif
			}
#ifdef _WIN32
			DWORD rc = WaitForMultipleObjects(threadNum, threads, true, INFINITE);
			for (unsigned long long c0(0); c0 < threadNum; ++c0)
				CloseHandle(threads[c0]);
#else
			for (unsigned long long c0(0); c0 < threadNum; ++c0)
				pthread_join(threads[c0], nullptr);
#endif
		}
		if (leftHeight > 1)
		{
			threadHeight = 2;
			threadNum = leftHeight / 2;
			for (unsigned long long c0(0); c0 < threadNum; ++c0)
			{
				paras[c0 * 9] = (void*)(source->data + source->width4d * (threadHeight * c0 + rowBeginning));			//A beginning
				paras[c0 * 9 + 1] = (void*)(source->data + source->width4d * (threadHeight * (c0 + 1) + rowBeginning));	//A ending
				paras[c0 * 9 + 2] = (void*)(a.data);																	//B
				paras[c0 * 9 + 3] = (void*)(b.data + a.width4d * (threadHeight * c0 + rowBeginning));					//C beginning
				paras[c0 * 9 + 4] = (void*)(source->width4d);															//width4d
				paras[c0 * 9 + 5] = (void*)(aWidthWarpFloor);															//aWidthWarpFloor
				paras[c0 * 9 + 6] = (void*)(minDim);																	//minDim
				paras[c0 * 9 + 7] = (void*)(aWidth256d);																//aWidth256d
				paras[c0 * 9 + 8] = (void*)(warpLeft);																	//warpLeft
#ifdef _WIN32
				threads[c0] = (HANDLE)_beginthread(lambda, 0, paras + c0 * 9);
#else
				pthread_create(threads + c0, 0, lambda, paras + c0 * 9);
				pthread_detach(threads[c0]);
#endif
			}
#ifdef _WIN32
			DWORD rc = WaitForMultipleObjects(threadNum, threads, true, INFINITE);
			for (unsigned long long c0(0); c0 < threadNum; ++c0)
				CloseHandle(threads[c0]);
#else
			for (unsigned long long c0(0); c0 < threadNum; ++c0)
				pthread_join(threads[c0], nullptr);
#endif
			rowBeginning += threadHeight * threadNum;
		}
		free(threads);
		free(paras);
		if (source->height & 1)
		{
			double* sourceBeginning(source->data + rowBeginning * source->width4d);
			bData += rowBeginning * aWidth256d;

			unsigned long long c1(0);
			for (; c1 < aWidthWarpFloor; c1 += warp)
			{
				__m256d ans0[warp] = { 0 };
				for (unsigned long long c2(0); c2 < minDim; ++c2)
				{
					__m256d tp0 = _mm256_set1_pd(sourceBeginning[c2]);
#pragma unroll(4)
					for (unsigned long long c3(0); c3 < warp; ++c3)
					{
						__m256d b = aData[aWidth256d * c2 + c1 + c3];
						ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
					}
				}
#pragma unroll(4)
				for (unsigned long long c3(0); c3 < warp; ++c3)
					bData[c1 + c3] = ans0[c3];
			}
			if (c1 < aWidth256d)
			{
				__m256d ans0[warp] = { 0 };
				for (unsigned long long c2(0); c2 < minDim; ++c2)
				{
					__m256d tp0 = _mm256_set1_pd(sourceBeginning[c2]);
					for (unsigned long long c3(0); c3 < warpLeft; ++c3)
					{
						__m256d b = aData[aWidth256d * c2 + c1 + c3];
						ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
					}
				}
				for (unsigned long long c3(0); c3 < warpLeft; ++c3)
					bData[c1 + c3] = ans0[c3];
			}
		}
	}
	return b;
}

BLAS::mat& matMultMT(BLAS::mat const& ts, BLAS::mat const& a, BLAS::mat& b)
{
	using namespace BLAS;

	unsigned long long threadNum;
#ifdef _WIN32
	SYSTEM_INFO systemInfo;
	GetSystemInfo(&systemInfo);
	threadNum = systemInfo.dwNumberOfProcessors;
#else
	threadNum = get_nprocs_conf();
#endif
	::printf("Number of processors: %llu\n", threadNum);

	unsigned long long minDim(ts.width > a.height ? a.height : ts.width);
	if (minDim)
	{
		bool overflow(ceiling4(a.width, ts.height) > b.width4d * b.height);
		if (overflow && b.type != Type::Native)return b;
		mat const* source(&ts);
		mat r;
		if (&b == &ts)
		{
			source = &r;
			r = ts;
		}
		if (overflow)b.reconstruct(a.width, source->height, false);
		__m256d* aData((__m256d*)a.data);
		__m256d* bData((__m256d*)b.data);
		constexpr unsigned long long warp = 16;
		unsigned long long aWidth256d(a.width4d / 4);
		unsigned long long aWidthWarpFloor((aWidth256d / warp) * warp);
		unsigned long long warpLeft(aWidth256d - aWidthWarpFloor);
		unsigned long long height2Floor(source->height & (-2));

		//ptr: A beginning
		//     A ending
		//     B
		//     C beginning
		//     width4d
		//     aWidthWarpFloor
		//     minDim
		//     aWidth256d
		//     warpLeft
#ifdef _WIN32
		void (*lambda)(void*) = [](void* ptr)
#else
		void* (*lambda)(void*) = [](void* ptr)->void*
#endif
		{
			void** ptrs((void**)ptr);
			double* source = ((double*)ptrs[0]);
			double* sourceEnding((double*)ptrs[1]);
			__m256d* aData((__m256d*)ptrs[2]);
			__m256d* bData((__m256d*)ptrs[3]);
			unsigned long long width4d((unsigned long long)ptrs[4]);
			unsigned long long aWidthWarpFloor((unsigned long long)ptrs[5]);
			unsigned long long minDim((unsigned long long)ptrs[6]);
			unsigned long long aWidth256d((unsigned long long)ptrs[7]);
			unsigned long long warpLeft((unsigned long long)ptrs[8]);

			constexpr unsigned long long warp = 16;
			for (; source < sourceEnding; source += width4d * 2, bData += aWidth256d * 2)
			{
				unsigned long long c1(0);
				for (; c1 < aWidthWarpFloor; c1 += warp)
				{
					__m256d ans0[warp] = { 0 };
					__m256d ans1[warp] = { 0 };
					for (unsigned long long c2(0); c2 < minDim; ++c2)
					{
						__m256d tp0 = _mm256_set1_pd(source[c2]);
						__m256d tp1 = _mm256_set1_pd(source[width4d + c2]);
#pragma unroll(4)
						for (unsigned long long c3(0); c3 < warp; ++c3)
						{
							__m256d b = aData[aWidth256d * c2 + c1 + c3];
							ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
							ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
						}
					}
#pragma unroll(4)
					for (unsigned long long c3(0); c3 < warp; ++c3)
					{
						bData[c1 + c3] = ans0[c3];
						bData[aWidth256d + c1 + c3] = ans1[c3];
					}
				}
				if (c1 < aWidth256d)
				{
					__m256d ans0[warp] = { 0 };
					__m256d ans1[warp] = { 0 };
					for (unsigned long long c2(0); c2 < minDim; ++c2)
					{
						__m256d tp0 = _mm256_set1_pd(source[c2]);
						__m256d tp1 = _mm256_set1_pd(source[width4d + c2]);
						for (unsigned long long c3(0); c3 < warpLeft; ++c3)
						{
							__m256d b = aData[aWidth256d * c2 + c1 + c3];
							ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
							ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
						}
					}
					for (unsigned long long c3(0); c3 < warpLeft; ++c3)
					{
						bData[c1 + c3] = ans0[c3];
						bData[aWidth256d + c1 + c3] = ans1[c3];
					}
				}
			}
#ifndef _WIN32
			return 0;
#endif
		};

		unsigned long long threadHeight((source->height / threadNum) & (-2));
		unsigned long long rowBeginning;
#ifdef _WIN32
		HANDLE* threads(nullptr);
#else
		pthread_t* threads(nullptr);
#endif
		void** paras(nullptr);
		if (threadHeight)
		{
			rowBeginning = threadHeight * (threadNum - 1);
#ifdef _WIN32
			threads = (HANDLE*)::malloc((threadNum - 1) * sizeof(HANDLE));
#else
			threads = (pthread_t*)::malloc((threadNum - 1) * sizeof(pthread_t));
#endif
			paras = (void**)::malloc((threadNum - 1) * 9 * sizeof(void*));

			for (unsigned long long c0(0); c0 < threadNum - 1; ++c0)
			{
				paras[c0 * 9] = (void*)(source->data + threadHeight * source->width4d * c0);			//A beginning
				paras[c0 * 9 + 1] = (void*)(source->data + threadHeight * source->width4d * (c0 + 1));	//A ending
				paras[c0 * 9 + 2] = (void*)(a.data);													//B
				paras[c0 * 9 + 3] = (void*)(b.data + threadHeight * a.width4d * c0);					//C beginning
				paras[c0 * 9 + 4] = (void*)(source->width4d);											//width4d
				paras[c0 * 9 + 5] = (void*)(aWidthWarpFloor);											//aWidthWarpFloor
				paras[c0 * 9 + 6] = (void*)(minDim);													//minDim
				paras[c0 * 9 + 7] = (void*)(aWidth256d);												//aWidth256d
				paras[c0 * 9 + 8] = (void*)(warpLeft);													//warpLeft

#ifdef _WIN32
				threads[c0] = (HANDLE)_beginthread(lambda, 0, paras + c0 * 9);
#else
				pthread_create(threads + c0, 0, lambda, paras + c0 * 9);
				pthread_detach(threads[c0]);
#endif
			}
		}
		else rowBeginning = 0;

		double* sourceBeginning(source->data + rowBeginning * source->width4d);
		double* sourceEnding(source->data + (source->height & -2) * source->width4d);
		unsigned long long width4d(source->width4d);
		bData += rowBeginning * aWidth256d;

		for (; sourceBeginning < sourceEnding; sourceBeginning += width4d * 2, bData += aWidth256d * 2)
		{
			unsigned long long c1(0);
			for (; c1 < aWidthWarpFloor; c1 += warp)
			{
				__m256d ans0[warp] = { 0 };
				__m256d ans1[warp] = { 0 };
				for (unsigned long long c2(0); c2 < minDim; ++c2)
				{
					__m256d tp0 = _mm256_set1_pd(sourceBeginning[c2]);
					__m256d tp1 = _mm256_set1_pd(sourceBeginning[width4d + c2]);
#pragma unroll(4)
					for (unsigned long long c3(0); c3 < warp; ++c3)
					{
						__m256d b = aData[aWidth256d * c2 + c1 + c3];
						ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
						ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
					}
				}
#pragma unroll(4)
				for (unsigned long long c3(0); c3 < warp; ++c3)
				{
					bData[c1 + c3] = ans0[c3];
					bData[aWidth256d + c1 + c3] = ans1[c3];
				}
			}
			if (c1 < aWidth256d)
			{
				__m256d ans0[warp] = { 0 };
				__m256d ans1[warp] = { 0 };
				for (unsigned long long c2(0); c2 < minDim; ++c2)
				{
					__m256d tp0 = _mm256_set1_pd(sourceBeginning[c2]);
					__m256d tp1 = _mm256_set1_pd(sourceBeginning[width4d + c2]);
					for (unsigned long long c3(0); c3 < warpLeft; ++c3)
					{
						__m256d b = aData[aWidth256d * c2 + c1 + c3];
						ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
						ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
					}
				}
				for (unsigned long long c3(0); c3 < warpLeft; ++c3)
				{
					bData[c1 + c3] = ans0[c3];
					bData[aWidth256d + c1 + c3] = ans1[c3];
				}
			}
		}
		if (source->height & 1)
		{
			unsigned long long c1(0);
			for (; c1 < aWidthWarpFloor; c1 += warp)
			{
				__m256d ans0[warp] = { 0 };
				for (unsigned long long c2(0); c2 < minDim; ++c2)
				{
					__m256d tp0 = _mm256_set1_pd(sourceBeginning[c2]);
#pragma unroll(4)
					for (unsigned long long c3(0); c3 < warp; ++c3)
					{
						__m256d b = aData[aWidth256d * c2 + c1 + c3];
						ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
					}
				}
#pragma unroll(4)
				for (unsigned long long c3(0); c3 < warp; ++c3)
					bData[c1 + c3] = ans0[c3];
			}
			if (c1 < aWidth256d)
			{
				__m256d ans0[warp] = { 0 };
				for (unsigned long long c2(0); c2 < minDim; ++c2)
				{
					__m256d tp0 = _mm256_set1_pd(sourceBeginning[c2]);
					for (unsigned long long c3(0); c3 < warpLeft; ++c3)
					{
						__m256d b = aData[aWidth256d * c2 + c1 + c3];
						ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
					}
				}
				for (unsigned long long c3(0); c3 < warpLeft; ++c3)
					bData[c1 + c3] = ans0[c3];
			}
		}
		if (threadHeight)
		{
#ifdef _WIN32
			DWORD rc = WaitForMultipleObjects(threadNum - 1, threads, true, INFINITE);
			int slot = rc - WAIT_OBJECT_0;
			if (slot >= 0 && slot < threadNum - 1)
				printf("All thread terminite\n");
			for (unsigned long long c0(0); c0 < threadNum - 1; ++c0)
				CloseHandle(threads[c0]);
#else
			for (unsigned long long c0(0); c0 < threadNum - 1; ++c0)
				pthread_join(threads[c0], nullptr);
#endif
			::free(threads);
			::free(paras);
		}
	}
	return b;
}

int main()
{
	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(-1.0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);

	remove("./matA.txt");
	remove("./matB.txt");
	remove("./matC.txt");

	Timer timer;
	using namespace BLAS;
	// 3950x avx2 multi thread test
	mat mA(1023, 1023, false);
	mat mB(1023, 1023, false);
	mat mC(1023, 1023, false);
	randomMat(mA, mt, rd);
	randomMat(mB, mt, rd);
	mA.printToTableTxt("./matA.txt");
	mB.printToTableTxt("./matB.txt");

	//timer.begin();
	//mA(mB, mC);
	//timer.end();
	//timer.print("Single Thread: ");

	timer.begin();
	for (unsigned long long c0(0); c0 < 1; ++c0)
		//mA(mB, mC);
		matMultMT_(mA, mB, mC);
	timer.end();
	timer.print("Multi Thread: ");
	mC.printToTableTxt("./matC.txt");
}