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
//still have bugs
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
		_beginthreadex_proc_type lambda = [](void* ptr)->unsigned int
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

			constexpr unsigned long long warp = 2;
			constexpr unsigned long long group = 4;
			for (; source < sourceEnding; source += width4d * group, bData += aWidth256d * group)
			{
				unsigned long long c1(0);
				for (; c1 < aWidthWarpFloor; c1 += warp)
				{
					__m256d ans0[warp] = { 0 };
					__m256d ans1[warp] = { 0 };
					__m256d ans2[warp] = { 0 };
					__m256d ans3[warp] = { 0 };
					for (unsigned long long c2(0); c2 < minDim; ++c2)
					{
						__m256d tp0 = _mm256_set1_pd(source[c2]);
						__m256d tp1 = _mm256_set1_pd(source[width4d + c2]);
						__m256d tp2 = _mm256_set1_pd(source[2 * width4d + c2]);
						__m256d tp3 = _mm256_set1_pd(source[3 * width4d + c2]);
						__m256d a = aData[aWidth256d * c2 + c1];
						__m256d b = aData[aWidth256d * c2 + c1 + 1];
						ans0[0] = _mm256_fmadd_pd(tp0, a, ans0[0]);
						ans0[1] = _mm256_fmadd_pd(tp0, b, ans0[1]);
						ans1[0] = _mm256_fmadd_pd(tp1, a, ans1[0]);
						ans1[1] = _mm256_fmadd_pd(tp1, b, ans1[1]);
						ans2[0] = _mm256_fmadd_pd(tp2, a, ans2[0]);
						ans2[1] = _mm256_fmadd_pd(tp2, b, ans2[1]);
						ans3[0] = _mm256_fmadd_pd(tp3, a, ans3[0]);
						ans3[1] = _mm256_fmadd_pd(tp3, b, ans3[1]);
					}
#pragma unroll(4)
					for (unsigned long long c3(0); c3 < warp; ++c3)
					{
						bData[c1 + c3] = ans0[c3];
						bData[aWidth256d + c1 + c3] = ans1[c3];
						bData[2 * aWidth256d + c1 + c3] = ans2[c3];
						bData[3 * aWidth256d + c1 + c3] = ans3[c3];
					}
				}
				/*if (c1 < aWidth256d)
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
				}*/
			}
			return 0;
		};

		unsigned long long threadHeight;
		unsigned long long rowBeginning;
		unsigned long long leftHeight;
		if (source->height / threadNum >= 4)
		{
			threadHeight = (source->height / threadNum) & (-4);
		}
		else
		{
			threadHeight = 4;
			threadNum = source->height / 4;
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
				threads[c0] = (HANDLE)_beginthreadex(nullptr, 0, lambda, paras + c0 * 9, 0, nullptr);
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
		/*if (leftHeight > 1)
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
				threads[c0] = (HANDLE)_beginthreadex(nullptr, 0, lambda, paras + c0 * 9, 0, nullptr);
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
		}*/
		free(threads);
		free(paras);
		/*if (source->height & 1)
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
		}*/
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

	/*for (unsigned int dim(1024); dim < 2050; ++dim)
	{
		mat mA(dim, dim, false);
		mat mB(dim, dim, false);
		mat mC(dim, dim, false);
		mat mD(dim, dim, false);
		randomMat(mA, mt, rd);
		randomMat(mB, mt, rd);
		//mA(mB, mC);
		//matMultMT(mA, mB, mD);
		//double eps(0);
		//mC -= mD;
		//for (unsigned int c0(0); c0 < dim; ++c0)
		//	eps += mC.row(c0).norm2();
		//if (eps > 0)
		//	printf("\r%4u: %.6e\n", dim, eps);
		//else
		printf("\r%4u", dim);
	}*/

	/*_beginthreadex_proc_type lambda = //(void(__stdcall*)(void*))
		[](void* ptr)->unsigned
	{
		printf("ahh\n");
		return 0;
	};

	void* ptr(nullptr);
	unsigned int threadID;

	HANDLE thread((HANDLE)_beginthreadex(nullptr, 0, lambda, ptr, 0, &threadID));
	WaitForSingleObjectEx(thread, INFINITE, false);
	CloseHandle(thread);*/

	mat mA(1024, 1024, false);
	mat mB(1024, 1024, false);
	mat mC(1024, 1024, false);
	randomMat(mA, mt, rd);
	randomMat(mB, mt, rd);

	//timer.begin();
	//mA(mB, mC);
	//timer.end();
	//timer.print("Single Thread: ");

	timer.begin();
	for (unsigned long long c0(0); c0 < 10; ++c0)
		//mA(mB, mC);
		matMultMT(mA, mB, mC);
	timer.end();
	timer.print("3950X AVX2 Multi Thread: ");
	//mA.printToTableTxt("./matA.txt");
	//mB.printToTableTxt("./matB.txt");
	//mC.printToTableTxt("./matC.txt");
}