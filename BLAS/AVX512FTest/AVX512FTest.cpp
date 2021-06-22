#include <_BLAS.h>
#include <_Time.h>
#include <random>
#include <Windows.h>
#include <process.h>

namespace BLAS
{
	inline unsigned long long ceiling8(unsigned long long length)
	{
		return (((length - 1) >> 3) + 1) << 3;
	}
	inline unsigned long long floor8(unsigned long long length)
	{
		return (length >> 3) << 3;
	}
	inline unsigned long long ceiling8(unsigned long long width, unsigned long long height)
	{
		return ((((width - 1) >> 3) + 1) << 3) * height;
	}
	inline size_t ceiling512dSize(unsigned long long length)
	{
		return (((size_t(length) - 1) >> 3) + 1) << 6;
	}
	inline size_t ceiling512dSize(unsigned long long width, unsigned long long height)
	{
		return ((((size_t(width) - 1) >> 3) + 1) << 6) * height;
	}
	inline double* malloc64d_512(unsigned long long length)
	{
		return (double*)_mm_malloc(length * sizeof(double), 64);
	}
	inline double* malloc512d(unsigned long long length)
	{
		return (double*)_mm_malloc(ceiling512dSize(length), 64);
	}
	inline double* malloc512d(unsigned long long width, unsigned long long height)
	{
		return (double*)_mm_malloc(ceiling512dSize(width) * height, 64);
	}
	inline void* memcpy512d(void* dst, void const* src, unsigned long long length)
	{
		return ::memcpy(dst, src, ceiling512dSize(length));
	}
	inline void* memcpy512d(void* dst, void const* src, unsigned long long width, unsigned long long height)
	{
		return ::memcpy(dst, src, ceiling512dSize(width, height));
	}
	inline void* memset512d(void* dst, int val, unsigned long long length)
	{
		return ::memset(dst, val, ceiling512dSize(length));
	}
	inline void* memset512d(void* dst, int val, unsigned long long width, unsigned long long height)
	{
		return ::memset(dst, val, ceiling512dSize(width, height));
	}

	inline unsigned long long getPtrOffset64d_512(double* ptr)
	{
		return ((unsigned long long)ptr >> 3) & 7;
	}
	inline double* getPtr512d(double* ptr)
	{
		return (double*)((unsigned long long)ptr & -64);
	}


	struct mat512
	{
		double* data;
		unsigned long long width;
		unsigned long long height;
		unsigned long long width8d;
		mat512() :data(nullptr), width(0), height(0), width8d(0) {}
		mat512(unsigned long long _width, unsigned long long _height, bool _clear = true)
			:
			data((_width&& _height) ? malloc512d(_width, _height) : nullptr),
			width(data ? _width : 0),
			height(data ? _height : 0),
			width8d(data ? ceiling8(_width) : 0)
		{
			if (_clear && data)
				memset512d(data, 0, width, height);
		}
		~mat512()
		{
			_mm_free(data);
			data = nullptr;
			width = height = 0;
		}
		inline double& operator() (unsigned long long a, unsigned long long b)
		{
			return data[a * width8d + b];
		}
		void clear()
		{
			if (data)memset64d(data, 0, height * width8d);
		}
		void reconstruct(unsigned long long _width, unsigned long long _height, bool _clear = true)
		{
			if (data)_mm_free(data);
			unsigned long long s(ceiling512dSize(_width, _height));
			if (s)
			{
				data = (double*)malloc64d_512(s);
				if (_clear)memset64d(data, 0, s);
				width = _width;
				height = _height;
			}
		}

		mat512& operator()(mat512 const& a, mat512& b)const
		{
			unsigned long long minDim(width > a.height ? a.height : width);
			if (minDim)
			{
				bool overflow(ceiling8(a.width, height) < b.width8d * b.height);
				if (overflow)return b;
				mat512 const* source(this);
				mat512 r;
				if (&b == this)
				{
					source = &r;
					r = *this;
				}
				if (overflow)b.reconstruct(a.width, height, false);
				__m512d* aData((__m512d*)a.data);
				__m512d* bData((__m512d*)b.data);
				constexpr unsigned long long warp = 2;
				constexpr unsigned long long group = 8;
				unsigned long long aWidth512d(a.width8d / 8);
				unsigned long long aWidthWarpFloor(aWidth512d / warp * warp);
				unsigned long long warpLeft(aWidth512d - aWidthWarpFloor);
				unsigned long long height4Floor(height & (-long long(group)));
				unsigned long long c0(0);
				for (; c0 < height4Floor; c0 += group)
				{
					unsigned long long c1(0);
					for (; c1 < aWidthWarpFloor; c1 += warp)
					{
						__m512d ans0[warp] = { 0 };
						__m512d ans1[warp] = { 0 };
						__m512d ans2[warp] = { 0 };
						__m512d ans3[warp] = { 0 };
						__m512d ans4[warp] = { 0 };
						__m512d ans5[warp] = { 0 };
						__m512d ans6[warp] = { 0 };
						__m512d ans7[warp] = { 0 };
						for (unsigned long long c2(0); c2 < minDim; ++c2)
						{
							__m512d tp0 = _mm512_set1_pd(source->data[c0 * width8d + c2]);
							__m512d tp1 = _mm512_set1_pd(source->data[(c0 + 1) * width8d + c2]);
							__m512d tp2 = _mm512_set1_pd(source->data[(c0 + 2) * width8d + c2]);
							__m512d tp3 = _mm512_set1_pd(source->data[(c0 + 3) * width8d + c2]);
							__m512d tp4 = _mm512_set1_pd(source->data[(c0 + 4) * width8d + c2]);
							__m512d tp5 = _mm512_set1_pd(source->data[(c0 + 5) * width8d + c2]);
							__m512d tp6 = _mm512_set1_pd(source->data[(c0 + 6) * width8d + c2]);
							__m512d tp7 = _mm512_set1_pd(source->data[(c0 + 7) * width8d + c2]);
							__m512d a = aData[aWidth512d * c2 + c1];
							__m512d b = aData[aWidth512d * c2 + c1 + 1];
							ans0[0] = _mm512_fmadd_pd(tp0, a, ans0[0]);
							ans0[1] = _mm512_fmadd_pd(tp0, b, ans0[1]);
							ans1[0] = _mm512_fmadd_pd(tp1, a, ans1[0]);
							ans1[1] = _mm512_fmadd_pd(tp1, b, ans1[1]);
							ans2[0] = _mm512_fmadd_pd(tp2, a, ans2[0]);
							ans2[1] = _mm512_fmadd_pd(tp2, b, ans2[1]);
							ans3[0] = _mm512_fmadd_pd(tp3, a, ans3[0]);
							ans3[1] = _mm512_fmadd_pd(tp3, b, ans3[1]);
							ans4[0] = _mm512_fmadd_pd(tp4, a, ans4[0]);
							ans4[1] = _mm512_fmadd_pd(tp4, b, ans4[1]);
							ans5[0] = _mm512_fmadd_pd(tp5, a, ans5[0]);
							ans5[1] = _mm512_fmadd_pd(tp5, b, ans5[1]);
							ans6[0] = _mm512_fmadd_pd(tp6, a, ans6[0]);
							ans6[1] = _mm512_fmadd_pd(tp6, b, ans6[1]);
							ans7[0] = _mm512_fmadd_pd(tp7, a, ans7[0]);
							ans7[1] = _mm512_fmadd_pd(tp7, b, ans7[1]);
						}
						for (unsigned long long c3(0); c3 < warp; ++c3)
						{
							bData[c0 * aWidth512d + c1 + c3] = ans0[c3];
							bData[(c0 + 1) * aWidth512d + c1 + c3] = ans1[c3];
							bData[(c0 + 2) * aWidth512d + c1 + c3] = ans2[c3];
							bData[(c0 + 3) * aWidth512d + c1 + c3] = ans3[c3];
							bData[(c0 + 4) * aWidth512d + c1 + c3] = ans4[c3];
							bData[(c0 + 5) * aWidth512d + c1 + c3] = ans5[c3];
							bData[(c0 + 6) * aWidth512d + c1 + c3] = ans6[c3];
							bData[(c0 + 7) * aWidth512d + c1 + c3] = ans7[c3];
						}
					}
					/*if (c1 < aWidth512d)
					{
						__m512d ans0[warp] = { 0 };
						__m512d ans1[warp] = { 0 };
						for (unsigned long long c2(0); c2 < minDim; ++c2)
						{
							__m512d tp0 = _mm512_set1_pd(source->data[c0 * width8d + c2]);
							__m512d tp1 = _mm512_set1_pd(source->data[(c0 + 1) * width8d + c2]);
							for (unsigned long long c3(0); c3 < warpLeft; ++c3)
							{
								__m512d b = aData[aWidth512d * c2 + c1 + c3];
								ans0[c3] = _mm512_fmadd_pd(tp0, b, ans0[c3]);
								ans1[c3] = _mm512_fmadd_pd(tp1, b, ans1[c3]);
							}
						}
						for (unsigned long long c3(0); c3 < warpLeft; ++c3)
						{
							bData[c0 * aWidth512d + c1 + c3] = ans0[c3];
							bData[(c0 + 1) * aWidth512d + c1 + c3] = ans1[c3];
						}
					}*/
				}
				/*if (c0 < height)
				{
					unsigned long long c1(0);
					for (; c1 < aWidthWarpFloor; c1 += warp)
					{
						__m512d ans0[warp] = { 0 };
						for (unsigned long long c2(0); c2 < minDim; ++c2)
						{
							__m512d tp0 = _mm512_set1_pd(source->data[c0 * width8d + c2]);
							for (unsigned long long c3(0); c3 < warp; ++c3)
							{
								__m512d b = aData[aWidth512d * c2 + c1 + c3];
								ans0[c3] = _mm512_fmadd_pd(tp0, b, ans0[c3]);
							}
						}
						for (unsigned long long c3(0); c3 < warp; ++c3)
							bData[c0 * aWidth512d + c1 + c3] = ans0[c3];
					}
					if (c1 < aWidth512d)
					{
						__m512d ans0[warp] = { 0 };
						for (unsigned long long c2(0); c2 < minDim; ++c2)
						{
							__m512d tp0 = _mm512_set1_pd(source->data[c0 * width8d + c2]);
							for (unsigned long long c3(0); c3 < warpLeft; ++c3)
							{
								__m512d b = aData[aWidth512d * c2 + c1 + c3];
								ans0[c3] = _mm512_fmadd_pd(tp0, b, ans0[c3]);
							}
						}
						for (unsigned long long c3(0); c3 < warpLeft; ++c3)
							bData[c0 * aWidth512d + c1 + c3] = ans0[c3];
					}
				}*/
			}
			return b;
		}

		void printToTableTxt(char const* name)const
		{
			//in the form of Mathematica matrix(for Import)
			if (data)
			{
				FILE* temp(::fopen(name, "w+"));
				for (unsigned long long c0(0); c0 < height; ++c0)
				{
					for (unsigned long long c1(0); c1 < width; ++c1)
						::fprintf(temp, "%.16e ", data[width8d * c0 + c1]);
					::fprintf(temp, "\n");
				}
				::fclose(temp);
			}
		}
	};

	template<class T>void randomMat(mat512& a, std::mt19937& mt, T& rd)
	{
		for (unsigned long long c0(0); c0 < a.height; ++c0)
			for (unsigned long long c1(0); c1 < a.width; ++c1)
				a(c0, c1) = rd(mt);
	}

	mat512& matMultMT(mat512 const& ts, mat512 const& a, mat512& b)
	{
		unsigned long long threadNum;
		SYSTEM_INFO systemInfo;
		GetSystemInfo(&systemInfo);
		threadNum = systemInfo.dwNumberOfProcessors;
		unsigned long long minDim(ts.width > a.height ? a.height : ts.width);
		if (minDim)
		{
			bool overflow(ceiling8(a.width, ts.height) < b.width8d * b.height);
			if (overflow)return b;
			mat512 const* source(&ts);
			mat512 r;
			if (&b == &ts)
			{
				source = &r;
				r = ts;
			}
			if (overflow)b.reconstruct(a.width, source->height, false);
			__m512d* aData((__m512d*)a.data);
			__m512d* bData((__m512d*)b.data);
			constexpr unsigned long long warp = 2;
			constexpr unsigned long long group = 8;
			unsigned long long aWidth512d(a.width8d / 8);
			unsigned long long aWidthWarpFloor(aWidth512d / warp * warp);
			unsigned long long warpLeft(aWidth512d - aWidthWarpFloor);
			unsigned long long height4Floor(source->height & (-long long(group)));

			//ptr: A beginning
			//     A ending
			//     B
			//     C beginning
			//     width4d
			//     aWidthWarpFloor
			//     minDim
			//     aWidth256d
			//     warpLeft
			_beginthreadex_proc_type lambda = [](void* ptr)->unsigned int
			{
				void** ptrs((void**)ptr);
				double* source = ((double*)ptrs[0]);
				double* sourceEnding((double*)ptrs[1]);
				__m512d* aData((__m512d*)ptrs[2]);
				__m512d* bData((__m512d*)ptrs[3]);
				unsigned long long width8d((unsigned long long)ptrs[4]);
				unsigned long long aWidthWarpFloor((unsigned long long)ptrs[5]);
				unsigned long long minDim((unsigned long long)ptrs[6]);
				unsigned long long aWidth512d((unsigned long long)ptrs[7]);
				unsigned long long warpLeft((unsigned long long)ptrs[8]);

				constexpr unsigned long long warp = 2;
				constexpr unsigned long long group = 8;
				for (; source < sourceEnding; source += width8d * group, bData += aWidth512d * group)
				{
					unsigned long long c1(0);
					for (; c1 < aWidthWarpFloor; c1 += warp)
					{
						__m512d ans0[warp] = { 0 };
						__m512d ans1[warp] = { 0 };
						__m512d ans2[warp] = { 0 };
						__m512d ans3[warp] = { 0 };
						__m512d ans4[warp] = { 0 };
						__m512d ans5[warp] = { 0 };
						__m512d ans6[warp] = { 0 };
						__m512d ans7[warp] = { 0 };
						for (unsigned long long c2(0); c2 < minDim; ++c2)
						{
							__m512d tp0 = _mm512_set1_pd(source[c2]);
							__m512d tp1 = _mm512_set1_pd(source[width8d + c2]);
							__m512d tp2 = _mm512_set1_pd(source[2 * width8d + c2]);
							__m512d tp3 = _mm512_set1_pd(source[3 * width8d + c2]);
							__m512d tp4 = _mm512_set1_pd(source[4 * width8d + c2]);
							__m512d tp5 = _mm512_set1_pd(source[5 * width8d + c2]);
							__m512d tp6 = _mm512_set1_pd(source[6 * width8d + c2]);
							__m512d tp7 = _mm512_set1_pd(source[7 * width8d + c2]);
							__m512d a = aData[aWidth512d * c2 + c1];
							__m512d b = aData[aWidth512d * c2 + c1 + 1];
							ans0[0] = _mm512_fmadd_pd(tp0, a, ans0[0]);
							ans0[1] = _mm512_fmadd_pd(tp0, b, ans0[1]);
							ans1[0] = _mm512_fmadd_pd(tp1, a, ans1[0]);
							ans1[1] = _mm512_fmadd_pd(tp1, b, ans1[1]);
							ans2[0] = _mm512_fmadd_pd(tp2, a, ans2[0]);
							ans2[1] = _mm512_fmadd_pd(tp2, b, ans2[1]);
							ans3[0] = _mm512_fmadd_pd(tp3, a, ans3[0]);
							ans3[1] = _mm512_fmadd_pd(tp3, b, ans3[1]);
							ans4[0] = _mm512_fmadd_pd(tp4, a, ans4[0]);
							ans4[1] = _mm512_fmadd_pd(tp4, b, ans4[1]);
							ans5[0] = _mm512_fmadd_pd(tp5, a, ans5[0]);
							ans5[1] = _mm512_fmadd_pd(tp5, b, ans5[1]);
							ans6[0] = _mm512_fmadd_pd(tp6, a, ans6[0]);
							ans6[1] = _mm512_fmadd_pd(tp6, b, ans6[1]);
							ans7[0] = _mm512_fmadd_pd(tp7, a, ans7[0]);
							ans7[1] = _mm512_fmadd_pd(tp7, b, ans7[1]);
						}
						for (unsigned long long c3(0); c3 < warp; ++c3)
						{
							bData[c1 + c3] = ans0[c3];
							bData[aWidth512d + c1 + c3] = ans1[c3];
							bData[2 * aWidth512d + c1 + c3] = ans2[c3];
							bData[3 * aWidth512d + c1 + c3] = ans3[c3];
							bData[4 * aWidth512d + c1 + c3] = ans4[c3];
							bData[5 * aWidth512d + c1 + c3] = ans5[c3];
							bData[6 * aWidth512d + c1 + c3] = ans6[c3];
							bData[7 * aWidth512d + c1 + c3] = ans7[c3];
						}
					}
				}
				return 0;
			};

			unsigned long long threadHeight;
			unsigned long long rowBeginning;
			unsigned long long leftHeight;
			if (source->height / threadNum >= 8)
			{
				threadHeight = (source->height / threadNum) & (-8);
			}
			else
			{
				threadHeight = 8;
				threadNum = source->height / 8;
			}
			rowBeginning = threadHeight * threadNum;
			leftHeight = source->height - rowBeginning;

			HANDLE* threads(nullptr);

			void** paras(nullptr);
			if (threadNum)
			{
				threads = (HANDLE*)::malloc(threadNum * sizeof(HANDLE));
				paras = (void**)::malloc(threadNum * 9 * sizeof(void*));
				for (unsigned long long c0(0); c0 < threadNum; ++c0)
				{
					paras[c0 * 9] = (void*)(source->data + threadHeight * source->width8d * c0);			//A beginning
					paras[c0 * 9 + 1] = (void*)(source->data + threadHeight * source->width8d * (c0 + 1));	//A ending
					paras[c0 * 9 + 2] = (void*)(a.data);													//B
					paras[c0 * 9 + 3] = (void*)(b.data + threadHeight * a.width8d * c0);					//C beginning
					paras[c0 * 9 + 4] = (void*)(source->width8d);											//width4d
					paras[c0 * 9 + 5] = (void*)(aWidthWarpFloor);											//aWidthWarpFloor
					paras[c0 * 9 + 6] = (void*)(minDim);													//minDim
					paras[c0 * 9 + 7] = (void*)(aWidth512d);												//aWidth256d
					paras[c0 * 9 + 8] = (void*)(warpLeft);													//warpLeft
					threads[c0] = (HANDLE)_beginthreadex(nullptr, 0, lambda, paras + c0 * 9, 0, nullptr);
				}
				DWORD rc = WaitForMultipleObjects(threadNum, threads, true, INFINITE);
				for (unsigned long long c0(0); c0 < threadNum; ++c0)
					CloseHandle(threads[c0]);
			}
			/*if (leftHeight > 1)
			{
				threadHeight = 2;
				threadNum = leftHeight / 2;
				for (unsigned long long c0(0); c0 < threadNum; ++c0)
				{
					paras[c0 * 9] = (void*)(source->data + source->width8d * (threadHeight * c0 + rowBeginning));			//A beginning
					paras[c0 * 9 + 1] = (void*)(source->data + source->width8d * (threadHeight * (c0 + 1) + rowBeginning));	//A ending
					paras[c0 * 9 + 2] = (void*)(a.data);																	//B
					paras[c0 * 9 + 3] = (void*)(b.data + a.width8d * (threadHeight * c0 + rowBeginning));					//C beginning
					paras[c0 * 9 + 4] = (void*)(source->width8d);															//width4d
					paras[c0 * 9 + 5] = (void*)(aWidthWarpFloor);															//aWidthWarpFloor
					paras[c0 * 9 + 6] = (void*)(minDim);																	//minDim
					paras[c0 * 9 + 7] = (void*)(aWidth512d);																//aWidth256d
					paras[c0 * 9 + 8] = (void*)(warpLeft);																	//warpLeft
					threads[c0] = (HANDLE)_beginthreadex(nullptr, 0, lambda, paras + c0 * 9, 0, nullptr);
				}
				DWORD rc = WaitForMultipleObjects(threadNum, threads, true, INFINITE);
				for (unsigned long long c0(0); c0 < threadNum; ++c0)
					CloseHandle(threads[c0]);
				rowBeginning += threadHeight * threadNum;
			}*/
			free(threads);
			free(paras);
		}
		return b;
	}
}

int main()
{
	using namespace BLAS;

	Timer timer;

	std::mt19937 mt(time(nullptr));
	std::uniform_real_distribution<double> rd(-1.0, 1.0);
	std::uniform_int_distribution<unsigned long long> rduint(1, 10);

	remove("./matA.txt");
	remove("./matB.txt");
	remove("./matC.txt");

	mat512 mA(1024, 1024, false);
	mat512 mB(1024, 1024, false);
	mat512 mC(1024, 1024, false);

	randomMat(mA, mt, rd);
	randomMat(mB, mt, rd);

	timer.begin();
	//mA(mB, mC);
	for (unsigned int c0(0); c0 < 10; ++c0)
		matMultMT(mA, mB, mC);
	timer.end();
	timer.print("11800H AVX512 Multi Thread: ");

	//mA.printToTableTxt("./matA.txt");
	//mB.printToTableTxt("./matB.txt");
	//mC.printToTableTxt("./matC.txt");
}