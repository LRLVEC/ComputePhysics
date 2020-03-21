#pragma once
#include <cmath>
#include <malloc.h>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <cstdio>
#include <immintrin.h>

//if you can change it, than change it
namespace BLAS
{
	static constexpr double Pi = 3.14159265358979323846264338327950288L;
	static constexpr double E = 2.71828182845904523536028747135266250L;

	enum class Type
	{
		Native = 0,
		Parasitic = 1,
		Non32Aligened = 2,//must be Parasitic as well...
	};
	enum class MatType
	{
		NormalMat,
		SquareMat,
		DiagonalMat,
		SymmetricMat,
		LMat,
		UMat,
		//for band mat, width == height
		BandMat,
		LBandMat,
		UBandMat,
	};

	//note: if a vec is Non32Aligened, then don't do any +-*/ operation to 
	// other vecs unless other ones has the same offset from 32-byte-aligened as this one!

	inline unsigned int ceiling4(unsigned int length)
	{
		return (((length - 1) >> 2) + 1) << 2;
	}
	inline unsigned int floor(unsigned int length)
	{
		return (length >> 2) << 2;
	}
	inline unsigned long long ceiling4(unsigned int width, unsigned int height)
	{
		return ((((unsigned long long(width) - 1) >> 2) + 1) << 2) * height;
	}
	inline size_t ceiling256dSize(unsigned int length)
	{
		return (((size_t(length) - 1) >> 2) + 1) << 5;
	}
	inline size_t ceiling256dSize(unsigned int width, unsigned int height)
	{
		return ((((size_t(width) - 1) >> 2) + 1) << 5) * height;
	}
	inline double* malloc64d(unsigned int length)
	{
		return (double*)_mm_malloc(length * sizeof(double), 32);
	}
	inline double* malloc256d(unsigned int length)
	{
		return (double*)_mm_malloc(ceiling256dSize(length), 32);
	}
	inline double* malloc256d(unsigned int width, unsigned int height)
	{
		return (double*)_mm_malloc(ceiling256dSize(width) * height, 32);
	}
	inline void* memcpy64d(void* dst, void const* src, unsigned int length)
	{
		return ::memcpy(dst, src, length * sizeof(double));
	}
	inline void* memcpy256d(void* dst, void const* src, unsigned int length)
	{
		return ::memcpy(dst, src, ceiling256dSize(length));
	}
	inline void* memcpy256d(void* dst, void const* src, unsigned int width, unsigned int height)
	{
		return ::memcpy(dst, src, ceiling256dSize(width, height));
	}
	inline void* memset64d(void* dst, int val, unsigned int length)
	{
		return ::memset(dst, val, length * sizeof(double));
	}
	inline void* memset256d(void* dst, int val, unsigned int length)
	{
		return ::memset(dst, val, ceiling256dSize(length));
	}
	inline void* memset256d(void* dst, int val, unsigned int width, unsigned int height)
	{
		return ::memset(dst, val, ceiling256dSize(width, height));
	}

	inline unsigned int getPtrOffset64d(double* ptr)
	{
		return (unsigned int(ptr) >> 3) & 3;
	}
	inline double* getPtr256d(double* ptr)
	{
		return (double*)(unsigned long long(ptr) & -32);
	}

	struct mat;
	struct vec
	{
		double* data;
		unsigned int dim;
		unsigned int begining;
		Type type;

		vec() :data(nullptr), dim(0), begining(0), type(Type::Native) {}
		vec(unsigned int _length, bool _clear = true)
			:
			data(_length ? malloc256d(_length) : nullptr),
			dim(_length),
			begining(0),
			type(Type::Native)
		{
			if (_clear && data)memset256d(data, 0, _length);
		}
		vec(vec const& a)
			:
			data(a.dim ? malloc256d(a.dim) : nullptr),
			dim(a.dim),
			begining(0),
			type(Type::Native)
		{
			if (dim)memcpy64d(data, a.data + a.begining, dim);
		}
		vec(vec&& a) :data(nullptr), dim(0), begining(0), type(Type::Native)
		{
			if (a.type == Type::Native)
			{
				data = a.data;
				dim = a.dim;
				a.data = nullptr;
				a.dim = 0;
			}
			else
			{
				if (a.dim)
				{
					data = malloc256d(a.dim);
					dim = a.dim;
					memcpy64d(data, a.data + a.begining, dim);
				}
			}
		}
		vec(double* _data, unsigned int _length, Type _type)
			:data(getPtr256d(_data)), dim(_length), begining(getPtrOffset64d(_data)), type(_type) {}
		vec(std::initializer_list<double>const& a)
			:
			data(a.size() ? malloc256d(a.size()) : nullptr),
			dim(a.size()),
			type(Type::Native)
		{
			if (dim)memcpy64d(data, a.begin(), dim);
		}
		~vec()
		{
			if (type == Type::Native)_mm_free(data);
			data = nullptr;
			dim = 0;
		}
		template<class T>inline double& operator[](T a)
		{
			return data[a + begining];
		}
		void realloc(unsigned int _dim, bool _clear = true)
		{
			if (type == Type::Native && _dim)
			{
				if (dim)
				{
					double* r((double*)malloc256d(_dim));
					memcpy64d(r, data, dim);
					if (_clear && _dim > dim)
						memset64d(r + dim, 0, _dim - dim);
					_mm_free(data);
					dim = _dim;
				}
				else
				{
					data = (double*)malloc256d(_dim);
					if (_clear)memset64d(data, 0, _dim);
					dim = _dim;
				}
			}
		}
		void reconstruct(unsigned int _dim, bool _clear = true)
		{
			if (type == Type::Native)
			{
				if (dim)_mm_free(data);
				if (_dim)
				{
					data = (double*)malloc256d(_dim);
					if (_clear)memset64d(data, 0, _dim);
					dim = _dim;
				}
			}
		}
		void clearTail()
		{
			if (dim & 3)
				for (unsigned int c0(dim); c0 < ceiling4(dim); ++c0)
					data[c0] = 0;
		}
		//moveTo
		//vec& moveTo(vec& a)
		//{
		//	if (type == Type::Native)
		//	{
		//		a.~vec();
		//		a.data = data;
		//		a.dim = dim;
		//		a.type = type;
		//		data = nullptr;
		//		dim = 0;
		//	}
		//	else a = *this;
		//	return a;
		//}

		//= += -= *= /=
		vec& operator =(vec&& a)
		{
			if (a.dim)
			{
				if (type == Type::Native)
				{
					if (a.type == Type::Native)
					{
						_mm_free(data);
						data = a.data;
						dim = a.dim;
						a.data = nullptr;
						a.dim = 0;
					}
					else
					{
						if (dim != a.dim)
						{
							_mm_free(data);
							data = malloc256d(a.dim);
							dim = a.dim;
						}
						memcpy64d(data, a.data + a.begining, dim);
					}
				}
				else
				{
					unsigned minDim(dim > a.dim ? a.dim : dim);
					memcpy64d(data, a.data + a.begining, minDim);
				}
			}
			else if (type == Type::Native)
			{
				_mm_free(data);
				data = nullptr;
				dim = 0;
			}
			return *this;
		}
		vec& operator =(vec const& a)
		{
			if (a.dim)
			{
				if (type == Type::Native)
				{
					if (dim != a.dim)
					{
						_mm_free(data);
						data = malloc256d(a.dim);
						dim = a.dim;
					}
					memcpy64d(data, a.data + a.begining, dim);
				}
				else
				{
					unsigned minDim(dim > a.dim ? a.dim : dim);
					memcpy64d(data, a.data + a.begining, minDim);
				}
			}
			else if (type == Type::Native)
			{
				_mm_free(data);
				data = nullptr;
				dim = 0;
			}
			return *this;
		}
		vec& operator+=(vec const& a)
		{
			if (a.dim && dim)
			{
				unsigned int minDim(dim > a.dim ? a.dim : dim);
				for (unsigned int c0(0); c0 < minDim; ++c0)
					data[c0] += a.data[c0];
				/*unsigned int minDim4(minDim >> 2);
				__m256d* aData((__m256d*)data);
				__m256d* bData((__m256d*)a.data);
				unsigned int c0(0);
				for (; c0 < minDim4; ++c0)
					aData[c0] = _mm256_add_pd(aData[c0], bData[c0]);
				if ((c0 << 2) < minDim)
					for (unsigned int c1(c0 << 2); c1 < minDim; ++c1)
						data[c1] += a.data[c1];*/
			}
			return *this;
		}
		vec& operator-=(vec const& a)
		{
			if (a.dim && dim)
			{
				unsigned int minDim(dim > a.dim ? a.dim : dim);
				for (unsigned int c0(0); c0 < minDim; ++c0)
					data[c0] -= a.data[c0];
			}
			return *this;
		}
		vec& operator*=(vec const& a)
		{
			if (a.dim && dim)
			{
				unsigned int minDim(dim > a.dim ? a.dim : dim);
				for (unsigned int c0(0); c0 < minDim; ++c0)
					data[c0] *= a.data[c0];
			}
			return *this;
		}
		vec& operator/=(vec const& a)
		{
			if (a.dim && dim)
			{
				unsigned int minDim(dim > a.dim ? a.dim : dim);
				for (unsigned int c0(0); c0 < minDim; ++c0)
					data[c0] /= a.data[c0];
			}
			return *this;
		}
		vec& operator =(double a)
		{
			if (dim)
			{
				for (unsigned int c0(0); c0 < dim; ++c0)
					data[c0] = a;
			}
			return *this;
		}
		vec& operator+=(double a)
		{
			if (dim)
			{
				for (unsigned int c0(0); c0 < dim; ++c0)
					data[c0] += a;
			}
			return *this;
		}
		vec& operator-=(double a)
		{
			if (dim)
			{
				for (unsigned int c0(0); c0 < dim; ++c0)
					data[c0] -= a;
			}
			return *this;
		}
		vec& operator*=(double a)
		{
			if (dim)
			{
				for (unsigned int c0(0); c0 < dim; ++c0)
					data[c0] *= a;
			}
			return *this;
		}
		vec& operator/=(double a)
		{
			if (dim)
			{
				for (unsigned int c0(0); c0 < dim; ++c0)
					data[c0] /= a;
			}
			return *this;
		}
		//vecA = a * vecB + vecA
		vec& fmadd(double a, vec const& b)
		{
			if (b.dim && dim)
			{
				unsigned int minDim(dim > b.dim ? b.dim : dim);
				/*for (unsigned int c0(0); c0 < minDim; ++c0)
					data[c0] += a * b.data[c0];*/
				unsigned int minDim4(minDim >> 2);
				__m256d* aData((__m256d*)data);
				__m256d* bData((__m256d*)b.data);
				unsigned int c0(0);
				__m256d tp = _mm256_broadcast_sd(&a);
				for (; c0 < minDim4; ++c0)
					aData[c0] = _mm256_fmadd_pd(tp, bData[c0], aData[c0]);
				if ((c0 << 2) < minDim)
					for (unsigned int c1(c0 << 2); c1 < minDim; ++c1)
						data[c1] += a * b.data[c1];
			}
			return *this;
		}
		//+-*/
		vec operator+(vec const& a)const
		{
			if (a.dim && dim)
			{
				unsigned int l(dim > a.dim ? a.dim : dim);
				double* d(malloc256d(l));
				for (unsigned int c0(0); c0 < l; ++c0)
					d[c0] = data[c0] + a.data[c0];
				return vec(d, l, Type::Native);
			}
			else return vec();
		}
		vec operator-(vec const& a)const
		{
			if (a.dim && dim)
			{
				unsigned int l(dim > a.dim ? a.dim : dim);
				double* d(malloc256d(l));
				for (unsigned int c0(0); c0 < l; ++c0)
					d[c0] = data[c0] - a.data[c0];
				return vec(d, l, Type::Native);
			}
			else return vec();
		}
		vec operator*(vec const& a)const
		{
			if (a.dim && dim)
			{
				unsigned int l(dim > a.dim ? a.dim : dim);
				double* d(malloc256d(l));
				for (unsigned int c0(0); c0 < l; ++c0)
					d[c0] = data[c0] * a.data[c0];
				return vec(d, l, Type::Native);
			}
			else return vec();
		}
		vec operator/(vec const& a)const
		{
			if (a.dim && dim)
			{
				unsigned int l(dim > a.dim ? a.dim : dim);
				double* d(malloc256d(l));
				for (unsigned int c0(0); c0 < l; ++c0)
					d[c0] = data[c0] / a.data[c0];
				return vec(d, l, Type::Native);
			}
			else return vec();
		}
		vec operator+(double a)const
		{
			if (dim)
			{
				unsigned int l(dim);
				double* d(malloc256d(l));
				for (unsigned int c0(0); c0 < l; ++c0)
					d[c0] = data[c0] + a;
				return vec(d, l, Type::Native);
			}
			else return vec();
		}
		vec operator-(double a)const
		{
			if (dim)
			{
				unsigned int l(dim);
				double* d(malloc256d(l));
				for (unsigned int c0(0); c0 < l; ++c0)
					d[c0] = data[c0] - a;
				return vec(d, l, Type::Native);
			}
			else return vec();
		}
		vec operator*(double a)const
		{
			if (dim)
			{
				unsigned int l(dim);
				double* d(malloc256d(l));
				for (unsigned int c0(0); c0 < l; ++c0)
					d[c0] = data[c0] * a;
				return vec(d, l, Type::Native);
			}
			else return vec();
		}
		vec operator/(double a)const
		{
			if (dim)
			{
				unsigned int l(dim);
				double* d(malloc256d(l));
				for (unsigned int c0(0); c0 < l; ++c0)
					d[c0] = data[c0] / a;
				return vec(d, l, Type::Native);
			}
			else return vec();
		}
		//dot
		double operator,(vec const& a)const
		{
			if (a.dim && dim)
			{
				unsigned int e0(dim + begining);
				unsigned int e1(a.dim + a.begining);
				unsigned int minE(e0 >= e1 ? e1 : e0);
				double s(0);
				/*for (unsigned int c0(0); c0 < minDim; ++c0)
					s += data[c0] * a.data[c0];*/
				unsigned int maxB(begining >= a.begining ? begining : a.begining);
				unsigned int minDim4(minE >> 2);
				__m256d* aData((__m256d*)data);
				__m256d* bData((__m256d*)a.data);
				unsigned int c0(0);
				__m256d tp = { 0 };
				if (begining)
				{
					tp = _mm256_fmadd_pd(aData[c0], bData[c0], tp);
					for (unsigned int c1(0); c1 < maxB; ++c1)
						tp.m256d_f64[c1] = 0;
					++c0;
				}
				for (unsigned int c1(minE); c1 < 4; ++c1)
					tp.m256d_f64[c1] = 0;
				for (; c0 < minDim4; ++c0)
					tp = _mm256_fmadd_pd(aData[c0], bData[c0], tp);
				for (unsigned int c1(c0 << 2); c1 < minE; ++c1)
					s += data[c1] * a.data[c1];
				s += tp.m256d_f64[0];
				s += tp.m256d_f64[1];
				s += tp.m256d_f64[2];
				s += tp.m256d_f64[3];
				return s;
			}
			else return 0;
		}
		//normalize
		vec& normalize()
		{
			double s(this->norm2());
			if (s)(*this) *= 1 / s;
			return *this;
		}
		//non-in-situ mult mat
		vec operator()(mat const& a)const;
		vec& operator()(mat const& a, vec& b)const;

		//norm
		double norm1()const
		{
			if (dim)
			{
				double s(0);
				/*for (unsigned int c0(0); c0 < dim; ++c0)
					s += abs(data[c0]);*/
				unsigned long long a((1llu << 63) - 1llu);
				__m256d gg = _mm256_broadcast_sd((double*)&a);
				unsigned int finalDim(dim + begining);
				unsigned int dim4(finalDim >> 2);
				__m256d* aData((__m256d*)data);
				unsigned int c0(0);
				__m256d tp = { 0 };
				if (begining)
				{
					tp = _mm256_add_pd(tp, _mm256_and_pd(gg, aData[c0++]));
					for (unsigned int c1(0); c1 < begining; ++c1)
						tp.m256d_f64[c1] = 0;
				}
				for (unsigned int c1(finalDim); c1 < 4; ++c1)
					tp.m256d_f64[c1] = 0;
				for (; c0 < dim4; ++c0)
					tp = _mm256_add_pd(tp, _mm256_and_pd(gg, aData[c0]));
				for (unsigned int c1(0); c1 < 4; ++c1)
					s += tp.m256d_f64[c1];
				for (unsigned int c1(c0 << 2); c1 < finalDim; ++c1)
					s += abs(data[c1]);
				return s;
			}
			return 0;
		}
		double norm2Square()const
		{
			if (dim)
			{
				double s(0);
				unsigned int finalDim(dim + begining);
				unsigned int dim4(finalDim >> 2);
				__m256d* aData((__m256d*)data);
				unsigned int c0(0);
				__m256d tp = { 0 };
				if (begining)
				{
					__m256d gg = aData[c0++];
					tp = _mm256_fmadd_pd(gg, gg, tp);
					for (unsigned int c1(0); c1 < begining; ++c1)
						tp.m256d_f64[c1] = 0;
				}
				for (unsigned int c1(finalDim); c1 < 4; ++c1)
					tp.m256d_f64[c1] = 0;
				for (; c0 < dim4; ++c0)
				{
					__m256d gg = aData[c0];
					tp = _mm256_fmadd_pd(gg, gg, tp);
				}
				for (unsigned int c1(c0 << 2); c1 < finalDim; ++c1)
					s += data[c1] * data[c1];
				s += tp.m256d_f64[0];
				s += tp.m256d_f64[1];
				s += tp.m256d_f64[2];
				s += tp.m256d_f64[3];
				return s;
			}
			else return 0;
		}
		double norm2()const
		{
			return sqrt(norm2Square());
		}
		double normInf()const
		{
			if (dim)
			{
				double s(0);
				/*for (unsigned int c0(0); c0 < dim; ++c0)
					if (s < abs(data[c0]))s = abs(data[c0]);*/
				unsigned long long a((1llu << 63) - 1llu);
				double g(*(double*)&a);
				__m256d gg = _mm256_broadcast_sd((double*)&a);
				unsigned int finalDim(dim + begining);
				unsigned int dim4(finalDim >> 2);
				__m256d* aData((__m256d*)data);
				unsigned int c0(0);
				__m256d tp = { 0 };
				if (begining)
				{
					tp = _mm256_and_pd(gg, aData[c0++]);
					for (unsigned int c1(0); c1 < begining; ++c1)
						tp.m256d_f64[c1] = 0;
				}
				for (unsigned int c1(finalDim); c1 < 4; ++c1)
					tp.m256d_f64[c1] = 0;
				for (; c0 < dim4; ++c0)
					tp = _mm256_max_pd(tp, _mm256_and_pd(gg, aData[c0]));
				for (unsigned int c1(0); c1 < 4; ++c1)
					if (s < tp.m256d_f64[c1])s = tp.m256d_f64[c1];
				for (unsigned int c1(c0 << 2); c1 < finalDim; ++c1)
					if (s < abs(data[c1]))s = abs(data[c1]);
				return s;
			}
			return 0;
		}
		double normP(double p)const
		{
			if (dim && p)
			{
				double s(0);
				/*for (unsigned int c0(0); c0 < dim; ++c0)
					s += pow(abs(data[c0]), p);*/
				unsigned long long a((1llu << 63) - 1llu);
				double g(*(double*)&a);
				__m256d gg = _mm256_broadcast_sd((double*)&a);
				__m256d pp = _mm256_broadcast_sd(&p);
				unsigned int finalDim(dim + begining);
				unsigned int dim4(finalDim >> 2);
				__m256d* aData((__m256d*)data);
				unsigned int c0(0);
				__m256d tp = { 0 };
				if (begining)
				{
					tp = _mm256_pow_pd(_mm256_and_pd(gg, aData[c0++]), pp);
					for (unsigned int c1(0); c1 < begining; ++c1)
						tp.m256d_f64[c1] = 0;
				}
				for (unsigned int c1(finalDim); c1 < 4; ++c1)
					tp.m256d_f64[c1] = 0;
				for (; c0 < dim4; ++c0)
					tp = _mm256_add_pd(
						_mm256_pow_pd(_mm256_and_pd(gg, aData[c0]), pp), tp);
				for (unsigned int c1(0); c1 < 4; ++c1)
					s += tp.m256d_f64[c1];
				for (unsigned int c1(c0 << 2); c1 < finalDim; ++c1)
					s += pow(abs(data[c1]), p);
				return pow(s, 1 / p);
			}
			return 0;
		}

		void print()const
		{
			::printf("[");
			unsigned int finalDim(dim + begining);
			if (dim)
			{
				for (unsigned int c0(begining); c0 < finalDim - 1; ++c0)
					::printf("%.4f, ", data[c0]);
				::printf("%.4f", data[finalDim - 1]);
			}
			::printf("]\n");
		}
		void printInfo()const
		{
			char const* str;
			switch (type)
			{
			case Type::Native:str = "Native"; break;
			case Type::Parasitic:str = "Parastic"; break;
			case Type::Non32Aligened:str = "Non32Aligened"; break;
			}
			::printf("{dim: %u, begining: %u, type: %s}\n", dim, begining, str);
		}
		void printToTxt(char const* name)const
		{
			//in the form of Mathematica matrix (paste)
			if (data)
			{
				FILE* temp(::fopen(name, "w+"));
				::fprintf(temp, "{");
				unsigned int finalDim(dim + begining);
				for (unsigned int c0(begining); c0 < finalDim - 1; ++c0)
					::fprintf(temp, "{%.14e}, ", data[c0]);
				::fprintf(temp, "{%.14e}}\n", data[finalDim - 1]);
				::fclose(temp);
			}
		}
		void printToTableTxt(char const* name, bool _inRow)const
		{
			//in the form of Mathematica matrix (table)
			if (data)
			{
				FILE* temp(::fopen(name, "w+"));
				for (unsigned int c0(begining); c0 < dim + begining; ++c0)
					::fprintf(temp, _inRow ? "%.14e " : "%.14e\n", data[c0]);
				::fclose(temp);
			}
		}
	};
	struct mat
	{
		double* data;
		union
		{
			unsigned int width;
			unsigned int halfBandWidth;//if is BandMat
		};
		unsigned int height;
		unsigned int width4d;
		Type type;
		MatType matType;

		mat() :data(nullptr), width(0), height(0), width4d(ceiling4(width)),
			type(Type::Native), matType(MatType::NormalMat) {}
		mat(unsigned int _width, unsigned int _height, bool _clear = true)
			:
			data((_width&& _height) ? malloc256d(_width, _height) : nullptr),
			width(data ? _width : 0),
			height(data ? _height : 0),
			width4d(data ? ceiling4(_width) : 0),
			type(Type::Native),
			matType(MatType::NormalMat)
		{
			if (_clear && data)
				memset256d(data, 0, width, height);
		}
		mat(unsigned int _halfBandWidth, unsigned int _height, MatType _type, bool _clear = true)
			:
			data(nullptr),
			halfBandWidth(0),
			height(0),
			width4d(0),
			type(Type::Native),
			matType(_type)
		{
			if (_type == MatType::NormalMat && width && height)
			{
				data = malloc256d(_halfBandWidth, _height);
				width = _halfBandWidth;
				height = _height;
				width4d = ceiling4(_halfBandWidth);
			}
			else if ((_type == MatType::LBandMat || _type == MatType::UBandMat) && _height)
			{
				data = malloc256d(_halfBandWidth + 4, _height);
				halfBandWidth = _halfBandWidth;
				height = _height;
				width4d = ceiling4(_halfBandWidth + 4);
			}
			if (_clear && data)
				memset64d(data, 0, width4d * height);
		}
		mat(mat const& a)
			:
			data((a.width&& a.height) ? malloc256d(a.width, a.height) : nullptr),
			width(a.width),
			height(a.height),
			width4d(a.width4d),
			type(Type::Native),
			matType(a.matType)
		{
			if (data)
				memcpy256d(data, a.data, width, height);
		}
		mat(mat&& a) :data(nullptr), width(0), height(0), type(Type::Native), matType(MatType::NormalMat)
		{
			if (a.data)
			{
				if (a.type == Type::Native)
				{
					data = a.data;
					width = a.width;
					height = a.height;
					matType = a.matType;
					a.data = nullptr;
					a.width = a.height = a.width4d = 0;
				}
				else
				{
					data = malloc256d(a.width, a.height);
					width = a.width;
					height = a.height;
					memcpy256d(data, a.data, a.width, a.height);
				}
				width4d = ceiling4(width);
			}
		}
		mat(double* _data, unsigned int _width, unsigned int _height, Type _type, MatType _matType)
			:
			data(_data), width(_width), height(_height), width4d(ceiling4(_width)),
			type(_type), matType(_matType)
		{
		}
		mat(std::initializer_list<std::initializer_list<double>>const& a)
			:
			data(nullptr), width(0), height(0), type(Type::Native), matType(MatType::NormalMat)
		{
			if (a.size())
			{
				unsigned int h(a.size());
				unsigned int w(0);
				for (unsigned int c0(0); c0 < h; ++c0)
					if ((a.begin() + c0)->size() > w)
						w = (a.begin() + c0)->size();
				if (w)
				{
					unsigned long long l(ceiling4(w, h));
					data = malloc64d(l);
					memset64d(data, 0, l);
					width = w;
					height = h;
					width4d = ceiling4(width);
					type = Type::Native;
					matType = MatType::NormalMat;
					for (unsigned int c0(0); c0 < height; ++c0)
						memcpy64d(data + width4d * c0, (a.begin() + c0)->begin(),
							(a.begin() + c0)->size());
				}
			}
		}
		~mat()
		{
			if (type == Type::Native)_mm_free(data);
			data = nullptr;
			width = height = 0;
		}
		template<class T>inline double& operator[](T a)
		{
			return data[a];
		}
		inline double& operator() (unsigned int a, unsigned int b)
		{
			return data[unsigned long long(a) * width4d + b];
		}
		inline double  LBandEle(unsigned int a, unsigned int b)const
		{
			if (a <= halfBandWidth)return data[a * width4d + b];
			else return data[a * width4d + b - (int(a - halfBandWidth) / 4) * 4];
		}
		inline double& LBandEleRef(unsigned int a, unsigned int b)
		{
			if (a <= halfBandWidth)return data[a * width4d + b];
			else return data[a * width4d + b - (int(a - halfBandWidth) / 4) * 4];
		}
		inline double  UBandEle(unsigned int a, unsigned int b)const
		{
			return data[a * width4d + b - (a & -4)];
		}
		inline double& UBandEleRef(unsigned int a, unsigned int b)
		{
			return data[a * width4d + b - (a & -4)];
		}
		inline unsigned int LBandBeginOffset(unsigned int a)const
		{
			if (a <= halfBandWidth)return 0;
			else return a - (int(a - halfBandWidth) / 4) * 4 - halfBandWidth;
		}
		inline unsigned int UBandBeginOffset(unsigned int a)const
		{
			return a % 4;
		}
		inline vec getLBandRow(unsigned int a)
		{
			return vec(data + a * width4d + LBandBeginOffset(a),
				a <= halfBandWidth ? a + 1 : halfBandWidth + 1, Type::Non32Aligened);
		}
		inline vec getUBandRow(unsigned int a)
		{
			return vec(data + a * width4d + a % 4,
				height - a <= halfBandWidth ? height - a : halfBandWidth + 1, Type::Non32Aligened);
		}
		inline vec const getLBandRow(unsigned int a)const
		{
			return vec(data + a * width4d + LBandBeginOffset(a),
				a <= halfBandWidth ? a + 1 : halfBandWidth + 1, Type::Non32Aligened);
		}
		inline vec const getUBandRow(unsigned int a)const
		{
			return vec(data + a * width4d + a % 4,
				height - a <= halfBandWidth ? height - a : halfBandWidth + 1, Type::Non32Aligened);
		}
		inline vec const getLBandRowL(unsigned int a)const
		{
			if (a)
				return vec(data + a * width4d + LBandBeginOffset(a),
					a <= halfBandWidth ? a : halfBandWidth, Type::Non32Aligened);
			return vec();
		}
		inline vec const getUBandRowU(unsigned int a)const
		{
			if (height - a)
				return vec(data + a * width4d + a % 4 + 1,
					height - a <= halfBandWidth ? height - a - 1 : halfBandWidth, Type::Non32Aligened);
			return vec();
		}
		void clear()
		{
			if (data)memset64d(data, 0, height * width4d);
		}
		void reconstruct(unsigned int _width, unsigned int _height, bool _clear = true)
		{
			if (type == Type::Native)
			{
				if (data)_mm_free(data);
				unsigned int s(ceiling256dSize(_width, _height));
				if (s)
				{
					data = (double*)malloc64d(s);
					if (_clear)memset64d(data, 0, s);
					width = _width;
					height = _height;
				}
			}
		}
		//= += -= *= /=
		mat& operator =(mat&& a)
		{
			if (type == Type::Native)
			{
				if (a.type == Type::Native)
				{
					_mm_free(data);
					data = a.data;
					width = a.width;
					height = a.height;
					width4d = a.width4d;
					matType = a.matType;
					a.data = nullptr;
					a.width = a.height = a.width4d = 0;
				}
				else
				{
					if (unsigned long long(a.width4d) * a.height !=
						unsigned long long(width4d) * height)
					{
						_mm_free(data);
						data = malloc256d(a.width, a.height);
					}
					width = a.width;
					height = a.height;
					width4d = a.width4d;
					memcpy256d(data, a.data, width, height);
				}
			}
			else
			{
				unsigned int minWidth(width > a.width ? a.width : width);
				unsigned int minHeight(height > a.height ? a.height : height);
				for (unsigned int c0(0); c0 < minHeight; ++c0)
					memcpy64d(data + width4d * c0, a.data + a.width4d * c0, minWidth);
			}
			return *this;
		}
		mat& operator =(mat const& a)
		{
			if (type == Type::Native)
			{
				if (unsigned long long(a.width4d) * a.height !=
					unsigned long long(width4d) * height)
				{
					_mm_free(data);
					data = malloc256d(a.width, a.height);
				}
				width = a.width;
				height = a.height;
				width4d = a.width4d;
				memcpy256d(data, a.data, width, height);
			}
			else
			{
				unsigned int minWidth(width > a.width ? a.width : width);
				unsigned int minHeight(height > a.height ? a.height : height);
				for (unsigned int c0(0); c0 < minHeight; ++c0)
					memcpy64d(data + width4d * c0, a.data + a.width4d * c0, minWidth);
			}
			return *this;
		}
		mat& operator+=(mat const& a)
		{
			if (data && a.data)
			{
				unsigned int minW(a.width > width ? width : a.width);
				unsigned int minH(a.height > height ? height : a.height);
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						data[c0 * width4d + c1] += a.data[c0 * a.width4d + c1];
			}
			return *this;
		}
		mat& operator-=(mat const& a)
		{
			if (data && a.data)
			{
				unsigned int minW(a.width > width ? width : a.width);
				unsigned int minH(a.height > height ? height : a.height);
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						data[c0 * width4d + c1] -= a.data[c0 * a.width4d + c1];
			}
			return *this;
		}
		mat& operator*=(mat const& a)
		{
			if (data && a.data)
			{
				unsigned int minW(a.width > width ? width : a.width);
				unsigned int minH(a.height > height ? height : a.height);
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						data[c0 * width4d + c1] *= a.data[c0 * a.width4d + c1];
			}
			return *this;
		}
		mat& operator/=(mat const& a)
		{
			if (data && a.data)
			{
				unsigned int minW(a.width > width ? width : a.width);
				unsigned int minH(a.height > height ? height : a.height);
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						data[c0 * width4d + c1] /= a.data[c0 * a.width4d + c1];
			}
			return *this;
		}
		mat& operator =(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width4d + c1] = a;
			}
			return *this;
		}
		mat& operator+=(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width4d + c1] += a;
			}
			return *this;
		}
		mat& operator-=(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width4d + c1] -= a;
			}
			return *this;
		}
		mat& operator*=(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width4d + c1] *= a;
			}
			return *this;
		}
		mat& operator/=(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width4d + c1] /= a;
			}
			return *this;
		}
		//+-*/ to do: mat types' operations?
		mat operator+(mat const& a)const
		{
			if (a.data && data)
			{
				unsigned int minW(a.width > width ? width : a.width);
				unsigned int minH(a.height > height ? height : a.height);
				unsigned int minW4d(ceiling4(minW));
				double* d(malloc256d(minW, minH));
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						d[c0 * minW4d + c1] = data[c0 * width4d + c1] + a.data[c0 * a.width4d + c1];
				return mat(d, minW, minH, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator-(mat const& a)const
		{
			if (a.data && data)
			{
				unsigned int minW(a.width > width ? width : a.width);
				unsigned int minH(a.height > height ? height : a.height);
				unsigned int minW4d(ceiling4(minW));
				double* d(malloc256d(minW, minH));
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						d[c0 * minW4d + c1] = data[c0 * width4d + c1] - a.data[c0 * a.width4d + c1];
				return mat(d, minW, minH, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator*(mat const& a)const
		{
			if (a.data && data)
			{
				unsigned int minW(a.width > width ? width : a.width);
				unsigned int minH(a.height > height ? height : a.height);
				unsigned int minW4d(ceiling4(minW));
				double* d(malloc256d(minW, minH));
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						d[c0 * minW4d + c1] = data[c0 * width4d + c1] * a.data[c0 * a.width4d + c1];
				return mat(d, minW, minH, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator/(mat const& a)const
		{
			if (a.data && data)
			{
				unsigned int minW(a.width > width ? width : a.width);
				unsigned int minH(a.height > height ? height : a.height);
				unsigned int minW4d(ceiling4(minW));
				double* d(malloc256d(minW, minH));
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						d[c0 * minW4d + c1] = data[c0 * width4d + c1] / a.data[c0 * a.width4d + c1];
				return mat(d, minW, minH, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator+(double a)const
		{
			if (data)
			{
				double* d(malloc256d(width, height));
				memcpy256d(d, data, width, height);
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						d[c0 * width4d + c1] += a;
				return mat(d, width, height, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator-(double a)const
		{
			if (data)
			{
				double* d(malloc256d(width, height));
				memcpy256d(d, data, width, height);
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						d[c0 * width4d + c1] -= a;
				return mat(d, width, height, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator*(double a)const
		{
			if (data)
			{
				double* d(malloc256d(width, height));
				memcpy256d(d, data, width, height);
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						d[c0 * width4d + c1] *= a;
				return mat(d, width, height, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator/(double a)const
		{
			if (data)
			{
				double* d(malloc256d(width, height));
				memcpy256d(d, data, width, height);
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						d[c0 * width4d + c1] /= a;
				return mat(d, width, height, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		//non-in-situ mult vec
		vec operator()(vec const& a)const
		{
			unsigned int minDim(width > a.dim ? a.dim : width);
			if (minDim && height)
			{
				vec r(height, false);
				return (*this)(a, r);
			}
			return vec();
		}
		vec& operator()(vec const& a, vec& b)const
		{
			if (matType < MatType::BandMat)
			{
				unsigned int minDim(width > a.dim ? a.dim : width);
				if (minDim)
				{
					bool overflow(ceiling4(minDim) > ceiling4(b.dim));
					if (overflow && b.type != Type::Native)return b;
					vec const* source(&a);
					vec r;
					if (&b == source)
					{
						source = &r;
						r = a;
					}
					if (overflow)b.reconstruct(minDim, false);
					constexpr unsigned int warp = 8;
					unsigned int minWidth4((minDim - 1) / 4 + 1);
					__m256d* aData((__m256d*)data);
					__m256d* bData((__m256d*)source->data);
					__m256d* rData((__m256d*)b.data);
					unsigned int heightFloor4((height >> 2) << 2);
					unsigned int widthWarp((minWidth4 / warp) * warp);
					unsigned int warpLeftFloor((minDim >> 2) - widthWarp);
					unsigned int warpLeftCeiling(minWidth4 - widthWarp);
					unsigned int c0(0);
					for (; c0 < heightFloor4; c0 += 4)
					{
						__m256d ans[4] = { 0 };
						__m256d tp[warp];
						unsigned int c1(0);
						for (; c1 < widthWarp; c1 += warp)
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
						if (c1 < minWidth4)
						{
							__m256d* s(aData + minWidth4 * c0 + c1);
#pragma unroll(4)
							for (unsigned int c2(0); c2 < warpLeftCeiling; ++c2)
								tp[c2] = bData[c1 + c2];
							unsigned int finalWidth(minDim - ((minDim >> 2) << 2));
							for (unsigned int c2(finalWidth); c2 < 4; ++c2)
								tp[warpLeftFloor].m256d_f64[c2] = 0;
							for (unsigned int c2(0); c2 < 4; ++c2, s += minWidth4)
							{
#pragma unroll(4)
								for (unsigned int c3(0); c3 < warpLeftCeiling; ++c3)
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
					if (c0 < height)
					{
						unsigned int heightLeft(height - heightFloor4);
						__m256d ans[4] = { 0 };
						__m256d tp[warp];
						unsigned int c1(0);
						for (; c1 < widthWarp; c1 += warp)
						{
							__m256d* s(aData + minWidth4 * c0 + c1);
#pragma unroll(4)
							for (unsigned int c2(0); c2 < warp; ++c2)
								tp[c2] = bData[c1 + c2];
							for (unsigned int c2(0); c2 < heightLeft; ++c2, s += minWidth4)
							{
#pragma unroll(4)
								for (unsigned int c3(0); c3 < warp; ++c3)
								{
									__m256d t = s[c3];
									ans[c2] = _mm256_fmadd_pd(t, tp[c3], ans[c2]);
								}
							}
						}
						if (c1 < minWidth4)
						{
							__m256d* s(aData + minWidth4 * c0 + c1);
#pragma unroll(4)
							for (unsigned int c2(0); c2 < warpLeftCeiling; ++c2)
								tp[c2] = bData[c1 + c2];
							unsigned int finalWidth(minDim - ((minDim >> 2) << 2));
							for (unsigned int c2(finalWidth); c2 < 4; ++c2)
								tp[warpLeftFloor].m256d_f64[c2] = 0;
							for (unsigned int c2(0); c2 < heightLeft; ++c2, s += minWidth4)
							{
#pragma unroll(4)
								for (unsigned int c3(0); c3 < warpLeftCeiling; ++c3)
								{
									__m256d t = s[c3];
									ans[c2] = _mm256_fmadd_pd(t, tp[c3], ans[c2]);
								}
							}
						}
						__m256d s;
						for (unsigned int c1(0); c1 < heightLeft; ++c1)
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
			//BandMat
			else if (matType == MatType::LBandMat)
			{
				unsigned int minDim(height > a.dim ? a.dim : height);
				if (minDim)
				{
					bool overflow(ceiling4(minDim) > ceiling4(b.dim));
					if (overflow && b.type != Type::Native)return b;
					vec const* source(&a);
					vec r;
					if (&b == source)
					{
						source = &r;
						r = a;
					}
					if (overflow)b.reconstruct(minDim, false);
					for (unsigned int c0(0); c0 < height; ++c0)
					{
						vec tp(getLBandRow(c0));
						unsigned int bgn(c0 <= halfBandWidth ? 0 : c0 - halfBandWidth);
						vec ta(a.data + bgn, tp.dim, Type::Non32Aligened);
						b.data[c0] = (tp, ta);
					}
				}
			}
			else if (matType == MatType::UBandMat)
			{
				unsigned int minDim(height > a.dim ? a.dim : height);
				if (minDim)
				{
					bool overflow(ceiling4(minDim) > ceiling4(b.dim));
					if (overflow && b.type != Type::Native)return b;
					vec const* source(&a);
					vec r;
					if (&b == source)
					{
						source = &r;
						r = a;
					}
					if (overflow)b.reconstruct(minDim, false);
					for (unsigned int c0(0); c0 < height; ++c0)
					{
						vec tp(getUBandRow(c0));
						vec ta(a.data + c0, tp.dim, Type::Non32Aligened);
						b.data[c0] = (tp, ta);
					}
				}
			}
			return b;
		}
		//non-in-situ mult mat
		mat operator()(mat const& a)const
		{
			unsigned int minDim(width > a.height ? a.height : width);
			if (minDim)
			{
				mat r(a.width, height, false);
				/*for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < minDim; ++c1)
						for (unsigned int c2(0); c2 < a.width; ++c2)
							r.data[c0 * a.width + c2] += data[c0 * width4d + c1] * a.data[c1 * a.width + c2];*/
				return (*this)(a, r);
			}
			return mat();
		}
		mat& operator()(mat const& a, mat& b)const
		{
			unsigned int minDim(width > a.height ? a.height : width);
			if (minDim)
			{
				bool overflow(ceiling4(a.width, height) > b.width4d * b.height);
				if (overflow && b.type != Type::Native)return b;
				mat const* source(this);
				mat r;
				if (&b == this)
				{
					source = &r;
					r = *this;
				}
				if (overflow)b.reconstruct(a.width, height, false);
				__m256d* aData((__m256d*)a.data);
				__m256d* bData((__m256d*)b.data);
				constexpr unsigned int warp = 16;
				unsigned int aWidth256d(a.width4d / 4);
				unsigned int aWidthWarpFloor(aWidth256d / warp * warp);
				unsigned int warpLeft(aWidth256d - aWidthWarpFloor);
				unsigned int height2Floor(height & (-2));
				unsigned int c0(0);
				for (; c0 < height2Floor; c0 += 2)
				{
					unsigned int c1(0);
					for (; c1 < aWidthWarpFloor; c1 += warp)
					{
						__m256d ans0[warp] = { 0 };
						__m256d ans1[warp] = { 0 };
						for (unsigned int c2(0); c2 < minDim; ++c2)
						{
							//__m256d t = _mm256_i32gather_pd(tempData, offset, 8);
							__m256d tp0 = _mm256_broadcast_sd(source->data + c0 * width4d + c2);
							__m256d tp1 = _mm256_broadcast_sd(source->data + (c0 + 1) * width4d + c2);
#pragma unroll(4)
							for (unsigned int c3(0); c3 < warp; ++c3)
							{
								__m256d b = aData[aWidth256d * c2 + c1 + c3];
								ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
								ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
							}
						}
#pragma unroll(4)
						for (unsigned int c3(0); c3 < warp; ++c3)
						{
							bData[c0 * aWidth256d + c1 + c3] = ans0[c3];
							bData[(c0 + 1) * aWidth256d + c1 + c3] = ans1[c3];
						}
					}
					if (c1 < aWidth256d)
					{
						__m256d ans0[warp] = { 0 };
						__m256d ans1[warp] = { 0 };
						for (unsigned int c2(0); c2 < minDim; ++c2)
						{
							__m256d tp0 = _mm256_broadcast_sd(source->data + c0 * width4d + c2);
							__m256d tp1 = _mm256_broadcast_sd(source->data + (c0 + 1) * width4d + c2);
							for (unsigned int c3(0); c3 < warpLeft; ++c3)
							{
								__m256d b = aData[aWidth256d * c2 + c1 + c3];
								ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
								ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
							}
						}
						for (unsigned int c3(0); c3 < warpLeft; ++c3)
						{
							bData[c0 * aWidth256d + c1 + c3] = ans0[c3];
							bData[(c0 + 1) * aWidth256d + c1 + c3] = ans1[c3];
						}
					}
				}
				if (c0 < height)
				{
					unsigned int c1(0);
					for (; c1 < aWidthWarpFloor; c1 += warp)
					{
						__m256d ans0[warp] = { 0 };
						for (unsigned int c2(0); c2 < minDim; ++c2)
						{
							__m256d tp0 = _mm256_broadcast_sd(source->data + c0 * width4d + c2);
#pragma unroll(4)
							for (unsigned int c3(0); c3 < warp; ++c3)
							{
								__m256d b = aData[aWidth256d * c2 + c1 + c3];
								ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
							}
						}
#pragma unroll(4)
						for (unsigned int c3(0); c3 < warp; ++c3)
							bData[c0 * aWidth256d + c1 + c3] = ans0[c3];
					}
					if (c1 < aWidth256d)
					{
						__m256d ans0[warp] = { 0 };
						for (unsigned int c2(0); c2 < minDim; ++c2)
						{
							__m256d tp0 = _mm256_broadcast_sd(source->data + c0 * width4d + c2);
							for (unsigned int c3(0); c3 < warpLeft; ++c3)
							{
								__m256d b = aData[aWidth256d * c2 + c1 + c3];
								ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
							}
						}
						for (unsigned int c3(0); c3 < warpLeft; ++c3)
							bData[c0 * aWidth256d + c1 + c3] = ans0[c3];
					}
				}
			}
			return b;
		}
		//true blas
		mat& schmidtOrtho()
		{
			if (data && width == height)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
				{
					vec tp(data + width4d * c0, width, Type::Parasitic);
					for (unsigned int c1(0); c1 < c0; ++c1)
					{
						vec ts(data + width4d * c1, width, Type::Parasitic);
						tp.fmadd(-(tp, ts), ts);
					}
					tp.normalize();
				}
			}
			return *this;
		}
		vec& solveL(vec const& a, vec& b)const
		{
			unsigned int minDim(height > a.dim ? a.dim : height);
			if (!minDim)return b;
			if (b.dim < minDim)
			{
				if (b.type == Type::Native)b.reconstruct(minDim, false);
				else return b;
			}
			double ll(data[0]);
			if (ll == 0.0)return b;
			b.data[0] = a.data[0] / ll;
			unsigned int c0(1);
			if (matType == MatType::LBandMat && a.dim >= height)
			{
				for (; c0 < height; ++c0)
				{
					ll = LBandEle(c0, c0);
					if (ll == 0.0)return b;
					vec tp(getLBandRowL(c0));
					unsigned int bgn(c0 <= halfBandWidth ? 0 : c0 - halfBandWidth);
					vec tb(b.data + bgn, tp.dim, Type::Non32Aligened);
					b.data[c0] = (a.data[c0] - (tp, tb)) / ll;
				}
			}
			else
			{
				for (; c0 < minDim; ++c0)
				{
					ll = data[c0 * width4d + c0];
					if (ll == 0.0)return b;
					vec tp(data + c0 * width4d, c0, Type::Parasitic);
					b.data[c0] = (a.data[c0] - (tp, b)) / ll;
				}
			}
			return b;
		}
		vec& solveU(vec const& a, vec& b)const
		{
			unsigned int minDim(height > a.dim ? a.dim : height);
			if (!minDim)return b;
			if (b.dim < minDim)
			{
				if (b.type == Type::Native)b.reconstruct(minDim, false);
				else return b;
			}
			double ll(data[(minDim - 1) * (width4d + 1)]);
			if (ll == 0.0)return b;
			int c0(minDim - 1);
			b.data[c0] = a.data[c0] / ll;
			--c0;
			if (matType == MatType::UBandMat && a.dim >= height)
			{
				for (; c0 >= 0; --c0)
				{
					ll = UBandEle(c0, c0);
					if (ll == 0.0)return b;
					vec tp(getUBandRowU(c0));
					vec tb(b.data + c0 + 1, tp.dim, Type::Non32Aligened);
					b.data[c0] = (a.data[c0] - (tp, tb)) / ll;
				}
			}
			else
			{
				for (; c0 >= 0; --c0)
				{
					ll = data[c0 * width4d + c0];
					if (ll == 0.0)return b;
					unsigned int c01(c0 + 1);
					vec tp(data + c0 * width4d + c01, minDim - c01, Type::Non32Aligened);
					vec tb(b.data + c01, minDim - c01, Type::Non32Aligened);
					b.data[c0] = (a.data[c0] - (tp, tb)) / ll;
				}
			}
			return b;
		}
		vec& solveUid(vec const& a, vec& b)const
		{
			//assuming that mat[i][i]==1
			unsigned int minDim(height > a.dim ? a.dim : height);
			if (!minDim)return b;
			if (b.dim < minDim)
			{
				if (b.type == Type::Native)b.reconstruct(minDim, false);
				else return b;
			}
			int c0(minDim - 1);
			b.data[c0] = a.data[c0];
			--c0;
			if (matType == MatType::UBandMat && a.dim >= height)
			{
				for (; c0 >= 0; --c0)
				{
					vec tp(getUBandRowU(c0));
					vec tb(b.data + c0 + 1, tp.dim, Type::Non32Aligened);
					b.data[c0] = (a.data[c0] - (tp, tb));
				}
			}
			else
			{
				for (; c0 >= 0; --c0)
				{
					unsigned int c01(c0 + 1);
					vec tp(data + c0 * width4d + c01, minDim - c01, Type::Non32Aligened);
					vec tb(b.data + c01, minDim - c01, Type::Non32Aligened);
					b.data[c0] = (a.data[c0] - (tp, tb));
				}
			}
			return b;
		}
		vec& solveGauss(vec& a, vec& b)
		{
			//no column principal
			unsigned int minDim(height > a.dim ? a.dim : height);
			if (!minDim)return b;
			if (b.dim < minDim)
			{
				if (b.type == Type::Native)b.reconstruct(minDim, false);
				else return b;
			}
			for (int c0(minDim - 1); c0 > 0; --c0)
			{
				double ll(data[c0 * width4d + c0]);
				if (ll == 0.0)return b;
				ll = -1 / ll;
				double aa(a.data[c0]);
				vec tp(data + c0 * width4d, c0, Type::Parasitic);
				for (int c1(c0 - 1); c1 >= 0; --c1)
				{
					double le(data[c1 * width4d + c0] * ll);
					if (le == 0.0)continue;
					vec ts(data + c1 * width4d, c0, Type::Parasitic);
					ts.fmadd(le, tp);
					a.data[c1] += le * aa;
					//data[c1 * width4d + c0] = 0;
				}
			}
			return solveL(a, b);
		}
		vec& solveJacobiIter(vec const& a, vec& b, double eps)const
		{
			unsigned int minDim(height > a.dim ? a.dim : height);
			vec irll(minDim, false);
			vec ll(minDim, false);
			vec b0(b.data, minDim, Type::Parasitic);
			vec b1(minDim, false);
			vec delta(minDim, false);
			b0 = a;
			for (unsigned int c0(0); c0 < minDim; ++c0)
				irll.data[c0] = -data[c0 * width4d + c0];
			for (unsigned int c0(0); c0 < 100; ++c0)
			{
				for (unsigned int c1(0); c1 < 5; ++c1)
				{
					(*this)(b0, b1);
					b1 -= a;
					b1 *= irll;
					b1 += b0;
					(*this)(b1, b0);
					b0 -= a;
					b0 *= irll;
					b0 += b1;
				}
				delta = b1; delta -= b0;
				if (delta.norm1() < eps)
				{
					::printf("iters:\t%d\n", c0 * 10);
					return b;
				}
			}
			return b;
		}
		//symmetric
		vec& solveCholesky(vec const& a, vec& b)
		{
			//only use L mat since it's a symmetric matrix...
			unsigned int minDim(height > a.dim ? a.dim : height);
			if (!minDim)return b;
			if (b.dim < minDim)
			{
				if (b.type == Type::Native)b.reconstruct(minDim, false);
				else return b;
			}
			double s(1 / data[0]);
			vec tp(minDim, false);
			tp[0] = s;
			for (unsigned int c0(1); c0 < minDim; ++c0)
				data[c0] = data[c0 * width4d] * s;
			for (unsigned int c0(1); c0 < minDim; ++c0)
			{
				vec tll(data + c0 * width4d, c0, Type::Parasitic);
				vec bll(b.data, c0, Type::Parasitic);
				bll = tll; bll *= tp;
				tp[c0] = 1 / (data[c0 * width4d + c0] -= (bll, tll));
				for (unsigned int c1(c0 + 1); c1 < minDim; ++c1)
				{
					double s(data[c1 * width4d + c0] -= (bll, vec(data + c1 * width4d, c0, Type::Parasitic)));
					data[c0 * width4d + c1] = s * tp[c0];
				}
			}
			solveL(a, tp);
			solveUid(tp, b);
		}
		vec& solveConjugateGradient(vec const& a, vec& b, double eps)const
		{
			//...
			return b;
		}
		//symmetric band (for LBandMat)
		vec& solveCholeskyBand(vec const& a, vec& b)//buggy...
		{
			unsigned int minDim(height > a.dim ? a.dim : height);
			if (!minDim)return b;
			if (b.dim < minDim)
			{
				if (b.type == Type::Native)b.reconstruct(minDim, false);
				else return b;
			}
			double s(1 / data[0]);
			vec tp(minDim, false);
			tp[0] = s;
			mat uM(halfBandWidth, height, MatType::UBandMat, true);
			for (unsigned int c0(1); c0 < 1 + halfBandWidth; ++c0)
				uM.data[c0] = data[c0 * width4d] * s;
			for (unsigned int c0(1); c0 < minDim; ++c0)
			{
				unsigned int len(c0 <= halfBandWidth ? c0 : halfBandWidth);
				unsigned int bgn(LBandBeginOffset(c0));
				unsigned int len4(ceiling4(len + bgn));
				vec tll(data + c0 * width4d, len4, Type::Parasitic);
				vec bll(b.data, len4, Type::Parasitic);
				unsigned int tbgn(c0 <= halfBandWidth ? 0 : c0 - (int(c0 - halfBandWidth) / 4) * 4);
				vec pll(tp.data + ((c0 - len) & -4), len4, Type::Parasitic);
				bll = tll; bll *= pll;
				vec bn(b.data + bgn, len, Type::Non32Aligened);
				tp[c0] = 1 / (LBandEleRef(c0, c0) -= (bn, tll));
				unsigned int c1(c0 + 1);
				for (; c1 < c0 + halfBandWidth && c1 < minDim; ++c1)
				{
					unsigned int bgn1(LBandBeginOffset(c1));
					unsigned int end1(c1 <= halfBandWidth ? c0 : c0 - (int(c1 - halfBandWidth) / 4) * 4);
					if (c1 > halfBandWidth && bgn1 == 0)
					{
						bll.data += 4;
						bll.dim -= 4;
					}
					double s(data[c1 * width4d + end1] -=
					(bll, vec(data + c1 * width4d + bgn1, end1 - bgn1, Type::Non32Aligened)));
					uM.UBandEleRef(c0, c1) = s * tp[c0];
				}
				if (c1 == c0 + halfBandWidth && c1 < minDim)
				{
					unsigned int end1(c1 <= halfBandWidth ? c0 : c0 - (int(c1 - halfBandWidth) / 4) * 4);
					uM.UBandEleRef(c0, c1) = data[c1 * width4d + end1] * tp[c0];
				}
			}
			solveL(a, tp);
			uM.solveUid(tp, b);
			return b;
		}

		void print()const
		{
			::printf("[\n");
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
				{
					::printf("\t[%.4e", data[width4d * c0]);
					unsigned int ed(width);
					if (matType == MatType::LBandMat || matType == MatType::UBandMat)
						ed = width4d;
					for (unsigned int c1(1); c1 < ed; ++c1)
						::printf(", %.4e", data[width4d * c0 + c1]);
					::printf("]\n");
				}
			}
			::printf("]\n");
		}
		void printInfo()const
		{
			char const* str = "";
			switch (matType)
			{
			case MatType::NormalMat:str = "NormalMat"; break;
			case MatType::SquareMat:str = "SquareMat"; break;
			case MatType::DiagonalMat:str = "DiagonalMat"; break;
			case MatType::LMat:str = "LMat"; break;
			case MatType::UMat:str = "UMat"; break;
			case MatType::BandMat:str = "BandMat"; break;
			case MatType::LBandMat:str = "LBandMat"; break;
			case MatType::UBandMat:str = "UBandMat"; break;
			}
			::printf("{width: %u, height: %u, type: %s, matType: %s}\n", width, height,
				type == Type::Native ? "Native" : "Parasitic", str);
		}
		void printToTxt(char const* name)const
		{
			//in the form of Mathematica matrix(for paste)
			if (data)
			{
				FILE* temp(::fopen(name, "w+"));
				::fprintf(temp, "{\n");
				for (unsigned int c0(0); c0 < height - 1; ++c0)
				{
					::fprintf(temp, "{%.14e", data[width4d * c0]);
					for (unsigned int c1(1); c1 < width; ++c1)
						::fprintf(temp, ", %.14e", data[width4d * c0 + c1]);
					::fprintf(temp, "},\n");
				}
				::fprintf(temp, "{%.14e", data[width4d * (height - 1)]);
				for (unsigned int c1(1); c1 < width; ++c1)
					::fprintf(temp, ", %.14e", data[width4d * (height - 1) + c1]);
				::fprintf(temp, "}\n}");
				::fclose(temp);
			}
		}
		void printToTableTxt(char const* name)const
		{
			//in the form of Mathematica matrix(for Import)
			if (data)
			{
				FILE* temp(::fopen(name, "w+"));
				for (unsigned int c0(0); c0 < height; ++c0)
				{
					for (unsigned int c1(0); c1 < width; ++c1)
						::fprintf(temp, "%.14e ", data[width4d * c0 + c1]);
					::fprintf(temp, "\n");
				}
			}
		}
	};

	//non-in-situ mult mat
	vec vec::operator()(mat const& a)const
	{
		unsigned int minDim(a.height > dim ? dim : a.height);
		if (minDim)
		{
			vec r(a.width, false);
			/*for (unsigned int c0(0); c0 < minDim; ++c0)
				for (unsigned int c1(0); c1 < a.width; ++c1)
					r.data[c1] += data[c0] * a.data[c0 * a.width4d + c1];*/
			return (*this)(a, r);
		}
		return vec();
	}
	vec& vec::operator()(mat const& a, vec& b)const
	{
		unsigned int minDim(a.height > dim ? dim : a.height);
		if (minDim)
		{
			bool overflow(ceiling4(minDim) > ceiling4(b.dim));
			if (overflow && b.type != Type::Native)return b;
			vec const* source(this);
			vec r;
			if (&b == this)
			{
				source = &r;
				r = *this;
			}
			if (overflow)b.reconstruct(minDim, false);
			constexpr unsigned int warp = 16;
			unsigned int width4(a.width4d >> 2);
			unsigned int widthWarpFloor((width4 / warp) * warp);
			unsigned int minDim4Floor(minDim & -4);
			__m256d* aData((__m256d*)a.data);
			__m256d* bData((__m256d*)b.data);
			unsigned int c0(0);
			for (; c0 < widthWarpFloor; c0 += warp)
			{
				__m256d ans[warp] = { 0 };
				unsigned int c1(0);
				__m256d tp[4];
				for (; c1 < minDim4Floor; c1 += 4)
				{
					tp[0] = _mm256_broadcast_sd(source->data + c1);
					tp[1] = _mm256_broadcast_sd(source->data + c1 + 1);
					tp[2] = _mm256_broadcast_sd(source->data + c1 + 2);
					tp[3] = _mm256_broadcast_sd(source->data + c1 + 3);
					__m256d* s(aData + width4 * c1 + c0);
#pragma unroll(4)
					for (unsigned int c2(0); c2 < 4; ++c2)
					{
#pragma unroll(4)
						for (unsigned int c3(0); c3 < warp; ++c3)
							ans[c3] = _mm256_fmadd_pd(s[c2 * width4 + c3], tp[c2], ans[c3]);
					}
				}
				if (c1 < minDim)
				{
					unsigned int deltaMinDim(minDim - c1);
					for (unsigned int c2(0); c2 < deltaMinDim; ++c2)
					{
						double b = source->data[c1 + c2];
						tp[c2] = { b,b,b,b };
					}
					__m256d* s(aData + width4 * c1 + c0);
#pragma unroll(4)
					for (unsigned int c2(0); c2 < deltaMinDim; ++c2)
					{
#pragma unroll(4)
						for (unsigned int c3(0); c3 < warp; ++c3)
							ans[c3] = _mm256_fmadd_pd(s[c2 * width4 + c3], tp[c2], ans[c3]);
					}
				}
#pragma unroll(4)
				for (unsigned int c3(0); c3 < warp; ++c3)
					bData[c0 + c3] = ans[c3];
			}
			if (c0 < a.width4d)
			{
				unsigned int warpLeft(width4 - widthWarpFloor);
				__m256d ans[warp] = { 0 };
				unsigned int c1(0);
				__m256d tp[4];
				for (; c1 < minDim4Floor; c1 += 4)
				{
					tp[0] = _mm256_broadcast_sd(source->data + c1);
					tp[1] = _mm256_broadcast_sd(source->data + c1 + 1);
					tp[2] = _mm256_broadcast_sd(source->data + c1 + 2);
					tp[3] = _mm256_broadcast_sd(source->data + c1 + 3);
					__m256d* s(aData + width4 * c1 + c0);
#pragma unroll(4)
					for (unsigned int c2(0); c2 < 4; ++c2)
					{
#pragma unroll(4)
						for (unsigned int c3(0); c3 < warpLeft; ++c3)
							ans[c3] = _mm256_fmadd_pd(s[c2 * width4 + c3], tp[c2], ans[c3]);
					}
				}
				if (c1 < minDim)
				{
					unsigned int deltaMinDim(minDim - c1);
					for (unsigned int c2(0); c2 < deltaMinDim; ++c2)
					{
						double b = source->data[c1 + c2];
						tp[c2] = { b,b,b,b };
					}
					__m256d* s(aData + width4 * c1 + c0);
#pragma unroll(4)
					for (unsigned int c2(0); c2 < deltaMinDim; ++c2)
					{
#pragma unroll(4)
						for (unsigned int c3(0); c3 < warpLeft; ++c3)
							ans[c3] = _mm256_fmadd_pd(s[c2 * width4 + c3], tp[c2], ans[c3]);
					}
				}
#pragma unroll(4)
				for (unsigned int c3(0); c3 < warpLeft; ++c3)
					bData[c0 + c3] = ans[c3];
			}
		}
		return b;
	}
}