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
		Non32Aligened = 2,
	};
	enum class MatType
	{
		NormalMat,
		SquareMat,
		DiagonalMat,
		LMat,
		UMat,
		BandMat,
	};

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
		return ((((unsigned long long(width) - 1) >> 2) + 1) << 2)* height;
	}
	inline size_t ceiling256dSize(unsigned int length)
	{
		return (((size_t(length) - 1) >> 2) + 1) << 5;
	}
	inline size_t ceiling256dSize(unsigned int width, unsigned int height)
	{
		return ((((size_t(width) - 1) >> 2) + 1) << 5)* height;
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

	struct mat;
	struct vec
	{
		double* data;
		unsigned int dim;
		Type type;

		vec() :data(nullptr), dim(0), type(Type::Native) {}
		vec(unsigned int _length, bool _clear = true)
			:
			data(_length ? malloc256d(_length) : nullptr),
			dim(_length),
			type(Type::Native)
		{
			if (_clear && data)memset256d(data, 0, _length);
		}
		vec(vec const& a)
			:
			data(a.dim ? malloc256d(a.dim) : nullptr),
			dim(a.dim),
			type(Type::Native)
		{
			if (dim)memcpy64d(data, a.data, dim);
		}
		vec(vec&& a) :data(nullptr), dim(0), type(Type::Native)
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
					memcpy64d(data, a.data, dim);
				}
			}
		}
		vec(double* _data, unsigned int _length, Type _type) :data(_data), dim(_length), type(_type) {}
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
			return data[a];
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
					_mm_free(data);
					if (a.data)data = malloc256d(a.dim);
					else data = nullptr;
					dim = a.dim;
				}
			}
			else
			{
				unsigned minDim(dim > a.dim ? a.dim : dim);
				if (minDim)memcpy64d(data, a.data, minDim);
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
					memcpy64d(data, a.data, dim);
				}
				else
				{
					unsigned minDim(dim > a.dim ? a.dim : dim);
					memcpy64d(data, a.data, minDim);
				}
			}
			return *this;
		}
		vec& operator+=(vec const& a)
		{
			if (a.dim && dim)
			{
				for (unsigned int c0(0); c0 < (dim > a.dim ? a.dim : dim); ++c0)
					data[c0] += a.data[c0];
			}
			return *this;
		}
		vec& operator-=(vec const& a)
		{
			if (a.dim && dim)
			{
				for (unsigned int c0(0); c0 < (dim > a.dim ? a.dim : dim); ++c0)
					data[c0] -= a.data[c0];
			}
			return *this;
		}
		vec& operator*=(vec const& a)
		{
			if (a.dim && dim)
			{
				for (unsigned int c0(0); c0 < (dim > a.dim ? a.dim : dim); ++c0)
					data[c0] *= a.data[c0];
			}
			return *this;
		}
		vec& operator/=(vec const& a)
		{
			if (a.dim && dim)
			{
				for (unsigned int c0(0); c0 < (dim > a.dim ? a.dim : dim); ++c0)
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
				double s(0);
				unsigned int l(dim > a.dim ? a.dim : dim);
				for (unsigned int c0(0); c0 < l; ++c0)
					s += data[c0] * a.data[c0];
				return s;
			}
			else return 0;
		}
		//non-in-situ mult mat
		vec operator()(mat const& a)const;

		//norm
		double norm1()const
		{
			if (dim)
			{
				double s(0);
				for (unsigned int c0(0); c0 < dim; ++c0)
					s += abs(data[c0]);
				return s;
			}
			return 0;
		}
		double norm2()const
		{
			if (dim)
			{
				double s(0);
				for (unsigned int c0(0); c0 < dim; ++c0)
					s += data[c0] * data[c0];
				return sqrt(s);
			}
			return 0;
		}
		double normInf()const
		{
			if (dim)
			{
				double s(0);
				for (unsigned int c0(0); c0 < dim; ++c0)
					if (s < abs(data[c0]))s = abs(data[c0]);
				return s;
			}
			return 0;
		}
		double normP(double p)const
		{
			if (dim && p)
			{
				double s(0);
				for (unsigned int c0(0); c0 < dim; ++c0)
				{
					s += pow(abs(data[c0]), p);
				}
				return pow(s, 1 / p);
			}
			return 0;
		}

		void print()const
		{
			::printf("[");
			if (dim)
			{
				for (unsigned int c0(0); c0 < dim - 1; ++c0)
					::printf("%.4f, ", data[c0]);
				::printf("%.4f", data[dim - 1]);
			}
			::printf("]\n");
		}
		void printInfo()const
		{
			::printf("{length: %u, type: %s}\n", dim,
				type == Type::Native ? "Native" : "Parasitic");
		}
		void printToTxt(char const* name)const
		{
			//in the form of Mathematica matrix
			if (data)
			{
				FILE* temp(::fopen(name, "w+"));
				::fprintf(temp, "{");
				for (unsigned int c0(0); c0 < dim - 1; ++c0)
					::fprintf(temp, "{%.8f}, ", data[c0]);
				::fprintf(temp, "{%.8f}}\n", data[dim - 1]);
				::fclose(temp);
			}
		}
	};
	struct mat
	{
		double* data;
		unsigned int width;
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
		template<class T, class R>inline double& operator() (T a, R b)
		{
			return data[a * width4d + b];
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
				/*for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < minDim; ++c1)
						r.data[c0] += a.data[c1] * data[c0 * width4d + c1];*/
						/*constexpr unsigned int warp = 4;
								int d(width4d);
								__m128i offset;
								offset.m128i_i32[0] = 0;
								offset.m128i_i32[1] = d;
								offset.m128i_i32[2] = 2 * d;
								offset.m128i_i32[3] = 3 * d;
								__m256d* rData((__m256d*)r.data);
								for (unsigned int c0(0); c0 < height; c0 += 4 * warp)
								{
									__m256d ans[warp] = { 0 };
									for (unsigned int c1(0); c1 < minDim; c1 += 4)
									{
										__m256d tp[4];
										double b = a.data[c1];
										tp[0] = { b,b,b,b };
										b = a.data[c1 + 1];
										tp[1] = { b,b,b,b };
										b = a.data[c1 + 2];
										tp[2] = { b,b,b,b };
										b = a.data[c1 + 3];
										tp[3] = { b,b,b,b };
										double* ad(data + width4d * c0 + c1);
				#pragma unroll(4)
										for (unsigned int c2(0); c2 < warp; ++c2)
										{
											double* adf(ad + c2 * 4 * width4d);
				#pragma unroll(4)
											for (unsigned int c3(0); c3 < 4; ++c3)
											{
												__m256d t = _mm256_i32gather_pd(adf + c3, offset, 8);
												ans[c2] = _mm256_fmadd_pd(t, tp[c3], ans[c2]);
											}
										}
									}
				#pragma unroll(4)
									for (unsigned int c2(0); c2 < warp; ++c2)
										rData[(c0 >> 2) + c2] = ans[c2];
								}*/
				constexpr unsigned int warp = 8;
				unsigned int minWidth4((minDim - 1) / 4 + 1);
				__m256d* aData((__m256d*)data);
				__m256d* bData((__m256d*)a.data);
				__m256d* rData((__m256d*)r.data);
				for (unsigned int c0(0); c0 < height; c0 += 4)
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
				return r;
			}
			return vec();
		}
		//non-in-situ mult mat
		mat operator()(mat const& a)const
		{
			if (height && a.width)
			{
				unsigned int minDim(width > a.height ? a.height : width);
				mat r(a.width, height, false);
				/*for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < minDim; ++c1)
						for (unsigned int c2(0); c2 < a.width; ++c2)
							r.data[c0 * a.width + c2] += data[c0 * width4d + c1] * a.data[c1 * a.width + c2];*/

				__m256d* aData((__m256d*)a.data);
				__m256d* rData((__m256d*)r.data);
				unsigned int aWidth256d(a.width4d / 4);
				constexpr unsigned int warp = 16;
				for (unsigned int c0(0); c0 < height; c0 += 2)
					for (unsigned int c1(0); c1 < aWidth256d; c1 += warp)
					{
						__m256d ans0[warp] = { 0 };
						__m256d ans1[warp] = { 0 };
						for (unsigned int c2(0); c2 < minDim; ++c2)
						{
							//__m256d t = _mm256_i32gather_pd(tempData, offset, 8);
							double s = data[c0 * width4d + c2];
							__m256d tp0 = { s,s,s,s };
							s = data[(c0 + 1) * width4d + c2];
							__m256d tp1 = { s,s,s,s };
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
							rData[c0 * aWidth256d + c1 + c3] = ans0[c3];
							rData[(c0 + 1) * aWidth256d + c1 + c3] = ans1[c3];
						}
					}
				return r;
			}
			return mat();
		}

		void print()const
		{
			::printf("[\n");
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
				{
					::printf("\t[%.4f", data[width4d * c0]);
					for (unsigned int c1(1); c1 < width; ++c1)
						::printf(", %.4f", data[width4d * c0 + c1]);
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
			}
			::printf("{width: %u, height: %u, type: %s, matType: %s}\n", width, height,
				type == Type::Native ? "Native" : "Parasitic", str);
		}
		void printToTxt(char const* name)const
		{
			//in the form of Mathematica matrix
			if (data)
			{
				FILE* temp(::fopen(name, "w+"));
				::fprintf(temp, "{\n");
				for (unsigned int c0(0); c0 < height - 1; ++c0)
				{
					::fprintf(temp, "{%.8f", data[width4d * c0]);
					for (unsigned int c1(1); c1 < width; ++c1)
						::fprintf(temp, ", %.8f", data[width4d * c0 + c1]);
					::fprintf(temp, "},\n");
				}
				::fprintf(temp, "{%.8f", data[width4d * (height - 1)]);
				for (unsigned int c1(1); c1 < width; ++c1)
					::fprintf(temp, ", %.8f", data[width4d * (height - 1) + c1]);
				::fprintf(temp, "}\n}");
				::fclose(temp);
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
			constexpr unsigned int warp = 16;
			unsigned int width4((a.width - 1) / 4 + 1);
			__m256d* aData((__m256d*)a.data);
			__m256d* rData((__m256d*)r.data);
			for (unsigned int c0(0); c0 < width4; c0 += warp)
			{
				__m256d ans[warp] = { 0 };
				for (unsigned int c1(0); c1 < minDim; c1 += 4)
				{
					__m256d tp[4];
					double b = data[c1];
					tp[0] = { b,b,b,b };
					b = data[c1 + 1];
					tp[1] = { b,b,b,b };
					b = data[c1 + 2];
					tp[2] = { b,b,b,b };
					b = data[c1 + 3];
					tp[3] = { b,b,b,b };
					__m256d* s(aData + width4 * c1 + c0);
#pragma unroll(4)
					for (unsigned int c2(0); c2 < warp; ++c2, ++s)
					{
#pragma unroll(4)
						for (unsigned int c3(0); c3 < 4; ++c3)
						{
							__m256d t = s[c3 * width4];
							ans[c2] = _mm256_fmadd_pd(t, tp[c3], ans[c2]);
						}
					}
				}
#pragma unroll(4)
				for (unsigned int c3(0); c3 < warp; ++c3)
					rData[c0 + c3] = ans[c3];
			}
			return r;
		}
		return vec();
	}
}