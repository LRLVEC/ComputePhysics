#pragma once
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <cstdio>
#include <immintrin.h>

namespace BLAS
{
	static constexpr double Pi = 3.14159265358979323846264338327950288L;
	static constexpr double E = 2.71828182845904523536028747135266250L;

	enum class Type
	{
		Native = 0,
		Parasitic = 1,
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

	struct mat;
	struct vec
	{
		double* data;
		unsigned int dim;
		Type type;

		vec() :data(nullptr), dim(0), type(Type::Native) {}
		vec(unsigned int _length, bool _clear = true)
			:
			data((double*)::malloc(_length * sizeof(double))),
			dim(_length),
			type(Type::Native)
		{
			if (_clear)memset(data, 0, _length * sizeof(double));
		}
		vec(vec const& a)
			:
			data(a.dim ? (double*)::malloc(a.dim * sizeof(double)) : nullptr),
			dim(a.dim),
			type(Type::Native)
		{
			if (dim)
				::memcpy(data, a.data, dim * sizeof(double));
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
					data = (double*)::malloc(a.dim * sizeof(double));
					dim = a.dim;
				}
			}
		}
		vec(double* _data, unsigned int _length, Type _type) :data(_data), dim(_length), type(_type) {}
		vec(std::initializer_list<double>const& a)
			:
			data(a.size() ? (double*)::malloc(a.size() * sizeof(double)) : nullptr),
			dim(a.size()),
			type(Type::Native)
		{
			if (dim)::memcpy(data, a.begin(), dim * sizeof(double));
		}
		~vec()
		{
			if (type == Type::Native)::free(data);
			data = nullptr;
			dim = 0;
		}
		template<class T>inline double& operator[](T a)
		{
			return data[a];
		}
		//moveTo
		vec& moveTo(vec& a)
		{
			if (type == Type::Native)
			{
				a.~vec();
				a.data = data;
				a.dim = dim;
				a.type = type;
				data = nullptr;
				dim = 0;
			}
			else a = *this;
			return a;
		}
		//= += -= *= /=
		vec& operator =(vec&& a)
		{
			if (a.type == Type::Native)
			{
				::free(data);
				data = a.data;
				dim = a.dim;
				a.data = nullptr;
				a.dim = 0;
			}
			else
			{
				if (dim >= a.dim)
					::memcpy(data, a.data, a.dim * sizeof(double));
				else
				{
					if (type == Type::Native)
					{
						::free(data);
						data = (double*)::malloc(a.dim * sizeof(double));
						dim = a.dim;
					}
					::memcpy(data, a.data, dim * sizeof(double));
				}
			}
		}
		vec& operator =(vec const& a)
		{
			if (a.dim)
			{
				if (dim >= a.dim)
					::memcpy(data, a.data, a.dim * sizeof(double));
				else
				{
					if (type == Type::Native)
					{
						::free(data);
						data = (double*)::malloc(a.dim * sizeof(double));
						dim = a.dim;
					}
					::memcpy(data, a.data, dim * sizeof(double));
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
				double* d((double*)::malloc(l * sizeof(double)));
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
				double* d((double*)::malloc(l * sizeof(double)));
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
				double* d((double*)::malloc(l * sizeof(double)));
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
				double* d((double*)::malloc(l * sizeof(double)));
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
				double* d((double*)::malloc(l * sizeof(double)));
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
				double* d((double*)::malloc(l * sizeof(double)));
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
				double* d((double*)::malloc(l * sizeof(double)));
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
				double* d((double*)::malloc(l * sizeof(double)));
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
	};
	struct mat
	{
		double* data;
		unsigned int width;
		unsigned int height;
		Type type;
		MatType matType;

		mat() :data(nullptr), width(0), height(0), type(Type::Native), matType(MatType::NormalMat) {}
		mat(unsigned int _width, unsigned int _height, bool _clear = true)
			:
			data((_width&& _height) ? (double*)::malloc(_width * sizeof(double) * _height) : nullptr),
			width(data ? _width : 0),
			height(data ? _height : 0),
			type(Type::Native),
			matType(MatType::NormalMat)
		{
			if (_clear && data)
				memset(data, 0, width * sizeof(double) * height);
		}
		mat(mat const& a)
			:
			data((a.width&& a.height) ? (double*)::malloc(a.width * sizeof(double) * a.height) : nullptr),
			width(a.width),
			height(a.height),
			type(Type::Native),
			matType(a.matType)
		{
			if (data)
				::memcpy(data, a.data, a.width * sizeof(double) * a.height);
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
					a.width = a.height = 0;
				}
				else
				{
					data = (double*)::malloc(a.width * sizeof(double) * a.height);
					width = a.width;
					height = a.height;
					memcpy(data, a.data, a.width * sizeof(double) * a.height);
				}
			}
		}
		mat(double* _data, unsigned int _width, unsigned int _height, Type _type, MatType _matType)
			:
			data(_data), width(_width), height(_height), type(_type), matType(_matType)
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
					size_t s(w * sizeof(double) * h);
					data = (double*)::malloc(s);
					memset(data, 0, s);
					width = w;
					height = h;
					type = Type::Native;
					matType = MatType::NormalMat;
					for (unsigned int c0(0); c0 < h; ++c0)
						memcpy(data + w * c0, (a.begin() + c0)->begin(),
						(a.begin() + c0)->size() * sizeof(double));
				}
			}
		}
		~mat()
		{
			if (type == Type::Native)::free(data);
			data = nullptr;
			width = height = 0;
		}
		template<class T>inline double& operator[](T a)
		{
			return data[a];
		}
		template<class T, class R>inline double& operator() (T a, R b)
		{
			return data[a * width + b];
		}
		//= += -= *= /=
		mat& operator =(mat&& a)
		{
			if (a.data)
			{
				if (type == Type::Native)
				{
					if (a.type == Type::Native)
					{
						::free(data);
						data = a.data;
						width = a.width;
						height = a.height;
						matType = a.matType;
						a.data = nullptr;
						a.width = a.height = 0;
					}
					else
					{
						if (unsigned long long(a.width) * a.height == unsigned long long(width) * height)
						{
							::free(data);
							data = (double*)::malloc(a.width * sizeof(double) * a.height);
						}
						width = a.width;
						height = a.height;
						memcpy(data, a.data, a.width * sizeof(double) * a.height);
					}
				}
				else
				{
					unsigned int minWidth(width <= a.width ? width : a.width);
					unsigned int minHeight(height <= a.height ? height : a.height);
					for (unsigned int c0(0); c0 < minHeight; ++c0)
						::memcpy(data + width * c0, a.data + a.width * c0, minWidth * sizeof(double));
				}
			}
			return *this;
		}
		mat& operator =(mat const& a)
		{
			if (a.data)
			{
				if (unsigned long long(a.width) * a.height == unsigned long long(width) * height)
					::memcpy(data, a.data, width * sizeof(double) * height);
				else
				{
					if (type == Type::Native)
					{
						::free(data);
						data = (double*)::malloc(a.width * sizeof(double) * a.height);
						width = a.width;
						height = a.height;
						matType = a.matType;
						::memcpy(data, a.data, width * sizeof(double) * height);
					}
					else
					{
						unsigned int minWidth(width <= a.width ? width : a.width);
						unsigned int minHeight(height <= a.height ? height : a.height);
						for (unsigned int c0(0); c0 < minHeight; ++c0)
							::memcpy(data + width * c0, a.data + a.width * c0, minWidth * sizeof(double));
					}
				}
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
						data[c0 * width + c1] += a.data[c0 * a.width + c1];
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
						data[c0 * width + c1] -= a.data[c0 * a.width + c1];
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
						data[c0 * width + c1] *= a.data[c0 * a.width + c1];
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
						data[c0 * width + c1] /= a.data[c0 * a.width + c1];
			}
			return *this;
		}
		mat& operator =(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width + c1] = a;
			}
			return *this;
		}
		mat& operator+=(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width + c1] += a;
			}
			return *this;
		}
		mat& operator-=(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width + c1] -= a;
			}
			return *this;
		}
		mat& operator*=(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width + c1] *= a;
			}
			return *this;
		}
		mat& operator/=(double a)
		{
			if (data)
			{
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
						data[c0 * width + c1] /= a;
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
				double* d((double*)::malloc(minW * sizeof(double) * minH));
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						d[c0 * minW + c1] = data[c0 * width + c1] + a.data[c0 * a.width + c1];
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
				double* d((double*)::malloc(minW * sizeof(double) * minH));
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						d[c0 * minW + c1] = data[c0 * width + c1] - a.data[c0 * a.width + c1];
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
				double* d((double*)::malloc(minW * sizeof(double) * minH));
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						d[c0 * minW + c1] = data[c0 * width + c1] * a.data[c0 * a.width + c1];
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
				double* d((double*)::malloc(minW * sizeof(double) * minH));
				for (unsigned int c0(0); c0 < minH; ++c0)
					for (unsigned int c1(0); c1 < minW; ++c1)
						d[c0 * minW + c1] = data[c0 * width + c1] / a.data[c0 * a.width + c1];
				return mat(d, minW, minH, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator+(double a)const
		{
			if (data)
			{
				double* d((double*)::malloc(width * sizeof(double) * height));
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
					{
						unsigned int id(c0 * width + c1);
						d[id] = data[id] + a;
					}
				return mat(d, width, height, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator-(double a)const
		{
			if (data)
			{
				double* d((double*)::malloc(width * sizeof(double) * height));
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
					{
						unsigned int id(c0 * width + c1);
						d[id] = data[id] - a;
					}
				return mat(d, width, height, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator*(double a)const
		{
			if (data)
			{
				double* d((double*)::malloc(width * sizeof(double) * height));
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
					{
						unsigned int id(c0 * width + c1);
						d[id] = data[id] * a;
					}
				return mat(d, width, height, Type::Native, MatType::NormalMat);
			}
			return mat();
		}
		mat operator/(double a)const
		{
			if (data)
			{
				double* d((double*)::malloc(width * sizeof(double) * height));
				for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < width; ++c1)
					{
						unsigned int id(c0 * width + c1);
						d[id] = data[id] / a;
					}
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
				vec r(height);
				for (unsigned int c0(0); c0 < minDim; ++c0)
					for (unsigned int c1(0); c1 < height; ++c1)
						r.data[c1] += a.data[c0] * data[c1 * width + c0];
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
				mat r(a.width, height);
				/*for (unsigned int c0(0); c0 < height; ++c0)
					for (unsigned int c1(0); c1 < minDim; ++c1)
						for (unsigned int c2(0); c2 < a.width; ++c2)
							r.data[c0 * a.width + c2] += data[c0 * width + c1] * a.data[c1 * a.width + c2];*/

				__m256d* aData((__m256d*)a.data);
				__m256d* rData((__m256d*)r.data);
				unsigned int aWidth4(a.width / 4);
				constexpr unsigned int warp = 16;
				for (unsigned int c0(0); c0 < height; c0 += 2)
					for (unsigned int c1(0); c1 < aWidth4; c1 += warp)
					{
						__m256d ans0[warp] = { 0 };
						__m256d ans1[warp] = { 0 };
						for (unsigned int c2(0); c2 < minDim; ++c2)
						{
							//__m256d t = _mm256_i32gather_pd(tempData, offset, 8);
							double s = data[c0 * width + c2];
							__m256d tp0 = { s,s,s,s };
							s = data[(c0 + 1) * width + c2];
							__m256d tp1 = { s,s,s,s };
#pragma unroll(4)
							for (unsigned int c3(0); c3 < warp; ++c3)
							{
								__m256d b = aData[aWidth4 * c2 + c1 + c3];
								ans0[c3] = _mm256_fmadd_pd(tp0, b, ans0[c3]);
								ans1[c3] = _mm256_fmadd_pd(tp1, b, ans1[c3]);
							}
						}
#pragma unroll(4)
						for (unsigned int c3(0); c3 < warp; ++c3)
						{
							rData[c0 * aWidth4 + c1 + c3] = ans0[c3];
							rData[(c0 + 1) * aWidth4 + c1 + c3] = ans1[c3];
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
					::printf("\t[%.4f", data[width * c0]);
					for (unsigned int c1(1); c1 < width; ++c1)
						::printf(", %.4f", data[width * c0 + c1]);
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
					::fprintf(temp, "{%.8f", data[width * c0]);
					for (unsigned int c1(1); c1 < width; ++c1)
						::fprintf(temp, ", %.8f", data[width * c0 + c1]);
					::fprintf(temp, "},\n");
				}
				::fprintf(temp, "{%.8f", data[width * (height - 1)]);
				for (unsigned int c1(1); c1 < width; ++c1)
					::fprintf(temp, ", %.8f", data[width * (height - 1) + c1]);
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
			vec r(a.width);
			for (unsigned int c0(0); c0 < minDim; ++c0)
				for (unsigned int c1(0); c1 < a.width; ++c1)
					r.data[c1] += data[c0] * a.data[c0 * a.width + c1];
			return r;
		}
		return vec();
	}
}