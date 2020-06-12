#pragma once
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <_TemplateMeta.h>

namespace Math
{
#define CheckNumType(T)	static_assert(NumType<T>::value, "Wrong NumType!")
#define CheckDim		static_assert(_dim, "Dimension cannot be 0!")
#define CheckRowDim		static_assert(_rowDim, "Dimension cannot be 0!")
#define CheckColDim		static_assert(_colDim, "Dimension cannot be 0!")
	static constexpr double Pi = 3.14159265358979323846264338327950288L;
	static constexpr double E = 2.71828182845904523536028747135266250L;

	template<class T, unsigned int _dim>struct vec;
	template<class T, unsigned int _row, unsigned int _col>struct mat;
	template<class T>using vec2 = vec<T, 2>;
	template<class T>using vec3 = vec<T, 3>;
	template<class T>using vec4 = vec<T, 4>;
	template<class T>using mat2 = mat<T, 2, 2>;
	template<class T>using mat3 = mat<T, 3, 3>;
	template<class T>using mat4 = mat<T, 4, 4>;

	template<class T, unsigned int _dim>struct vec
	{
		CheckNumType(T);
		CheckDim;
		static constexpr unsigned int dim = _dim;

		T data[_dim];

		vec();
		vec(vec<T, _dim>const&) = default;
		template<class R>vec(R const&);
		template<class R, unsigned int _dim1>vec(vec<R, _dim1>const&);
		vec(std::initializer_list<T>const&);
		~vec() = default;
		T& operator[](int);
		vec<T, _dim>& operator=(vec<T, _dim>const&) = default;
		template<class R>vec<T, _dim>& operator= (R const&);
		template<class R>vec<T, _dim>& operator+=(R const&);
		template<class R>vec<T, _dim>& operator-=(R const&);
		template<class R>vec<T, _dim>& operator*=(R const&);
		template<class R>vec<T, _dim>& operator/=(R const&);
		bool operator==(vec<T, _dim>const&)const;
		template<class R>bool operator==(R const&)const;
		template<class R>bool operator!=(R const&)const;
		template<class R, unsigned int _dim1>vec<T, _dim>& operator =(vec<R, _dim1>const&);
		template<class R, unsigned int _dim1>vec<T, _dim>& operator+=(vec<R, _dim1>const&);
		template<class R, unsigned int _dim1>vec<T, _dim>& operator-=(vec<R, _dim1>const&);
		template<class R, unsigned int _dim1>vec<T, _dim>& operator*=(vec<R, _dim1>const&);
		template<class R, unsigned int _dim1>vec<T, _dim>& operator/=(vec<R, _dim1>const&);

		vec<T, _dim> operator-()const;
		template<class R>auto operator+(R const&)const;
		template<class R>auto operator-(R const&)const;
		template<class R>auto operator*(R const&)const;
		template<class R>auto operator/(R const&)const;
		template<class R, unsigned int _dim1>auto operator+(vec<R, _dim1>const&)const;
		template<class R, unsigned int _dim1>auto operator-(vec<R, _dim1>const&)const;
		template<class R, unsigned int _dim1>auto operator*(vec<R, _dim1>const&)const;
		template<class R, unsigned int _dim1>auto operator/(vec<R, _dim1>const&)const;
		//dot
		template<class R, unsigned int _dim1>auto operator,(vec<R, _dim1>const&)const;
		//cross
		template<class R, unsigned int _dim1>auto operator|(vec<R, _dim1>const&);
		mat3<T> crossMat();
		//tensor product
		template<class R, unsigned int _dim1>auto operator^(vec<R, _dim1>const&);
		//transfer
		template<class R, unsigned int _rowDim, unsigned int _colDim>auto operator()(mat<R, _rowDim, _colDim>const&);
		//rotate
		template<class R, class Y, unsigned int _dim1>auto operator()(vec<R, _dim1>const&, Y);
		template<class Y>mat3<double> rotMat(Y);
		template<class R, class Y, unsigned int _rowDim1, unsigned int _colDim1>auto& operator()(mat<R, _rowDim1, _colDim1>*, Y);
		template<class R, class Y, unsigned int _rowDim1, unsigned int _colDim1>auto operator()(mat<R, _rowDim1, _colDim1>const&, Y);
		//square
		T square()const;
		T square(int)const;
		//length
		T length()const;
		T length(int)const;
		//max
		T max()const;
		T max(int)const;
		//min
		T min()const;
		T min(int)const;
		//normaliaze
		vec<T, _dim>& normaliaze();
		vec<T, _dim>& normaliaze(int);
		//print
		void print()const;
		void printInfo(char const*)const;
	};
	template<class T, unsigned int _rowDim, unsigned int _colDim>struct mat
	{
		CheckNumType(T);
		CheckRowDim;
		CheckColDim;
		static constexpr unsigned int rowDim = _rowDim;
		static constexpr unsigned int colDIm = _colDim;
		using row = vec<T, _colDim>;
		using col = vec<T, _rowDim>;

		union
		{
			T array[_rowDim][_colDim];
			row rowVec[_rowDim];
		};

		mat();
		mat(mat<T, _rowDim, _colDim>const&) = default;
		mat(std::initializer_list<row>const&);
		template<class R>mat(R const&);
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>mat(mat<R, _rowDim1, _colDim1>const&);
		~mat() = default;
		row& operator[](int);
		col column(int)const;
		template<class R, unsigned int _dim>mat<T, _rowDim, _colDim>& setCol(vec<R, _dim>&, unsigned int);
		mat<T, _rowDim, _colDim>& operator=(mat<T, _rowDim, _colDim>const&) = default;
		bool operator==(mat<T, _rowDim, _colDim>const&)const;

		template<class R>mat<T, _rowDim, _colDim>& operator= (R const&);
		template<class R>mat<T, _rowDim, _colDim>& operator+=(R const&);
		template<class R>mat<T, _rowDim, _colDim>& operator-=(R const&);
		template<class R>mat<T, _rowDim, _colDim>& operator*=(R const&);
		template<class R>mat<T, _rowDim, _colDim>& operator/=(R const&);
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>mat<T, _rowDim, _colDim>& operator =(mat<R, _rowDim1, _colDim1>const&);
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>mat<T, _rowDim, _colDim>& operator+=(mat<R, _rowDim1, _colDim1>const&);
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>mat<T, _rowDim, _colDim>& operator-=(mat<R, _rowDim1, _colDim1>const&);
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>mat<T, _rowDim, _colDim>& operator*=(mat<R, _rowDim1, _colDim1>const&);
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>mat<T, _rowDim, _colDim>& operator/=(mat<R, _rowDim1, _colDim1>const&);
		template<class R>auto operator+(R const&)const;
		template<class R>auto operator-(R const&)const;
		template<class R>auto operator*(R const&)const;
		template<class R>auto operator/(R const&)const;
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>auto operator+(mat<R, _rowDim1, _colDim1>const&)const;
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>auto operator-(mat<R, _rowDim1, _colDim1>const&)const;
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>auto operator*(mat<R, _rowDim1, _colDim1>const&)const;
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>auto operator/(mat<R, _rowDim1, _colDim1>const&)const;
		//vec transfer
		template<class R, unsigned int _dim>auto operator,(vec<R, _dim>const&)const;
		//mat transfer
		template<class R, unsigned int _rowDim1, unsigned int _colDim1>auto operator,(mat<R, _rowDim1, _colDim1>const&)const;
		//id
		static auto id();
		static auto id(T);
		//det
		T det();
		T det(int, int);
		//inverse
		auto operator~();
		//transposition
		mat<T, _colDim, _rowDim> operator!()const;
		//print
		void print()const;
		void printInfo(char const*)const;
	};

	template<class T>vec3<T>eulerAngle(vec3<T>const& a)
	{
		CheckNumType(T);
		T c_phi(cos(a.data[0])), s_phi(sin(a.data[0]));
		T c_theta(cos(a.data[1])), s_theta(sin(a.data[1]));
		T c_psi(cos(a.data[2])), s_psi(sin(a.data[2]));
		return vec3<T>
		{
			c_phi* c_psi - s_phi * c_theta * s_psi,
				c_psi* s_phi + c_phi * c_theta * s_psi,
				s_theta* s_psi
		};
	}
	//==============================================vec====================================
	template<class T, unsigned int _dim>										inline vec<T, _dim>::vec()
		:
		data{ 0 }
	{
	}
	template<class T, unsigned int _dim>template<class R>						inline vec<T, _dim>::vec(R const& a)
	{
		CheckNumType(R);
		for (auto& d : data)d = (T)a;
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline vec<T, _dim>::vec(vec<R, _dim1>const& a)
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		static constexpr unsigned int const DifferDim = Differ<unsigned int, _dim, _dim1>::value;
		if constexpr (IsSameType<T, R>::value)
		{
			memcpy(data, a.data, sizeof(T) * MinDim);
			if constexpr (_dim > _dim1)
			{
				memset(data + MinDim, 0, sizeof(T) * DifferDim);
			}
		}
		else
		{
			if constexpr (_dim <= _dim1)
			{
				for (int c1 = 0; c1 < _dim; c1++)data[c1] = a.data[c1];
			}
			else
			{
				for (int c1 = 0; c1 < MinDim; c1++)data[c1] = a.data[c1];
				memset(data + MinDim, 0, sizeof(T) * DifferDim);
			}
		}
	}
	template<class T, unsigned int _dim>										inline vec<T, _dim>::vec(std::initializer_list<T>const& a)
	{
		if (_dim > a.size())
		{
			memcpy(data, a.begin(), sizeof(T) * _dim);
			memset(data + a.size(), 0, sizeof(T) * (_dim - a.size()));
		}
		else
			memcpy(data, a.begin(), sizeof(T) * _dim);

	}

	template<class T, unsigned int _dim>										inline T& vec<T, _dim>::operator[](int a)
	{
		return data[a];
	}



	template<class T, unsigned int _dim>template<class R>						inline vec<T, _dim>& vec<T, _dim>::operator= (R const& a)
	{
		CheckNumType(T);
		CheckNumType(R);
		for (T& d : data)d = (T)a;
		return *this;
	}
	template<class T, unsigned int _dim>template<class R>						inline vec<T, _dim>& vec<T, _dim>::operator+=(R const& a)
	{
		CheckNumType(T);
		CheckNumType(R);
		for (T& d : data)d += (T)a;
		return *this;
	}
	template<class T, unsigned int _dim>template<class R>						inline vec<T, _dim>& vec<T, _dim>::operator-=(R const& a)
	{
		CheckNumType(T);
		CheckNumType(R);
		for (T& d : data)d -= (T)a;
		return *this;
	}
	template<class T, unsigned int _dim>template<class R>						inline vec<T, _dim>& vec<T, _dim>::operator*=(R const& a)
	{
		CheckNumType(T);
		CheckNumType(R);
		for (T& d : data)d *= (T)a;
		return *this;
	}
	template<class T, unsigned int _dim>template<class R>						inline vec<T, _dim>& vec<T, _dim>::operator/=(R const& a)
	{
		CheckNumType(T);
		CheckNumType(R);
		for (T& d : data)d /= (T)a;
		return *this;
	}
	template<class T, unsigned int _dim>										inline bool vec<T, _dim>::operator==(vec<T, _dim> const& a)const
	{
		return !memcmp(data, a.data, sizeof(data));
	}
	template<class T, unsigned int _dim>template<class R>						inline bool vec<T, _dim>::operator==(R const& a)const
	{
		CheckNumType(T);
		CheckNumType(R);
		for (T& d : data)
			if (d != T(a))
				return false;
		return true;
	}
	template<class T, unsigned int _dim>template<class R>						inline bool vec<T, _dim>::operator!=(R const& a)const
	{
		CheckNumType(T);
		CheckNumType(R);
		for (T const& d : data)
			if (d != T(a))
				return true;
		return false;
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline vec<T, _dim>& vec<T, _dim>::operator =(vec<R, _dim1>const& a)
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		static constexpr unsigned int const DifferDim = Differ<unsigned int, _dim, _dim1>::value;
		if constexpr (IsSameType<T, R>::value)
		{
			memcpy(data, a.data, sizeof(T) * MinDim);
		}
		else
		{
			for (int c1 = 0; c1 < MinDim; c1++)data[c1] = (T)a.data[c1];
		}
		return *this;
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline vec<T, _dim>& vec<T, _dim>::operator+=(vec<R, _dim1>const& a)
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		for (int c1 = 0; c1 < MinDim; c1++)data[c1] += (T)a.data[c1];
		return *this;
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline vec<T, _dim>& vec<T, _dim>::operator-=(vec<R, _dim1>const& a)
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		for (int c1 = 0; c1 < MinDim; c1++)data[c1] -= (T)a.data[c1];
		return *this;
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline vec<T, _dim>& vec<T, _dim>::operator*=(vec<R, _dim1>const& a)
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		for (int c1 = 0; c1 < MinDim; c1++)data[c1] *= (T)a.data[c1];
		return *this;
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline vec<T, _dim>& vec<T, _dim>::operator/=(vec<R, _dim1>const& a)
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		for (int c1 = 0; c1 < MinDim; c1++)data[c1] /= (T)a.data[c1];
		return *this;
	}
	template<class T, unsigned int _dim>inline vec<T, _dim>vec<T, _dim>::operator-()const
	{
		vec<T, _dim> temp;
		for (int c0(0); c0 < _dim; ++c0)temp.data[c0] = -data[c0];
		return temp;
	}
	template<class T, unsigned int _dim>template<class R>						inline auto vec<T, _dim>::operator+(R const& a)const
	{
		CheckNumType(R);
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _dim>temp{ *this };
		for (HigherType& d : temp.data)d += HigherType(a);
		return temp;
	}
	template<class T, unsigned int _dim>template<class R>						inline auto vec<T, _dim>::operator-(R const& a)const
	{
		CheckNumType(R);
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _dim>temp{ *this };
		for (HigherType& d : temp.data)d -= HigherType(a);
		return temp;
	}
	template<class T, unsigned int _dim>template<class R>						inline auto vec<T, _dim>::operator*(R const& a)const
	{
		CheckNumType(R);
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _dim>temp{ *this };
		for (HigherType& d : temp.data)d *= HigherType(a);
		return temp;
	}
	template<class T, unsigned int _dim>template<class R>						inline auto vec<T, _dim>::operator/(R const& a)const
	{
		CheckNumType(R);
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _dim>temp{ *this };
		for (HigherType& d : temp.data)d /= HigherType(a);
		return temp;
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline auto vec<T, _dim>::operator+(vec<R, _dim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MaxDim = Max<unsigned int, _dim, _dim1>::value;
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		if constexpr (MaxDim == _dim)
		{
			vec<HigherType, MaxDim>temp{ *this };
			for (int c1 = 0; c1 < MinDim; c1++)temp.data[c1] += (HigherType)a.data[c1];
			return temp;
		}
		else
		{
			vec<HigherType, MaxDim>temp{ a };
			for (int c1 = 0; c1 < MinDim; c1++)temp.data[c1] += (HigherType)data[c1];
			return temp;
		}
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline auto vec<T, _dim>::operator-(vec<R, _dim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MaxDim = Max<unsigned int, _dim, _dim1>::value;
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		if constexpr (MaxDim == _dim)
		{
			vec<HigherType, MaxDim>temp{ *this };
			for (int c1 = 0; c1 < MinDim; c1++)temp.data[c1] -= (HigherType)a.data[c1];
			return temp;
		}
		else
		{
			vec<HigherType, MaxDim>temp{ a };
			for (int c1 = 0; c1 < MinDim; c1++)temp.data[c1] -= (HigherType)data[c1];
			return temp;
		}
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline auto vec<T, _dim>::operator*(vec<R, _dim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MaxDim = Max<unsigned int, _dim, _dim1>::value;
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		if constexpr (MaxDim == _dim)
		{
			vec<HigherType, MaxDim>temp{ *this };
			for (int c1 = 0; c1 < MinDim; c1++)temp.data[c1] *= (HigherType)a.data[c1];
			return temp;
		}
		else
		{
			vec<HigherType, MaxDim>temp{ a };
			for (int c1 = 0; c1 < MinDim; c1++)temp.data[c1] *= (HigherType)data[c1];
			return temp;
		}
	}
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline auto vec<T, _dim>::operator/(vec<R, _dim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MaxDim = Max<unsigned int, _dim, _dim1>::value;
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		if constexpr (MaxDim == _dim)
		{
			vec<HigherType, MaxDim>temp{ *this };
			for (int c1 = 0; c1 < MinDim; c1++)temp.data[c1] /= (HigherType)a.data[c1];
			return temp;
		}
		else
		{
			vec<HigherType, MaxDim>temp{ a };
			for (int c1 = 0; c1 < MinDim; c1++)temp.data[c1] /= (HigherType)data[c1];
			return temp;
		}
	}
	//append
	template<class T, unsigned int _dim, class R>inline auto operator+(R const& a, vec<T, _dim>const& b)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _dim>temp{ b };
		for (HigherType& d : temp.data)d += HigherType(a);
		return temp;
	}
	template<class T, unsigned int _dim, class R>inline auto operator-(R const& a, vec<T, _dim>const& b)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _dim>temp{ b };
		for (HigherType& d : temp.data)d = HigherType(a) - d;
		return temp;
	}
	template<class T, unsigned int _dim, class R>inline auto operator*(R const& a, vec<T, _dim>const& b)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _dim>temp{ b };
		for (HigherType& d : temp.data)d *= HigherType(a);
		return temp;
	}
	template<class T, unsigned int _dim, class R>inline auto operator/(R const& a, vec<T, _dim>const& b)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _dim>temp{ b };
		for (HigherType& d : temp.data)d = HigherType(a) / d;
		return temp;
	}
	//dot
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline auto vec<T, _dim>::operator,(vec<R, _dim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MinDim = Min<unsigned int, _dim, _dim1>::value;
		HigherType temp(0);
		for (int c1 = 0; c1 < MinDim; ++c1)temp += (HigherType)data[c1] * (HigherType)a.data[c1];
		return temp;
	}
	//cross
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>	inline auto vec<T, _dim>::operator|(vec<R, _dim1>const& a)
	{
		static_assert(_dim >= 3 && _dim1 >= 3, "Cannot cross vec whose dimension is lower than 3!");
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		return vec3<HigherType>
		{
			(HigherType)data[1] * (HigherType)a.data[2] - (HigherType)data[2] * (HigherType)a.data[1],
				(HigherType)data[2] * (HigherType)a.data[0] - (HigherType)data[0] * (HigherType)a.data[2],
				(HigherType)data[0] * (HigherType)a.data[1] - (HigherType)data[1] * (HigherType)a.data[0]
		};
	}
	template<class T, unsigned int _dim>inline mat3<T> vec<T, _dim>::crossMat()
	{
		static_assert(_dim >= 3, "Cannot build with vec whose dimsion is less than 3!");
		return mat3<T>
		{
			{ 0, -data[2], data[1] },
			{ data[2], 0, -data[0] },
			{ -data[1], data[0], 0 },
		};
	}
	//tensor product
	template<class T, unsigned int _dim>template<class R, unsigned int _dim1>inline auto vec<T, _dim>::operator^(vec<R, _dim1>const& a)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _dim, _dim1>temp;
		int c1{ 0 };
		for (auto& d : temp.rowVec)
			d = HigherType(data[c1++]) * a;
		return temp;
	}
	//transfer
	template<class T, unsigned int _dim>template<class R, unsigned int _rowDim, unsigned int _colDim>auto vec<T, _dim>::operator()(mat<R, _rowDim, _colDim>const& a)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		constexpr unsigned int MinDim = Min<unsigned int, _dim, _rowDim>::value;
		vec<HigherType, _colDim>temp;
		for (int c1(0); c1 < MinDim; ++c1)
			temp += a.rowVec[c1] * HigherType(data[c1]);
		return temp;
	}
	//rotate
	template<class T, unsigned int _dim>template<class R, class Y, unsigned int _dim1>inline auto vec<T, _dim>::operator()(vec<R, _dim1>const& a, Y b)
	{
		static_assert(NumType<Y>::value, "Wrong NumType!");
		static_assert(_dim >= 3, "Axis vec dimsion must be more than 2!");
		if (this->square(3) == 0)return vec3<double>();
		vec3<double>n{ *this };
		vec3<double>r{ a };
		n.normaliaze();
		double c{ cos(b) };
		return ((1.0 - c) * (n, r))* n + c * r + sin(b) * (n | r);
	}
	template<class T, unsigned int _dim>template<class Y>inline mat3<double> vec<T, _dim>::rotMat(Y a)
	{
		static_assert(NumType<Y>::value, "Wrong NumType!");
		static_assert(_dim >= 3, "Axis vec dimsion must be more than 2!");
		if (this->square(3) == 0)return mat3<double>();
		vec3<double>n{ *this };
		n.normaliaze();
		double c{ cos(a) };
		return (1 - c)* (n ^ n) + mat3<double>::id(c) + sin(a) * n.crossMat();
	}
	//note: rotate the mat(other coordinate to this) in this coordinate.
	template<class T, unsigned int _dim>template<class R, class Y, unsigned int _rowDim1, unsigned int _colDim1>auto& vec<T, _dim>::operator()(mat<R, _rowDim1, _colDim1>* a, Y b)
	{
		CheckNumType(Y);
		static_assert(_dim >= 3, "Axis vec dimsion must be more than 2!");
		static constexpr const unsigned int MinDim = Min<unsigned int, _rowDim1, _colDim1>::value;
		static_assert(MinDim >= 3, "Cannot rotate mat under 3 order!");
		auto r{ rotMat(b) };
		return *a = (r, *a);
	}
	template<class T, unsigned int _dim>template<class R, class Y, unsigned int _rowDim1, unsigned int _colDim1>auto vec<T, _dim>::operator()(mat<R, _rowDim1, _colDim1>const& a, Y b)
	{
		CheckNumType(Y);
		static_assert(_dim >= 3, "Axis vec dimsion must be more than 2!");
		static constexpr const unsigned int MinDim = Min<unsigned int, _rowDim1, _colDim1>::value;
		static_assert(MinDim >= 3, "Cannot rotate mat under 3 order!");
		auto r{ rotMat(b) };
		return  (r, a);
	}
	//square
	template<class T, unsigned int _dim>inline T vec<T, _dim>::square()const
	{
		T t{ 0 };
		for (auto& d : data)t += d * d;
		return t;
	}
	template<class T, unsigned int _dim>inline T vec<T, _dim>::square(int a)const
	{
		T t{ 0 };
		a = a < _dim ? a : _dim;
		for (int c0{ 0 }; c0 < a; ++c0)t += data[c0] * data[c0];
		return t;
	}
	//length
	template<class T, unsigned int _dim>inline T vec<T, _dim>::length()const
	{
		return T(sqrt(square()));
	}
	template<class T, unsigned int _dim>inline T vec<T, _dim>::length(int a)const
	{
		return T(square(a));
	}
	//max
	template<class T, unsigned int _dim>inline T vec<T, _dim>::max()const
	{
		T t{ data[0] };
		for (auto& d : data)if (d > t)t = d;
		return t;
	}
	template<class T, unsigned int _dim>inline T vec<T, _dim>::max(int a)const
	{
		T t{ 0 };
		a = a < _dim ? a : _dim;
		for (int c0{ 0 }; c0 < a; ++c0)if (data[c0] > t)t = data[c0];
		return t;
	}
	//min
	template<class T, unsigned int _dim>inline T vec<T, _dim>::min()const
	{
		T t{ data[0] };
		for (auto& d : data)if (d < t)t = d;
		return t;
	}
	template<class T, unsigned int _dim>inline T vec<T, _dim>::min(int a)const
	{
		T t{ 0 };
		a = a < _dim ? a : _dim;
		for (int c0{ 0 }; c0 < a; ++c0)if (data[c0] < t)t = data[c0];
		return t;
	}
	//normaliaze
	template<class T, unsigned int _dim>inline vec<T, _dim>& vec<T, _dim>::normaliaze()
	{
		double temp(0);
		if constexpr (NumType<T>::serial < 11)
			for (T const& d : data)temp += (unsigned long long)d * (unsigned long long)d;
		else
			for (T const& d : data)temp += d * d;
		temp = sqrt(temp);
		for (T& d : data)d /= temp;
		return *this;
	}
	template<class T, unsigned int _dim>inline vec<T, _dim>& vec<T, _dim>::normaliaze(int a)
	{
		double temp(0);
		if constexpr (NumType<T>::serial < 11)
			for (int c0{ 0 }; c0 < a; ++c0)temp += (unsigned long long)data[c0] * (unsigned long long)data[c0];
		else
			for (int c0{ 0 }; c0 < a; ++c0)temp += data[c0] * data[c0];
		temp = sqrt(temp);
		for (int c0{ 0 }; c0 < a; ++c0)data[c0] /= temp;
		return *this;
	}
	//print
	template<class T, unsigned int _dim>inline void vec<T, _dim>::print()const
	{
		printf("[");
		char str[8] = { ", " };
		strcat(str, NumType<T>::printInfo);
		printf(NumType<T>::printInfo, data[0]);
		for (int c0(1); c0 < _dim; ++c0)
			printf(str, data[c0]);
		printf("]");
	}
	template<class T, unsigned int _dim>inline void vec<T, _dim>::printInfo(char const* a) const
	{
		::printf("%s", a);
		print();
	}
	//==============================================mat====================================
	template<class T, unsigned int _rowDim, unsigned int _colDim>	inline mat<T, _rowDim, _colDim>::mat()
		:array{ 0 }
	{
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>	inline mat<T, _rowDim, _colDim>::mat(std::initializer_list<row>const& a)
	{
		memcpy(array, a.begin(), sizeof(row) * a.size());
		if (_rowDim > a.size())
			memset(rowVec + a.size(), 0, sizeof(row) * (_rowDim - a.size()));
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>
	template<class R>												inline mat<T, _rowDim, _colDim>::mat(R const& a)
	{
		CheckNumType(R);
		for (auto& d0 : array)for (auto& d1 : d0)d1 = (T)a;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>
	template<class R, unsigned int _rowDim1, unsigned int _colDim1>	inline mat<T, _rowDim, _colDim>::mat(mat<R, _rowDim1, _colDim1>const& a)
	{
		static constexpr unsigned int const MaxRow = Max<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MaxCol = Max<unsigned int, _colDim, _colDim1>::value;
		static constexpr unsigned int const MinCol = Min<unsigned int, _colDim, _colDim1>::value;
		static constexpr unsigned int const DifferRow = Differ<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const DifferCol = Differ<unsigned int, _colDim, _colDim1>::value;
		if constexpr (IsSameType<T, R>::value)
		{
			if constexpr (_colDim == _colDim1)
			{
				memcpy(array, a.array, MinRow * sizeof(row));
				if constexpr (MinRow != _rowDim)
				{
					memset(rowVec + MinRow, 0, DifferRow * sizeof(row));
				}
			}
			else
			{
				for (int c1 = 0; c1 < MinRow; c1++)
				{
					rowVec[c1] = a.rowVec[c1];
					if constexpr (MinCol != _colDim)
					{
						memset(array[c1] + MinCol, 0, DifferCol * sizeof(T));
					}
				}
				memset(rowVec + MinRow, 0, DifferRow * sizeof(row));
			}
		}
		else
		{
			for (int c1 = 0; c1 < MinRow; c1++)
			{
				rowVec[c1] = a.rowVec[c1];
				if constexpr (MinCol != _colDim)
				{
					memset(&array[c1][MinCol], 0, DifferCol * sizeof(T));
				}
			}
			if constexpr (_colDim < MaxRow)
				memset(array[MinRow], 0, DifferRow * sizeof(row));
		}
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>	inline vec<T, _colDim>& mat<T, _rowDim, _colDim>::operator[](int a)
	{
		return rowVec[a];
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>	inline vec<T, _rowDim> mat<T, _rowDim, _colDim>::column(int a)const
	{
		vec<T, _rowDim>temp;
		int c0{ 0 };
		for (T& d : temp.data)d = array[c0++][a];
		return temp;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>	inline bool mat<T, _rowDim, _colDim>::operator==(mat<T, _rowDim, _colDim> const& a) const
	{
		return !memcmp(array, a.array, sizeof(array));
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>
	template<class R, unsigned int _dim>							inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::setCol(vec<R, _dim>& a, unsigned int b)
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _rowDim, _dim>::value;
		for (int c0{ 0 }; c0 < MinDim; c0++)array[c0][b] = a.data[c0];
		return *this;
	}

	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator= (R const& a)
	{
		CheckNumType(R);
		for (auto& d0 : array)
			for (T& d1 : d0)d1 = (T)a;
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator+=(R const& a)
	{
		CheckNumType(R);
		for (auto& d0 : array)
			for (T& d1 : d0)d1 += (T)a;
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator-=(R const& a)
	{
		CheckNumType(R);
		for (auto& d0 : array)
			for (T& d1 : d0)d1 -= (T)a;
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator*=(R const& a)
	{
		CheckNumType(R);
		for (auto& d0 : array)
			for (T& d1 : d0)d1 *= (T)a;
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator/=(R const& a)
	{
		CheckNumType(R);
		for (auto& d0 : array)
			for (T& d1 : d0)d1 /= (T)a;
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator= (mat<R, _rowDim1, _colDim1>const& a)
	{
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		if constexpr (IsSameType<T, R>::value)
		{
			if constexpr (_colDim == _colDim1)
				memcpy(array, a.array, MinRow * sizeof(row));
			else
				for (int c1 = 0; c1 < MinRow; c1++)
					rowVec[c1] = a.rowVec[c1];
		}
		else
			for (int c1 = 0; c1 < MinRow; c1++)
				rowVec[c1] = a.rowVec[c1];
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator+=(mat<R, _rowDim1, _colDim1>const& a)
	{
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		if constexpr (IsSameType<T, R>::value)
		{
			if constexpr (_colDim == _colDim1)
				memcpy(array, a.array, MinRow * sizeof(row));
			else
				for (int c1 = 0; c1 < MinRow; c1++)
					rowVec[c1] += a.rowVec[c1];
		}
		else
			for (int c1 = 0; c1 < MinRow; c1++)
				rowVec[c1] += a.rowVec[c1];
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator-=(mat<R, _rowDim1, _colDim1>const& a)
	{
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		if constexpr (IsSameType<T, R>::value)
		{
			if constexpr (_colDim == _colDim1)
				memcpy(array, a.array, MinRow * sizeof(row));
			else
				for (int c1 = 0; c1 < MinRow; c1++)
					rowVec[c1] -= a.rowVec[c1];

		}
		else
			for (int c1 = 0; c1 < MinRow; c1++)
				rowVec[c1] -= a.rowVec[c1];
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator*=(mat<R, _rowDim1, _colDim1>const& a)
	{
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		if constexpr (IsSameType<T, R>::value)
		{
			if constexpr (_colDim == _colDim1)
				memcpy(array, a.array, MinRow * sizeof(row));
			else
				for (int c1 = 0; c1 < MinRow; c1++)
					rowVec[c1] *= a.rowVec[c1];
		}
		else
			for (int c1 = 0; c1 < MinRow; c1++)
				rowVec[c1] *= a.rowVec[c1];
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline mat<T, _rowDim, _colDim>& mat<T, _rowDim, _colDim>::operator/=(mat<R, _rowDim1, _colDim1>const& a)
	{
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		if constexpr (IsSameType<T, R>::value)
		{
			if constexpr (_colDim == _colDim1)
				memcpy(array, a.array, MinRow * sizeof(row));
			else
				for (int c1 = 0; c1 < MinRow; c1++)
					rowVec[c1] /= a.rowVec[c1];
		}
		else
			for (int c1 = 0; c1 < MinRow; c1++)
				rowVec[c1] /= a.rowVec[c1];
		return *this;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline auto mat<T, _rowDim, _colDim>::operator+(R const& a)const
	{
		CheckNumType(R);
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim>temp{ *this };
		for (auto& d0 : temp.array)
			for (auto& d1 : d0)
				d1 += (HigherType)a;
		return temp;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline auto mat<T, _rowDim, _colDim>::operator-(R const& a)const
	{
		CheckNumType(R);
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim>temp{ *this };
		for (auto& d0 : temp.array)
			for (auto& d1 : d0)
				d1 -= (HigherType)a;
		return temp;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline auto mat<T, _rowDim, _colDim>::operator*(R const& a)const
	{
		CheckNumType(R);
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim>temp{ *this };
		for (auto& d0 : temp.array)
			for (auto& d1 : d0)
				d1 *= (HigherType)a;
		return temp;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R>inline auto mat<T, _rowDim, _colDim>::operator/(R const& a)const
	{
		CheckNumType(R);
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim>temp{ *this };
		for (auto& d0 : temp.array)
			for (auto& d1 : d0)
				d1 /= (HigherType)a;
		return temp;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline auto mat<T, _rowDim, _colDim>::operator+(mat<R, _rowDim1, _colDim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MaxRow = Max<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MaxCol = Max<unsigned int, _colDim, _colDim1>::value;
		if constexpr (MaxRow == _rowDim)
		{
			mat<HigherType, MaxRow, MaxCol>temp{ *this };
			for (int c1 = 0; c1 < MinRow; c1++)temp.rowVec[c1] += a.rowVec[c1];
			return temp;
		}
		else if constexpr (MaxRow == _rowDim1)
		{
			mat<HigherType, MaxRow, MaxCol>temp{ a };
			for (int c1 = 0; c1 < MinRow; c1++)temp.rowVec[c1] += rowVec[c1];
			return temp;
		}
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline auto mat<T, _rowDim, _colDim>::operator-(mat<R, _rowDim1, _colDim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MaxRow = Max<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MaxCol = Max<unsigned int, _colDim, _colDim1>::value;
		if constexpr (MaxRow == _rowDim)
		{
			mat<HigherType, MaxRow, MaxCol>temp{ *this };
			for (int c1 = 0; c1 < MinRow; c1++)temp.rowVec[c1] -= a.rowVec[c1];
			return temp;
		}
		else if constexpr (MaxRow == _rowDim1)
		{
			mat<HigherType, MaxRow, MaxCol>temp{ a };
			for (int c1 = 0; c1 < MinRow; c1++)temp.rowVec[c1] -= rowVec[c1];
			return temp;
		}
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline auto mat<T, _rowDim, _colDim>::operator*(mat<R, _rowDim1, _colDim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MaxRow = Max<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MaxCol = Max<unsigned int, _colDim, _colDim1>::value;
		if constexpr (MaxRow == _rowDim)
		{
			mat<HigherType, MaxRow, MaxCol>temp{ *this };
			for (int c1 = 0; c1 < MinRow; c1++)temp.rowVec[c1] *= a.rowVec[c1];
			return temp;
		}
		else if constexpr (MaxRow == _rowDim1)
		{
			mat<HigherType, MaxRow, MaxCol>temp{ a };
			for (int c1 = 0; c1 < MinRow; c1++)temp.rowVec[c1] *= rowVec[c1];
			return temp;
		}
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline auto mat<T, _rowDim, _colDim>::operator/(mat<R, _rowDim1, _colDim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		static constexpr unsigned int const MaxRow = Max<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MinRow = Min<unsigned int, _rowDim, _rowDim1>::value;
		static constexpr unsigned int const MaxCol = Max<unsigned int, _colDim, _colDim1>::value;
		if constexpr (MaxRow == _rowDim)
		{
			mat<HigherType, MaxRow, MaxCol>temp{ *this };
			for (int c1 = 0; c1 < MinRow; c1++)temp.rowVec[c1] /= a.rowVec[c1];
			return temp;
		}
		else if constexpr (MaxRow == _rowDim1)
		{
			mat<HigherType, MaxRow, MaxCol>temp{ a };
			for (int c1 = 0; c1 < MinRow; c1++)temp.rowVec[c1] /= rowVec[c1];
			return temp;
		}
	}
	//append
	template<class T, unsigned int _rowDim, unsigned int _colDim, class R>inline auto operator+(R const& a, mat<T, _rowDim, _colDim>const& b)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim>temp{ b };
		for (auto& d0 : temp.array)
			for (auto& d1 : d0)
				d1 += (HigherType)a;
		return temp;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim, class R>inline auto operator-(R const& a, mat<T, _rowDim, _colDim>const& b)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim>temp{ b };
		for (auto& d0 : temp.array)
			for (auto& d1 : d0)
				d1 = (HigherType)a - d1;
		return temp;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim, class R>inline auto operator*(R const& a, mat<T, _rowDim, _colDim>const& b)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim>temp{ b };
		for (auto& d0 : temp.array)
			for (auto& d1 : d0)
				d1 *= (HigherType)a;
		return temp;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim, class R>inline auto operator/(R const& a, mat<T, _rowDim, _colDim>const& b)
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim>temp{ b };
		for (auto& d0 : temp.array)
			for (auto& d1 : d0)
				d1 = (HigherType)a / d1;
		return temp;
	}
	//vec transfer
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _dim>inline auto mat<T, _rowDim, _colDim>::operator,(vec<R, _dim> const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		vec<HigherType, _rowDim>temp;
		int c1{ 0 };
		for (auto& d : temp.data)d = (rowVec[c1++], a);
		return temp;
	}
	//mat transfer
	template<class T, unsigned int _rowDim, unsigned int _colDim>template<class R, unsigned int _rowDim1, unsigned int _colDim1>inline auto mat<T, _rowDim, _colDim>::operator,(mat<R, _rowDim1, _colDim1>const& a)const
	{
		using HigherType = typename GetNumType<HigherNumTypeTable[NumType<T>::serial][NumType<R>::serial]>::Result;
		mat<HigherType, _rowDim, _colDim1>temp;
		int c1{ 0 };
		for (auto& d0 : temp.rowVec)
		{
			for (int c2(0); c2 < _rowDim; ++c2)
				d0 += HigherType(array[c1][c2]) * a.rowVec[c2];
			++c1;
		}
		return temp;
	}
	//id
	template<class T, unsigned int _rowDim, unsigned int _colDim>inline auto mat<T, _rowDim, _colDim>::id()
	{
		static_assert(_rowDim == _colDim, "Invalid operation!");
		mat<T, _rowDim, _colDim>t;
		int c0{ 0 };
		for (auto& d : t.array)d[c0++] = 1;
		return t;
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>inline auto mat<T, _rowDim, _colDim>::id(T a)
	{
		static_assert(_rowDim == _colDim, "Invalid operation!");
		mat<T, _rowDim, _colDim>t;
		int c0{ 0 };
		for (auto& d : t.array)d[c0++] = a;
		return t;
	}
	//det
	template<class T, unsigned int _rowDim, unsigned int _colDim>inline T mat<T, _rowDim, _colDim>::det()
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _rowDim, _colDim>::value;
		static_assert(MinDim <= 4, "Support det only for mat under 5 order!");
		if constexpr (MinDim == 1)
		{
			return array[0][0];
		}
		else if constexpr (MinDim == 2)
		{
			return array[0][0] * array[1][1] - array[1][0] * array[0][1];
		}
		else if constexpr (MinDim == 3)
		{
			return
				array[0][0] * (array[1][1] * array[2][2] - array[2][1] * array[1][2]) +
				array[1][0] * (array[2][1] * array[0][2] - array[0][1] * array[2][2]) +
				array[2][0] * (array[0][1] * array[1][2] - array[0][2] * array[1][1]);
		}
		else
		{
			return
				array[0][0] * (
					array[1][1] * (array[2][2] * array[3][3] - array[3][2] * array[2][3]) +
					array[2][1] * (array[3][2] * array[1][3] - array[1][2] * array[3][3]) +
					array[3][1] * (array[1][2] * array[2][3] - array[1][3] * array[2][2])) -
				array[1][0] * (
					array[0][1] * (array[2][2] * array[3][3] - array[3][2] * array[2][3]) +
					array[2][1] * (array[3][2] * array[0][3] - array[0][2] * array[3][3]) +
					array[3][1] * (array[0][2] * array[2][3] - array[0][3] * array[2][2])) +
				array[2][0] * (
					array[0][1] * (array[1][2] * array[3][3] - array[3][2] * array[1][3]) +
					array[1][1] * (array[3][2] * array[0][3] - array[0][2] * array[3][3]) +
					array[3][1] * (array[0][2] * array[1][3] - array[0][3] * array[1][2])) -
				array[3][0] * (
					array[0][1] * (array[1][2] * array[2][3] - array[2][2] * array[1][3]) +
					array[1][1] * (array[2][2] * array[0][3] - array[0][2] * array[2][3]) +
					array[2][1] * (array[0][2] * array[1][3] - array[0][3] * array[1][2]));
		}
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>inline T mat<T, _rowDim, _colDim>::det(int i, int j)
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _rowDim, _colDim>::value;
		static_assert(MinDim >= 1, "Cannot get the algebraic complement!");
		static_assert(MinDim <= 4, "Support only mat under 5 order!");
#define _move( q , p ) ( ( q < p ) ? q : q + 1 )
		if constexpr (MinDim == 2)
		{
			return pow(-1, i + j)* array[_move(i, 0)][_move(j, 0)];
		}
		else if constexpr (MinDim == 3)
		{
			return
				pow(-1, i + j)* (
					array[_move(0, i)][_move(0, j)] * array[_move(1, i)][_move(1, j)] -
					array[_move(1, i)][_move(0, j)] * array[_move(0, i)][_move(1, j)]
					);
		}
		else
		{
			return
				pow(-1.0, i + j)* (
					array[_move(0, i)][_move(0, j)] * (array[_move(1, i)][_move(1, j)] * array[_move(2, i)][_move(2, j)] - array[_move(2, i)][_move(1, j)] * array[_move(1, i)][_move(2, j)]) +
					array[_move(1, i)][_move(0, j)] * (array[_move(2, i)][_move(1, j)] * array[_move(0, i)][_move(2, j)] - array[_move(0, i)][_move(1, j)] * array[_move(2, i)][_move(2, j)]) +
					array[_move(2, i)][_move(0, j)] * (array[_move(0, i)][_move(1, j)] * array[_move(1, i)][_move(2, j)] - array[_move(0, i)][_move(2, j)] * array[_move(1, i)][_move(1, j)])
					);
		}
#undef _move( q, p )
	}
	//inverse
	template<class T, unsigned int _rowDim, unsigned int _colDim>inline auto mat<T, _rowDim, _colDim>::operator~()
	{
		static constexpr unsigned int const MinDim = Min<unsigned int, _rowDim, _colDim>::value;
		static_assert(MinDim <= 4, "Support inverse only for mat under 5 order!");
		if constexpr (MinDim == 1)
		{
			return mat<T, MinDim, MinDim>{1 / array[0][0]};
		}
		else if constexpr (MinDim == 2)
		{
			T d{ this->det() };
			return mat<T, MinDim, MinDim>
			{
				{ array[1][1] / d, array[0][1] / d },
				{ array[1][0] / d, array[0][0] / d }
			};
		}
		else if constexpr (MinDim == 3)
		{
			T d{ this->det() };
			return mat<T, MinDim, MinDim>
			{
				{ this->det(0, 0) / d, this->det(1, 0) / d, this->det(2, 0) / d },
				{ this->det(0, 1) / d, this->det(1, 1) / d, this->det(2, 1) / d },
				{ this->det(0, 2) / d, this->det(1, 2) / d, this->det(2, 2) / d }
			};
		}
		else if constexpr (MinDim == 4)
		{
			T d{ this->det() };
			return mat<T, MinDim, MinDim>
			{
				{ this->det(0, 0) / d, this->det(1, 0) / d, this->det(2, 0) / d, this->det(3, 0) / d },
				{ this->det(0, 1) / d, this->det(1, 1) / d, this->det(2, 1) / d, this->det(3, 1) / d },
				{ this->det(0, 2) / d, this->det(1, 2) / d, this->det(2, 2) / d, this->det(3, 2) / d },
				{ this->det(0, 3) / d, this->det(1, 3) / d, this->det(2, 3) / d, this->det(3, 3) / d }
			};
		}
	}
	//transposition
	template<class T, unsigned int _rowDim, unsigned int _colDim>inline mat<T, _colDim, _rowDim> mat<T, _rowDim, _colDim>::operator!()const
	{
		int c0{ 0 };
		mat<T, _colDim, _rowDim>temp;
		for (auto& d0 : temp.array)
		{
			int c1{ 0 };
			for (auto& d1 : d0)d1 = array[c1++][c0];
			++c0;
		}
		return temp;
	}
	//print
	template<class T, unsigned int _rowDim, unsigned int _colDim>inline void mat<T, _rowDim, _colDim>::print()const
	{
		printf("[\n");
		for (row const& d : rowVec)
		{
			printf("\t");
			d.print();
			printf("\n");
		}
		printf("]\n");
	}
	template<class T, unsigned int _rowDim, unsigned int _colDim>inline void mat<T, _rowDim, _colDim>::printInfo(char const* a)const
	{
		::printf("%s", a);
		print();
	}
	//==============================================test====================================
	void testVecMat()
	{
		vec<int, 3>va{ 1,2,3 };
		va.printInfo("va:		");
		vec<double, 4>vb(va);
		vb.printInfo("vb:		");
		(vb = { 1,2,3,4 }).printInfo("vb = {1,2,3,4}:	");
		vec<double, 3>vc(va);
		vc.printInfo("vc:		");
		::printf("vb[1]:		%lf\n", vb[1]);
		(va += 3).printInfo("va += 3:	");
		(va -= 3).printInfo("va -= 3:	");
		(va *= 3).printInfo("va *= 3:	");
		(va /= 3).printInfo("va /= 3:	");
		(vb += 3).printInfo("vb += 3:	");
		(vb -= 3).printInfo("vb -= 3:	");
		(vb *= 3).printInfo("vb *= 3:	");
		(vb /= 3).printInfo("vb /= 3:	");
		(va = vb).printInfo("va = vb:	");
		(vb = va).printInfo("vb = va:	");
		(vb += va).printInfo("vb += va:	");
		(vb -= va).printInfo("vb -= va:	");
		(vb *= va).printInfo("vb *= va:	");
		(vb /= va).printInfo("vb /= va:	");
		(va += vb).printInfo("va += vb:	");
		(va -= vb).printInfo("va -= vb:	");
		(va *= vb).printInfo("va *= vb:	");
		(va /= vb).printInfo("va /= vb:	");
		(va + vb).printInfo("va + vb:	");
		(va - vb).printInfo("va - vb:	");
		(va * vb).printInfo("va * vb:	");
		(va / vb).printInfo("va / vb:	");
		(va + 3).printInfo("va + 3:		");
		(va - 3).printInfo("va - 3:		");
		(va * 3).printInfo("va * 3:		");
		(va / 3).printInfo("va / 3:		");
		(3 + va).printInfo("3 + va:		");
		(3 - va).printInfo("3 - va:		");
		(3 * va).printInfo("3 * va:		");
		(3 / va).printInfo("3 / va:		");
		(vb + 3).printInfo("vb + 3:		");
		(vb - 3).printInfo("vb - 3:		");
		(vb * 3).printInfo("vb * 3:		");
		(vb / 3).printInfo("vb / 3:		");
		(3 + vb).printInfo("3 + vb:		");
		(3 - vb).printInfo("3 - vb:		");
		(3 * vb).printInfo("3 * vb:		");
		(3 / vb).printInfo("3 / vb:		");
		::printf("(va, vb):	%lf\n", (va, vb));
		(va = { 3,2,1 }).printInfo("va = {3,2,1}	");
		(va | vb).printInfo("va | vb:	");
		::printf("va.square():	%d\n", va.square());
		::printf("vb.square():	%lf\n", vb.square());
		::printf("va.square(2):	%d\n", va.square(2));
		::printf("vb.square(2):	%lf\n", vb.square(2));
		::printf("va.length():	%d\n", va.length());
		::printf("vb.length():	%lf\n", vb.length());
		(vb = { 3,2,1,0 }).printInfo("vb = {3,2,1,0}:	");
		va.normaliaze().printInfo("va.normaliaze:	");
		vb.normaliaze().printInfo("vb.normaliaze:	");
		(vc -= (vb, vc) * vb).printInfo("vc-=(vb,vc)*vb:	");
		vb(vc, Pi).printInfo("vb(vc, Pi):	");
		vb.rotMat(Pi).printInfo("vb.rotMat(Pi):\n");
		(vb.rotMat(Pi), vc).printInfo("(vb.rotMat(Pi),vc):	");
		(va = { 1,2,3 }).printInfo("va = {1,2,3}:	");
		(vc = { 1,2,3 }).printInfo("vc = {1,2,3}:	");
		(va ^ vc).printInfo("vb^vc:\n");
		mat3<double>ma
		{
			{1,2,3},
			{2,4,6},
			{3,6,9}
		};
		ma.printInfo("ma:\n");
		vec3<double>{1, 0, 0}(&ma, Pi).printInfo("vec3<double>{1,0,0}(&ma,Pi):\n");
		vec3<double>{1, 0, 0}(mat3<double>{ { 1, 2, 3 }, { 4,5,6 }, { 7,8,9 }}, Pi).
			printInfo("vec3<double>{1,0,0}(mat3<double>{{1,2,3},{4,5,6},{7,8,9}}, Pi):\n");
		mat3<double>::id(2).printInfo("mat3<double>::id(2):\n");
	}
}


