#pragma once
#include <_Vector.h>

//TODO: 

template<class T>inline bool judgeUp(T* const a, int p, int q)
{
	for (int c0(p); c0 < q - 1; ++c0)
		if (a[c0] > a[c0 + 1])return false;
	return true;
}

template<class T>inline void maxHeap(T* const a, int p, int q)
{
	T l((p << 1) + 1);
	T r((p + 1) << 1);
	if (a[p] < a[l])
	{
		if (r < q && a[l] < a[r])
		{
			T t(a[p]);
			a[p] = a[r];
			a[r] = t;
			if ((r + 1) << 1 < q)
				maxHeap(a, r, q);
		}
		else
		{
			T t(a[p]);
			a[p] = a[l];
			a[l] = t;
			if ((l << 1) + 1 < q)
				maxHeap(a, l, q);
		}
	}
	else
	{
		if (r < q && a[p] < a[r])
		{
			T t(a[p]);
			a[p] = a[r];
			a[r] = t;
			if ((r + 1) << 1 < q)
				maxHeap(a, r, q);
		}
	}
}
template<class T>inline void maxTopHeap(T* const a, int q)
{
	if (1 < q && *a < a[1])
	{
		if (2 < q && a[1] < a[2])
		{
			T t(*a);
			*a = a[2];
			a[2] = t;
			if (5 < q)
				maxHeap(a, 2, q);
		}
		else
		{
			T t(*a);
			*a = a[1];
			a[1] = t;
			if (3 < q)
				maxHeap(a, 1, q);
		}
	}
	else
	{
		if (2 < q && *a < a[2])
		{
			T t(*a);
			*a = a[2];
			a[2] = t;
			if (5 < q)
				maxHeap(a, 2, q);
		}
	}
}
template<class T>inline void buildMaxHeap(T* const a, int q)
{
	for (int c((q >> 1) - 1); c; --c)maxHeap(a, c, q);
	maxTopHeap(a, q);
}
template<class T>inline void heapSort(T* const a, int q)
{
	for (int c((q >> 1) - 1); c; --c)maxHeap(a, c, q);
	maxTopHeap(a, q);
	for (int c(q - 1); c; --c)
	{
		T t(a[c]);
		a[c] = a[0];
		a[0] = t;
		maxTopHeap(a, c);
	}
}

// T must has operator<, operator<=
// All sections are closed intervals.
template<class T>inline void qsort(T* const a, int p, int q)
{
	if (p + 1 < q)
	{
		T& const k(a[p]);
		int m(p + 1), n(p);
		while (++n != q)
			if (a[n] < k) { T t = a[m]; a[m++] = a[n]; a[n] = t; }
		T t = a[m - 1]; a[m - 1] = a[p]; a[p] = t;
		if (p + 2 < m)qsort(a, p, m - 1);
		if (m + 1 < n)qsort(a, m, n);
	}
}

template<class T>struct Interval
{
	static_assert(NumType<T>::value == true, "Non-numeric type not supported yet!");
	union
	{
		T data[2];
		struct
		{
			T a, b;
		};
	};

	Interval()
		:
		a(0),
		b(0)
	{
	}
	template<class R, class S>Interval(R const& _a, S const& _b)
		:
		a(_a),
		b(_b)
	{
		static_assert(NumType<R>::value == true, "Non-numeric type not supported yet!");
		static_assert(NumType<S>::value == true, "Non-numeric type not supported yet!");
	}
	template<class R, class S>Interval(R&& _a, S&& _b)
		:
		a(_a),
		b(_b)
	{
		static_assert(NumType<R>::value == true, "Non-numeric type not supported yet!");
		static_assert(NumType<S>::value == true, "Non-numeric type not supported yet!");
	}
	template<class R, class S>Interval(R const& _a, S&& _b)
		:
		a(_a),
		b(_b)
	{
		static_assert(NumType<R>::value == true, "Non-numeric type not supported yet!");
		static_assert(NumType<S>::value == true, "Non-numeric type not supported yet!");
	}
	template<class R, class S>Interval(R&& _a, S const& _b)
		:
		a(_a),
		b(_b)
	{
		static_assert(NumType<R>::value == true, "Non-numeric type not supported yet!");
		static_assert(NumType<S>::value == true, "Non-numeric type not supported yet!");
	}
	bool operator<(Interval<T> const& s)const
	{
		return a < s.a;
	}
	//If A^B == null, then do nothing and return this.
	Interval<T>	operator+(Interval<T>const& s)const
	{
		if ((s.a < a ? a : s.a) <= (b < s.b ? b : s.b))
		{
			return { a < s.a ? a : s.a, s.b < b ? b : s.b };
		}
		return *this;
	}
	Interval<T>& operator+=(Interval<T>const& s)
	{

		if ((s.a < a ? a : s.a) <= (b < s.b ? b : s.b))
		{
			if (s.a < a) a = s.a;
			if (b < s.b) b = s.b;
		}
		return *this;
	}
	//If A^B == null, return a invalid Interval.
	Interval<T>  operator^(Interval<T>const& s)const
	{
		if ((b < s.b ? b : s.b) < (s.a < a ? a : s.a) || !s.valid())
			return { 1,0 };
		return { s.a < a ? a : s.a, b < s.b ? b : s.b };
	}
	Interval<T>& operator^=(Interval<T>const& s)
	{
		if ((b < s.b ? b : s.b) < (s.a < a ? a : s.a) || !s.valid())
		{
			a = 1;
			b = 0;
		}
		else
		{
			if (a < s.a) a = s.a;
			if (s.b < b) b = s.b;
		}
		return *this;
	}
	IntervalSet<T>  operator^(IntervalSet<T>const&)const;
	Interval<T>  move(T x)const
	{
		return { a + x,b + x };
	}
	Interval<T>& moveSelf(T x)
	{
		a += x; b += x; return *this;
	}
	bool contains(Interval<T>const& s)const
	{
		if (valid())
			return !((s.a < a) || (b < s.b));
		else return false;
	}
	bool hasIntersectionWith(Interval<T>const& s)const
	{
		if (valid())
			return (s.a < a ? a : s.a) <= (b < s.b ? b : s.b);
		else return false;
	}
	bool isContainedIn(Interval<T>const& s)const
	{
		if (valid())
			return !((a < s.a) || (s.b < b));
		else return false;
	}
	bool valid()const
	{
		return a <= b;
	}
	T area()const
	{
		return b - a;
	}
	void print()const;
	void print(char const*)const;
	void print(char const*, char const*)const;
};
template<class T>struct IntervalSet :Vector<Interval<T>>
{
	//b is the section width, which must be positive
	using B = Vector<Interval<T>>;
	IntervalSet() :Vector<Interval<T>>()
	{
	}
	IntervalSet(Interval<T>const& a) :Vector<Interval<T>>(a) {}
	IntervalSet(IntervalSet<T>const& a) :Vector<Interval<T>>(a) {}
	IntervalSet(Vector<Interval<T>>const& a) :Vector<Interval<T>>(a) {}
	template<class R>IntervalSet(Vector<R>const& a, T const& b, bool withBorder) : Vector<Interval<T>>()
	{
		static_assert(NumType<R>::value == true, "Non-numeric type not supported yet!");
		B::malloc(a.length);
		withBorder &= (NumType<T>::numType == IsInteger);
		for (int c0(0); c0 < a.length; c0++)
			B::pushBack(Interval<T>(T(a.data[c0]), T(a.data[c0]) + b - withBorder));
	}
	void print()const
	{
		::printf("{");
		for (int c0(0); c0 < B::length - 1; ++c0)
			B::data[c0].print("", ", ");
		if (B::length)B::data[B::length - 1].print();
		::printf("}\n");
	}
	void print(char const* a)const
	{
		::printf("%s", a);
		print();
	}
	bool checkOrder()const
	{
		int c0(0);
		while (c0 + 1 < B::length)
		{
			if (B::data[c0 + 1] < B::data[c0])return false;
			++c0;
		}
		return true;
	}
	bool checkSimplified()const
	{
		if (B::length)
		{
			if (checkOrder())
			{
				int c0(0);
				while (c0 + 1 < B::length)
				{
					if (B::data[c0].hasIntersectionWith(B::data[c0 + 1]))
						return false;
					++c0;
				}
			}
			else
			{
				IntervalSet<T>tp(*this);
				tp.sort();
				int c0(0);
				while (c0 + 1 < B::length)
				{
					if (tp.data[c0].hasIntersectionWith(tp.data[c0 + 1]))
						return false;
					++c0;
				}
			}
			return true;
		}
		return true;
	}
	T area(bool withBorder)const
	{
		if (B::length)
		{
			withBorder &= (NumType<T>::numType == IsInteger);
			T a(0);
			unsigned int n;
			if (checkSimplified())
			{
				n = B::length;
				for (int c0(0); c0 < B::length; ++c0)
					a += B::data[c0].area();
			}
			else
			{
				IntervalSet<T>tp(*this);
				tp.simplify();
				n = tp.length;
				for (int c0(0); c0 < tp.length; ++c0)
					a += tp.data[c0].area();
			}
			return a + withBorder * n;
		}
		return 0;
	}
	IntervalSet<T>& simplify()
	{
		if (B::length)
		{
			sort();
			Vector<Interval<T>>vp;
			vp.pushBack(B::data[0]);
			Interval<T>* ip(&vp.end());
			int c0(1);
			while (c0 < B::length)
			{
				if ((*ip).hasIntersectionWith(B::data[c0]))
				{
					(*ip) += B::data[c0++];
				}
				else
				{
					vp.pushBack(B::data[c0++]);
					ip = &vp.end();
				}
			}
			vp.moveTo(*this);
		}
		return *this;
	}
	IntervalSet<T>& sort()
	{
		if (!checkOrder())
			qsort(B::data, 0, B::length);
		return *this;
	}
	IntervalSet<T>  operator^(Interval<T>const& s)const
	{
		if (!s.valid() || B::length == 0)return IntervalSet<T>();
		IntervalSet<T>tp(*this);
		tp.simplify();
		if (s.a > tp.end().b || s.b < tp.data[0].a)return IntervalSet<T>();
		int p, q;//[p+1, q-1] is the interval which has intersection with s
		if (s.a <= tp.data[0].b)p = -1;
		else
		{
			int n0(0), n1(tp.length - 1);
			while (n1 - n0 > 1)
				if (tp.data[(n0 + n1) / 2].b < s.a)n0 = (n0 + n1) / 2;
				else n1 = (n0 + n1) / 2;
			p = n0;
		}
		if (tp.end().a <= s.b)q = tp.length;
		else
		{
			int n0(0), n1(tp.length - 1);
			while (n1 - n0 > 1)
				if (s.b < tp.data[(n0 + n1) / 2].a)n1 = (n0 + n1) / 2;
				else n0 = (n0 + n1) / 2;
			q = n1;
		}
		switch (q - p)
		{
			case 1:return IntervalSet<T>();
			case 2:return tp.data[p + 1] ^ s;
			default:
			{
				tp.truncateSelf(p + 1, q - p - 1);
				tp.begin() ^= s;
				tp.end() ^= s;
				return tp;
			}
		}
	}
	IntervalSet<T>  operator^(IntervalSet<T>const& s)const
	{
		IntervalSet<T>t0(*this);
		IntervalSet<T>t1(s);
		t0.simplify();
		t1.simplify();
		if (t0.length == 0 || t1.length == 0)return IntervalSet<T>();
		IntervalSet<T>as;
		as.malloc(t0.length + t1.length);
		int p(0), q(0);
		while (p < t0.length && q < t1.length)
		{
			bool flag(false);
			if (t0.data[p].hasIntersectionWith(t1.data[q]))
			{
				as.pushBack(t0.data[p] ^ t1.data[q]);
				flag = true;
			}
			if (t0.data[p] < t1.data[q]) { p++; if (flag) q++; }
			else { q++; if (flag) p++; }
		}
		if (p == t0.length)
			as.concat(t1.data + q, t1.length - q);
		else
			as.concat(t0.data + p, t0.length - p);
		return as;
	}
	IntervalSet<T>& operator^=(Interval<T>const& s)
	{
		if (!s.valid() || B::length == 0)
		{
			this->B::~Vector();
			return *this;
		}
		simplify();
		if (s.a > B::end().b || s.b < B::data[0].a)
		{
			this->B::~Vector();
			return *this;
		}
		int p, q;//[p+1, q-1] is the interval which has intersection with s
		if (s.a <= B::data[0].b)p = -1;
		else
		{
			int n0(0), n1(B::length - 1);
			while (n1 - n0 > 1)
				if (B::data[(n0 + n1) / 2].b < s.a)n0 = (n0 + n1) / 2;
				else n1 = (n0 + n1) / 2;
			p = n0;
		}
		if (B::end().a <= s.b)q = B::length;
		else
		{
			int n0(0), n1(B::length - 1);
			while (n1 - n0 > 1)
				if (s.b < B::data[(n0 + n1) / 2].a)n1 = (n0 + n1) / 2;
				else n0 = (n0 + n1) / 2;
			q = n1;
		}
		switch (q - p)
		{
			case 1:this->B::~Vector(); break;
			case 2:
			{
				Interval<T>tp(B::data[p + 1] ^ s);
				this->B::~Vector();
				B::pushBack(tp);
				break;
			}
			default:
			{
				B::truncateSelf(p + 1, q - p - 1);
				B::begin() ^= s;
				B::end() ^= s;
			}
		}
		return *this;
	}
	IntervalSet<T>& operator^=(IntervalSet<T>const& s)
	{
		IntervalSet<T>t1(s);
		simplify();
		t1.simplify();
		if (B::length == 0 || t1.length == 0)
		{
			this->B::~Vector();
			return *this;
		}
		IntervalSet<T>as;
		as.malloc(B::length + t1.length);
		int p(0), q(0);
		while (p < B::length && q < t1.length)
		{
			bool flag(false);
			if (B::data[p].hasIntersectionWith(t1.data[q]))
			{
				as.pushBack(B::data[p] ^ t1.data[q]);
				flag = true;
			}
			if (B::data[p] < t1.data[q]) { p++; if (flag) q++; }
			else { q++; if (flag) p++; }
		}
		if (p == B::length)
			as.concat(t1.data + q, t1.length - q);
		else
			as.concat(B::data + p, B::length - p);
		as.moveTo(*this);
		return *this;
	}
	IntervalSet<T>  move(T x)const
	{
		IntervalSet<T>tp(*this);
		return tp.moveSelf(x);
	}
	IntervalSet<T>& moveSelf(T x)
	{
		for (int c0(0); c0 < B::length; ++c0)B::data[c0].moveSelf(x);
		return *this;
	}
};


//Interval
template<class T>IntervalSet<T> Interval<T>::operator^(IntervalSet<T>const& s)const
{
	return s ^ (*this);
}

//_Vector.h
template<class T>inline Vector<T>  Vector<T>::truncate(Interval<int> const& s)const
{
	Interval<int>itvl(0, length - 1);
	itvl ^= s;
	if (itvl.valid())return truncate(itvl.a, itvl.b - itvl.a + 1);
	return Vector<T>();
}
template<class T>inline Vector<T>  Vector<T>::truncate(IntervalSet<int> const& s)const
{
	IntervalSet<int>tp(s ^ Interval<int>(0, length - 1));
	int _lengthAll(tp.area(true));
	if (!_lengthAll)return Vector<T>();
	Vector<T>as;
	as.malloc(_lengthAll);
	as.length = _lengthAll;
	int n(0);
	for (int c0(0); c0 < tp.length; ++c0)
		for (int c1(tp.data[c0].a); c1 <= tp.data[c0].b; ++c1)
			new(as.data + n++)T(data[c1]);
	return as;
}
template<class T>inline Vector<T>& Vector<T>::truncateSelf(Interval<int> const& s)
{
	Interval<int>itvl(0, length - 1);
	itvl ^= s;
	if (itvl.valid())return truncateSelf(itvl.a, itvl.b - itvl.a + 1);
	this->~Vector();
	return *this;
}
template<class T>inline Vector<T>& Vector<T>::truncateSelf(IntervalSet<int> const& s)
{
	IntervalSet<int>tp(s ^ Interval<int>(0, length - 1));
	length = tp.area(true);
	if (!length)
	{
		this->~Vector();
		return *this;
	}
	int n(0);
	for (int c0(0); c0 < tp.length; ++c0)
		for (int c1(tp.data[c0].a); c1 <= tp.data[c0].b; ++c1)
			if (n != c1)
			{
				(data + n)->~T();
				new(data + n++)T(data[c1]);
				(data + c1)->~T();
			}
	return *this;
}
