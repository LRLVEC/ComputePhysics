#pragma once
#include <new>
#include <typeinfo>
#include <initializer_list>
#include <_Vector.h>
//Just a array, you need to implement a control system yourself...  -.-
template<class T, unsigned int _length>struct Array
{
	using elementType = T;
	static_assert(_length, "Array length cannot be 0!");
	static constexpr unsigned int length = _length;

	T data[_length];

	//Construction
	Array();
	Array(std::initializer_list<T>const&);
	template<unsigned int _length1>Array(Array<T, _length1>const&);
	Array(Array<T, _length>const&);
	Array(T const&);
	//Destrucion
	~Array();
	//operator=
	Array<T, _length>& operator=	(Array<T, _length>&&);
	Array<T, _length>& operator=	(Array<T, _length>const&);
	//operator+
	template<unsigned int _length1>auto operator+(Array<T, _length1>&&);
	template<unsigned int _length1>auto operator+(Array<T, _length1>const&);
	//element
	T& begin();
	T& end();
	T* endptr();
	T& operator[](unsigned int);
	//Find position
	int posFirst(const T&)const;
	Vector<unsigned int>posAll(T const&)const;
	//Find element
	T& findFirst(T const&);
	T& findFirst(bool(*cmp)(const T&, const T&), T const&);
	T& findFirst(bool(*cmp)(const T&, const T&), T&&);
	Vector<T*> find(T const&);
	//traverse
	bool traverse(bool(*p)(T&));
	bool traverse(bool(*p)(T const&))const;
	//printInfo
	void printInfo()const;
};





//Construction
template<class T, unsigned int _length>inline Array<T, _length>::Array()
{
}
template<class T, unsigned int _length>inline Array<T, _length>::Array(std::initializer_list<T>const& a)
{
	unsigned int tempLength = _length <= a.size() ? length : a.size();
	T const* p(a.begin());
	for (unsigned int c0 = 0; c0 < tempLength; ++c0)
		new(data + c0)T(*(p + c0));
}

template<class T, unsigned int _length>template<unsigned int _length1>inline Array<T, _length>::Array(Array<T, _length1>const& a)
{
	if constexpr (_length <= _length1)
	{
		for (unsigned int c0 = 0; c0 < _length; ++c0)
			new(data + c0)T(a.data[c0]);
	}
	else
	{
		for (unsigned int c0 = 0; c0 < _length1; ++c0)
			new(data + c0)T(a.data[c0]);
	}
}
template<class T, unsigned int _length>inline Array<T, _length>::Array(Array<T, _length>const& a)
{
	for (unsigned int c0 = 0; c0 < _length; ++c0)
		new(data + c0)T(a.data[c0]);
}
template<class T, unsigned int _length>inline Array<T, _length>::Array(T const& a)
{
	*data = a;
}
//Destruction
template<class T, unsigned int _length>inline Array<T, _length>::~Array()
{
	for (T& d : data)d.~T();
}
//operator=
template<class T, unsigned int _length>inline Array<T, _length>& Array<T, _length>::operator=(Array<T, _length> && a)
{
	for (T& d : data)
	{
		d.~T();
		new(&d)T(*(a.data + &d - data));
	}
}
template<class T, unsigned int _length>inline Array<T, _length>& Array<T, _length>::operator=(Array<T, _length>const& a)
{
	for (T& d : data)
	{
		d.~T();
		new(&d)T(*(a.data + &d - data));
	}
}
//operator+
template<class T, unsigned int _length>template<unsigned int _length1>inline auto Array<T, _length>::operator+(Array<T, _length1> && a)
{
	Array<T, _length + _length1>r;
	int c0 = 0;
	for (; c0 < _length; ++c0)new(r.data + c0)T(data[c0]);
	for (; c0 < r.length; ++c0)new(r.data + c0)T(a.data[c0 - _length]);
	return r;
}
template<class T, unsigned int _length>template<unsigned int _length1>inline auto Array<T, _length>::operator+(Array<T, _length1>const& a)
{
	Array<T, _length + _length1>r;
	int c0 = 0;
	for (; c0 < _length; ++c0)new(r.data + c0)T(data[c0]);
	for (; c0 < r.length; ++c0)new(r.data + c0)T(a.data[c0 - _length]);
	return r;
}
//element
template<class T, unsigned int _length>inline T & Array<T, _length>::begin()
{
	return *data;
}
template<class T, unsigned int _length>inline T& Array<T, _length>::end()
{
	return data[_length - 1];
}
template<class T, unsigned int _length>inline T* Array<T, _length>::endptr()
{
	return data + _length - 1;
}
template<class T, unsigned int _length>inline T& Array<T, _length>::operator[](unsigned int a)
{
	return data[a];
}
//Find pos
template<class T, unsigned int _length>inline int Array<T, _length>::posFirst(const T & a) const
{
	for (T const& d : data)
		if (d == a)return &d - data;
	return -1;
}
template<class T, unsigned int _length>inline Vector<unsigned int> Array<T, _length>::posAll(T const& a) const
{
	Vector<unsigned int>r;
	for (T const& d : data)
		if (d == a)
			r.pushBack(&d - data);
	return r;
}
//Find element
template<class T, unsigned int _length>inline T& Array<T, _length>::findFirst(T const& a)
{
	for (T const& d : data)
		if (d == a)return d;
	return *(T*)nullptr;
}
template<class T, unsigned int _length>inline T & Array<T, _length>::findFirst(bool(*cmp)(const T&, const T&), T const& a)
{
	for (T const& d : data)
		if (cmp(a, d))return d;
	return *(T*)nullptr;
}
template<class T, unsigned int _length>inline T & Array<T, _length>::findFirst(bool(*cmp)(const T&, const T&), T && a)
{
	for (T const& d : data)
		if (cmp(a, d))return d;
	return *(T*)nullptr;
}
template<class T, unsigned int _length>inline Vector<T*> Array<T, _length>::find(T const& a)
{
	Vector<T*>r;
	for (T const& d : data)
		if (d == a)
			r.pushBack(&d);
	return r;
}
//traverse
template<class T, unsigned int _length>inline bool Array<T, _length>::traverse(bool(*p)(T&))
{
	for (T& d : data)
		if (!p(d))return false;
	return true;
}
template<class T, unsigned int _length>inline bool Array<T, _length>::traverse(bool(*p)(T const&)) const
{
	for (T const& d : data)
		if (!p(d))return false;
	return true;
}
//printInfo
template<class T, unsigned int _length>inline void Array<T, _length>::printInfo() const
{
	::printf("%s:", typeid(*this).name());
}

