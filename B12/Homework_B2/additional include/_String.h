#pragma once
#include <clocale>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <_Vector.h>
#include <_Algorithm.h>

/*
	To add:

	1. Concat some strings to construct
	2. omit
	3. if String<T>::data == nullptr, then operator+... do nothing
*/


static bool wchar_tInited(false);
static void wchar_tInit();
//Support char and wchar_t in Windows only!
//Note: String<wchar_t> stores Chinese charactors in unicode!
template<class T>struct String
{
	static_assert(CharType<T>::value, "Wrong CharType!");
	template<class R>static T* transfer(R const*);

	T* data;
	unsigned int length;
	unsigned int lengthAll;

	//Construction
	String();
	template<class R>String(String<R>const&);
	String(String<T>const&);
	String(String<T>&&);
	template<class R>String(R const*);
	String(T const*);
	String(T*&, unsigned int, unsigned int);
	String(T*&&, unsigned int, unsigned int);
	//Destruction
	~String();
	//operator T*()
	operator T* ();
	operator T const* ()const;
	//operator=
	template<class R>String<T>& operator=(String<R>const&);
	String<T>& operator=(String<T>const&);
	//operator==
	template<class R>bool operator==(String<R>const&)const;
	template<class R>bool operator==(R const*)const;
	//operator+
	template<class R>auto operator+(String<R>const&)const;
	template<class R>auto operator+(R const*)const;
	//operator+=
	template<class R>auto& operator+=(String<R>const&);
	template<class R>auto& operator+=(R const*);
	template<class R>auto& operator+=(R*);
	template<class R>auto& operator+=(R);
	//exchange
	String<T>& exchange(String<T>&);
	//find
	template<class R>int findFirst(String<R>const&);
	int* getKmpNext()const;
	int KMP(String<T>const&)const;
	template<class R>int findFirstKMP(String<R>const&)const;
	template<class R>int findFirst(R const*)const;
	template<class R>int findNext(String<R>const&, int)const;
	template<class R>int findNext(R const*,int)const;
	template<class R>int findAhead(String<R>const&,int)const;
	template<class R>int findAhead(R const*,int)const;
	template<class R>Vector<int>find(String<R>const&)const;
	template<class R>Vector<int>find(R const*)const;
	//malloc
	String<T>& malloc(unsigned int);
	//omit
	String<T>  omit(int)const;
	String<T>  omit(int, int)const;
	String<T>  omit(Interval<int>const&)const;
	String<T>  omit(IntervalSet<int>const&)const;
	String<T>& omitSelf(int);
	String<T>& omitSelf(int, int);
	String<T>& omitSelf(Interval<int>const&);
	String<T>& omitSelf(IntervalSet<int>const&);
	template<class R>String<T>  omitSelf(String<R>const&)const;
	template<class R>String<T>& omitSelf(R const*);
	template<class R>String<T>& omitSelf(String<R>const&);
	//truncate: if _length < 0: truncate the left
	String<T>  truncate(int, int)const;
	String<T>  truncate(Interval<int>const&)const;
	String<T>  truncate(IntervalSet<int>const&)const;
	String<T>& truncateSelf(int, int);
	String<T>& truncateSelf(Interval<int>const&);
	String<T>& truncateSelf(IntervalSet<int>const&);
	//print
	void print()const;
	void print(char const*)const;
	void print(char const*, char const*)const;
	void printInfo()const;
};


inline void wchar_tInit()
{
	setlocale(LC_CTYPE, "chs");
	wchar_tInited = true;
}
//Must add %n at the ending...
inline unsigned int sweepStr(char const* _str, char const* _mode)
{
	unsigned int t(0);
	::sscanf(_str, _mode, &t);
	return t;
}
//transfer
template<class T>template<class R>	inline T* String<T>::transfer(R const* a)
{
	static_assert(!IsSameType<T, R>::value, "Cannot tranfer same CharType!");
	static_assert(CharType<R>::value, "Wrong CharType!");
	if (!a)return (T*)nullptr;
	if constexpr (IsSameType<R, char>::value)
	{
		wchar_t* temp;
		unsigned int tempLength((unsigned int)::mbstowcs(nullptr, a, 0) + 1);
		::mbstowcs(temp = (wchar_t*)::malloc(tempLength << 1), a, tempLength);
		return temp;
	}
	else
	{
		char* temp;
		unsigned int tempLength((unsigned int)::wcstombs(nullptr, a, 0) + 1);
		::wcstombs(temp = (char*)::malloc(tempLength), a, tempLength);
		return temp;
	}
}
//Construction
template<class T>					inline String<T>::String() :data(nullptr), length(0), lengthAll(0)
{
}
template<class T>template<class R>	inline String<T>::String(String<R>const& a)
{
	if (!a.data)
	{
		data = nullptr;
		length = lengthAll = 0;
		return;
	}
	if constexpr (IsSameType<T, char>::value)
	{
		length = (unsigned int)::wcstombs(nullptr, a.data, 0);
		lengthAll = 1;
		while (lengthAll < length + 1)lengthAll <<= 1;
		::wcstombs(data = (T*)::malloc(lengthAll), a.data, length + 1);
	}
	else
	{
		length = (unsigned int)::mbstowcs(nullptr, a.data, 0);
		lengthAll = 1;
		while (lengthAll < length + 1)lengthAll <<= 1;
		::mbstowcs(data = (T*)::malloc(lengthAll << 1), a.data, length + 1);
	}
}
template<class T>					inline String<T>::String(String<T>const& a)
{
	if (&a == this)return;
	if (!a.data)
	{
		data = nullptr;
		length = lengthAll = 0;
		return;
	}
	if constexpr (IsSameType<T, char>::value)
	{
		length = a.length;
		lengthAll = a.lengthAll;
		data = (T*)::malloc(lengthAll);
		memcpy(data, a.data, length + 1);
	}
	else
	{
		length = a.length;
		lengthAll = a.lengthAll;
		data = (T*)::malloc(lengthAll << 1);
		memcpy(data, a.data, (length + 1) << 1);
	}
}
template<class T>					inline String<T>::String(String<T>&& a)
	:
	data(a.data),
	length(a.length),
	lengthAll(a.lengthAll)
{
	a.data = nullptr;
}
template<class T>template<class R>	inline String<T>::String(R const* a)
{
	static_assert(CharType<T>::value, "Wrong CharType!");
	if (!a)
	{
		data = nullptr;
		length = lengthAll = 0;
		return;
	}
	if constexpr (IsSameType<T, char>::value)
	{
		length = (unsigned int)::wcstombs(nullptr, a, 0);
		lengthAll = 1;
		while (lengthAll < length + 1)lengthAll <<= 1;
		::wcstombs(data = (T*)::malloc(lengthAll), a, length + 1);
	}
	else
	{
		length = (unsigned int)::mbstowcs(nullptr, a, 0);
		lengthAll = 1;
		while (lengthAll < length + 1)lengthAll <<= 1;
		::mbstowcs(data = (T*)::malloc(lengthAll << 1), a, length + 1);
	}
}
template<class T>					inline String<T>::String(T const* a)
{
	if constexpr (IsSameType<T, char>::value)
	{
		length = (unsigned int)::strlen(a);
		lengthAll = 1;
		while (lengthAll < length + 1)lengthAll <<= 1;
		data = (T*)::malloc(lengthAll);
		memcpy(data, a, length + 1);
	}
	else
	{
		length = (unsigned int)::wcslen(a);
		lengthAll = 1;
		while (lengthAll < length + 1)lengthAll <<= 1;
		data = (T*)::malloc(lengthAll << 1);
		memcpy(data, a, (length + 1) << 1);
	}
}
template<class T>					inline String<T>::String(T*& a, unsigned int _length, unsigned int _lengthAll)
{
	if (_length)
	{
		data = a;
		a = nullptr;
		length = _length;
		if (_lengthAll)lengthAll = _lengthAll;
		else
		{
			lengthAll = 1;
			while (lengthAll < length + 1)lengthAll <<= 1;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
			data[length] = 0;
		}
	}
	else
	{
		data = a;
		a = nullptr;
		if constexpr (IsSameType<T, char>::value)length = (unsigned int)::strlen(data);
		else length = (unsigned int)::wcslen(data);
		if (_lengthAll)lengthAll = _lengthAll;
		else
		{
			lengthAll = 1;
			while (lengthAll < length + 1)lengthAll <<= 1;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
			data[length] = 0;
		}
	}
}
template<class T>					inline String<T>::String(T*&& a, unsigned int _length, unsigned int _lengthAll)
{
	if (_length)
	{
		data = a;
		a = nullptr;
		length = _length;
		if (_lengthAll)lengthAll = _lengthAll;
		else
		{
			lengthAll = 1;
			while (lengthAll < length + 1)lengthAll <<= 1;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
			data[length] = 0;
		}
	}
	else
	{
		data = a;
		a = nullptr;
		if constexpr (IsSameType<T, char>::value)length = (unsigned int)::strlen(data);
		else length = (unsigned int)::wcslen(data);
		if (_lengthAll)lengthAll = _lengthAll;
		else
		{
			lengthAll = 1;
			while (lengthAll < length + 1)lengthAll <<= 1;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
			data[length] = 0;
		}
	}
}
//Destrucion
template<class T>					inline String<T>::~String()
{
	if (data)::free(data);
	data = nullptr;
}
//operator T*()
template<class T>					inline String<T>::operator T* ()
{
	return data;
}
template<class T>inline String<T>::operator T const* ()const
{
	return (T const*)data;
}
//operator=
template<class T>template<class R>	inline String<T>& String<T>::operator=(String<R>const& a)
{
	if (data)::free(data);
	if constexpr (IsSameType<T, R>::value)
	{
		if constexpr (IsSameType<T, char>::value)
		{
			length = a.length;
			lengthAll = a.lengthAll;
			data = (T*)::malloc(lengthAll);
			memcpy(data, a.data, length + 1);
		}
		else
		{
			length = a.length;
			lengthAll = a.lengthAll;
			data = (T*)::malloc(lengthAll << 1);
			memcpy(data, a.data, (length + 1) << 1);
		}
		return *this;
	}
	else
	{
		if constexpr (IsSameType<T, char>::value)
		{
			length = (unsigned int)::wcstombs(nullptr, a.data, 0);
			lengthAll = 1;
			while (lengthAll < length + 1)lengthAll <<= 1;
			::wcstombs(data = (T*)::malloc(lengthAll), a.data, length + 1);
		}
		else
		{
			length = (unsigned int)::mbstowcs(nullptr, a.data, 0);
			lengthAll = 1;
			while (lengthAll < length + 1)lengthAll <<= 1;
			::mbstowcs(data = (T*)::malloc(lengthAll << 1), a.data, length + 1);
		}
	}
	return *this;
}
template<class T>					inline String<T>& String<T>::operator=(String<T>const& a)
{
	if (data)::free(data);
	if constexpr (IsSameType<T, char>::value)
	{
		length = a.length;
		lengthAll = a.lengthAll;
		data = (T*)::malloc(lengthAll);
		memcpy(data, a.data, length + 1);
	}
	else
	{
		length = a.length;
		lengthAll = a.lengthAll;
		data = (T*)::malloc(lengthAll << 1);
		memcpy(data, a.data, (length + 1) << 1);
	}
	return *this;
}
//operator==
template<class T>template<class R>	inline bool String<T>::operator==(String<R>const& a)const
{
	if constexpr (IsSameType<T, R>::value)
	{
		if constexpr (IsSameType<T, char>::value)return !strcmp(data, a.data);
		else return !wcscmp(data, a.data);
	}
	else
	{
		if constexpr (IsSameType<T, char>::value)
		{
			char* temp(String<char>::transfer(a.data));
			bool tb(!strcmp(data, temp));
			::free(temp);
			return tb;
		}
		else
		{
			char* temp(String<char>::transfer(data));
			bool tb(!strcmp(a.data, temp));
			::free(temp);
			return tb;
		}
	}
}
template<class T>template<class R>	inline bool String<T>::operator==(R const* a)const
{
	if constexpr (IsSameType<T, R>::value)
	{
		if constexpr (IsSameType<T, char>::value)return !strcmp(data, a);
		else return !wcscmp(data, a);
	}
	else
	{
		static_assert(CharType<R>::value, "Wrong CharType!");
		if constexpr (IsSameType<T, char>::value)
		{
			char* temp(String<char>::transfer(a));
			bool tb(!strcmp(data, temp));
			::free(temp);
			return tb;
		}
		else
		{
			char* temp(String<char>::transfer(data));
			bool tb(!strcmp(a, temp));
			::free(temp);
			return tb;
		}
	}
}
template<class T, class R>			inline bool operator==(R const* a, String<T>const& b)
{
	if constexpr (IsSameType<T, R>::value)
	{
		if constexpr (IsSameType<T, char>::value)return !strcmp(b.data, a);
		else return !wcscmp(b.data, a);
	}
	else
	{
		static_assert(CharType<R>::value, "Wrong CharType!");
		if constexpr (IsSameType<T, char>::value)
		{
			char* temp(String<char>::transfer(a));
			bool tb(!strcmp(b.data, temp));
			::free(temp);
			return tb;
		}
		else
		{
			char* temp(String<char>::transfer(b.data));
			bool tb(!strcmp(a, temp));
			::free(temp);
			return tb;
		}
	}
}
//operator+
template<class T>template<class R>	inline auto String<T>::operator+(String<R>const& a)const
{
	unsigned int _length, _lengthAll(1);
	if constexpr (CharType<T>::serial + CharType<R>::serial > 2)
	{
		wchar_t* temp;
		if constexpr (IsSameType<T, char>::value)
		{
			unsigned int tempLength((unsigned int)::mbstowcs(nullptr, data, 0));
			_length = (tempLength + a.length);
			while (_lengthAll < _length + 1)_lengthAll <<= 1;
			temp = (wchar_t*)::malloc(_lengthAll << 1);
			::mbstowcs(temp, data, tempLength);
			::memcpy(temp + tempLength, a.data, a.length << 1);
			temp[_length] = 0;
		}
		else
		{
			if constexpr (IsSameType<R, char>::value)
			{
				unsigned int tempLength((unsigned int)::mbstowcs(nullptr, a.data, 0));
				_length = length + tempLength;
				while (_lengthAll < _length + 1)_lengthAll <<= 1;
				temp = (wchar_t*)::malloc(_lengthAll << 1);
				::memcpy(temp, data, lengthAll << 1);
				::mbstowcs(temp + length, a.data, tempLength);
				temp[_length] = 0;
			}
			else
			{

				_length = length + a.length;
				while (_lengthAll < _length + 1)_lengthAll <<= 1;
				temp = (wchar_t*)::malloc(_lengthAll << 1);
				::memcpy(temp, data, length << 1);
				::memcpy(temp + length, a.data, a.length << 1);
				temp[_length] = 0;
			}
		}
		return String<wchar_t>(temp, _length, _lengthAll);
	}
	else
	{
		char* temp;
		_length = length + a.length;
		while (_lengthAll < _length + 1)_lengthAll <<= 1;
		temp = (char*)::malloc(_lengthAll);
		::memcpy(temp, data, length);
		::memcpy(temp + length, a.data, a.length);
		temp[_length] = 0;
		return String<char>(temp, _length, _lengthAll);
	}
}
template<class T>template<class R>	inline auto String<T>::operator+(R const* a)const
{
	unsigned int _length, _lengthAll(1);
	if constexpr (CharType<T>::serial + CharType<R>::serial > 2)
	{
		wchar_t* temp;
		if constexpr (IsSameType<T, char>::value)
		{
			unsigned int tempLength0((unsigned int)::mbstowcs(nullptr, data, 0));
			unsigned int tempLength1((unsigned int)::wcslen(a));
			_length = (tempLength0 + tempLength1);
			while (_lengthAll < _length + 1)_lengthAll <<= 1;
			temp = (wchar_t*)::malloc(_lengthAll << 1);
			::mbstowcs(temp, data, tempLength0);
			::memcpy(temp + tempLength0, a, tempLength1 << 1);
			temp[_length] = 0;
		}
		else
		{
			if constexpr (IsSameType<R, char>::value)
			{
				unsigned int tempLength((unsigned int)::mbstowcs(nullptr, a, 0));
				_length = length + tempLength;
				while (_lengthAll < _length + 1)_lengthAll <<= 1;
				temp = (wchar_t*)::malloc(_lengthAll << 1);
				::memcpy(temp, data, lengthAll << 1);
				::mbstowcs(temp + length, a, tempLength);
				temp[_length] = 0;
			}
			else
			{
				unsigned int tempLength((unsigned int)::wcslen(a));
				_length = length + tempLength;
				while (_lengthAll < _length + 1)_lengthAll <<= 1;
				temp = (wchar_t*)::malloc(_lengthAll << 1);
				::memcpy(temp, data, length << 1);
				::memcpy(temp + length, a, tempLength << 1);
				temp[_length] = 0;
			}
		}
		return String<wchar_t>(temp, _length, _lengthAll);
	}
	else
	{
		char* temp;
		unsigned int tempLength((unsigned int)::strlen(a));
		_length = length + tempLength;
		while (_lengthAll < _length + 1)_lengthAll <<= 1;
		temp = (char*)::malloc(_lengthAll);
		::memcpy(temp, data, length);
		::memcpy(temp + length, a, tempLength);
		temp[_length] = 0;
		return String<char>(temp, _length, _lengthAll);
	}
}
template<class T, class R>			inline auto operator+(R const* a, String<T>const& b)
{
	static_assert(CharType<T>::value, "Wrong CharType!");
	unsigned int length, lengthAll(1);
	if constexpr (CharType<T>::serial + CharType<R>::serial > 2)
	{
		wchar_t* temp;
		if constexpr (IsSameType<R, char>::value)
		{
			unsigned int tempLength((unsigned int)::mbstowcs(nullptr, a, 0));
			length = (tempLength + b.length);
			while (lengthAll < length + 1)lengthAll <<= 1;
			temp = (wchar_t*)::malloc(lengthAll << 1);
			::mbstowcs(temp, a, tempLength);
			::memcpy(temp + tempLength, b.data, b.length << 1);
			temp[length] = 0;
		}
		else
		{
			if constexpr (IsSameType<T, char>::value)
			{
				unsigned int tempLength0((unsigned int)::wcslen(a));
				unsigned int tempLength1((unsigned int)::mbstowcs(nullptr, b.data, 0));
				length = tempLength0 + tempLength1;
				while (lengthAll < length + 1)lengthAll <<= 1;
				temp = (wchar_t*)::malloc(lengthAll << 1);
				::memcpy(temp, a, tempLength0 << 1);
				::mbstowcs(temp + tempLength0, b.data, tempLength1);
				temp[length] = 0;
			}
			else
			{
				unsigned int tempLength((unsigned int)::wcslen(a));
				length = tempLength + b.length;
				while (lengthAll < length + 1)lengthAll <<= 1;
				temp = (wchar_t*)::malloc(lengthAll << 1);
				::memcpy(temp, a, tempLength << 1);
				::memcpy(temp + tempLength, b.data, b.length << 1);
				temp[length] = 0;
			}
		}
		return String<wchar_t>(temp, length, lengthAll);
	}
	else
	{
		char* temp;
		unsigned int tempLength((unsigned int)::strlen(a));
		length = tempLength + b.length;
		while (lengthAll < length + 1)lengthAll <<= 1;
		temp = (char*)::malloc(lengthAll);
		::memcpy(temp, a, tempLength);
		::memcpy(temp + tempLength, b.data, b.length);
		temp[length] = 0;
		return String<char>(temp, length, lengthAll);
	}
}
//operator+=
template<class T>template<class R>	inline auto& String<T>::operator+=(String<R>const& a)
{
	if constexpr (IsSameType<T, R>::value)
	{
		lengthAll = lengthAll ? lengthAll : 1;
		if (lengthAll <= length + a.length)
		{
			while (lengthAll <= length + a.length)lengthAll <<= 1;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
		}
		::memcpy(data + length, a.data, a.length * sizeof(T));
		length += a.length;
		data[length] = 0;
		return *this;
	}
	else
	{
		unsigned int tempLength;
		if constexpr (IsSameType<T, char>::value)
		{
			tempLength = (unsigned int)::wcstombs(nullptr, a.data, 0);
			lengthAll = lengthAll ? lengthAll : 1;
			if (lengthAll <= length + tempLength)
			{
				while (lengthAll <= length + tempLength)lengthAll <<= 1;
				data = (T*)::realloc(data, lengthAll);
			}
			::wcstombs(data + length, a.data, tempLength + 1);
			length += tempLength;
		}
		else
		{
			tempLength = (unsigned int)::mbstowcs(nullptr, a.data, 0);
			lengthAll = lengthAll ? lengthAll : 1;
			if (lengthAll <= length + tempLength)
			{
				while (lengthAll <= length + tempLength)lengthAll <<= 1;
				data = (T*)::realloc(data, lengthAll << 1);
			}
			::mbstowcs(data + length, a.data, tempLength + 1);
			length += tempLength;
		}
		return *this;
	}
}
template<class T>template<class R>	inline auto& String<T>::operator+=(R const* a)
{
	static_assert(CharType<R>::value, "Wrong CharType!");
	if constexpr (IsSameType<T, R>::value)
	{
		unsigned int tempLength;
		if constexpr (IsSameType<T, char>::value)tempLength = (unsigned int)::strlen(a);
		else tempLength = (unsigned int)::wcslen(a);
		if (lengthAll <= length + tempLength)
		{
			lengthAll = lengthAll ? lengthAll : 1;
			while (lengthAll <= length + tempLength)lengthAll <<= 1;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
		}
		::memcpy(data + length, a, tempLength * sizeof(T));
		length += tempLength;
		data[length] = 0;
		return *this;
	}
	else
	{
		unsigned int tempLength;
		if constexpr (IsSameType<T, char>::value)
		{
			tempLength = (unsigned int)::wcstombs(nullptr, a, 0);
			lengthAll = lengthAll ? lengthAll : 1;
			if (lengthAll <= length + tempLength)
			{
				while (lengthAll <= length + tempLength)lengthAll <<= 1;
				data = (T*)::realloc(data, lengthAll);
			}
			::wcstombs(data + length, a, tempLength + 1);
			length += tempLength;
		}
		else
		{
			tempLength = (unsigned int)::mbstowcs(nullptr, a, 0);
			lengthAll = lengthAll ? lengthAll : 1;
			if (lengthAll <= length + tempLength)
			{
				while (lengthAll <= length + tempLength)lengthAll <<= 1;
				data = (T*)::realloc(data, lengthAll << 1);
			}
			::mbstowcs(data + length, a, tempLength + 1);
			length += tempLength;
		}
		return *this;
	}
}
template<class T>template<class R>	inline auto& String<T>::operator+=(R* a)
{
	static_assert(CharType<R>::value, "Wrong CharType!");
	if constexpr (IsSameType<T, R>::value)
	{
		unsigned int tempLength;
		if constexpr (IsSameType<T, char>::value)tempLength = (unsigned int)::strlen(a);
		else tempLength = (unsigned int)::wcslen(a);
		if (lengthAll <= length + tempLength)
		{
			lengthAll = lengthAll ? lengthAll : 1;
			while (lengthAll <= length + tempLength)lengthAll <<= 1;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
		}
		::memcpy(data + length, a, tempLength * sizeof(T));
		length += tempLength;
		data[length] = 0;
		return *this;
	}
	else
	{
		unsigned int tempLength;
		if constexpr (IsSameType<T, char>::value)
		{
			tempLength = (unsigned int)::wcstombs(nullptr, a, 0);
			lengthAll = lengthAll ? lengthAll : 1;
			if (lengthAll <= length + tempLength)
			{
				while (lengthAll <= length + tempLength)lengthAll <<= 1;
				data = (T*)::realloc(data, lengthAll);
			}
			::wcstombs(data + length, a, tempLength + 1);
			length += tempLength;
		}
		else
		{
			tempLength = (unsigned int)::mbstowcs(nullptr, a, 0);
			lengthAll = lengthAll ? lengthAll : 1;
			if (lengthAll <= length + tempLength)
			{
				while (lengthAll <= length + tempLength)lengthAll <<= 1;
				data = (T*)::realloc(data, lengthAll << 1);
			}
			::mbstowcs(data + length, a, tempLength + 1);
			length += tempLength;
		}
		return *this;
	}
}
template<class T>template<class R>	inline auto& String<T>::operator+=(R a)
{
	static_assert(CharType<R>::value, "Wrong CharType!");
	if constexpr (IsSameType<T, R>::value)
	{
		lengthAll = lengthAll ? lengthAll : 1;
		if (lengthAll <= length + 2)
		{
			lengthAll <<= 1;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
		}
		data[length++] = a;
		data[length] = 0;
		return *this;
	}
	else if constexpr (IsSameType<R, char>::value)
	{
		if (lengthAll <= length + 2)
		{
			lengthAll = lengthAll > 1 ? lengthAll <<= 1 : 2;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
		}
		data[length++] = a;
		data[length] = 0;
		return *this;
	}
	else
	{
		char32_t temp(a);
		unsigned int tempLength((unsigned int)::wcstombs(nullptr, (wchar_t*)&temp, 0));
		if (lengthAll <= length + tempLength)
		{
			lengthAll = lengthAll > 2 ? lengthAll <<= 1 : 4;
			data = (T*)::realloc(data, lengthAll * sizeof(T));
		}
		::wcstombs(data + length, (wchar_t*)&temp, tempLength + 1);
		data[length += tempLength] = 0;
		return *this;
	}
}
//exchange
template<class T>					inline String<T>& String<T>::exchange(String<T>& a)
{
	T* tp(a.data);
	a.data = data;
	data = tp;
	unsigned int n(a.length);
	a.length = length;
	length = n;
	n = a.lengthAll;
	a.lengthAll = lengthAll;
	lengthAll = n;
	return *this;
}
//find
template<class T>template<class R>	inline int String<T>::findFirst(String<R>const& a)
{
	T* temp;
	if constexpr (IsSameType<T, R>::value)
	{
		if constexpr (IsSameType<T, char>::value) temp = strstr(data, a.data);
		else temp = wcsstr(data, a.data);
	}
	else
	{
		void* tpp;
		if constexpr (IsSameType<T, char>::value) temp = strstr(data, tpp = String<char>::transfer(a.data));
		else temp = wcsstr(data, tpp = String<wchar_t>::transfer(a.data));
		::free(tpp);
	}
	if (temp)return int(temp - data);
	return -1;
}
template<class T>					inline int* String<T>::getKmpNext()const
{
	int* r((int*)::malloc(4 * length));
	int c0(0), c1(-1);
	r[0] = -1;
	while (c0 < length - 1)
	{
		while (c1 >= 0 && data[c0] != data[c1])
			c1 = r[c1];
		++c0; ++c1;
		if (data[c0] == data[c1])r[c0] = r[c1];
		else r[c0] = c1;
	}
	return r;
}
template<class T>					inline int String<T>::KMP(String<T>const& a)const
{
	int* next(a.getKmpNext());
	int c0(0), c1(0);
	while ((c0 < length) && (c1 < int(a.length)))
		if (c1 == -1 || data[c0] == a.data[c1])
		{
			++c0; ++c1;
		}
		else c1 = next[c1];
	::free(next);
	if (c1 >= a.length)return c0 - a.length;
	else return -1;
}
template<class T>template<class R>	inline int String<T>::findFirstKMP(String<R>const& a)const
{
	int temp;
	if constexpr (IsSameType<T, R>::value)
		temp = KMP(a);
	else
	{
		void* tpp;
		if constexpr (IsSameType<T, char>::value) temp = KMP(String<char>(tpp = String<char>::transfer(a.data)));
		else temp = KMP(String<wchar_t>(tpp = String<wchar_t>::transfer(a.data)));
		::free(tpp);
	}
	return temp;
}
template<class T>template<class R>	inline int String<T>::findFirst(R const* a)const
{
	static_assert(CharType<R>::value, "Wrong CharType!");
	T* temp;
	if constexpr (IsSameType<T, R>::value)
	{
		if constexpr (IsSameType<T, char>::value) temp = strstr(data, a);
		else temp = wcsstr(data, a);
	}
	else
	{

		void* tpp;
		if constexpr (IsSameType<T, char>::value) temp = strstr(data, tpp = String<char>::transfer(a));
		else temp = wcsstr(data, tpp = String<wchar_t>::transfer(a));
		::free(tpp);
	}
	if (temp)return int(temp - data);
	return -1;
}
template<class T>template<class R>	inline Vector<int> String<T>::find(String<R>const& a)const
{
	Vector<int>temp;
	if constexpr (IsSameType<T, R>::value)
	{
		if constexpr (IsSameType<T, char>::value)
		{
			T* flag(strstr(data, a.data));
			int tempAddress(flag - data);
			while (flag)
			{
				temp.pushBack(tempAddress);
				tempAddress = (flag = strstr(data + tempAddress + 1, a.data)) - data;
			}
			return temp;
		}
		else
		{
			T* flag(wcsstr(data, a.data));
			int tempAddress(flag - data);
			while (flag)
			{
				temp.pushBack(tempAddress);
				tempAddress = (flag = wcsstr(data + tempAddress + 1, a.data)) - data;
			}
			return temp;
		}
	}
	else
	{
		T* ts(String<T>::transfer(a.data));
		if constexpr (IsSameType<T, char>::value)
		{
			T* flag(strstr(data, ts));
			int tempAddress(flag - data);
			while (flag)
			{
				temp.pushBack(tempAddress);
				tempAddress = (flag = strstr(data + tempAddress + 1, ts)) - data;
			}
			::free(ts);
			return temp;
		}
		else
		{
			T* flag(wcsstr(data, ts));
			int tempAddress(flag - data);
			while (flag)
			{
				temp.pushBack(tempAddress);
				tempAddress = (flag = wcsstr(data + tempAddress + 1, ts)) - data;
			}
			::free(ts);
			return temp;
		}
	}
}
template<class T>template<class R>	inline Vector<int> String<T>::find(R const* a)const
{
	static_assert(CharType<R>::value, "Wrong CharType!");
	Vector<int>temp;
	if constexpr (IsSameType<T, R>::value)
	{
		if constexpr (IsSameType<T, char>::value)
		{
			T* flag(strstr(data, a));
			int tempAddress(int(flag - data));
			while (flag)
			{
				temp.pushBack(tempAddress);
				tempAddress = int((flag = strstr(data + tempAddress + 1, a)) - data);
			}
			return temp;
		}
		else
		{
			T* flag(wcsstr(data, a));
			int tempAddress(int(flag - data));
			while (flag)
			{
				temp.pushBack(tempAddress);
				tempAddress = int((flag = wcsstr(data + tempAddress + 1, a)) - data);
			}
			return temp;
		}
	}
	else
	{
		T* ts(String<T>::transfer(a));
		if constexpr (IsSameType<T, char>::value)
		{
			T* flag(strstr(data, ts));
			int tempAddress(int(flag - data));
			while (flag)
			{
				temp.pushBack(tempAddress);
				tempAddress = int((flag = strstr(data + tempAddress + 1, ts)) - data);
			}
			::free(ts);
			return temp;
		}
		else
		{
			T* flag(wcsstr(data, ts));
			int tempAddress(int(flag - data));
			while (flag)
			{
				temp.pushBack(tempAddress);
				tempAddress = int((flag = wcsstr(data + tempAddress + 1, ts)) - data);
			}
			::free(ts);
			return temp;
		}
	}
}
//malloc
template<class T>inline String<T>& String<T>::malloc(unsigned int a)
{
	if (!a)return *this;
	if (length + a + 1 > lengthAll)
	{
		if (!lengthAll)lengthAll = 1;
		while (length + a + 1 > lengthAll)lengthAll <<= 1;
		data = (T*)::realloc(data, lengthAll * sizeof(T));
		data[length] = 0;
	}
	return *this;
}
//omit
template<class T>					inline String<T>& String<T>::omitSelf(IntervalSet<int>const& a)
{
	IntervalSet<int>tp(a ^ Interval<int>(0, length - 1));
	if (!tp.length)return*this;
	unsigned int _length = length - tp.area(true);
	if (!_length)
	{
		data = (T*)::realloc(data, sizeof(T));
		data[0] = 0;
		lengthAll = 1;
		return *this;
	}
	int n(0);
	int m(0);
	for (int c0(0); c0 < tp.length; ++c0)
	{
		int dn(tp.data[c0].a - m);
		if (!dn)
		{
			m = tp.data[c0].b + 1;
			continue;
		}
		::memmove(data + n, data + m, dn * sizeof(T));
		m = tp.data[c0].b + 1;
		n += dn;
	}
	::memmove(data + n, data + m, (length - m) * sizeof(T));
	data[length = _length] = 0;//bug here
	lengthAll = length + 1;
	data = (T*)::realloc(data, lengthAll * sizeof(T));
	return *this;
}
template<class T>template<class R>  inline String<T>& String<T>::omitSelf(R const* a)
{
	static_assert(CharType<R>::value, "Wrong CharType!");
	String<T> ar(a);
	Vector<int>pos(find(ar));
	IntervalSet<int>tp(pos, ar.length, true);
	tp.simplify();
	return omitSelf(tp);
}
template<class T>template<class R>  inline String<T>& String<T>::omitSelf(String<R> const& a)
{
	Vector<int>pos(find(a));
	IntervalSet<int>tp(pos, a.length, true);
	tp.simplify();
	return omitSelf(tp);
}
//truncate: if _length < 0 then truncate the part in the left;
//			if _head is out side the interval [0, length - 1], then do nothing
//				(include ~String<T> and return String<T>(*this))
template<class T>					inline String<T>  String<T>::truncate(int _head, int _length)const
{
	if (_head < 0 || _head >= (int)length)return String<T>();
	unsigned int _lengthAll(1);
	if (_head + _length > (int)length || _length <= 0) _length = length - _head;
	while ((int)_lengthAll < _length + 1)_lengthAll <<= 1;
	T* temp((T*)::malloc(_lengthAll * sizeof(T)));
	::memcpy(temp, data + _head, _length * sizeof(T));
	temp[_length] = 0;
	return String<T>(temp, _length, _lengthAll);
}
template<class T>					inline String<T>  String<T>::truncate(Interval<int> const& a)const
{
	Interval<int>itvl(0, length - 1);
	itvl ^= a;
	if (itvl.valid())return truncate(itvl.a, itvl.b - itvl.a + 1);
	return String<T>();
}
template<class T>					inline String<T>  String<T>::truncate(IntervalSet<int> const& a)const
{
	IntervalSet<int>tp(a ^ Interval<int>(0, length - 1));
	int _lengthAll(tp.area(true));
	if (!_lengthAll)return String<T>();
	String<T>as;
	as.malloc(_lengthAll);
	as.length = _lengthAll;
	int n(0);
	for (int c0(0); c0 < tp.length; ++c0)
	{
		int dn(tp.data[c0].b - tp.data[c0].a + 1);
		::memcpy(as.data + n, data + tp.data[c0].a, dn * sizeof(T));
		n += dn;
	}
	as.data[as.length] = 0;
	return as;
}
template<class T>					inline String<T>& String<T>::truncateSelf(int _head, int _length)
{
	if (_head < 0 || _head >= (int)length)return *this;
	if (_head + _length > (int)length || _length <= 0) length -= _head;
	while (lengthAll < length + 1)lengthAll <<= 1;
	::memmove(data, data + _head, length * sizeof(T));
	data[length] = 0;
	::realloc(data, sizeof(T) * lengthAll);
	return *this;
}
template<class T>					inline String<T>& String<T>::truncateSelf(Interval<int> const& a)
{
	Interval<int>itvl(0, length - 1);
	itvl ^= a;
	if (itvl.valid())truncateSelf(itvl.a, itvl.b - itvl.a + 1);
	return *this;
}
template<class T>					inline String<T>& String<T>::truncateSelf(IntervalSet<int> const& a)
{
	IntervalSet<int>tp(a ^ Interval<int>(0, length - 1));
	length = tp.area(true);
	if (!length)
	{
		this->~String();
		return *this;
	}
	int n(0);
	for (int c0(0); c0 < tp.length; ++c0)
	{
		int dn(tp.data[c0].b - tp.data[c0].a + 1);
		::memmove(data + n, data + tp.data[c0].a, dn * sizeof(T));
		n += dn;
	}
	data[length] = 0;
	lengthAll = length + 1;
	data = (T*)::realloc(data, lengthAll * sizeof(T));
	return *this;
}
//print
template<class T>					inline void String<T>::print()const
{
	if constexpr (IsSameType<T, char>::value)::printf("%s", data);
	else ::wprintf(L"%ls", data);
}
template<class T>					inline void String<T>::print(char const* a) const
{
	::printf("%s", a);
	print();
}
template<class T>					inline void String<T>::print(char const* a, char const* b) const
{
	::printf("%s", a);
	print();
	::printf("%s", b);
}
template<class T>					inline void String<T>::printInfo()const
{
	if constexpr (IsSameType<T, char>::value)::printf("[\"%s\", %u ,%u]\n", data, length, lengthAll);
	else ::wprintf(L"[\"%ls\", %u ,%u]\n", data, length, lengthAll);
}


//_Vector.h
template<>inline Vector<String<char>>& Vector<String<char>>::inverse()
{
	for (int c0(0); c0 < length / 2; ++c0)
		data[c0].exchange(data[length - c0 - 1]);
	return *this;
}
template<>inline Vector<String<wchar_t>>& Vector<String<wchar_t>>::inverse()
{
	for (int c0(0); c0 < length / 2; ++c0)
		data[c0].exchange(data[length - c0 - 1]);
	return *this;
}

//_Algorithm.h
template<class T>void Interval<T>::print()const
{
	String<char>tp("[");
	tp += NumType<T>::printInfo;
	tp += ", ";
	tp += NumType<T>::printInfo;
	tp += "]";
	::printf(tp.data, a, b);
}
template<class T>void Interval<T>::print(char const* p)const
{
	::printf("%s", p);
	print();
}
template<class T>inline void Interval<T>::print(char const* p, char const* q) const
{
	::printf("%s", p);
	print();
	::printf("%s", q);
}
