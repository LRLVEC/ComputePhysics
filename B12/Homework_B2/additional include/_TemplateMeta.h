#pragma once


template<class T, T a, T b>struct Min
{
	static constexpr const T value = b > a ? a : b;
};
template<class T, T a, T b>struct Max
{
	static constexpr const T value = b < a ? a : b;
};
template<class T, T a, T b>struct Differ
{
	static constexpr const T value = a > b ? a - b : b - a;
};


template<class T, T _value>struct Constant
{
	static constexpr T value = _value;
};

template<bool _value>using BoolConstant = Constant<bool, _value>;

using True = Constant<bool, true>;
using False = Constant<bool, false>;

template<class T, class R>struct IsSameType : False {};
template<class T>struct IsSameType<T, T> :True {};

template<class T>struct BoolCondition
{
	using Result = T;
};
using TrueCondition = BoolCondition<True>;
using FalseCondition = BoolCondition<False>;

template<class Condition, class Then, class Else>struct If
{
};
template<class Then, class Else>struct If<True, Then, Else>
{
	using Result = Then;
};
template<class Then, class Else>struct If<False, Then, Else>
{
	using Result = Else;
};

template<unsigned int n>
struct F
{
	static constexpr const unsigned int a = F<n - 1>::a + F<n - 2>::a;
};
template<>struct F<0>
{
	static constexpr const unsigned int a = 1;
};
template<>struct F<1>
{
	static constexpr const unsigned int a = 1;
};

template<class Condition>struct Not
{
	using Result = If<Condition, FalseCondition, TrueCondition>;
};

//==============================================
enum NumericalType
{
	IsNonNumerical,
	IsInteger,
	IsFloat,
};
template<class T>struct NumType
{
	static constexpr bool value = false;
	static constexpr NumericalType numType = IsNonNumerical;
	static constexpr unsigned char serial = 0;
};
template<>struct NumType<char>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 1;
	static constexpr const char* printInfo = "%d";
};
template<>struct NumType<short>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 2;
	static constexpr const char* printInfo = "%d";
};
template<>struct NumType<int>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 3;
	static constexpr const char* printInfo = "%d";
};
template<>struct NumType<long>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 4;
	static constexpr const char* printInfo = "%d";
};
template<>struct NumType<long long>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 5;
	static constexpr const char* printInfo = "%ld";
};
template<>struct NumType<unsigned char>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 6;
	static constexpr const char* printInfo = "%u";
};
template<>struct NumType<unsigned short>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 7;
	static constexpr const char* printInfo = "%u";
};
template<>struct NumType<unsigned int>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 8;
	static constexpr const char* printInfo = "%u";
};
template<>struct NumType<unsigned long>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 9;
	static constexpr const char* printInfo = "%lu";
};
template<>struct NumType<unsigned long long>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsInteger;
	static constexpr unsigned char serial = 10;
	static constexpr const char* printInfo = "%llu";
};
template<>struct NumType<float>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsFloat;
	static constexpr unsigned char serial = 11;
	static constexpr const char* printInfo = "%f";
};
template<>struct NumType<double>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsFloat;
	static constexpr unsigned char serial = 12;
	static constexpr const char* printInfo = "%lf";
};
template<>struct NumType<long double>
{
	static constexpr bool value = true;
	static constexpr NumericalType numType = IsFloat;
	static constexpr unsigned char serial = 13;
	static constexpr const char* printInfo = "%lf";
};

template<unsigned char>struct GetNumType {};
template<>struct GetNumType< 1>
{
	using Result = char;
};
template<>struct GetNumType< 2>
{
	using Result = short;
};
template<>struct GetNumType< 3>
{
	using Result = int;
};
template<>struct GetNumType< 4>
{
	using Result = long;
};
template<>struct GetNumType< 5>
{
	using Result = long long;
};
template<>struct GetNumType< 6>
{
	using Result = unsigned char;
};
template<>struct GetNumType< 7>
{
	using Result = unsigned short;
};
template<>struct GetNumType< 8>
{
	using Result = unsigned int;
};
template<>struct GetNumType< 9>
{
	using Result = unsigned long;
};
template<>struct GetNumType<10>
{
	using Result = unsigned long long;
};
template<>struct GetNumType<11>
{
	using Result = float;
};
template<>struct GetNumType<12>
{
	using Result = double;
};
template<>struct GetNumType<13>
{
	using Result = long double;
};

static constexpr const unsigned char HigherNumTypeTable[14][14] =
{
	{0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
	{0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13},// 1: char
	{0,  2,  2,  3,  4,  5,  1,  7,  8,  9, 10, 11, 12, 13},// 2: short
	{0,  3,  3,  3,  4,  5,  3,  3,  8,  9, 10, 11, 12, 13},// 3: int
	{0,  4,  4,  4,  4,  5,  4,  4,  4,  9, 10, 11, 12, 13},// 4: long
	{0,  5,  5,  5,  5,  5,  5,  5,  5,  5, 10, 11, 12, 13},// 5: long long
	{0,  6,  1,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13},// 6: unsigned char
	{0,  7,  7,  3,  4,  5,  7,  7,  8,  9, 10, 11, 12, 13},// 7: unsigned short
	{0,  8,  8,  8,  4,  5,  8,  8,  8,  9, 10, 11, 12, 13},// 8: unsigned int
	{0,  9,  9,  9,  9,  5,  9,  9,  9,  9, 10, 11, 12, 13},// 9: unsigned long
	{0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 12, 13},//10: unsigned long long
	{0, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 13},//11: float
	{0, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13},//12: double
	{0, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13},//13: long souble
};
//==============================================
template<class T>struct CharType
{
	static constexpr bool value = false;
	static constexpr unsigned char serial = 0;
};

template<>struct CharType<char>
{
	static constexpr bool value = true;
	static constexpr unsigned char serial = 1;
	static constexpr const char* printInfo = "%s";
};
template<>struct CharType<wchar_t>
{
	static constexpr bool value = true;
	static constexpr unsigned char serial = 2;
	static constexpr const wchar_t* printInfo = L"%ls";
};


