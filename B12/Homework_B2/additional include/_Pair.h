#pragma once

template<class T, class R>struct Pair
{
	T data0;
	R data1;

	Pair() = default;
	Pair(T const&);
	Pair(T const&, R const&);
};


template<class T, class R>inline Pair<T, R>::Pair(T const& _data0)
	:
	data0(_data0),
	data1()
{
}
template<class T, class R>inline Pair<T, R>::Pair(T const& _data0, R const& _data1)
	:
	data0(_data0),
	data1(_data1)
{
}
