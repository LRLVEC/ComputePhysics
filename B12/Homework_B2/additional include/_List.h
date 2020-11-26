#pragma once
#include <_Vector.h>
#include <typeinfo>
#include <initializer_list>
#include <new>

template<class T>struct List
{
	struct ListNode
	{
		T data;
		ListNode* pre;
		ListNode* suc;

		ListNode();
		ListNode(T const&);
		ListNode(ListNode*, T const&);
		ListNode(T const&, ListNode*);
		~ListNode();
	};

	ListNode* begin;
	ListNode* end;
	unsigned int length;

	List();
	List(T const&);
	List(List<T> const&);
	List(std::initializer_list<T>const&);


	~List();

	ListNode& operator[] (unsigned int);

	List<T>& pushBack(T const&);
	List<T>& popBack();
	List<T>& insert(T const&, unsigned int);
	List<T>& omit(unsigned int);

	template<class R>ListNode& find(R const&);
	ListNode& find(T const&);
	bool traverse(bool(*p)(T const&))const;
	bool check(bool(*p)(T const&));

	void printInfo()const;
	void printInfo(char const*, bool(*p)(T const&))const;
};

//ListNode
template<class T>inline List<T>::ListNode::ListNode()
	:
	data(),
	pre(nullptr),
	suc(nullptr)
{

}
template<class T>inline List<T>::ListNode::ListNode(T const& _data)
	:
	data(_data),
	pre(nullptr),
	suc(nullptr)
{
}
template<class T>inline List<T>::ListNode::ListNode(ListNode* _pre, T const& _data)
	:
	data(_data),
	pre(_pre)
{
	if (pre->suc)
	{
		suc = pre->suc;
		suc->pre = this;
		pre->suc = this;
	}
	else
	{
		pre->suc = this;
		suc = nullptr;
	}
}
template<class T>inline List<T>::ListNode::ListNode(T const& _data, ListNode* _suc)
	:
	data(_data),
	suc(_suc)
{
	if (suc->pre)
	{
		pre = suc->pre;
		pre->suc = this;
		suc->pre = this;
	}
	else
	{
		suc->pre = this;
		pre = nullptr;
	}
}

template<class T>inline List<T>::ListNode::~ListNode()
{
	if (pre)
	{
		if (suc)
		{
			pre->suc = suc;
			suc->pre = pre;
			suc = nullptr;
		}
		else
		{
			pre->suc = nullptr;
		}
		pre = nullptr;
	}
	else if (suc)
		suc->pre = nullptr;
}

//List
template<class T>inline List<T>::List()
	:
	begin(nullptr),
	end(nullptr),
	length(0)
{
}
template<class T>inline List<T>::List(T const& _data)
	:
	begin(new((ListNode*)::malloc(sizeof(ListNode)))ListNode(_data)),
	end(begin),
	length(1)
{
}
template<class T>inline List<T>::List(List<T> const& a)
	:
	length(a.length)
{
	ListNode* t(a.begin);
	if (t)begin = new((ListNode*)::malloc(sizeof(ListNode)))ListNode(a.begin->data);
	else
	{
		begin = nullptr;
		end = nullptr;
		return;
	};
	ListNode* k(begin);
	while (t = t->suc)
		k = new((ListNode*)::malloc(sizeof(ListNode)))ListNode(k, t->data);
	end = k;
}
template<class T>inline List<T>::List(std::initializer_list<T> const& a)
	:
	length(a.size())
{
	if (length)
	{
		unsigned int c0(0);
		T const* p(a.begin());
		begin = new((ListNode*)::malloc(sizeof(ListNode)))ListNode(*p);
		ListNode* t(begin);
		while (++c0 < length)
			t = new((ListNode*)::malloc(sizeof(ListNode)))ListNode(t, p[c0]);
		end = t;
	}
	else
	{
		begin = nullptr;
		end = nullptr;
	}
}

template<class T>inline List<T>::~List()
{
	while (begin)
	{
		(&(begin->data))->~T();
		ListNode* t(begin->suc);
		::free(begin);
		begin = t;
	}
}

template<class T>inline typename List<T>::ListNode& List<T>::operator[](unsigned int n)
{
	if (n >= length)
		return *(ListNode*)nullptr;
	ListNode * r;
	if (n <= (length << 1))
	{
		r = begin;
		while (n--)
			r = r->suc;
		return *r;
	}
	else
	{
		r = end;
		n = length - n;
		while (--n)
			r = r->pre;
		return *r;
	}
}

template<class T>inline List<T>& List<T>::pushBack(T const& _data)
{
	if (begin)
		end = new((ListNode*)::malloc(sizeof(ListNode)))ListNode(end, _data);
	else
		end = begin = new((ListNode*)::malloc(sizeof(ListNode)))ListNode(_data);
	++length;
	return *this;
}
template<class T>inline List<T>& List<T>::popBack()
{
	if (begin != end)
	{
		end = end->pre;
		end->suc->~ListNode();
		::free(end->suc);
		end->suc = nullptr;
		--length;
	}
	else
	{
		if (begin)
		{
			begin->~ListNode();
			::free(begin);
			begin = end = nullptr;
			--length;
		}
	}
	return *this;
}
template<class T>inline List<T>& List<T>::insert(T const& a, unsigned int n)
{
	if (length)
	{
		if (n)
		{
			if (length >= n)
				new((ListNode*)::malloc(sizeof(ListNode)))ListNode(&operator[](n - 1), a);
			else
				end = new((ListNode*)::malloc(sizeof(ListNode)))ListNode(end, a);
		}
		else
			new((ListNode*)::malloc(sizeof(ListNode)))ListNode(a, begin);
	}
	else
		begin = end = new((ListNode*)::malloc(sizeof(ListNode)))ListNode(a);
	++length;
	return *this;
}
template<class T>inline List<T>& List<T>::omit(unsigned int n)
{
	if (length)
	{
		if (length > 1)
		{
			if (n + 1 >= length)
			{
				end = end->pre;
				end->suc->~ListNode();
				::free(end->suc);
				end->suc = nullptr;
			}
			else
			{
				if (n)
				{
					ListNode& t(operator[](n));
					t.pre->suc = t.suc;
					t.suc->pre = t.pre;
					(&t)->~ListNode();
					::free(&t);
				}
				else
				{
					begin = begin->suc;
					begin->pre->~ListNode();
					::free(begin->pre);
				}
			}
		}
		else
		{
			begin->~ListNode();
			::free(begin);
			begin = end = nullptr;
		}
		--length;
	}
	return *this;
}
template<class T>template<class R>	inline typename List<T>::ListNode& List<T>::find(R const& a)
{
	ListNode* t(begin);
	while (t)
	{
		if (t->data == a)
			return *t;
		t = t->suc;
	}
	return *(ListNode*)nullptr;
}
template<class T>					inline typename List<T>::ListNode& List<T>::find(T const& a)
{
	ListNode* t(begin);
	while (t)
	{
		if (a == t->data)
			return *t;
		t = t->suc;
	}
	return *(ListNode*)nullptr;
}
template<class T>inline bool List<T>::traverse(bool(*p)(T const&)) const
{
	ListNode* t(begin);
	while (t)
	{
		if (!p(t->data))return false;
		t = t->suc;
	}
	return true;
}
template<class T>inline bool List<T>::check(bool(*p)(T const&))
{
	ListNode* t(begin);
	while (t)
	{
		if (!p(t->data))
		{
			ListNode* k(t->suc);
			if (t == begin)
			{
				if (t != end)begin = k;
				else begin = end = nullptr;
			}
			else if (t == end)
			{
				end = t->pre;
			}
			--length;
			t->~ListNode();
			::free(t);
			t = k;
			continue;
		}
		t = t->suc;
	}
	return false;
}

template<class T>inline void List<T>::printInfo() const
{
	::printf("[List<%s>, %u]\n", typeid(T).name(), length);
}
template<class T>inline void List<T>::printInfo(char const* a, bool(*p)(T const&)) const
{
	::printf("[%s: List<%s>, [%u]", a, typeid(T).name(), length);
	traverse(p);
	::printf("]\n");
}


