#pragma once
#include <cstdlib>
#include <new>
#include <_String.h>

namespace Expression
{
	enum NodeType
	{
		IsOprand,
		IsOperator
	};
	template<class T>struct Oprand
	{
		static constexpr const NodeType type = IsOprand;
		virtual T operator()()const = 0;
		virtual T operator=(T const&) = 0;
	};
	template<class T>struct Operator
	{
		static constexpr const NodeType type = IsOperator;
		virtual T operator()(Oprand<T>const*, Oprand<T>const*)const = 0;
	};

	// All memories are controlled inside.
	// R must has default construct function.
	template<class T, class R, class S>struct Node
	{
		Oprand<T>* oprd;
		Operator<T>* oprt;
		Node<T, R, S>* left, * right;
		static_assert(R::type == IsOprand, "R is not a correct Oprand<> struct!");
		static_assert(S::type == IsOperator, "S is not a correct Operator<> struct!");
		Node(T const& a)
			:
			oprd((R*)::malloc(sizeof(R))),
			oprt(nullptr),
			left(nullptr),
			right(nullptr)
		{
			new(oprd)R(a);
		}
		Node(T* a)
			:
			oprd((R*)::malloc(sizeof(R))),
			oprt(nullptr),
			left(nullptr),
			right(nullptr)
		{
			new(oprd)R(a);
		}
		template<class M>Node(M const& a)
			:
			oprd((R*)::malloc(sizeof(R))),
			oprt((S*)::malloc(sizeof(S))),
			left(nullptr),
			right(nullptr)
		{
			new(oprd)R;
			new(oprt)S(a);
		}
		~Node()
		{
			::free(oprd);
			::free(oprt);
			oprd = nullptr;
			oprt = nullptr;
			if (left) { left->~Node(); left = nullptr; }
			if (right) { right->~Node(); left = nullptr; }
		}
		Oprand<T> const& operator()()
		{
			if (left || right)(*oprd) = (*oprt)(left ? &(*left)() : nullptr, right ? &(*right)() : nullptr);
			return *oprd;
		}
	};
	//some assumptions: "()" 's meaning is just "()"; " " 's position doesn't change it's meaning
	struct Expr
	{
		String<char> expr;
		Vector<Expr> subExprs;
		Expr() {}
		Expr(String<char>const& a) :expr(a) {}
		bool checkParentheses()const
		{
			if (!expr.length)return false;
			if (expr.findFirst("(") == -1)return true;
			else return false;
		}
		Vector<IntervalSet<int>> getParentheses()const
		{
			if (!expr.length)return Vector<IntervalSet<int>>();
			Vector<int>a(expr.find("("));
			Vector<int>b(expr.find(")"));
			if (a.length != b.length || a.length == 0)return Vector<IntervalSet<int>>();
			int major(-1);
			int p(0), q(0);
			Vector<IntervalSet<int>>tp;
			bool last(true);//true: (, false: )
			while (p < a.length || q < a.length)
			{
				if (a.data[p] < b.data[q] && p < a.length)
				{
					if (last)++major; last = true;
				}
				else { if (!last)--major; last = false; }
				if (major < 0)return Vector<IntervalSet<int>>();
				if (major >= tp.length)tp.pushBack();
				if (last)tp.data[major].pushBack({ int(a.data[p++]), 0 });
				else tp.data[major].end().b = b.data[q++];
			}
			return tp;
		}
		/*Expr& preSimplify()
		{
			if (checkParentheses())
			{
				Vector<IntervalSet<int>>vi(getParentheses());
				int bg(0);
				for (int c0(0); c0 < vi[0].length; ++c0)
				{
					if (vi[0][c0].a > bg)subExprs.pushBack(expr.truncate(bg, vi[0][c0].a - bg));
					subExprs.pushBack(expr.truncate(vi[0][c0].a + 1, vi[0][c0].b - vi[0][c0].a - 1));
					subExprs.end().preSimplify(vi, 0, c0);
					bg = vi[0][c0].b + 1;
				}
				if (bg < expr.length)subExprs.pushBack(expr.truncate(bg, expr.length - bg));
			}
		}
		Expr& preSimplify(Vector<IntervalSet<int>>& const vi, int depth, int pos)
		{

		}*/
	};
}