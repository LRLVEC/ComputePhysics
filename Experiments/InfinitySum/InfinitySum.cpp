#include <cstdio>

int main()
{
	float sum(0.0f);
	for (int c0(0); c0 < 1000; ++c0)
	{
		for (int c1(1); c1 < 10001; ++c1)
		{
			unsigned int id(c1 + c0 * 10000);
			float delta(1.0f / id);
			sum += delta;
			if (c0 == 209)
				::printf("\t%4d: %.8f\n", c1, sum);
		}
		::printf("%4d: %.8f\n", c0, sum);
	}
	//result: 2097151: 15.40368271
}