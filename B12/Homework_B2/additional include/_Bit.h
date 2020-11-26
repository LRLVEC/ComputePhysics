#pragma once
#include <cstdio>



void printBit(void const* a, int bits, int num)
{
	char*ptr((char*)a);
	while (num--)
	{
		for (int c0((bits >> 3) - 1); c0 >= 0; --c0)
		{
			for (int c1(7); c1 >= 0; --c1)
				::printf("%d", (*(ptr + c0)&(1 << c1)) != 0);
		}
		::printf(" ");
		ptr += bits >> 3;
	}
	::printf("\n");
}
void printFloatBit(float const* a, int num)
{
	unsigned int*ptr((unsigned int*)a);
	while (num--)
	{
		::printf("%d ", (*ptr & 0x80000000U) != 0);
		for (int c1(30); c1 > 22; --c1)::printf("%d", (*ptr & (1U << c1)) != 0);
		::printf(" ");
		for (int c1(22); c1 >= 0; --c1)::printf("%d", (*ptr & (1U << c1)) != 0);
		::printf("\n");
		ptr++;
	}
}

void printFloatInfo(float const* a, int num)
{
	unsigned int*ptr((unsigned int*)a);
	while (num--)
	{
		if (*ptr & 0x80000000U)::printf("[-, ");
		else ::printf("[+, ");
		printf("%d, ", (char)((unsigned char)(*ptr >> 23) - 0x7fU));
		for (int c1(22); c1 >= 0; --c1)::printf("%d", (*ptr & (1U << c1)) != 0);
		::printf("]\n");
		ptr++;
	}
}
unsigned short floatToUnsignedShort(float a)
{
	unsigned int*ptr((unsigned int*)&a);
	return (unsigned short)((*ptr) >> 7);
}