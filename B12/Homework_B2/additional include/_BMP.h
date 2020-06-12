#pragma once
#include <_File.h>

//Only support 24bit bmp!
struct BMP
{
#pragma pack(1)
	struct Header
	{
		char identifier[2];
		unsigned int fileSize;
		unsigned int reserved;
		unsigned int dataOffset;
		unsigned int headerSize;
		unsigned int width;
		unsigned int height;
		unsigned short planeNum;
		unsigned short bitsPerPixel;
		unsigned int compressionMethod;
		unsigned int dataSize;
		unsigned int verticalPixelsPerMeter;
		unsigned int horizontalPixelsPerMeter;
		unsigned int colorNum;
		unsigned int importantColorNum;
		void printInfo()const
		{
			::printf("BMP file info:\n");
			::printf("\tFile Size:\t%u\n", fileSize);
			::printf("\tBits Per Pixel:\t%u\n", bitsPerPixel);
			::printf("\tWidth:\t\t%u\n", width);
			::printf("\tHeight:\t\t%u\n", height);
		}
	};
	struct Pixel_24
	{
		unsigned char b;
		unsigned char g;
		unsigned char r;
		unsigned int bgr()
		{
			return b | (g << 8) | (r << 16) | (255 << 24);
		}
		unsigned int rgb()
		{
			return r | (g << 8) | (b << 16) | (255 << 24);
		}
	};
	struct Pixel_32
	{
		unsigned char r;
		unsigned char g;
		unsigned char b;
		unsigned char a;

		Pixel_32(Pixel_24 tp) : r(tp.r), g(tp.g), b(tp.b), a(0) {}
	};
#pragma pack()

	Header header;
	Pixel_24* data_24;
	Pixel_32* data_32;
	unsigned char* textureData;

	BMP()
		:
		header(),
		data_24(nullptr),
		data_32(nullptr),
		textureData(nullptr)
	{
	}
	BMP(String<char>const& _path)
		:
		data_24(nullptr),
		data_32(nullptr),
		textureData(nullptr)
	{
		FILE* temp(::fopen(_path.data, "rb+"));
		::fseek(temp, 0, SEEK_SET);
		::fread(&header, 1, 54, temp);
		::fseek(temp, header.dataOffset, SEEK_SET);
		if (header.width % 4)
		{
			data_24 = (BMP::Pixel_24*)::malloc(3u * header.width * header.height + 4);
			for (int c0(0); c0 < header.height; ++c0)
				::fread((data_24 + header.width * c0), 4, 1 + header.width * 3 / 4, temp);
		}
		else
		{
			data_24 = (BMP::Pixel_24*)::malloc(3u * header.width * header.height + 4);
			for (int c0(0); c0 < header.height; ++c0)
				::fread((data_24 + header.width * c0), 4, header.width * 3 / 4, temp);
		}
		::fclose(temp);
	}
	BMP(String<char>const& _path, bool is32bit)//!!!!!!!!!!!!
		:
		data_24(nullptr),
		data_32(nullptr),
		textureData(nullptr)
	{
		FILE* temp(::fopen(_path.data, "rb+"));
		::fseek(temp, 0, SEEK_SET);
		::fread(&header, 1, 54, temp);
		::fseek(temp, header.dataOffset, SEEK_SET);
		data_32 = (BMP::Pixel_32*)::malloc(sizeof(BMP::Pixel_32) * header.width * header.height);
		unsigned int rowSize((header.width * 3 + 3) / 4);
		for (int c0(0); c0 < header.height; ++c0)
		{
			Pixel_32* ptr(data_32 + c0 * header.width);
			::fread(ptr, 4, rowSize, temp);
			for (int c1(header.width - 1); c1 >= 0; --c1)
			{
				Pixel_24 ahh(*(((Pixel_24*)ptr) + c1));
				new(ptr + c1) Pixel_32(ahh);
			}
		}
		::fclose(temp);
	}
	~BMP()
	{
		::free(data_24);
		::free(data_32);
		data_24 = nullptr;
		data_32 = nullptr;
	}

	bool checkType()const
	{
		return header.identifier[0] == 'B' && header.identifier[1] == 'M';
	}
	void printInfo()const
	{
		header.printInfo();
	}
};
struct BMPCube
{
	BMP bmp[6];
	BMPCube(String<char>const& _path) :bmp{ _path + "front.bmp",_path + "back.bmp",_path + "down.bmp",_path + "up.bmp",_path + "right.bmp",_path + "left.bmp" } {}
};


//File...
inline BMP File::readBMP()const
{
	if (!this)return BMP();
	return BMP(property.path + property.file.name);
}
inline BMP File::readBMP(String<char> const& _name)const
{
	if (!this)return BMP();
	return BMP(property.path + _name);
}
inline BMP File::readBMP32bit()const
{
	if (!this)return BMP();
	return BMP(property.path + property.file.name, true);
}
inline BMP File::readBMP32bit(String<char> const& _name)const
{
	if (!this)return BMP();
	return BMP(property.path + _name, true);
}