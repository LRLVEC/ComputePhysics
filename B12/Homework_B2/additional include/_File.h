#pragma once
#include <_Vector.h>
#include <_String.h>
#include <cstdlib>
#include <direct.h>
#include <io.h>
//#include <direct.h>

struct STL;
struct BMP;
namespace RayTracing
{
	struct Model;
}
//Support path format: [C:/] [./abc/] [../abc/]
struct File
{
	struct Property
	{
		bool isFolder;			//if this is a folder, true
		String<char> path;		//folder path
		__finddata64_t file;	//file infomation
		FILE* filePtr;			//FILE* pointer
		/*======================================================================*/
		Property();
		Property(String<char>const&);
		Property(__finddata64_t const&, String<char>const&);
		/*--------------------------------------------------------*/
		~Property();
	};

	//Read path into Vector<String<char>>
	static Vector<String<char>>readPath(String<char>const&);
	//Path Simplify
	static Vector<String<char>>pathSimplify(Vector<String<char>>const&);
	static Vector<String<char>>pathSimplify(String<char>const&);
	//Path translation: Write Vector<String<char>> into String<char>
	static String<char> pathTranslate(Vector<String<char>>const&);

	//Menbers
	bool valid;					//if this file valid, true
	Property property;			//file attribs
	Vector<File>childs;			//childs
	File* father;				//pointer to father

	//Construction
	File();
	File(String<char>const&);
	File(__finddata64_t const&, String<char>const&, File*);
	//Destruction
	~File();
	//Build a File
	void build();
	void build(String<char>const&);
	//Find in this folder, allow only name like [abc]
	File& findInThis(String<char>const&);
	//Find such a File
	File& find(String<char>const&);
	File& find(Vector<String<char>>const&);
	File& find(Vector<String<char>>const&, int);
	//Create file
	File& createDirectory(String<char>const&);
	File& createText(String<char>const&, String<char>const&);
	File& createBinary(String<char>const&, void*, unsigned int);
	File& createSTL(String<char>const&, STL const&);
	//Read file
	String<char>readText()const;
	String<char>readText(String<char>const&)const;
	STL readSTL()const;
	STL readSTL(String<char>const&)const;
	BMP readBMP()const;
	BMP readBMP(String<char>const&)const;
	BMP readBMP32bit()const;
	BMP readBMP32bit(String<char>const&)const;
	RayTracing::Model readModel()const;
	Vector<unsigned char> readBinary()const;
	Vector<unsigned char> readBinary(String<char>const&)const;
	//Print info
	void print()const;
};

/*==================================================static funtions================================================*/
//Read path into Vector<String<char>>
inline Vector<String<char>>File::readPath(String<char>const& _path)
{
	Vector<int>offset(_path.find("/"));
	Vector<String<char>>r;
	int tempPos(0);
	if (offset.length == 0)return Vector<String<char>>(_path);
	int n(0);
	while (n < offset.length)
	{
		r.pushBack(_path.truncate(tempPos, offset.data[n] - tempPos));
		if (offset.data[n] == -1)break;
		tempPos = offset.data[n++] + 1;
	}
	if (offset.data[n - 1] < (int)_path.length - 1)r.pushBack(_path.truncate(tempPos, -1));
	return r;
}
//Simplify path like ,/a/../
inline Vector<String<char>>File::pathSimplify(Vector<String<char>>const& _path)
{
	Vector<unsigned int>upper(_path.posAll(".."));
	Vector<unsigned int>current(_path.posAll("."));
	Vector<int>left(-1);
	int upperFlag(0), currentFlag(0);
	Vector<String<char>>r;
	int c1(0);
	while (_path.length - c1)
	{
		if (currentFlag < current.length && current.data[currentFlag] == c1)
		{
			++currentFlag;
			c1++;
			continue;
		}
		if (upperFlag < upper.length && upper.data[upperFlag] == c1)
		{
			++upperFlag;
			if (left.end() != -1 && !(r.end() == ".."))
			{
				r.popBack();
				left.popBack();
				c1++;
				continue;
			}
		}
		r.pushBack(_path.data[c1]);
		left.pushBack(c1++);
	}
	return r;
}
inline Vector<String<char>>File::pathSimplify(String<char>const& _path)
{
	Vector<String<char>>temp(File::readPath(_path));
	return File::pathSimplify(temp);
}
//Path translation: Write Vector<String<char>> into String<char>
inline String<char> File::pathTranslate(Vector<String<char>>const& a)
{
	if (!a.length)return String<char>();
	unsigned int _length(0), _lengthAll(1);
	for (int c1(0); c1 < a.length; c1++)_length += a.data[c1].length;
	_length += a.length;
	while (_lengthAll < _length + 1)_lengthAll <<= 1;
	char* temp((char*)::malloc(_lengthAll));
	temp[0] = 0;
	unsigned int posTemp(0);
	for (int c1(0); c1 < a.length; c1++)
	{
		::strcat(temp, a.data[c1].data);
		posTemp += a.data[c1].length;
		temp[posTemp] = '/';
		temp[++posTemp] = 0;
	}
	return String<char>(temp, _length, _lengthAll);
}
/*==================================================Property=========================================================*/
//Construction
inline File::Property::Property()
	:
	isFolder(false),
	path(),
	filePtr(nullptr)
{
}
inline File::Property::Property(String<char> const& _path)
	:
	isFolder(true),
	path(_path),
	filePtr(nullptr)
{
}
inline File::Property::Property(__finddata64_t const& _file, String<char> const& _path)
	:
	isFolder(_file.attrib& _A_SUBDIR),
	path(_path),
	file(_file),
	filePtr(nullptr)
{
}
//Destruction
inline File::Property::~Property()
{
	if (filePtr)fclose(filePtr);
	filePtr = nullptr;
}
/*==================================================File==============================================================*/
//Construction
inline File::File()
	:
	valid(false),
	property(),
	childs(),
	father(nullptr)
{
}
inline File::File(String<char>const& _path)
	:
	valid(false),
	property(_path),
	childs(),
	father(nullptr)
{
	build(_path);
}
inline File::File(__finddata64_t const& _file, String<char> const& _path, File* _father)
	:
	valid(true),
	property(_file, _path),
	childs(),
	father(_father)
{
	if (property.isFolder)build(_path);
}
//Destruction
inline File::~File()
{
}
//Build a File
inline void File::build()
{
	if (!property.path.data)return;
	intptr_t handle;
	__finddata64_t tempFileInfo;
	String<char>tempPath(property.path);
	handle = _findfirst64(property.path + "*.*", &tempFileInfo);
	if (handle == -1)return;
	bool flag(true);
	while (flag)
	{
		if (tempFileInfo.attrib & _A_SUBDIR)
			childs.pushBack(File(tempFileInfo, property.path + tempFileInfo.name + "/", this));
		else
			childs.pushBack(File(tempFileInfo, tempPath, this));
		flag = !_findnext64(handle, &tempFileInfo);
	}
	_findclose(handle);
}
inline void File::build(String<char>const& _path)
{
	intptr_t handle;
	__finddata64_t tempFileInfo;
	String<char>tempPath(_path);
	handle = _findfirst64(_path + "*.*", &tempFileInfo);
	if (handle == -1)return;
	bool flag(true);
	while (flag)
	{
		if (strcmp(tempFileInfo.name, ".") && strcmp(tempFileInfo.name, ".."))
		{
			if (tempFileInfo.attrib & _A_SUBDIR)
				childs.pushBack(File(tempFileInfo, _path + tempFileInfo.name + "/", this));
			else
				childs.pushBack(File(tempFileInfo, tempPath, this));
		}
		flag = !_findnext64(handle, &tempFileInfo);
	}
	_findclose(handle);
}
//Find in this folder, allow only name like [abc]
inline File& File::findInThis(String<char>const& a)
{
	for (int c1(0); c1 < childs.length; c1++)
	{
		if (a == childs.data[c1].property.file.name)return childs.data[c1];
	}
	return *(File*)nullptr;
}
//Find such a File
inline File& File::find(String<char>const& _path)
{
	return find(File::pathSimplify(File::readPath(_path)));
}
inline File& File::find(Vector<String<char>>const& _path)
{
	return find(_path, 0);
}
inline File& File::find(Vector<String<char>>const& _path, int n)
{
	if (n < _path.length)
	{
		if (_path.data[n] == ".")return find(_path, n + 1);
		if (_path.data[n] == "..")
		{
			if (father)return father->find(_path, n + 1);
			else return *(File*)nullptr;
		}
		File& temp(findInThis(_path.data[n]));
		if (&temp)return temp.find(_path, n + 1);
		else return *(File*)nullptr;
	}
	else if (n == _path.length)return *this;
	return *(File*)nullptr;
}
//Create file
inline File& File::createDirectory(String<char>const& _path)
{
	::_mkdir((property.path + _path).data);
	build();
	return *this;
}
inline File& File::createText(String<char>const& _name, String<char>const& _text)
{
	FILE* temp(::fopen((property.path + _name).data, "w+"));
	::fprintf(temp, "%s", _text.data);
	::fclose(temp);
	build();
	return *this;
}
inline File& File::createBinary(String<char>const& _name, void* _data, unsigned int _length)
{
	FILE* temp(::fopen((property.path + _name).data, "wb+"));
	::fwrite(_data, 1, _length, temp);
	::fclose(temp);
	build();
	return *this;
}
//Read file
inline String<char>File::readText()const
{
	if (!this)return String<char>();
	FILE* temp(::fopen((property.path + property.file.name).data, "r+"));
	::fseek(temp, 0, SEEK_END);
	unsigned int _length((unsigned int)::ftell(temp) + 1);
	char* r((char*)::malloc(_length + 1));
	::fseek(temp, 0, SEEK_SET);
	unsigned int __length = (unsigned int)fread(r, 1, _length, temp);
	r[__length] = 0;
	::fclose(temp);
	return String<char>(r, __length, _length + 1);
}
inline String<char>File::readText(String<char>const& _name)const
{
	if (!this)return String<char>();
	FILE* temp(::fopen((property.path + _name).data, "r+"));
	::fseek(temp, 0, SEEK_END);
	unsigned int _length((unsigned int)::ftell(temp) + 1);
	char* r((char*)::malloc(_length + 1));
	::fseek(temp, 0, SEEK_SET);
	_length = (unsigned int)fread(r, 1, _length, temp);
	r[_length] = 0;
	::fclose(temp);
	return String<char>(r, _length, 0);
}
inline Vector<unsigned char> File::readBinary()const
{
	if (!this)return Vector<unsigned char>();
	FILE* temp(::fopen((property.path + property.file.name).data, "r+"));
	::fseek(temp, 0, SEEK_END);
	unsigned int _length((unsigned int)::ftell(temp) + 1);
	unsigned char* r((unsigned char*)::malloc(_length + 1));
	::fseek(temp, 0, SEEK_SET);
	unsigned int __length = (unsigned int)fread(r, 1, _length, temp);
	::fclose(temp);
	return Vector<unsigned char>(r, __length, _length);
}
inline Vector<unsigned char> File::readBinary(String<char>const& _name)const
{
	if (!this)return Vector<unsigned char>();
	FILE* temp(::fopen((property.path + _name).data, "r+"));
	::fseek(temp, 0, SEEK_END);
	unsigned int _length((unsigned int)::ftell(temp) + 1);
	unsigned char* r((unsigned char*)::malloc(_length + 1));
	::fseek(temp, 0, SEEK_SET);
	unsigned int __length = (unsigned int)fread(r, 1, _length, temp);
	::fclose(temp);
	return Vector<unsigned char>(r, __length, _length);
}

//Print info
inline void File::print()const
{
	if (property.path.data)property.path.print();
	if (!property.isFolder)printf("%s", property.file.name);
	::printf("\n");
	for (int c1(0); c1 < childs.length; c1++)childs.data[c1].print();
}