#pragma once
#include <cstdio>
#include <_Vector.h>
#include <_String.h>
#include <_File.h>
#include <_Math.h>
#include <GL/_OpenGL.h>
#include <RayTracing/_RayTracing.h>
#include <unordered_map>
#ifdef __OptiX__
#include <OptiX/_OptiX.h>
#endif
struct STL
{
#pragma pack(2)
	struct Triangle
	{
		Math::vec3<float>normal;
		Math::mat3<float>vertices;
		unsigned short attrib;
		bool operator==(Triangle const&);
		double getMinEdgeLength()const;
		void print()const;
	};
#pragma pack()

	String<char>name;
	Vector<Triangle>triangles;

	bool verticesUpdate;
	Vector<Math::vec3<float>>vertices;
	bool verticesRepeatedUpdate;
	Vector<Math::vec3<float>>verticesRepeated;
	Vector<Math::vec3<float>>normals;


	STL();
	STL(String<char>const&);
	STL(STL const&);

	Vector<Math::vec3<float>>& getVertices();
	Vector<Math::vec3<float>>& getVerticesRepeated();
	Vector<Math::vec3<float>>& getNormals();
	void removeUseless();
	double getMinTriangleScale();
	void printInfo(bool)const;
};
struct _Vertex
{
	Math::vec3<float>p;
	unsigned int id;
	_Vertex(Math::vec3<float>const& _p) :p(_p)
	{
	}
	_Vertex(Math::vec3<float>const& _p, unsigned int _id) :p(_p), id(_id)
	{
	}
	bool operator==(_Vertex const& a)const
	{
		return p == a.p;
	}
};
namespace std
{
	template<>struct hash<_Vertex>
	{
		std::size_t operator()(_Vertex const& p) const
		{
			return ((hash<float>()(p.p.data[0]) ^
				(hash<float>()(p.p.data[1]) << 1)) >> 1) ^
				(hash<float>()(p.p.data[2]) << 1);
		}
	};
}
struct STLIndex
{
	using Triangle = Math::vec3<unsigned int>;
	struct Vertex
	{
		struct TriangleId
		{
			unsigned int id;
			unsigned char p;
		};
		Math::vec3<float>position;
		Vector<TriangleId>ids;
		Math::vec3<float>normal;
		void avergeNormal()
		{
			normal = (normal / ids.length).normaliaze();
		}
	};
	Vector<Triangle>triangles;
	Vector<Vertex>vertices;
	STLIndex(STL const& stl)
	{
		build(stl);
	}
	void build(STL const& stl)
	{
		triangles.malloc(stl.triangles.length);
		std::unordered_map<_Vertex, unsigned int>map;
		for (unsigned int c0(0); c0 < stl.triangles.length * 3; ++c0)
		{
			Math::vec3<float> const& p(stl.triangles.data[c0 / 3].vertices.rowVec[c0 % 3]);
			auto it(map.find(p));
			if (it == map.end())
			{
				map.insert(std::make_pair(_Vertex(p, vertices.length), c0));
				triangles[c0 / 3][c0 % 3] = vertices.length;
				vertices.pushBack({ p, {{ c0 / 3,c0 % 3 }}, stl.triangles.data[c0 / 3].normal });
			}
			else
			{
				triangles[c0 / 3][c0 % 3] = it->first.id;
				vertices[it->first.id].ids.pushBack({ c0 / 3,c0 % 3 });
				vertices[it->first.id].normal += stl.triangles.data[c0 / 3].normal;
			}
		}
		vertices.traverse([](Vertex& a) { a.avergeNormal(); return true; });
	}
};

namespace OpenGL
{
	struct STLVertices :Buffer::Data
	{
		STL* stl;
		virtual void* pointer()override
		{
			return stl->verticesRepeated.data;
		}
		virtual unsigned int size()override
		{
			return sizeof(Math::vec3<float>) * stl->verticesRepeated.length;
		}
		STLVertices() = default;
		STLVertices(STL*);
		STLVertices(STLVertices const&) = default;
	};
	struct STLNormals :Buffer::Data
	{
		STL* stl;
		virtual void* pointer()override
		{
			return stl->normals.data;
		}
		virtual unsigned int size()override
		{
			return sizeof(Math::vec4<float>) * stl->triangles.length;
		}
		STLNormals() = default;
		STLNormals(STL*);
		STLNormals(STLNormals const&) = default;
	};


	inline STLVertices::STLVertices(STL* _stl)
		:
		stl(_stl)
	{
	}
	inline STLNormals::STLNormals(STL* _stl)
		:
		stl(_stl)
	{
	}
}

namespace RayTracing
{
	inline void Model::addSTL(STL const& _stl, Color const& _color, unsigned int num)
	{
		num = num < _stl.triangles.length ? num : _stl.triangles.length;
		for (int c0(0); c0 < num; ++c0)
		{
			vec3 d0(_stl.triangles.data[c0].vertices.rowVec[1] - _stl.triangles.data[c0].vertices.rowVec[0]);
			vec3 d1(_stl.triangles.data[c0].vertices.rowVec[2] - _stl.triangles.data[c0].vertices.rowVec[0]);
			vec3 n(d0 | d1);
			if ((n, _stl.triangles.data[c0].normal) > 0)
			{
				triangles.trianglesOrigin.trianglesOrigin.pushBack
				({
					_stl.triangles.data[c0].vertices,
					{ 0,0 },
					{ 1,0 },
					{ 0,1 },
					-1,
					_color
					});
			}
			else
			{
				triangles.trianglesOrigin.trianglesOrigin.pushBack
				({
					{
						_stl.triangles.data[c0].vertices.rowVec[0],
						_stl.triangles.data[c0].vertices.rowVec[2],
						_stl.triangles.data[c0].vertices.rowVec[1],
					},
					{ 0,0 },
					{ 1,0 },
					{ 0,1 },
					-1,
					_color
					});
			}
		}
	}
}
#ifdef __OptiX__
namespace OpenGL
{
	namespace OptiX
	{
		//Assuming that the format is linear compact layout...
		void GeometryTriangles::addSTL(String<char>const& name, STL const& stl, unsigned int num)
		{
			num = num < stl.triangles.length ? num : stl.triangles.length;
			Buffer* vertexBuffer(&buffers[name].buffer);
			unsigned long long size(vertexBuffer->getSize1());
			vertexBuffer->setSize(size + num * 3);
			Math::mat3<float>* a((Math::mat3<float>*)vertexBuffer->map() + size / 3);
			for (int c0(0); c0 < num; ++c0)
			{
				Math::vec3<float> d0(stl.triangles.data[c0].vertices.rowVec[1] - stl.triangles.data[c0].vertices.rowVec[0]);
				Math::vec3<float> d1(stl.triangles.data[c0].vertices.rowVec[2] - stl.triangles.data[c0].vertices.rowVec[0]);
				Math::vec3<float> n(d0 | d1);
				if ((n, stl.triangles.data[c0].normal) > 0)
					a[c0] = stl.triangles.data[c0].vertices;
				else
					a[c0] = Math::mat3<float>{
						stl.triangles.data[c0].vertices.rowVec[0],
						stl.triangles.data[c0].vertices.rowVec[2],
						stl.triangles.data[c0].vertices.rowVec[1] };
			}
			vertexBuffer->unmap();
			setCount(count + num);
			setVertices(vertexBuffer, size + num * 3, 0, 12);
		}
		void GeometryTriangles::addSTL(String<char>const& vertices, String<char>const& normals, String<char>const& indices, STL const& stl)
		{
			STLIndex ahh(stl);
			Buffer* vertexBuffer(&buffers[vertices].buffer);//float3
			Buffer* normalBuffer(&buffers[normals].buffer);//float3
			Buffer* indexBuffer(&buffers[indices].buffer);//uint3
			vertexBuffer->setSize(ahh.vertices.length);
			normalBuffer->setSize(ahh.vertices.length);
			indexBuffer->setSize(ahh.triangles.length);
			Math::vec3<float>* v((Math::vec3<float>*)vertexBuffer->map());
			Math::vec3<float>* n((Math::vec3<float>*)normalBuffer->map());
			Math::vec3<unsigned int>* i((Math::vec3<unsigned int>*)indexBuffer->map());
			for (unsigned int c0(0); c0 < ahh.vertices.length; ++c0)
			{
				v[c0] = ahh.vertices[c0].position;
				n[c0] = ahh.vertices[c0].normal;
			}
			for (unsigned int c0(0); c0 < ahh.triangles.length; ++c0)
				i[c0] = ahh.triangles[c0];
			vertexBuffer->unmap();
			normalBuffer->unmap();
			indexBuffer->unmap();
			setCount(ahh.triangles.length);
			setVertices(vertexBuffer, ahh.vertices.length, 0, 12);
			setIndices(indexBuffer, 0, 12);
		}
	}
}
#endif

inline bool STL::Triangle::operator==(Triangle const& a)
{
	return vertices == a.vertices;
}
inline double STL::Triangle::getMinEdgeLength() const
{
	double t((vertices.rowVec[0] - vertices.rowVec[1]).square());
	double s((vertices.rowVec[1] - vertices.rowVec[2]).square());
	double k((vertices.rowVec[2] - vertices.rowVec[0]).square());
	double n(t <= k ? t : k);
	return (n <= s ? n : s);
}
inline void STL::Triangle::print() const
{
	::printf("[");
	normal.print();
	::printf(", [");
	vertices.rowVec[0].print();
	::printf(", ");
	vertices.rowVec[1].print();
	::printf(", ");
	vertices.rowVec[2].print();
	::printf("], %u]\n", attrib);
}

inline STL::STL()
	:
	name(),
	triangles(),
	verticesUpdate(false),
	vertices(),
	verticesRepeatedUpdate(false),
	verticesRepeated()
{
}
inline STL::STL(String<char> const& _name)
	:
	name(_name),
	triangles(),
	verticesUpdate(false),
	vertices(),
	verticesRepeatedUpdate(false),
	verticesRepeated()
{
}
inline STL::STL(STL const& a)
	:
	name(a.name),
	triangles(a.triangles),
	verticesUpdate(a.verticesUpdate),
	vertices(a.vertices),
	verticesRepeatedUpdate(a.verticesRepeatedUpdate),
	verticesRepeated(a.verticesRepeated)
{
}



inline Vector<Math::vec3<float>>& STL::getVertices()
{

}
inline Vector<Math::vec3<float>>& STL::getVerticesRepeated()
{
	verticesRepeated.malloc(triangles.length * 3);
	verticesRepeated.length = 3 * triangles.length;
	for (int c0(0); c0 < triangles.length; ++c0)
	{
		verticesRepeated.data[3 * c0] = triangles.data[c0].vertices.rowVec[0];
		verticesRepeated.data[3 * c0 + 1] = triangles.data[c0].vertices.rowVec[1];
		verticesRepeated.data[3 * c0 + 2] = triangles.data[c0].vertices.rowVec[2];
	}
	return verticesRepeated;
}
inline Vector<Math::vec3<float>>& STL::getNormals()
{
	normals.malloc(triangles.length);
	normals.length = triangles.length;
	for (int c0(0); c0 < triangles.length; ++c0)
		normals.data[c0] = triangles.data[c0].normal;
	return normals;
}
inline void STL::removeUseless()
{
	for (int c0(0); c0 < triangles.length; ++c0)
		if (triangles[c0].getMinEdgeLength() == 0.0)
			triangles.omit(c0--);
}
inline double STL::getMinTriangleScale()
{
	if (!triangles.length)return 0.0;
	double t(triangles[0].getMinEdgeLength());
	for (int c0(1); c0 < triangles.length; ++c0)
	{
		double s(triangles[c0].getMinEdgeLength());
		if (t > s)t = s;
	}
	return sqrt(t);
}
inline void STL::printInfo(bool a) const
{
	::printf("[");
	name.print();
	::printf(": num: %d, vertices(repeated) num: %d, normals num: %d]\n",
		triangles.length, verticesRepeated.length, normals.length);
	if (a)
		triangles.traverse([](Triangle const& a)
			{
				a.print();
				return true;
			});
}

//File...
inline File& File::createSTL(String<char>const& _name, STL const& _stl)
{
	FILE* temp(::fopen((property.path + _name).data, "wb+"));
	::fwrite(_stl.name, 1, _stl.name.length + 1, temp);
	::fseek(temp, 80, SEEK_SET);
	::fwrite(&_stl.triangles.length, 4, 1, temp);
	::fwrite(_stl.triangles.data, 1, _stl.triangles.length * 50, temp);
	::fclose(temp);
	build();
	return *this;
}
inline STL File::readSTL() const
{
	if (!this)return String<char>();
	FILE* temp(::fopen((property.path + property.file.name).data, "rb+"));
	::fseek(temp, 0, SEEK_SET);
	char t[100];
	::fread(t, 1, 80, temp);
	STL r(property.file.name);
	unsigned int _num;
	::fread(&_num, 4, 1, temp);
	int c(1);
	while (c < _num)c <<= 1;
	r.triangles.data = (STL::Triangle*)::malloc(sizeof(STL::Triangle) * c);
	r.triangles.length = _num;
	r.triangles.lengthAll = c;
	::fread(r.triangles.data, sizeof(STL::Triangle), _num, temp);
	::fclose(temp);
	return r;
}
inline STL File::readSTL(String<char> const& _name) const
{
	if (!this)return String<char>();
	FILE* temp(::fopen((property.path + _name).data, "rb+"));
	::fseek(temp, 0, SEEK_SET);
	char t[100];
	::fread(t, 1, 80, temp);
	STL r(t);
	unsigned int _num;
	::fread(&_num, 4, 1, temp);
	int c(1);
	while (c < _num)c <<= 1;
	r.triangles.data = (STL::Triangle*)::malloc(sizeof(STL::Triangle) * c);
	r.triangles.length = _num;
	r.triangles.lengthAll = c;
	::printf("%d\n", ::fread(r.triangles.data, sizeof(STL::Triangle), _num, temp));
	::printf("%d\n", ::ftell(temp));
	::fclose(temp);
	return r;
}

