#ifndef _MESH_H_
#define _MESH_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>

typedef unsigned int uint;

typedef struct {
	//3d coordinates
	float x;
	float y;
	float z;
} Vertex;

typedef struct{
	//three vertex ids
	uint a,b,c;
} Face;

typedef struct{
	//2 vertex ids
	uint a,b;
} Edge;

class Mesh{
private:
	std::vector<Vertex> V;
	std::vector<Face> F;
	std::vector<Edge> E;
	std::set<std::pair<uint, uint>> Es;			//set of edges
	std::vector<std::vector<int>> VFMapping;	//vertex to faces mapping, 1 vector for each vertex
	std::set<int> FToRemove;					//faces to remove after mesh simplification
	std::vector<Face> FTemp;					//faces connected to both vertices in an edge to be collapsed
	std::set<int> PassedVertices;				//vertices part of edges selected to be collapsed
	std::set<int> PassedFaces;					//faces connected to edges selected to be collapsed

public:
	Mesh() {};
	Mesh(const char*);
	//load a Mesh from .mesh file
	void loadMF(const char*);
	//write a Mesh to .mesh file (no header)
	void writeMF(const char*);
	//simplify a mesh
	void simplifyMesh(const char* input, const char* output, int faceCnt);
	//return vertex count
	int Vcnt();
	//return face count
	int Fcnt();
	//replaces the E vector with the correct Edge values after 1 pass of simplifying the mesh
	void recomputeEdgesVFMapping();
};
#endif