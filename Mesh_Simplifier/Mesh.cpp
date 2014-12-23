#include "Mesh.h"
#include <set>

// Indexes of vectors = number-id of the object - 1
// Value inside the vectors = number-id of the object
// Objects refer to Vertex, Face, Edge, VFMapping

Mesh::Mesh(const char* filename){	
	loadMF(filename);
}

void Mesh::loadMF(const char* filename){
	if(V.size()>0) V.clear();
	if(F.size()>0) F.clear();
	std::ifstream infile;
	infile.open(filename, std::ios::in);
	std::string strbuff;
	while(std::getline(infile,strbuff)){
		std::stringstream ss;
		ss<<strbuff;
		char type;
		ss>>type;
		if(type=='v'){
			Vertex v;
			ss>>v.x>>v.y>>v.z;
			V.push_back(v);

			/* New code begin */
			std::vector<int> VFMappingSingle;
			VFMapping.push_back(VFMappingSingle);
			/* New code end */
		}
		else if(type=='f'){
			Face f;
			ss>>f.a>>f.b>>f.c;
			F.push_back(f);

			/* New code begin */
			VFMapping.at(f.a-1).push_back(F.size());	// For each vertex, track the faces which use it
			VFMapping.at(f.b-1).push_back(F.size());
			VFMapping.at(f.c-1).push_back(F.size());

			Edge e, e1, e2;								// 3 edges which make up the current face
			e.a = f.a; e.b = f.b;						// An edge comprise of 2 vertices
			e1.a = f.a; e1.b = f.c;
			e2.a = f.b; e2.b = f.c;

			std::pair<uint, uint> newP, newP1, newP2, newPOpp, newPOpp1, newPOpp2;
			newP.first = e.a;   newPOpp.first = e.b;
			newP.second= e.b;   newPOpp.second = e.a;
			newP1.first = e1.a; newPOpp1.first = e1.b;
			newP1.second= e1.b; newPOpp1.second = e1.a;
			newP2.first = e2.a; newPOpp2.first = e2.b;
			newP2.second= e2.b; newPOpp2.second = e2.a;

			if(Es.find(newP) == Es.end() && Es.find(newPOpp) == Es.end()) {    // Check if edge(a, b) has already been added. Also edge(b, a) is checked since it is the same
				Es.insert(newP);
				E.push_back(e);												   // Track all the edges being used
			}
			if(Es.find(newP1) == Es.end() && Es.find(newPOpp1) == Es.end()) {
				Es.insert(newP1);
				E.push_back(e1);
			}			
			if(Es.find(newP2) == Es.end() && Es.find(newPOpp2) == Es.end()) {
				Es.insert(newP2);
				E.push_back(e2);
			}
			/* New code end */
		}
	}
	infile.close();
}

void Mesh::writeMF(const char* filename){
	std::ofstream outfile;
	outfile.open(filename, std::ios::out);
	std::string strbuff;
	for(uint i=0;i<V.size();i++){
		outfile<<"v "<<V[i].x<<" "<<V[i].y<<" "<<V[i].z<<std::endl;
	}
	for(uint i=0;i<F.size();i++){
		outfile<<"f "<<F[i].a<<" "<<F[i].b<<" "<<F[i].c<<std::endl;
	}
	outfile.close();
}

// Used during 2nd onwards pass through of the mesh during simplification
void Mesh::recomputeEdgesVFMapping(){
	for(int i=0; i<V.size(); ++i)
	{
		std::vector<int> emptyVect;
		VFMapping.push_back(emptyVect);
	}
	//Recompute vertex to face mapping
	for(int i=0; i<F.size(); ++i)
	{
		Face f = F.at(i);
		VFMapping.at(f.a-1).push_back(i+1);
		VFMapping.at(f.b-1).push_back(i+1);
		VFMapping.at(f.c-1).push_back(i+1);

		//For each face add its 3 edges. Don't add duplicates
		Edge e, e1, e2;
		e.a = f.a; e.b = f.b;
		e1.a = f.a; e1.b = f.c;
		e2.a = f.b; e2.b = f.c;

		std::pair<uint, uint> newP, newP1, newP2, newPOpp, newPOpp1, newPOpp2;
		newP.first = e.a;   newPOpp.first = e.b;
		newP.second= e.b;   newPOpp.second = e.a;
		newP1.first = e1.a; newPOpp1.first = e1.b;
		newP1.second= e1.b; newPOpp1.second = e1.a;
		newP2.first = e2.a; newPOpp2.first = e2.b;
		newP2.second= e2.b; newPOpp2.second = e2.a;

		if(Es.find(newP) == Es.end() && Es.find(newPOpp) == Es.end()) {
			Es.insert(newP);
			E.push_back(e);
		}
		if(Es.find(newP1) == Es.end() && Es.find(newPOpp1) == Es.end()) {
			Es.insert(newP1);
			E.push_back(e1);
		}			
		if(Es.find(newP2) == Es.end() && Es.find(newPOpp2) == Es.end()) {
			Es.insert(newP2);
			E.push_back(e2);
		}
	}
}

// This function uses Edge Collapse to simplify a mesh until it has faceCnt number of faces
void Mesh::simplifyMesh(const char* input, const char* output, int faceCnt){
	// you may assume inputs are always valid
	loadMF(input);
	bool firstPass = true;
	while(Fcnt() > faceCnt)		// Loops until the desired faceCnt is achieved. Loops required to prevent collapsing edges which have overlapping connected faces
	{
		if(!firstPass)
			recomputeEdgesVFMapping();

		int edgeCnt = (Fcnt() - faceCnt)/2;		// Take faceCnt divided by 2 number of edges that are geographically seperated. edgeCnt is number of edges to remove
		int iIncrement;							// Used to determine number of edges to skip when traversing edge vector
		double iIncrementDouble = (double)E.size()/edgeCnt;
		if(iIncrementDouble - (int)iIncrementDouble >= 0.5)
			iIncrement = (int)iIncrementDouble + 1;
		else
			iIncrement = (int)iIncrementDouble;

		// For each of those edges take its midpoint and create a new vertex
		for(int i=0; i<E.size(); i += iIncrement)
		{
			if(edgeCnt == 0)	// Guard to ensure correct final face count
				break;
			Edge e = E.at(i);
			if(PassedVertices.find(e.a) != PassedVertices.end() || PassedVertices.find(e.b) != PassedVertices.end())	// If the edge has already been collapsed skip it. Referring to vertices which no longer exist
				continue;

			// Get all faces connected to this edge and add it into passed faces, to avoid choosing 2 edges that have overlapping connected faces
			std::vector<int> faces1 = VFMapping.at(e.a-1);
			std::vector<int> faces2 = VFMapping.at(e.b-1);
			bool facePassed = false;
			for(int i=0; i<faces1.size(); ++i)
			{
				if(PassedFaces.find(faces1.at(i)) != PassedFaces.end())
				{
					facePassed = true;
					break;
				}
			}
			if(facePassed)
				continue;

			for(int i=0; i<faces2.size(); ++i)
			{
				if(PassedFaces.find(faces2.at(i)) != PassedFaces.end())
				{
					facePassed = true;
					break;
				}
			}
			if(facePassed)
				continue;

			PassedVertices.insert(e.a);
			PassedVertices.insert(e.b);
			for (int i = 0; i < faces1.size(); i++)
			{
				PassedFaces.insert(faces1.at(i));
			}
			for (int i = 0; i < faces2.size(); i++)
			{
				PassedFaces.insert(faces2.at(i));
			}

			Vertex v1 = V.at(e.a-1);
			Vertex v2 = V.at(e.b-1);
			Vertex midpoint;
			midpoint.x = (v1.x + v2.x)/2;
			midpoint.y = (v1.y + v2.y)/2;
			midpoint.z = (v1.z + v2.z)/2;
			// The vertices have to remain with the same indices, since you don't want to affect the indices of the other vertices which are referenced by Face in F
			// Therefore we add to the end of the Vertex vector
			V.push_back(midpoint);

			// Get the 2 faces which contain both vertices from the edge to-be-deleted, for deletion in the future and for ignoring during the vertex replacement step
			std::vector<int> FToRemovePartial;
			for (int i = 0; i < faces1.size(); i++)
			{
				if(F.at(faces1.at(i)-1).a == e.a && ( F.at(faces1.at(i)-1).b == e.b || F.at(faces1.at(i)-1).c == e.b ))
				{
					FToRemovePartial.push_back(faces1.at(i));
					FToRemove.insert(faces1.at(i));
				}
				else if(F.at(faces1.at(i)-1).b == e.a && ( F.at(faces1.at(i)-1).a == e.b || F.at(faces1.at(i)-1).c == e.b ))
				{
					FToRemovePartial.push_back(faces1.at(i));
					FToRemove.insert(faces1.at(i));
				}
				else if(F.at(faces1.at(i)-1).c == e.a && ( F.at(faces1.at(i)-1).a == e.b || F.at(faces1.at(i)-1).b == e.b ))
				{
					FToRemovePartial.push_back(faces1.at(i));
					FToRemove.insert(faces1.at(i));
				}
			}

			// Connect the vertices connected to the to-be-deleted edge to this new vertex instead
			// Do so by retrieving the faces containing the to-be-deleted vertices and replacing the old vertices with the new vertex
			for (int i = 0; i < faces1.size(); i++)
			{
				bool skip = false;
				for(int j=0; j<FToRemovePartial.size(); ++j)
				{
					if(FToRemovePartial.at(j) == faces1.at(i))
					{
						skip = true;
						break;
					}
				}
				if(skip)
					continue;

				// Replacement
				if(F.at(faces1.at(i)-1).a == e.a)
					F.at(faces1.at(i)-1).a = V.size();	
				else if(F.at(faces1.at(i)-1).b == e.a)
					F.at(faces1.at(i)-1).b = V.size();
				else if(F.at(faces1.at(i)-1).c == e.a)
					F.at(faces1.at(i)-1).c = V.size();
			}

			for (int i = 0; i < faces2.size(); i++)
			{
				bool skip = false;
				for(int j=0; j<FToRemovePartial.size(); ++j)
				{
					if(FToRemovePartial.at(j) == faces2.at(i))
					{
						skip = true;
						break;
					}
				}
				if(skip)
					continue;

				// Replacement
				if(F.at(faces2.at(i)-1).a == e.b)
					F.at(faces2.at(i)-1).a = V.size();
				else if(F.at(faces2.at(i)-1).b == e.b)
					F.at(faces2.at(i)-1).b = V.size();
				else if(F.at(faces2.at(i)-1).c == e.b)
					F.at(faces2.at(i)-1).c = V.size();
			}
			--edgeCnt;
		}

		// The faces directly connected to the edges that were were collapsed will need to be deleted
		for (int i = 0; i < F.size(); i++)
		{
			if(FToRemove.find(i+1) == FToRemove.end())
				FTemp.push_back( F.at(i) );
		}
		F.clear();
		F = FTemp;

		// Cleanup
		FTemp.clear();
		PassedVertices.clear();
		PassedFaces.clear();
		FToRemove.clear();
		E.clear();
		Es.clear();
		VFMapping.clear();
		firstPass = false;
	}
	writeMF(output);

	std::cout<<"Original face count: "<<F.size()<<std::endl;
}

int Mesh::Vcnt(){
	return V.size();
}

int Mesh::Fcnt(){
	return F.size();
}