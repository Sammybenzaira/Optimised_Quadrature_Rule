#ifndef polygon_header
#define polygon_header

//////////////////////////////////////////////////////////////
//Module that implements a triangulation of various 
//polygons
//////////////////////////////////////////////////////////////

#include "linear_algebra.hpp"

//////////////////////////////////////////////////////////////
//Struct Definitions
//////////////////////////////////////////////////////////////

//General Triangle - just stores the vertices of 
//a triangle
;

struct Triangle
{
	double vertices[3][2];	
};

//General Polygon - composed of individual triangles
struct GeneralPolygon
{
	int no_triangles;
	Triangle* component_triangles;
};

//Triangle Tree - Used for refinement of a triangle
struct TriangleTree
{
	Triangle parent;
	int no_children = 0;
	TriangleTree* children;
};

//////////////////////////////////////////////////////////////
//Function Prototypes
//////////////////////////////////////////////////////////////
void DeallocateGeneralPolygon(GeneralPolygon polygon);
void PrintGeneralPolygon(GeneralPolygon polygon);
GeneralPolygon GenerateSquareMesh(int noRefinements);
GeneralPolygon GenerateLShapedMesh(int noRefinements);
GeneralPolygon GenerateHexagonalMesh(int noRefinements);
void RefineTree(TriangleTree& node,int no_refinements);
void DeleteTree(TriangleTree& node);
void WalkTree(TriangleTree& node,int& nodeNo,GeneralPolygon& polygon);


#endif
