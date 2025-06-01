//////////////////////////////////////////////////////////////
//Module that implements a triangulation of various 
//polygons
//////////////////////////////////////////////////////////////

#include <cassert>
#include <cmath>

#include "linear_algebra.hpp"
#include "polygon.hpp"

//////////////////////////////////////////////////////////////
void DeallocateGeneralPolygon(GeneralPolygon polygon)
//Deallocates a previously allocated polygon
{
    delete[] polygon.component_triangles;
    polygon.no_triangles = 0;
}

//////////////////////////////////////////////////////////////
void PrintGeneralPolygon(GeneralPolygon polygon)
//Prints a polygonal domain to the screen
{
    std::cout << "No triangles = " << polygon.no_triangles << std::endl;

    for (int k=0;k<polygon.no_triangles;k++)
    {
        std::cout << "Triangle " << k << std::endl;

        for (int i=0;i<3;i++)
	    {
		    for (int j=0;j<2;j++)
		    {
			    std::cout << polygon.component_triangles[k].vertices[i][j] << " ";
		    }
		    std::cout << std::endl;
	    }
    }
}

//////////////////////////////////////////////////////////////
GeneralPolygon GenerateSquareMesh(int noRefinements)
//Creates a triangular decomposition of the square [0,1]^2
//0 refinements produces a triangulation comprising two right angled
//triangles. 1 refinement splits each triangle into 4 equally
//shaped triangles and so on.
//Inputs:
//  noRefinements - Number of refinements (>= 0) to carry out
//Ouput:
// GeneralPolygon
{
    assert(noRefinements >= 0);

    GeneralPolygon polygon;

    polygon.no_triangles = 2*pow(4,noRefinements);
    polygon.component_triangles = new Triangle[polygon.no_triangles];

    TriangleTree current_master_element;
    int current_triangle = 0;

    //First master triangle
    current_master_element.parent.vertices[0][0] = 0.0;
    current_master_element.parent.vertices[0][1] = 0.0;
    current_master_element.parent.vertices[1][0] = 1.0;
    current_master_element.parent.vertices[1][1] = 0.0;
    current_master_element.parent.vertices[2][0] = 0.0;
    current_master_element.parent.vertices[2][1] = 1.0;

    RefineTree(current_master_element,noRefinements);

    WalkTree(current_master_element,current_triangle,polygon);

    DeleteTree(current_master_element);

    //Second master triangle
    current_master_element.parent.vertices[0][0] = 1.0;
    current_master_element.parent.vertices[0][1] = 0.0;
    current_master_element.parent.vertices[1][0] = 1.0;
    current_master_element.parent.vertices[1][1] = 1.0;
    current_master_element.parent.vertices[2][0] = 0.0;
    current_master_element.parent.vertices[2][1] = 1.0;

    RefineTree(current_master_element,noRefinements);

    WalkTree(current_master_element,current_triangle,polygon);

    DeleteTree(current_master_element);


    return polygon;
}

/////////////////////////////////////////////////////////////
GeneralPolygon GenerateLShapedMesh(int noRefinements)
//Creates a triangular decomposition of the L-shaped domain
//[0,1]^2\[0.5,1]^2.
//0 refinements produces a triangulation comprising 6 right angled
//triangles. 1 refinement splits each triangle into 4 equally
//shaped triangles and so on.
//Inputs:
//  noRefinements - Number of refinements (>= 0) to carry out
//Ouput:
// GeneralPolygon
{
    assert(noRefinements >= 0);

    GeneralPolygon polygon;

    polygon.no_triangles = 6*pow(4,noRefinements);
    polygon.component_triangles = new Triangle[polygon.no_triangles];

    double coordinates[8][2] = {{0.0,0.0},{0.5,0.0},{1.0,0.0},
        {0.0,0.5},{0.5,0.5},{1.0,0.5},
        {0.0,1.0},{0.5,1.0}};

    int master_tris[6][3] = {{0,1,3},{1,4,3},{1,2,4},{2,5,4},
        {3,4,6},{4,7,6}};

    TriangleTree current_master_element;
    int current_triangle = 0;

    for (int k=0;k<6;k++)
    {
        current_master_element.parent.vertices[0][0] = coordinates[master_tris[k][0]][0];
        current_master_element.parent.vertices[0][1] = coordinates[master_tris[k][0]][1];
        current_master_element.parent.vertices[1][0] = coordinates[master_tris[k][1]][0];
        current_master_element.parent.vertices[1][1] = coordinates[master_tris[k][1]][1];
        current_master_element.parent.vertices[2][0] = coordinates[master_tris[k][2]][0];
        current_master_element.parent.vertices[2][1] = coordinates[master_tris[k][2]][1];

        RefineTree(current_master_element,noRefinements);

        WalkTree(current_master_element,current_triangle,polygon);

        DeleteTree(current_master_element);
    }

    return polygon;
}

//////////////////////////////////////////////////////////////
GeneralPolygon GenerateHexagonalMesh(int noRefinements)
//Creates a triangular decomposition of the regular hexagon
//contained within the unit circle.
//Inputs:
//  noRefinements - Number of refinements (>= 0) to carry out
//Ouput:
// GeneralPolygon
{
    assert(noRefinements >= 0);

    GeneralPolygon polygon;

    polygon.no_triangles = 6*pow(4,noRefinements);
    polygon.component_triangles = new Triangle[polygon.no_triangles];

    DoubleMatrix hexagon_points = AllocateDoubleMatrix(7,2);
    for (int k=0;k<7;k++)
    {
        hexagon_points.matrix_entries[k][0] = cos((double)k*M_PI/3.0);
        hexagon_points.matrix_entries[k][1] = sin((double)k*M_PI/3.0);
    }

    TriangleTree current_sector;
    int current_triangle = 0;
    //Loop over each sector
    for (int k=0;k<6;k++)
    {
       
        current_sector.parent.vertices[0][0] = 0.0;
        current_sector.parent.vertices[0][1] = 0.0;
        current_sector.parent.vertices[1][0] = hexagon_points.matrix_entries[k][0];
        current_sector.parent.vertices[1][1] = hexagon_points.matrix_entries[k][1];
        current_sector.parent.vertices[2][0] = hexagon_points.matrix_entries[k+1][0];
        current_sector.parent.vertices[2][1] = hexagon_points.matrix_entries[k+1][1];

        RefineTree(current_sector,noRefinements);

        WalkTree(current_sector,current_triangle,polygon);

        DeleteTree(current_sector);

    }

    DeallocateMatrix(hexagon_points);

    return polygon;
}

//////////////////////////////////////////////////////////////
void RefineTree(TriangleTree& node,int no_refinements)
//Function to recursively refine a quad tree data structure and
//subvide triangles into 4 equally shaped children
{
    if (no_refinements > 0)
    {
        //Set up refined elements
        node.no_children = 4;
        node.children = new TriangleTree[4];

        //Child 0
        node.children[0].parent.vertices[0][0] = node.parent.vertices[0][0];
        node.children[0].parent.vertices[0][1] = node.parent.vertices[0][1];
        node.children[0].parent.vertices[1][0] = (node.parent.vertices[0][0]+node.parent.vertices[1][0])/2.0;
        node.children[0].parent.vertices[1][1] = (node.parent.vertices[0][1]+node.parent.vertices[1][1])/2.0;
        node.children[0].parent.vertices[2][0] = (node.parent.vertices[0][0]+node.parent.vertices[2][0])/2.0;
        node.children[0].parent.vertices[2][1] = (node.parent.vertices[0][1]+node.parent.vertices[2][1])/2.0;

        //Child 1
        node.children[1].parent.vertices[0][0] = (node.parent.vertices[0][0]+node.parent.vertices[1][0])/2.0;
        node.children[1].parent.vertices[0][1] = (node.parent.vertices[0][1]+node.parent.vertices[1][1])/2.0;
        node.children[1].parent.vertices[1][0] = (node.parent.vertices[2][0]+node.parent.vertices[1][0])/2.0;
        node.children[1].parent.vertices[1][1] = (node.parent.vertices[2][1]+node.parent.vertices[1][1])/2.0;
        node.children[1].parent.vertices[2][0] = (node.parent.vertices[0][0]+node.parent.vertices[2][0])/2.0;
        node.children[1].parent.vertices[2][1] = (node.parent.vertices[0][1]+node.parent.vertices[2][1])/2.0;

        //Child 2
        node.children[2].parent.vertices[0][0] = (node.parent.vertices[0][0]+node.parent.vertices[1][0])/2.0;
        node.children[2].parent.vertices[0][1] = (node.parent.vertices[0][1]+node.parent.vertices[1][1])/2.0;
        node.children[2].parent.vertices[1][0] = node.parent.vertices[1][0];
        node.children[2].parent.vertices[1][1] = node.parent.vertices[1][1];
        node.children[2].parent.vertices[2][0] = (node.parent.vertices[2][0]+node.parent.vertices[1][0])/2.0;
        node.children[2].parent.vertices[2][1] = (node.parent.vertices[2][1]+node.parent.vertices[1][1])/2.0;

        //Child 3
        node.children[3].parent.vertices[0][0] = (node.parent.vertices[0][0]+node.parent.vertices[2][0])/2.0;
        node.children[3].parent.vertices[0][1] = (node.parent.vertices[0][1]+node.parent.vertices[2][1])/2.0;
        node.children[3].parent.vertices[1][0] = (node.parent.vertices[2][0]+node.parent.vertices[1][0])/2.0;
        node.children[3].parent.vertices[1][1] = (node.parent.vertices[2][1]+node.parent.vertices[1][1])/2.0;
        node.children[3].parent.vertices[2][0] = node.parent.vertices[2][0];
        node.children[3].parent.vertices[2][1] = node.parent.vertices[2][1];

        for (int j=0;j<4;j++)
        {
            RefineTree(node.children[j],no_refinements-1);
        }
    }
}

//////////////////////////////////////////////////////////////
void DeleteTree(TriangleTree& node)
//Function to delete a quad tree data structure recursively
{
    if (node.no_children > 0)
    {
        for (int k=0;k<4;k++)
        {
            DeleteTree(node.children[k]);
        }
        node.no_children = 0;
        delete[] node.children;
    }
}

//////////////////////////////////////////////////////////////
void WalkTree(TriangleTree& node,int& nodeNo,GeneralPolygon& polygon)
//Function to walk a quad tree and create a triangulation of 
//a polygon
{
    if (node.no_children > 0)
    {
        for (int k=0;k<4;k++)
        {
            WalkTree(node.children[k],nodeNo,polygon);
        }
    }
    else
    {
        //Add triangle to polygon
        polygon.component_triangles[nodeNo].vertices[0][0] = node.parent.vertices[0][0];
        polygon.component_triangles[nodeNo].vertices[0][1] = node.parent.vertices[0][1];
        polygon.component_triangles[nodeNo].vertices[1][0] = node.parent.vertices[1][0];
        polygon.component_triangles[nodeNo].vertices[1][1] = node.parent.vertices[1][1];
        polygon.component_triangles[nodeNo].vertices[2][0] = node.parent.vertices[2][0];
        polygon.component_triangles[nodeNo].vertices[2][1] = node.parent.vertices[2][1];

        nodeNo++;
    }
}