/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef  MARCHINGCUBE_HPP
#define  MARCHINGCUBE_HPP

#include <cmath>
#include <string>
#include <vector>
#include <yade/yade-lib-wm3-math/Vector3.hpp>

using namespace std;

class MarchingCube 
{

/// ATTRIBUTES

	private : vector<Vector3f> triangles;
	public  : const vector<Vector3f>& getTriangles() { return triangles; }
	
	private : vector<Vector3f> normals;
	public  : const vector<Vector3f>& getNormals() { return normals; }

	private : int nbTriangles;
	public  : int getNbTriangles() { return nbTriangles; }
  
	private : Vector3f beta;
	private : Vector3f alpha;
	private : int sizeX,sizeY,sizeZ;
	private : float isoValue;
  
	private : vector<vector<vector<Vector3f> > > positions;
	private : static const int edgeArray[256];
	private : static const int triTable[256][16];
	Vector3f aNormal;	
	
/// PRIVATE METHOD
  
	/** triangulate cell  (x,y,z) **/
	private : void polygonize (const vector<vector<vector<float> > >& scalarField, int x, int y, int z);

	/** compute normals of the triangles previously found with polygonizecalcule les normales des triangles trouver dans la case (x,y,z)
		@param n : indice of the first triangle to process
	**/
  	private : void computeNormal(const vector<vector<vector<float> > >& scalarField, int x, int y, int z,int offset, int triangleNum);
	
	/** interpolate coordinates of point vect (that is on isosurface) from coordinates of points vect1 et vect2 **/
	private : void interpolate (const Vector3f&  vect1, const Vector3f& vect2, float val1, float val2,Vector3f& vect);

	/** Same as interpolate but in 1D **/
	private : float interpolateValue(float val1, float val2, float val_cible1, float val_cible2);

	/** Compute normal to vertice or triangle inside cell (x,y,z) **/
	private : const Vector3f& computeNormalX(const vector<vector<vector<float> > >& scalarField, int x, int y, int z);
	private : const Vector3f& computeNormalY(const vector<vector<vector<float> > >& scalarField, int x, int y, int z );
	private : const Vector3f& computeNormalZ(const vector<vector<vector<float> > >& scalarField, int x, int y, int z); 

/// CONSTRUCTOR/DESTRUCTOR

	public  : MarchingCube ();
	public  : ~MarchingCube ();

/// PULIC METHODS

	public  : void computeTriangulation(const vector<vector<vector<float> > >& scalarField, float iso);

	public  : void init(int sx, int sy, int sz, const Vector3f& min, const Vector3f& max);

	public  : void resizeScalarField(vector<vector<vector<float> > >& scalarField, int sx, int sy, int sz);
};

#endif // MARCHINGCUBE_HPP

