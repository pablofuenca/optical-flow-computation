
#ifndef OPTICALFLOW_H
#define OPTICALFLOW_H

#include "matrix.h"
#include "data.h"

using namespace std;

class opticalFlow
{
public:
	opticalFlow();
	void computeImputDerivatives(data& data, matrix& _Ix, matrix& _Iy, matrix& _Iz, matrix& _It,  matrix& _m, matrix& _m2);
	void computeMotionTensor(data& data, matrix& _Ix, matrix& _Iy, matrix& _Iz, matrix& _It, 
		matrix& _J00, matrix& _J01, matrix& _J02, matrix& _J03,
		matrix& _J11, matrix& _J12, matrix& _J13, matrix& _J22, matrix& _J23);
	void computeSolver(data& data, int _n_iter, matrix& _u, matrix& _v, matrix& _w, matrix& _J00, matrix& _J01, matrix& _J02, matrix& _J03,
		matrix& _J11, matrix& _J12, matrix& _J13, matrix& _J22, matrix& _J23);


	void saveData(matrix &u, matrix& v, matrix& w, data& data);
	void opticalFlow::saveDataGround(float x_dat, float y_dat, float z_dat, int cols_x, int rows_y, int depth_z, int x_pos, int y_pos, int z_pos);
	void deallocateOptFlow();
	void deallocateGround();
	~opticalFlow();

	friend class VolumeDisplayer;

	int index, rows, cols, depth;
	float divx, divy, divz;
	float aux;
	float *OptFlow;
	float *ground;
	float Ixx, Iyy, Izz, Ixy, Ixz, Iyz, Ixt, Iyt, Izt;
	float temp00, temp01, temp02, temp03, temp11, temp12, temp13, temp22, temp23;
	

private:
	
	float J[4][4], Jgrad[4][4];
};

#endif


