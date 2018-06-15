#define _USE_MATH_DEFINES
#include <cmath>

#include "matrix.h"
#include "data.h"
#include "VolumeDisplayer.h"
#include "opticalFlow.h"
#include "stdio.h"
#include <math.h>
#include <fstream>
#include <iostream>



using namespace std;

void save3DVectorField(string _url, matrix& _x, matrix& _y, matrix& _z){
	fstream file(_url, fstream::out | fstream::binary);

	if (!file.is_open())
		return;

	float t_x, t_y, t_z;

	for (int z = _z.bound_z; z<(_z.depth + _z.bound_z); z++)
	{
		for (int y = _y.bound_y; y<(_y.rows + _y.bound_y); y++)
		{
			for (int x = _x.bound_x; x<(_x.cols + _x.bound_x); x++){
				t_x = (float)_x(x, y, z);
				t_y = (float)_y(x, y, z);
				t_z = (float)_z(x, y, z);
				file.write((char*)&t_x, sizeof(float));
				file.write((char*)&t_y, sizeof(float));
				file.write((char*)&t_z, sizeof(float));
			}
		}
	}

	file.close();
}


void save3DGroundTruth(string _url, opticalFlow& OptFlowGround, int cols_x, int rows_y, int depth_z){
	fstream file(_url, fstream::out | fstream::binary);

	if (!file.is_open())
		return;

	float t_x, t_y, t_z;
	int index;

	for (int z = 0; z<depth_z; z++)
	{
		for (int y = 0; y<rows_y; y++)
		{
			for (int x = 0; x<cols_x; x++){
				index = x + (y*cols_x) + (z*rows_y*cols_x);
				t_x = (float)OptFlowGround.ground[(index * 3) + 0];
				t_y = (float)OptFlowGround.ground[(index * 3) + 1];
				t_z = (float)OptFlowGround.ground[(index * 3) + 2];
				file.write((char*)&t_x, sizeof(float));
				file.write((char*)&t_y, sizeof(float));
				file.write((char*)&t_z, sizeof(float));
			}
		}
	}

	file.close();
}

void load3DVectorField(string _url, float** _target, data& _data){
	fstream file(_url, fstream::in | fstream::binary);

	int t_size = _data.x_dim * _data.y_dim * _data.z_dim * 3;

	(*_target) = new float[t_size];

	if (!file.is_open())
		return;

	float t_x, t_y, t_z;

	float max_x = 0;
	float max_y = 0;
	float max_z = 0;

	float min_x = 0;
	float min_y = 0;
	float min_z = 0;



	int index;
	for (int z = 0; z < _data.z_dim; z++){ //  each pixel
		for (int y = 0; y < _data.y_dim; y++){
			for (int x = 0; x < _data.x_dim; x++){

				file.read((char*)&t_x, sizeof(float));
				file.read((char*)&t_y, sizeof(float));
				file.read((char*)&t_z, sizeof(float));

				if (t_x < min_x)	min_x = t_x;
				if (t_x > max_x)	max_x = t_x;

				if (t_y < min_y)	min_y = t_y;
				if (t_y > min_y)	max_y = t_y;

				if (t_z < min_z)	min_z = t_z;
				if (t_z > min_z)	max_z = t_z;

				index = x + (y*_data.x_dim) + (z*_data.y_dim*_data.x_dim);
				(*_target)[(index * 3) + 0] = (float)t_x;
				(*_target)[(index * 3) + 1] = (float)t_y;
				(*_target)[(index * 3) + 2] = (float)t_z;
			}
		}
	}

	file.close();
}

void load3DVectorFieldInMatrices(string _url, matrix& _x, matrix& _y, matrix& _z){
	fstream file(_url, fstream::in | fstream::binary);

	if (!file.is_open())
		return;

	float t_x, t_y, t_z;

	for (int z = _z.bound_z; z<(_z.depth + _z.bound_z); z++)
	{
		for (int y = _y.bound_y; y<(_y.rows + _y.bound_y); y++)
		{
			for (int x = _x.bound_x; x<(_x.cols + _x.bound_x); x++)
			{
				file.read((char*)&t_x, sizeof(float));
				file.read((char*)&t_y, sizeof(float));
				file.read((char*)&t_z, sizeof(float));

				_x(x, y, z, t_x*10);
				_y(x, y, z, t_y*10);
				_z(x, y, z, t_z*10);
			}
		}
	}

	file.close();
}

void test(){

	int n_rows = 256;
	int n_cols = 256;
	int n_depth = 256;

	/*data* tdata = new data();

	string img1 = "skull.raw";
	string img2 = "skull_5.raw";


	tdata->loadVolumeScalarData(img1.c_str(), n_rows, n_cols, n_depth, 1);
	tdata->loadVolumeScalarData(img2.c_str(), n_rows, n_cols, n_depth, 2);

	delete tdata;*/

	matrix Ix(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix Iy(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix Iz(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix It(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix m(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix m2(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix i1(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix i2(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix u(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix v(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix w(n_rows, n_cols, n_depth, 0, 0, 0);
}


//void onlyCalc(string _name, string _img1, string _img2, int n_iter, int n_rows, int n_cols, int n_depth, int constancy){
//	data data;
//	opticalFlow OF;
//
//	data.x_dim = n_cols;
//	data.y_dim = n_rows;
//	data.z_dim = n_depth;
//
//	matrix m(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
//	matrix m2(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
//
//	m.loadVolumeScalarDataInMatrix(_img1.c_str());
//	m2.loadVolumeScalarDataInMatrix(_img2.c_str());
//
//	matrix Ix(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
//	matrix Iy(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
//	matrix Iz(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
//	matrix It(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
//
//	matrix u(n_rows, n_cols, n_depth, 0, 0, 0);
//	matrix v(n_rows, n_cols, n_depth, 0, 0, 0);
//	matrix w(n_rows, n_cols, n_depth, 0, 0, 0);
//
//	OF.computeImputDerivatives(data, Ix, Iy, Iz, It, m, m2, constancy);
//	OF.computeMotionTensor(data, n_iter, Ix, Iy, Iz, It, u, v, w);
//
//	m2.freeMatrix3D();
//	m.freeMatrix3D();
//
//	It.freeMatrix3D();
//	Iz.freeMatrix3D();
//	Iy.freeMatrix3D();
//	Ix.freeMatrix3D();
//
//	data.deallocate();
//
//	save3DVectorField(_name, u, v, w);
//}

//void doVisualization(string _name, int _dx, int _dy, int _dz){
//	data data;
//	opticalFlow OF;
//
//	data.x_dim = _dx;
//	data.y_dim = _dy;
//	data.z_dim = _dz;
//
//	load3DVectorField(_name, &(OF.OptFlow), data);
//
//	VolumeDisplayer VolumeDisplayer;
//	VolumeDisplayer.visualizeVectorVolume(data, OF);
//}


float  NormCoord(float _dim_coord, int _dim){
	return _dim_coord / ((float)(_dim - 1));
}

float DimCoord(float _norm_coord, int _dim){
	return _norm_coord*((float)(_dim - 1));
}


vector<float> NormCoords(float _x, float _y, float _z, int colx, int rowy, int depz){
	vector<float> ret;

	ret.push_back(NormCoord(_x, colx));
	ret.push_back(NormCoord(_y, rowy));
	ret.push_back(NormCoord(_z, depz));
	return ret;
}

vector<float> DimCoords(float _norm_x, float _norm_y, float _norm_z, int cols_x, int rows_y, int depth_z){
	vector<float> ret;

	ret.push_back(DimCoord(_norm_x, cols_x));
	ret.push_back(DimCoord(_norm_y, rows_y));
	ret.push_back(DimCoord(_norm_z, depth_z));

	return ret;
}


//linear interpolation between _fist and _second with alpha = _a
float linearInterpolation(float _first, float _second, float _a){
	return _first*(1.0f - _a) + _second*_a;
}

//bilinear interpolation in a quad
float bilinearInterpolation(float _tl, float _tr, float _bl, float _br, float _a, float _b){

	float top = linearInterpolation(_tl, _tr, _a);
	float bottom = linearInterpolation(_bl, _br, _a);

	return linearInterpolation(top, bottom, _b);
}

//trilinear interpolation in a quad
float trilinearInterpolation(matrix& _src, float _x, float _y, float _z){

	if ((_x >= (float)_src.cols) || (_x<0.0f))
		return 0.0f;

	if ((_y >= (float)_src.rows) || (_y<0.0f))
		return 0.0f;

	if ((_z >= (float)_src.depth) || (_z<0.0f))
		return 0.0f;


	int x = (int)floor(_x);
	int y = (int)floor(_y);
	int z = (int)floor(_z);

	float a = _x - (float)x;
	float b = _y - (float)y;
	float c = _z - (float)z;


	float front = bilinearInterpolation(_src(x, y, z), _src(x + 1, y, z), _src(x, y + 1, z), _src(x + 1, y + 1, z), a, b);

	float back = bilinearInterpolation(_src(x, y, z + 1), _src(x + 1, y, z + 1), _src(x, y + 1, z + 1), _src(x + 1, y + 1, z + 1), a, b);

	return  linearInterpolation(front, back, c);
}



//vector field integrator
vector<float> RungeKutta4(matrix& _vx, matrix& _vy, matrix& _vz, float _delta_t, float _norm_x, float _norm_y, float _norm_z, int cols_x, int rows_y, int depth_z){

	float h1 = 1.0f / 6.0f;
	float h2 = 1.0f / 3.0f;

	float t_value;
	vector<float> k1, k2, k3, k4;


	vector<float> dimCoords = DimCoords(_norm_x, _norm_y, _norm_z, cols_x, rows_y, depth_z);

	//-----------------------------------------------------------------------------------------------

	t_value = trilinearInterpolation(_vx, dimCoords[0], dimCoords[1], dimCoords[2]);
	k1.push_back(t_value*_delta_t);
	t_value = trilinearInterpolation(_vy, dimCoords[0], dimCoords[1], dimCoords[2]);
	k1.push_back(t_value*_delta_t);
	t_value = trilinearInterpolation(_vz, dimCoords[0], dimCoords[1], dimCoords[2]);
	k1.push_back(t_value*_delta_t);

	//-----------------------------------------------------------------------------------------------

	t_value = trilinearInterpolation(_vx, dimCoords[0] + k1[0] * 0.5f, dimCoords[1] + k1[1] * 0.5f, dimCoords[2] + k1[2] * 0.5f);
	k2.push_back(t_value*_delta_t);
	t_value = trilinearInterpolation(_vy, dimCoords[0] + k1[0] * 0.5f, dimCoords[1] + k1[1] * 0.5f, dimCoords[2] + k1[2] * 0.5f);
	k2.push_back(t_value*_delta_t);
	t_value = trilinearInterpolation(_vz, dimCoords[0] + k1[0] * 0.5f, dimCoords[1] + k1[1] * 0.5f, dimCoords[2] + k1[2] * 0.5f);
	k2.push_back(t_value*_delta_t);

	//-----------------------------------------------------------------------------------------------

	t_value = trilinearInterpolation(_vx, dimCoords[0] + k2[0] * 0.5f, dimCoords[1] + k2[1] * 0.5f, dimCoords[2] + k2[2] * 0.5f);
	k3.push_back(t_value*_delta_t);
	t_value = trilinearInterpolation(_vy, dimCoords[0] + k2[0] * 0.5f, dimCoords[1] + k2[1] * 0.5f, dimCoords[2] + k2[2] * 0.5f);
	k3.push_back(t_value*_delta_t);
	t_value = trilinearInterpolation(_vz, dimCoords[0] + k2[0] * 0.5f, dimCoords[1] + k2[1] * 0.5f, dimCoords[2] + k2[2] * 0.5f);
	k3.push_back(t_value*_delta_t);

	//-----------------------------------------------------------------------------------------------

	t_value = trilinearInterpolation(_vx, dimCoords[0] + k3[0], dimCoords[1] + k3[1], dimCoords[2] + k3[2]);
	k4.push_back(t_value*_delta_t);
	t_value = trilinearInterpolation(_vy, dimCoords[0] + k3[0], dimCoords[1] + k3[1], dimCoords[2] + k3[2]);
	k4.push_back(t_value*_delta_t);
	t_value = trilinearInterpolation(_vz, dimCoords[0] + k3[0], dimCoords[1] + k3[1], dimCoords[2] + k3[2]);
	k4.push_back(t_value*_delta_t);

	k1[0] *= h1;
	k1[1] *= h1;
	k1[2] *= h1;

	k2[0] *= h2;
	k2[1] *= h2;
	k2[2] *= h2;

	k3[0] *= h2;
	k3[1] *= h2;
	k3[2] *= h2;

	k4[0] *= h1;
	k4[1] *= h1;
	k4[2] *= h1;

	vector<float> ret;

	ret.push_back(k1[0] + k2[0] + k3[0] + k4[0]);
	ret.push_back(k1[1] + k2[1] + k3[1] + k4[1]);
	ret.push_back(k1[2] + k2[2] + k3[2] + k4[2]);

	return ret;
}


void interpolation(matrix& _x, matrix& _y, matrix& _z, int cols_x, int rows_y, int depth_z, float _stepsize, string _salida){

	opticalFlow OptFlowGround;

	//vector which saves the interpolated values
	vector<float> dir;
	//vector with normalised positions
	vector<float> position_norm;


	for (int z = 0; z < depth_z; z++){
		for (int y = 0; y < rows_y; y++){
			for (int x = 0; x < cols_x; x++){

				//normaised positions
				position_norm = NormCoords(x, y, z, cols_x, rows_y, depth_z);

				//RungeKutta method to obtain the interpolated values
				dir = RungeKutta4(_x, _y, _z, -_stepsize, position_norm[0], position_norm[1], position_norm[2], cols_x, rows_y, depth_z);

				//save this data in one vector
				OptFlowGround.saveDataGround(dir[0], dir[1], dir[2], cols_x, rows_y, depth_z, x, y, z);

			}
		}
	}

	//free memory
	_x.freeMatrix3D();
	_y.freeMatrix3D();
	_z.freeMatrix3D();

	//save the interpolated grund truth data in an output file
	save3DGroundTruth(_salida, OptFlowGround, cols_x, rows_y, depth_z);
	OptFlowGround.deallocateGround();

}

float vectorLength(vector<float>& _v){
	return sqrt(_v[0] * _v[0] + _v[1] * _v[1] + _v[2] * _v[2]);
}

float calcAngleBetweenVectors(vector<float>& i_v, vector<float> &r_v, int x, int y, int z){

	float ret = 0.0;

	float top = i_v[0] * r_v[0] + i_v[1] * r_v[1] + i_v[2] * r_v[2];

	float lv1 = vectorLength(i_v);
	float lv2 = vectorLength(r_v);

	float down = lv1*lv2;

	if (down == 0.0){
		return 0.0;
	}

	if (top >= down){
		return 0.0;
	}
	else{
		ret = acos(top/down);
	}

	return ret;
}

vector<float> getInterpolatedVector(matrix & _c_x, matrix & _c_y, matrix & _c_z, float _rel_x, float _rel_y, float _rel_z){

	vector<float> ret;

	float n_x = (_c_x.cols - 1)*_rel_x;
	float n_y = (_c_x.rows - 1)*_rel_y;
	float n_z = (_c_x.depth - 1)*_rel_z;

	float t_value;

	t_value = trilinearInterpolation(_c_x, n_x, n_y, n_z);
	ret.push_back(t_value);

	t_value = trilinearInterpolation(_c_y, n_x, n_y, n_z);
	ret.push_back(t_value);

	t_value = trilinearInterpolation(_c_z, n_x, n_y, n_z);
	ret.push_back(t_value);

	return ret;
}

float compareVectorLength(vector<float>& _gt_v, vector<float>& _v){
	float gt_length = vectorLength(_gt_v);
	float v_length = vectorLength(_v);

	return abs(v_length - gt_length);
}

void doError(string _outfile_length, string _outfile_angle, string _groundTruth, int _g_cols_x, int _g_rows_y, int _g_depth_z, float _gt_vector_prolong, matrix & _r_x, matrix & _r_y, matrix & _r_z){

	cout << "\nCalculating Error";

	fstream file_angle(_outfile_angle, fstream::out | fstream::binary);
	fstream file_length(_outfile_length, fstream::out | fstream::binary);

	int depth_z = _r_x.depth - 1;
	int rows_y = _r_x.rows - 1;
	int cols_x = _r_x.cols - 1;

	//create matrices to save the ground truth data
	matrix g_x(_g_cols_x, _g_rows_y, _g_depth_z, 0, 0, 0);
	matrix g_y(_g_cols_x, _g_rows_y, _g_depth_z, 0, 0, 0);
	matrix g_z(_g_cols_x, _g_rows_y, _g_depth_z, 0, 0, 0);

	//load the ground truth data in the created matrices
	load3DVectorFieldInMatrices(_groundTruth, g_x, g_y, g_z);

	double A = 0.0;
	double angle_degree = 0;

	double avg_aae = 0;

	vector<float> interpolated_vector;
	vector<float> result_vector = vector<float>(3, 0.0);

	float compared_length = 0;

	int counter = 0;

	for (int z = 0; z <= depth_z; z++){
		for (int y = 0; y <= rows_y; y++){
			for (int x = 0; x <= cols_x; x++){

				result_vector[0] = _r_x(x, y, z);
				result_vector[1] = _r_y(x, y, z);
				result_vector[2] = _r_z(x, y, z);

				interpolated_vector = getInterpolatedVector(g_x, g_y, g_z, (float)x / (float)cols_x, (float)y / (float)rows_y, (float)z / (float)depth_z);
				interpolated_vector[0] *= _gt_vector_prolong;
				interpolated_vector[1] *= _gt_vector_prolong;
				interpolated_vector[2] *= _gt_vector_prolong;

				A = calcAngleBetweenVectors(interpolated_vector, result_vector, x, y, z);
				
				/*if (isnan(A)){
					cout << "\nInd en (x,y,z): " << x << ", " << y << ", " << z << endl;
					system("pause");
				}*/

				angle_degree = double((180.0 / M_PI)*A);

				avg_aae += angle_degree;
				counter++;

				file_angle.write((char*)&angle_degree, sizeof(float));

				compared_length = compareVectorLength(interpolated_vector, result_vector);
				file_length.write((char*)&compared_length, sizeof(float));
			}
		}
	}

	file_angle.close();
	file_length.close();

	cout << "\ncounter= " << counter << "\navg_aae= " << avg_aae << "\n";
	cout << "\n\nAAE in degree: " << avg_aae / (float)(counter);

}


void doOF(string _name, string _img1, string _img2, string _gt, int _g_cols_x, int _g_rows_y, int _g_depth_z, 
	float _gt_vector_prolong, int n_iter, int n_rows, int n_cols, int n_depth, int _n_channels){
	
	data data;
	opticalFlow OF;

	data.x_dim = n_cols;
	data.y_dim = n_rows;
	data.z_dim = n_depth;

	cout << "\n\nCalculating motion tensor";

	matrix J00(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix J01(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix J02(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix J03(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix J11(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix J12(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix J13(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix J22(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix J23(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);

	//allocate memory for the matrices which store the data(one per channel of each image)
	matrix m(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix m2(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);

	//allocate momory for the image derivatives
	matrix Ix(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix Iy(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix Iz(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);
	matrix It(data.x_dim, data.y_dim, data.z_dim, 0, 0, 0);

	//implementamos el solver
	for (int n = 0; n < _n_channels; n++){
		
		//store values of the image in the matrices created before with the number of the channel that we are in
		m.loadVolumeScalarDataInMatrix(_img1.c_str(), n, _n_channels);
		m2.loadVolumeScalarDataInMatrix(_img2.c_str(), n, _n_channels);

		//create the derivatives, the motion tensor and compute the solver
		OF.computeImputDerivatives(data, Ix, Iy, Iz, It, m, m2);
		//allocate memory for the optical flow matrices
		OF.computeMotionTensor(data, Ix, Iy, Iz, It, J00, J01, J02, J03, J11, J12, J13, J22, J23);

	}

	//free memory of the matrices with the data and with the image derivatives
	m2.freeMatrix3D();
	m.freeMatrix3D();

	It.freeMatrix3D();
	Iz.freeMatrix3D();
	Iy.freeMatrix3D();
	Ix.freeMatrix3D();

	cout << "\n\nCalculating Optical Flow";

	matrix u(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix v(n_rows, n_cols, n_depth, 0, 0, 0);
	matrix w(n_rows, n_cols, n_depth, 0, 0, 0);

	//compute the solver
	OF.computeSolver(data, n_iter, u, v, w, J00, J01, J02, J03, J11, J12, J13, J22, J23);

	//free memory of the motion tensor matrices
	J00.freeMatrix3D();
	J01.freeMatrix3D();
	J02.freeMatrix3D();
	J03.freeMatrix3D();
	J11.freeMatrix3D();
	J12.freeMatrix3D();
	J13.freeMatrix3D();
	J22.freeMatrix3D();
	J23.freeMatrix3D();

	//save the optical flow in an optput file
	save3DVectorField(_name, u, v, w);

	//////open the file for reading and type binary
	//fstream file(_name, fstream::in | fstream::binary);

	//////check if the file is correctly opened
	////if (!file.is_open())
	////	return;

	//////auxiliar variable
	//unsigned char value;
	////int index;

	//////fill matrix whit the data
	//for (int z = 0; z < 256; z++){
	//	for (int y = 0; y < 256; y++){
	//		for (int x = 0; x < 256; x++){

	//			file.read((char*)&value, sizeof(unsigned char));

	//			cout << (float)value << " ";
	//		}
	//		system("pause");

	//	}
	//}
	//system("pause");

	//save the three components of thevector field 
	//in one single string
	OF.saveData(u, v, w, data);

	string length_name = "length_";
	string angle_name = "angle_";

	length_name.append(_name);
	angle_name.append(_name);

	if (!_gt.empty())
		doError(length_name, angle_name, _gt, _g_cols_x, _g_rows_y, _g_depth_z, _gt_vector_prolong, u, v, w);

	//free memory once we have saved the optical flow data
	u.freeMatrix3D();
	v.freeMatrix3D();
	w.freeMatrix3D();

	//represent the data
	VolumeDisplayer vol;
	vol.visualizeVectorVolume(data, OF);
	OF.deallocateOptFlow();
	data.deallocate();
}

void writeFile(string filename, matrix& _x, matrix& _y, matrix& _z, int _w, int _h, int _d){
	fstream file(filename, fstream::out | fstream::binary);

	if (!file.is_open())
		return;

	float t_x, t_y, t_z;

	for (int z = 0; z < _d; z++){
		for (int y = 0; y < _h; y++){
			for (int x = 0; x < _w; x++){

				t_x = (float)_x(x, y, z);
				t_y = (float)_y(x, y, z);
				t_z = (float)_z(x, y, z);

				file.write((char*)&t_x, sizeof(float));
				file.write((char*)&t_y, sizeof(float));
				file.write((char*)&t_z, sizeof(float));
				//cout << " " << t_x << " " << t_y << " " << t_z;
			}
		}
	}

	file.close();


}

void generate3DVec2TestImage(){
	int w = 81; //width
	int h = 81; //height
	int d = 81; //depth

	int tx = 3; //translation of each pixel in x direction
	int ty = 3; //translation of each pixel in y direction
	int tz = 3; //translation of each pixel in z direction

	matrix xx(w, h, d, 0, 0, 0); // channel 1 of first image
	matrix yy(w, h, d, 0, 0, 0); // channel 2 of first image
	matrix zz(w, h, d, 0, 0, 0); // channel 3 of first image

	matrix t_xx(w, h, d, 0, 0, 0); // channel 1 of second image
	matrix t_yy(w, h, d, 0, 0, 0); // channel 2 of second image
	matrix t_zz(w, h, d, 0, 0, 0); // channel 3 of second image

	//-------------------------------
	// generate unique vector according to the position
	for (int z = 0; z < d; z++){
		for (int y = 0; y < h; y++){
			for (int x = 0; x < w; x++){

				xx(x, y, z, x);
				yy(x, y, z, y);
				zz(x, y, z, z);

			}
		}
	}

	//shift first image and save it in second image
	for (int z = 0; z < d - tz; z++){
		for (int y = 0; y < h - ty; y++){
			for (int x = 0; x < w - tx; x++){

				t_xx(x + tx, y + ty, z + tz, xx(x, y, z));
				t_yy(x + tx, y + ty, z + tz, yy(x, y, z));
				t_zz(x + tx, y + ty, z + tz, zz(x, y, z));
				
			}
		}
	}

	// save images
	writeFile("0.dat", xx, yy, zz, w, h, d);
	writeFile("1.dat", t_xx, t_yy, t_zz, w, h, d);

	//-------------------------------
	xx.freeMatrix3D();
	yy.freeMatrix3D();
	zz.freeMatrix3D();
	t_xx.freeMatrix3D();
	t_yy.freeMatrix3D();
	t_zz.freeMatrix3D();
	
}

int main(){

	//generate3DVec2TestImage();

	int n_iter = 100;
	int n_channel;

	string img1 = "skull.raw";
	string img2 = "skull_5.raw";
	float gt_prolonged = 5.0;

	int n_cols_x = 256;
	int n_rows_y = 256;
	int n_depth_z = 256;

	

	string outfilename = "sk5.raw";
	string groundTruth = "27.024731.raw";

	int g_cols_x = 61;
	int g_rows_y = 31;
	int g_depth_z = 61;

	cout << "\nInsert number of channels: ";
	cin >> n_channel;

	doOF(outfilename, img1, img2, "", g_cols_x, g_rows_y, g_depth_z, gt_prolonged, n_iter, n_rows_y, n_cols_x, n_depth_z, n_channel);

	system("pause");


	return 0;
}