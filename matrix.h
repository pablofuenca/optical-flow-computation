#ifndef MATRIX_H
#define MATRIX_H
 
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

class matrix{
public:
	matrix(void);
	matrix(int _rows, int _cols, int _depth, int _bound_x, int _bound_y, int _bound_z);
	~matrix(void);

	float operator()(int _x, int _y, int _z);
	void operator()(int _x, int _y, int _z, float _value);
	float operator()(int _x, int _y, int _z, int _derivat, int _derivative);

	void create3Dmatrix(int _rows, int _cols, int _depth, int _bound_x, int _bound_y, int _bound_z);
	void fillWithScalar(float _value);

	void createVectorMatrix(int rows, int cols, int depth, int boundaries);
	void loadVolumeScalarDataInMatrix(const char* fileName, int channel, int number_channels);
	void loadVolumeVectorDataInMatrix(const char* fileName);
	void loadVolumeScalarDataInMatrixPrueba(const char* fileName);
	void loadVolumeScalarDataInMatrixPrueba2(const char* fileName);


	void print_matrix();
	

	int rows, cols, depth;
	int bound_x, bound_y, bound_z;
	
	int first, second, third;

	float ***m_vecmatrixX;
	float ***m_vecmatrixY;
	float ***m_vecmatrixZ;

	int boundaries;
	int aux1;
	int aux2;
	int aux3;

	void freeMatrix3D();

private:
	
	//vector<vector<vector<double>>> * tarray3D;
	float *array3D;
};


#endif // MATRIX_H