#include "opticalFlow.h"
#include "matrix.h"
#include <cmath>


opticalFlow::opticalFlow()
{
	this->OptFlow = NULL;
}



//we compute the input derivatives with the formula we develope
void opticalFlow::computeImputDerivatives(data& data, matrix& _Ix, matrix& _Iy, matrix& _Iz, matrix& _It, matrix& _m, matrix& _m2){

	//image derivatives
	for (int z = 0; z < data.z_dim; z++){
		for (int y = 0; y < data.y_dim; y++){
			for (int x = 0; x < data.x_dim; x++){
				_Iz(x, y, z, (_m2(x, y, z + 1, 0, 0) - _m2(x, y, z - 1, 0, 0) + _m(x, y, z + 1, 0, 0) - _m(x, y, z - 1, 0, 0)) / 4);
				_Iy(x, y, z, (_m2(x, y + 1, z, 0, 0) - _m2(x, y - 1, z, 0, 0) + _m(x, y + 1, z, 0, 0) - _m(x, y - 1, z, 0, 0)) / 4);
				_Ix(x, y, z, (_m2(x + 1, y, z, 0, 0) - _m2(x - 1, y, z, 0, 0) + _m(x + 1, y, z, 0, 0) - _m(x - 1, y, z, 0, 0)) / 4);
				_It(x, y, z, _m2(x, y, z, 0, 0) - _m(x, y, z, 0, 0));
				//cout << "\nm: " << _m(x, y, z) << endl;
			}
		}
	}
}



void opticalFlow::computeMotionTensor(data& data, matrix& _Ix, matrix& _Iy, matrix& _Iz, matrix& _It, 
	matrix& _J00, matrix& _J01, matrix& _J02, matrix& _J03,
	matrix& _J11, matrix& _J12, matrix& _J13, matrix& _J22, matrix& _J23){

	float gamma1 = 0.5;
	float gamma2 = 2.0;

	/**************************************round over each pixel****************************************/
		for (int z = 0; z < data.z_dim; z++){ //  each pixel
			for (int y = 0; y < data.y_dim; y++){
				for (int x = 0; x < data.x_dim; x++){

					Ixx = (_Ix(x + 1, y, z, 0, 0) - _Ix(x - 1, y, z, 0, 0)) / 2.0;
					Iyy = (_Iy(x, y + 1, z, 0, 0) - _Iy(x, y - 1, z, 0, 0)) / 2.0;
					Izz = (_Iz(x, y, z + 1, 0, 0) - _Iz(x, y, z - 1, 0, 0)) / 2.0;

					Ixy = (_Ix(x, y + 1, z, 0, 0) - _Ix(x, y - 1, z, 0, 0)) / 2.0;
					Ixz = (_Ix(x, y, z + 1, 0, 0) - _Ix(x, y, z - 1, 0, 0)) / 2.0;

					Iyz = (_Iy(x, y, z + 1, 0, 0) - _Iy(x, y, z - 1, 0, 0)) / 2.0;

					Ixt = (_It(x + 1, y, z, 0, 0) - _It(x - 1, y, z, 0, 0)) / 2.0;
					Iyt = (_It(x, y + 1, z, 0, 0) - _It(x, y - 1, z, 0, 0)) / 2.0;
					Izt = (_It(x, y, z + 1, 0, 0) - _It(x, y, z - 1, 0, 0)) / 2.0;

					/*********************************BRIGHTNESS CONSTANCY***********************************/
					J[0][0] = _Ix(x, y, z)*_Ix(x, y, z)*gamma1;
					J[0][1] = _Ix(x, y, z)*_Iy(x, y, z)*gamma1;
					J[0][2] = _Ix(x, y, z)*_Iz(x, y, z)*gamma1;
					J[0][3] = _Ix(x, y, z)*_It(x, y, z)*gamma1;
					J[1][1] = _Iy(x, y, z)*_Iy(x, y, z)*gamma1;
					J[1][2] = _Iy(x, y, z)*_Iz(x, y, z)*gamma1;
					J[1][3] = _Iy(x, y, z)*_It(x, y, z)*gamma1;
					J[2][2] = _Iz(x, y, z)*_Iz(x, y, z)*gamma1;
					J[2][3] = _Iz(x, y, z)*_It(x, y, z)*gamma1;
					/*****************************************************************************************/
					
					/***********************************GRADIENT CONSTANCY************************************/
					Jgrad[0][0] = Ixx  * Ixx + Ixy  * Ixy + Ixz  * Ixz;
					Jgrad[0][1] = Ixx  * Ixy + Ixy  * Iyy + Ixz  * Iyz;
					Jgrad[0][2] = Ixx  * Ixz + Ixy  * Iyz + Ixz  * Izz;
					Jgrad[0][3] = Ixx  * Ixt + Ixy  * Iyt + Ixz  * Izt;
					Jgrad[1][1] = Ixy  * Ixy + Iyy  * Iyy + Iyz  * Iyz;
					Jgrad[1][2] = Ixz  * Ixy + Iyz  * Iyy + Izz  * Iyz;
					Jgrad[1][3] = Ixt  * Ixy + Iyt  * Iyy + Izt  * Iyz;
					Jgrad[2][2] = Ixz  * Ixz + Iyz  * Iyz + Izz  * Izz;
					Jgrad[2][3] = Ixt  * Ixz + Iyt  * Iyz + Izz  * Izt;
					/******************************************************************************************/
					
					/*************************mixed brightness and gradient constancy**************************/
					J[0][0] += gamma2 * Jgrad[0][0];
					J[0][1] += gamma2 * Jgrad[0][1];
					J[0][2] += gamma2 * Jgrad[0][2];
					J[0][3] += gamma2 * Jgrad[0][3];
					J[1][1] += gamma2 * Jgrad[1][1];
					J[1][2] += gamma2 * Jgrad[1][2];
					J[1][3] += gamma2 * Jgrad[1][3];
					J[2][2] += gamma2 * Jgrad[2][2];
					J[2][3] += gamma2 * Jgrad[2][3];
					/******************************************************************************************/

					/********************************multichannel approach*************************************/
					temp00 = J[0][0] + _J00(x, y, z);
					temp01 = J[0][1] + _J01(x, y, z);
					temp02 = J[0][2] + _J02(x, y, z);
					temp03 = J[0][3] + _J03(x, y, z);
					temp11 = J[1][1] + _J11(x, y, z);
					temp12 = J[1][2] + _J12(x, y, z);
					temp13 = J[1][3] + _J13(x, y, z);
					temp22 = J[2][2] + _J22(x, y, z);
					temp23 = J[2][3] + _J23(x, y, z);
					_J00(x, y, z, temp00);
					_J01(x, y, z, temp01);
					_J02(x, y, z, temp02);
					_J03(x, y, z, temp03);
					_J11(x, y, z, temp11);
					_J12(x, y, z, temp12);
					_J13(x, y, z, temp13);
					_J22(x, y, z, temp22);
					_J23(x, y, z, temp23);
					/******************************************************************************************/
				}
			}
		}		
	/**********************************************************************************************************/
}




void opticalFlow::computeSolver(data& data, int iter, matrix& _u, matrix& _v, matrix& _w, matrix& _J00, matrix& _J01, matrix& _J02, matrix& _J03,
	matrix& _J11, matrix& _J12, matrix& _J13, matrix& _J22, matrix& _J23){

	float alpha = 500.0;
	float relax = 1.99;

	for (int n = 0; n < iter; n++){ // each iteration

		cout << "\nIteration: " << n + 1;

		for (int z = 0; z < data.z_dim; z++){ //  each pixel
			for (int y = 0; y < data.y_dim; y++){
				for (int x = 0; x < data.x_dim; x++){

					/***********************solver (GAUSS-SEIDEL METHOD)**************************/

					aux = 6;	//to compute the neighbours

					if (x == 0 || x == data.x_dim - 1)
						aux = aux - 1;

					if (y == 0 || y == data.y_dim - 1)
						aux = aux - 1;

					if (z == 0 || z == data.z_dim - 1)
						aux = aux - 1;

					divx = 1.0 / (aux * alpha + _J00(x,y,z));
					divy = 1.0 / (aux * alpha + _J11(x,y,z));
					divz = 1.0 / (aux * alpha + _J22(x,y,z));

					float tempx = -_J03(x, y, z) - (_J01(x, y, z) * _v(x, y, z) + _J02(x, y, z) * _w(x, y, z) -
						alpha * (_u(x - 1, y, z) + _u(x, y - 1, z) + _u(x, y, z - 1)) -
						alpha * (_u(x + 1, y, z) + _u(x, y + 1, z) + _u(x, y, z + 1)));
					tempx = tempx *divx;
					_u(x, y, z, (1 - relax)*_u(x, y, z) + relax * tempx);

					/*if (x == 247 && y == 233 && z == 88){
						cout << "\nu = " << _u(150, 130, 70) << endl;
						system("pause");
					}*/



					float tempy = -_J13(x, y, z) - (_J01(x, y, z) * _u(x, y, z) + _J12(x, y, z) * _w(x, y, z) -
						alpha * (_v(x - 1, y, z) + _v(x, y - 1, z) + _v(x, y, z - 1)) -
						alpha * (_v(x + 1, y, z) + _v(x, y + 1, z) + _v(x, y, z + 1)));
					tempy = tempy*divy;
					_v(x, y, z, (1 - relax) * _v(x, y, z) + relax * tempy);


					float tempz = -_J23(x, y, z) - (_J02(x, y, z) * _u(x, y, z) + _J12(x, y, z) * _v(x, y, z) -
						alpha * (_w(x - 1, y, z) + _w(x, y - 1, z) + _w(x, y, z - 1)) -
						alpha * (_w(x + 1, y, z) + _w(x, y + 1, z) + _w(x, y, z + 1)));
					tempz = tempz*divz;
					_w(x, y, z, (1 - relax) * _w(x, y, z) + relax * tempz);

					/******************************************************************************************/
				}
			}
		}
	}

}



void opticalFlow::saveData(matrix& _u, matrix& _v, matrix& _w, data& data){
	if (this->OptFlow != NULL)
		delete[] this->OptFlow;

	//allocation of memory of the final array
	int t_size = data.x_dim * data.y_dim * data.z_dim * 3;

	this->OptFlow = new float[t_size];


	cout << "\n\nSaving results..";

	//we save here the generated data in an array to represent it 
	for (int z = 0; z < data.z_dim; z++){ //  each pixel
		for (int y = 0; y < data.y_dim; y++){
			for (int x = 0; x < data.x_dim; x++){
				index = x + (y*data.x_dim) + (z*data.y_dim*data.x_dim);
					OptFlow[(index * 3) + 0] = _u(x, y, z);
					OptFlow[(index * 3) + 1] = _v(x, y, z);
					OptFlow[(index * 3) + 2] = _w(x, y, z);	
			}
		}
	}

	//cout << "\n\nOptFlow[0] = " << OptFlow[0];
}


void opticalFlow::saveDataGround(float x_dat, float y_dat, float z_dat, int cols_x, int rows_y, int depth_z, int x_pos, int y_pos, int z_pos){
	
	if (x_pos == 0 && y_pos == 0 && z_pos == 0){
		if (this->ground != NULL)
			delete[] this->ground;

		//allocation of memory of the final array
		int t_size = cols_x * rows_y * depth_z * 3;

		this->ground = new float[t_size];
	}


	index = x_pos + (y_pos*cols_x) + (z_pos*rows_y*cols_x);
	ground[(index * 3) + 0] = x_dat;
	ground[(index * 3) + 1] = y_dat;
	ground[(index * 3) + 2] = z_dat;
			
}

void opticalFlow::deallocateOptFlow(){
	if (this->OptFlow != NULL){
		delete[] this->OptFlow;
		this->OptFlow = NULL;
	}
}

void opticalFlow::deallocateGround(){
	if (this->ground != NULL){
		delete[] this->ground;
		this->ground = NULL;
	}
}


opticalFlow::~opticalFlow()
{
}