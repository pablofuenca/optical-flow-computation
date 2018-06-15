#include "data.h"
#include "VolumeDisplayer.h"
#include "opticalFlow.h"

//constructor
data::data(){
	this->m = NULL;
	this->m2 = NULL;
}

/*************************************/
//load scalar data from our 3D matrix
/*************************************/
void data::loadVolumeScalarDataFromMatrix(matrix& _m){
	x_dim = _m.rows;
	y_dim = _m.cols;
	z_dim = _m.depth;

	if(this->m != NULL)
		delete [] this->m;

	this->m = new float[x_dim*y_dim*z_dim];

//	double value;
	//fill array whit the data of the 3D matrix
	for(int z=0; z<z_dim; z++){
		for(int y=0; y<y_dim; y++){
			for(int x=0; x<x_dim; x++){
				//value = _m.matrix(x,y,z);
				index = x + (y*x_dim) + (z*y_dim*x_dim);
				//m[index] = (double)value;
			}
		}
	}
}


/*************************************/
//loads vector data from our 3D matrix
/*************************************/
void data::loadVolumeVecorDataFromMatrix(matrix& _m){
	x_dim = _m.rows;
	y_dim = _m.cols;
	z_dim = _m.depth;

	if(this->m != NULL)
		delete [] this->m;

	this->m = new float[x_dim*y_dim*z_dim * 3];

	float xx, yy, zz;
	//fill vector whit the data of the 3D matrix that store vector field

	//**********make 3 matrix
	for(int z=0; z<z_dim; z++){
		for(int y=0; y<y_dim; y++){
			for(int x=0; x<x_dim; x++){
				xx = _m.m_vecmatrixX[z][y][x];
				yy = _m.m_vecmatrixY[z][y][x];
				zz = _m.m_vecmatrixZ[z][y][x];

				index = x + (y*x_dim) + (z*y_dim*x_dim);

				m[(index*3)+0] = xx;
				m[(index*3)+1] = yy;
				m[(index*3)+2] = zz;
			}
		}
	}
}


/*****************************************************************/
//this function open the file and save the scalar data in the array
/*****************************************************************/
void data::loadVolumeScalarData(const char* fileName, int rows, int cols, int depth, int contador){
	x_dim = cols;
	y_dim = rows;
	z_dim = depth;


	//auxiliar variable
	unsigned char value;
	float temp;

	if (contador == 0){			//first frame
		if (this->m != NULL)
			delete[] this->m;
		
		//allocation of memory
		this->m = new float[x_dim*y_dim*z_dim];

		//open the file for reading and type binary
		fstream file(fileName, fstream::in | fstream::binary);

		//check if the file is correctly opened
		if (!file.is_open()){
			cout << "\n\nNo image was read\n\n";
			system("pause");
			return;
		}

		//fill array whit the data of img1
		for (int z = 0; z < z_dim; z++){
			for (int y = 0; y < y_dim; y++){
				for (int x = 0; x < x_dim; x++){
					file.read((char*)&value, sizeof(unsigned char));
					index = x + (y*x_dim) + (z*y_dim*x_dim);
					temp = (float)value;
					m[index] = temp;
				}
			}
		}
	}
	else{						//second frame

		if (this->m2 != NULL)
			delete[] this->m2;

		//allocation of memory
		this->m2 = new float[x_dim*y_dim*z_dim];

		//open the file for reading and type binary
		fstream file(fileName, fstream::in | fstream::binary);

		//check if the file is correctly opened
		if (!file.is_open()){
			cout << "\n\nNo se ha leido ninguna imagen\n\n";
			system("pause");
			return;
		}

		//fill array whit the data of img2
		for (int z = 0; z < z_dim; z++){
			for (int y = 0; y < y_dim; y++){
				for (int x = 0; x < x_dim; x++){
					file.read((char*)&value, sizeof(unsigned char));
					index = x + (y*x_dim) + (z*y_dim*x_dim); 
					m2[index] = (float)value;
				}
			}
		}
	}
}


/*****************************************************************/
//this function open the file and save the vector data in the array
/*****************************************************************/
void data::loadVolumeVectorData(const char* fileName, int rows, int cols, int depth){
	x_dim = rows;
	y_dim = cols;
	z_dim = depth;

	if(this->m != NULL)
		delete [] this->m;

	//allocation of memory
	this->m = new float[x_dim*y_dim*z_dim * 3];

	//open the file for reading and type binary
	fstream file(fileName, fstream::in|fstream::binary);

	//check if the file is correctly opened
	if (!file.is_open())
		return;

	//auxiliar variable
	float valuex;
	float valuey;
	float valuez;
	int contador = 0;

	//fill array whit the data
	for(int z=0; z<z_dim; z++){
		for(int y=0; y<y_dim; y++){
			for(int x=0; x<x_dim; x++){
				file.read((char*)&valuez, sizeof(float));
				file.read((char*)&valuey, sizeof(float));
				file.read((char*)&valuex, sizeof(float));
				index = x + (y*x_dim) + (z*y_dim*x_dim);
				m[(index * 3) + 0] = (float)valuez;
				m[(index * 3) + 1] = (float)valuey;
				m[(index * 3) + 2] = (float)valuex;
				contador = contador + 1;
			}
		}
	}
	system("pause");
}

void data::deallocate(){
	if (this->m != NULL){
		delete[] this->m;
		this->m = NULL;
	}

	if(this->m2 != NULL){
		delete[] this->m2;
		this->m2 = NULL;
	}
}

data::~data(){
	this->deallocate();
}