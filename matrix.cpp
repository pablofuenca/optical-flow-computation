#include "matrix.h"

/*****************************/
//creates an scalar 3D matrix//
/*****************************/
matrix::matrix(void){
	this->array3D = nullptr;

	this->rows    = -1;
	this->cols    = -1;
	this->depth   = -1;

	this->bound_x = -1;
	this->bound_y = -1;
	this->bound_z = -1;

	this->first  = -1;
	this->second = -1;
	this->third  = -1;
}

matrix::matrix(int _cols, int _rows, int _depth, int _bound_x, int _bound_y, int _bound_z){
	this->array3D = nullptr;

	create3Dmatrix(_cols, _rows, _depth, _bound_x, _bound_y, _bound_z);
}

matrix::~matrix(void){
	freeMatrix3D();
}

float matrix::operator()(int _x, int _y, int _z, int _derivat, int derivative){
	if ((_y >= this->rows))
		_y = this->rows - 1;
	if ((_x >= this->cols))
		_x = this->cols - 1;
	if ((_z >= this->depth))
		_z = this->depth - 1;
	if ((_x < 0))
		_x = 0;
	if ((_y < 0))
		_y = 0;
	if ((_z < 0))
		_z = 0;

	int index = _x + _y*this->third + _z*this->third*this->second;
	return (this->array3D)[index];
}



float matrix::operator()(int _x, int _y, int _z){
	if ((_y >= this->rows))
		return 0;
		//_y = this->rows - 1;
	if ((_x >= this->cols))
		return 0;
		//_x = this->cols - 1;
	if ((_z >= this->depth))
		return 0;
		//_z = this->depth - 1;
	if ((_x < 0))
		return 0;
		//_x = 0;
	if ( (_y < 0))
		return 0;
		//_y = 0;
	if ((_z < 0))
		return 0;
		//_z = 0;

	int index = _x + _y*this->third + _z*this->third*this->second;
	return (this->array3D)[index];

	//return (this->array3D)[_z][_y][_x];
}

void matrix::operator()(int _x, int _y, int _z, float _value){
	if ((_y >= this->rows) || (_y < 0))
		return;
	if ((_x >= this->cols) || (_x < 0))
		return;
	if ((_z >= this->depth) || (_z < 0))
		return;

	int index = _x + _y*this->third + _z*this->third*this->second;
	(this->array3D)[index] = _value;
	//(this->array3D)[_z][_y][_x] = _value;
}

void matrix::freeMatrix3D(){	//deallocate memory
	if (this->array3D == nullptr)
		return;

	delete[] this->array3D;

	this->array3D = nullptr;

	this->rows    = -1;
	this->cols    = -1;
	this->depth   = -1;

	this->bound_x = -1;
	this->bound_y = -1;
	this->bound_z = -1;

	this->first  = -1;
	this->second = -1;
	this->third  = -1;


}


void matrix::create3Dmatrix(int _cols, int _rows, int _depth, int _bound_x, int _bound_y, int _bound_z){
	freeMatrix3D();

	this->rows = _rows;
	this->cols = _cols;
	this->depth = _depth;

	this->bound_x = _bound_x;
	this->bound_y = _bound_y;
	this->bound_z = _bound_z;

	this->first  = this->depth + 2*this->bound_z;
	this->second = this->rows  + 2*this->bound_y;
	this->third  = this->cols  + 2*this->bound_x;

	this->array3D = new float[this->first*this->second*this->third];

	
	fillWithScalar(0.0f);
}


void matrix::fillWithScalar(float _value){
	
	int index = 0;

	for (int i = 0; i < this->first; i++){
		for (int j = 0; j < this->second; j++){
			for (int k = 0; k < this->third; k++){

				index = k + j*this->third + i*this->third*this->second;
				(this->array3D)[index] = _value;
			}
		}
	}
}

//*********************************************************************/
//this function open the file and save the scalar data in the 3D matrix
//*********************************************************************/
// we read the data of each channel, to do that, we need to read every value of the
//input image, therefore we should take only the value of the selected channel, to do 
//that we read the value of every channel, but save in array3D only the value of the 
//channel we are computing.
void matrix::loadVolumeScalarDataInMatrix(const char* fileName, int channel_act, int number_channels){
	//open the file for reading and type binary
	fstream file(fileName, fstream::in|fstream::binary);

	//check if the file is correctly opened
	if (!file.is_open())
		return;

	//auxiliar variable
	unsigned char value;
	int index;

	//fill matrix whit the data
	for(int z = this->bound_z; z < (this->depth+this->bound_z); z++){
		for(int y = this->bound_y; y < (this->rows+this->bound_y); y++){
			for(int x = this->bound_x; x < (this->cols+this->bound_x); x++){

				
				for (int nc = 0; nc < number_channels; nc++){
					file.read((char*)&value, sizeof(unsigned char));

					if (nc == channel_act){
						index = x + y*this->third + z*this->third*this->second;
						array3D[index] = value;
					}
				}
				
			}
		}
	}
}


///****************************/
////creates a vector 3D matrix//
///****************************/
void matrix::createVectorMatrix(int rows, int cols, int depth, int boundaries){


	//matrix for x component
	aux1 = rows;
	aux2 = cols;
	aux3 = depth;

	rows = rows   + (2 * boundaries);
	cols = cols   + (2 * boundaries);
	depth = depth + (2 * boundaries);

	m_vecmatrixX = new float**[depth];

	for (int i = 0; i < depth; i++){
		m_vecmatrixX[i] = new float*[rows];
	}

	for (int i = 0; i < depth; i++){
		for (int j = 0; j < cols; j++){
			m_vecmatrixX[i][j] = new float[cols];
		}
	}

	//matrix for y component

	m_vecmatrixY = new float**[depth];

	for (int i = 0; i < depth; i++){
		m_vecmatrixY[i] = new float*[rows];
	}

	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			m_vecmatrixY[i][j] = new float[cols];
		}
	}


	//matrix for z component

	m_vecmatrixZ = new float**[depth];

	for (int i = 0; i < depth; i++){
		m_vecmatrixZ[i] = new float*[rows];
	}

	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			m_vecmatrixZ[i][j] = new float[cols];
		}
	}
}





//*********************************************************************/
//this function open the file and save the vector data in the 3D matrix
//*********************************************************************/
//maybe I should call this function from createVectorMatrix()
//in order to open directly the folder, or can I ask the client
//for the file's name, after I create the matrix, and call  
//from the main
void matrix::loadVolumeVectorDataInMatrix(const char* fileName){
	//open the file for reading and type binary
	fstream file(fileName, fstream::in|fstream::binary);

	//check if the file is correctly opened
	if (!file.is_open())
		return;

	//auxiliar variable
	float valuex, valuey, valuez;

	//fill matrix whit the data
	for(int z=boundaries; z<(depth-boundaries); z++){
		for(int y=boundaries; y<(rows-boundaries); y++){
			for(int x=boundaries; x<(cols-boundaries); x++){
				file.read((char*)&valuex, sizeof(float));
				m_vecmatrixX[z][y][x] = (float)valuex;
				file.read((char*)&valuey, sizeof(float));
				m_vecmatrixY[z][y][x] = (float)valuey;
				file.read((char*)&valuez, sizeof(float));
				m_vecmatrixZ[z][y][x] = (float)valuez;
			}
		}
	}
//after create the matrix, and fill it with the data I should call
//the function that loads the scalar data from this matrix, but
//since the matrix.h is contained in data.h maybe I should call it from
//the main
}


