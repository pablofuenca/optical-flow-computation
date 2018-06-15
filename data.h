#ifndef DATA_H
#define DATA_H
 
#include <vtkObject.h>
#include "vtkImageData.h"
#include "vtkPointData.h"
#include <vtkDataArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkImageReader.h>
#include <vtkGlyphSource2D.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkMaskPoints.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkArrayCalculator.h>
#include <vtkContourFilter.h>
#include <vtkImageExtractComponents.h>
#include <vtkXMLImageDataReader.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkImageActor.h>
#include <vtkImageViewer2.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkPolyData.h>
#include <vtkImageGradientMagnitude.h>
//#include <vtkImageSliceMapper.h>
//#include <vtkImageSlice.h>
#include <vtkPolyDataMapper.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageResample.h>

#include "matrix.h"

using namespace std;

class data{
public:
	data();						//constructor
	void loadVolumeScalarData(const char* fileName, int, int, int, int);
	void loadVolumeScalarDataFromMatrix(matrix& _m);
	void loadVolumeVectorData(const char* fileName, int, int, int);
	void loadVolumeVecorDataFromMatrix(matrix& _m);
	void deallocate();
	~data();					//destructor

	//make VolumeDisplay a friend class to pass the 
	friend class VolumeDisplayer;

	float *m;
	float *m2;
	int x_dim;
	int y_dim;
	int z_dim;


private:
	int index;
};

#endif 