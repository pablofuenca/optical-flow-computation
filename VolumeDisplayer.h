#pragma once
#include "data.h"
#include "opticalFlow.h"

//#include <vtkObject.h>
//#include "vtkImageData.h"
//#include "vtkPointData.h"
//#include <vtkDataArray.h>
//#include <vtkUnsignedCharArray.h>
//#include <vtkDoubleArray.h>
//#include <vtkArrayCalculator.h>
//#include <vtkImageExtractComponents.h>
//#include <vtkXMLImageDataReader.h>
//
//#include <vtkVersion.h>
//#include <vtkSmartPointer.h>
//#include <vtkXMLImageDataWriter.h>
//#include <vtkImageData.h>
//#include <vtkPointData.h>
//#include <vtkDoubleArray.h>
//#include <vtkProperty.h>
//#include <vtkDataSetMapper.h>
//#include <vtkImageActor.h>
//#include <vtkImageViewer2.h>
//#include <vtkXMLImageDataReader.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>

class VolumeDisplayer
{
public:
	VolumeDisplayer(void);
	~VolumeDisplayer(void);

	//this is going to work as long as we create an object named data of the class data
	void visualizeScalarVolume(data& data);
	void visualizeVectorVolume(data& data, opticalFlow& OF);

	//void gaussianSmooth(data& data);

};

