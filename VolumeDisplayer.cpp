#include "VolumeDisplayer.h"

//constructor by defect
VolumeDisplayer::VolumeDisplayer(void)
{
}

//destructor by defect
VolumeDisplayer::~VolumeDisplayer(void)
{
}

/*****************************************/
//function that represents the Scalar Data
/*****************************************/
void VolumeDisplayer::visualizeScalarVolume(data& data)
{

	vtkFloatArray *ucPointer = vtkFloatArray::New();
	ucPointer->SetNumberOfComponents(1);
	ucPointer->SetArray(data.m, data.x_dim*data.y_dim*data.z_dim, 1);
	ucPointer->SetName("array");

	vtkImageData *image = vtkImageData::New();
	image->Initialize();
	image->SetDimensions(data.x_dim, data.y_dim, data.z_dim);
	image->SetExtent(0, data.x_dim - 1, 0, data.y_dim - 1, 0, data.z_dim - 1);
	//image->SetScalarTypeToDouble();
	//image->SetSpacing(0.00166667, 0.00166667, 0.00166667);
	image->SetNumberOfScalarComponents(1);
	image->GetPointData()->SetScalars(ucPointer);
	image->SetOrigin(0, 0, 0);
	image->Update();

	//gaussian smooth filter
	vtkSmartPointer<vtkImageGaussianSmooth> gaussianSmoothFilter = vtkSmartPointer<vtkImageGaussianSmooth>::New();
	gaussianSmoothFilter->SetInput(image);
	gaussianSmoothFilter->SetDimensionality(3);
	gaussianSmoothFilter->SetRadiusFactors(20, 20, 20);
	gaussianSmoothFilter->SetStandardDeviation(3, 3, 3);
	gaussianSmoothFilter->Update();

	//create iso surface object
	vtkSmartPointer<vtkContourFilter> iso_surface = vtkSmartPointer<vtkContourFilter>::New();
	//connect with calculator
	iso_surface->SetInputConnection(gaussianSmoothFilter->GetOutputPort());
	//set iso value
	iso_surface->SetValue(1, 100);
	//update module
	iso_surface->Update();


	//create the mapper
	vtkSmartPointer<vtkOpenGLPolyDataMapper> mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
	mapper->SetInputConnection( iso_surface->GetOutputPort() );


	//create the actor
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);


	//create the renderer
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);
	renderer->ResetCamera();
	renderer->SetBackground(0.8,0.8,0.8);


	//create the render window
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindow->SetSize( 600, 600 );


	//create the interactor
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

}


/*****************************************/
//function that represents the Vector Data
/*****************************************/
void VolumeDisplayer::visualizeVectorVolume(data& _data, opticalFlow& _OF){
	
	cout << "\n\nCreating visualization...";

	
	//read data from our vector, which stores the data from the input image
	//in this case the data of a vector field
	vtkFloatArray *ucPointer = vtkFloatArray::New();
	ucPointer->SetNumberOfComponents(3);
	ucPointer->SetArray(_OF.OptFlow, _data.x_dim*_data.y_dim*_data.z_dim*3, 1);
	ucPointer->SetName("Field");


	//create the ImageData, where we define the caracteristic of the image
	//that we want to represent. In this case I introduce the spacing
	//and the origin because I am trying to read the 27.024731.raw file
	//and this way I save some time
	vtkImageData *image = vtkImageData::New();
	image->Initialize();
	image->SetDimensions(_data.x_dim, _data.y_dim, _data.z_dim);
	image->SetExtent(0, _data.x_dim-1, 0, _data.y_dim-1, 0, _data.z_dim-1);
	image->SetNumberOfScalarComponents(3);
	image->SetScalarTypeToFloat();
	image->AllocateScalars();
	//image->SetSpacing(0.00166667, 0.00166667, 0.00166667);
	image->SetOrigin(0, 0, 0);
	//insert the data in the IMageData and activate the vectors
	image->GetPointData()->SetScalars(ucPointer);
	image->GetPointData()->SetActiveVectors("Field");
	image->Update();
	
	//create object to reduce number of data samples
	//without this module, all data samples (64x64x64 here) would be used
	//this can result in a slow and cluttered visualization
	vtkSmartPointer<vtkMaskPoints> mask = vtkSmartPointer<vtkMaskPoints>::New();
	//connect with reader object
	mask->SetInput(image);
	//set number of sampling points
	mask->SetMaximumNumberOfPoints(20000);
	//activate random sampling
	mask->RandomModeOn();
	//update module
	//**********same size
	mask->Update();

	//pitfall of VTK!!!
	//active vector field has to be set after vtkMaskPoints
	mask->GetOutput()->GetPointData()->SetActiveVectors("Field");


	//create object for generating the used glyph - a 2D arrow here
	vtkSmartPointer<vtkGlyphSource2D> arrow = vtkSmartPointer<vtkGlyphSource2D>::New();
	//set glyph type
	arrow->SetGlyphTypeToArrow();
	//set size
	arrow->SetScale(4);
	//arrow->SetColor(0.4, 0.5, 0.6);
	//set filled glyphs
	arrow->FilledOn();
	//update module
	arrow->UpdateWholeExtent();

	//create object for visualizing the vector field with glyphs
	vtkSmartPointer<vtkGlyph3D> glyph = vtkSmartPointer<vtkGlyph3D>::New();
	//connect to data source
	glyph->SetInput(mask->GetOutput());
	//connect to glyph source
	glyph->SetSource(arrow->GetOutput());
	//scale glyph with vector magnitude
	//glyph->ScalingOn();

	//to represent the arrows with the same size
	glyph->SetScaleModeToDataScalingOff();
	//scale with vector magnitude
	//glyph->SetScaleModeToScaleByVector();	
	//orient glyphs along vectors
	glyph->OrientOn();
	//glyph->SetRange(1.4, 1.5);
	//set scale factor
	//glyph->SetScaling(0.3);
	////color glyphs by vector magnitude 
    //glyph->SetColorModeToColorByVector();
    //update module
    glyph->Update();

	//create object to map the polygons of the glyph	
	vtkSmartPointer<vtkOpenGLPolyDataMapper> mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
	//connect with the glyph generator
	mapper->SetInputConnection( glyph->GetOutputPort() );

	//create actor object - required by VTK
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	//connect to mapper module
  	actor->SetMapper( mapper );

	//create render object
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	//connect with actor object
  	renderer->AddActor( actor );
	//set background color, default is black
	renderer->SetBackground(0.5, 0.5, 0.5);


	//create render window object
	vtkSmartPointer<vtkRenderWindow> render_window = vtkSmartPointer<vtkRenderWindow>::New();
	//connect with renderer
	render_window->AddRenderer( renderer );
	//set size of the output window
	render_window->SetSize( 1000, 1000 );


	//create interactor object handling user input
	vtkSmartPointer<vtkRenderWindowInteractor> render_window_interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	//connect with render window
	render_window_interactor->SetRenderWindow(render_window);

	//create object defining the input style
	vtkSmartPointer<vtkInteractorStyleTrackballCamera> interactor_style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	//connect with interactor object
	render_window_interactor->SetInteractorStyle( interactor_style );


	//init interactor
	render_window_interactor->Initialize();
	//start interaction and event loop
	render_window_interactor->Start();

}