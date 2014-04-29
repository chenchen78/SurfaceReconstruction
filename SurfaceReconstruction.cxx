//#include "itkImageFileWriter.h"

//#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

//#include "CLIModuleTemplateCLP.h"

//add 0814
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkPasteImageFilter.h"



#include <vtkFloatArray.h>


#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkMergePoints.h>
#include <vtkPointSource.h>
#include <vtkPolyDataNormals.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkMath.h>

#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>


#include "vtkPoissonReconstruction.h"
#include "SurfaceReconstructionCLP.h"


// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

//typedef itk::Image< PixelType,  3> InputImageType;
//typedef itk::Image< unsigned char,  3> OutputImageType;

typedef itk::Image< unsigned char, 3> MidImageType;
typedef itk::Image< unsigned char, 2> SliceImageType;

MidImageType::Pointer getVolumeContour (const MidImageType::Pointer midImage)
{

	
	typedef itk::ExtractImageFilter< MidImageType,SliceImageType > 
		extractFilterType;
	typedef itk::BinaryContourImageFilter <SliceImageType,SliceImageType >
		binaryContourImageFilterType;
	typedef itk::CastImageFilter< SliceImageType, MidImageType > 
		CastFilterType;
	typedef itk::PasteImageFilter<MidImageType,MidImageType > 
		PasteFilterType;

		
	//Get Resion, size and start of input midImage	
	MidImageType::RegionType inputRegion =
		midImage->GetLargestPossibleRegion();
	MidImageType::SizeType inputSize =
		inputRegion.GetSize();
	MidImageType::IndexType start = inputRegion.GetIndex();
	MidImageType::RegionType desiredRegion;

	//New a output volume used to contain contour volume
	MidImageType::Pointer  outputContour = MidImageType::New();
	outputContour->SetRegions(inputRegion);
	outputContour->Allocate();
	outputContour = midImage;
	

	extractFilterType::Pointer extractFilter = extractFilterType::New();
	extractFilter->InPlaceOn();
	extractFilter->SetDirectionCollapseToSubmatrix();

	binaryContourImageFilterType::Pointer binaryContourFilter = 
		binaryContourImageFilterType::New ();
	CastFilterType::Pointer castFilter = CastFilterType::New();
	PasteFilterType::Pointer pasteFilter = PasteFilterType::New();

	unsigned int zSize=inputSize[2];
	inputSize[2]=0;
	for (unsigned int k=0; k<zSize; k++)
	{
	//Extract one slice from midImage
	start[2] = k;
	desiredRegion.SetSize(  inputSize  );
	desiredRegion.SetIndex( start );
	extractFilter->SetExtractionRegion( desiredRegion );
	extractFilter->SetInput(midImage );

	//Extract a contour of one slice
	binaryContourFilter->SetInput(extractFilter->GetOutput() );

	//Case one contour slice of 2D to a slice of 3D which can be a input of pasteFilter
	castFilter->SetInput(binaryContourFilter->GetOutput());
	castFilter->UpdateLargestPossibleRegion();
	
	//Paste a contour slice to outputContour volume
	pasteFilter->SetSourceImage(castFilter->GetOutput() );
	pasteFilter->SetDestinationImage( outputContour );
	pasteFilter->SetDestinationIndex( start );

	//sliceImage3D = castFilter->GetOutput();
	pasteFilter->SetSourceRegion( castFilter->GetOutput()->GetBufferedRegion() );
	pasteFilter->Update();
	
	outputContour=pasteFilter->GetOutput();
	

	}
	//outputContour=pasteFilter->GetOutput();

	return outputContour;
}


template <class T>
int DoIt( int argc, char * argv[], T )
{
	PARSE_ARGS;

	typedef    T InputPixelType;
	typedef    T OutputPixelType;

	typedef itk::Image<InputPixelType,  3> InputImageType;
	typedef itk::Image<OutputPixelType, 3> OutputImageType;
	//typedef itk::Image< unsigned char, 3> MidImageType;

	typedef itk::ImageFileReader<InputImageType>  ReaderType;
	typedef itk::ImageFileWriter<OutputImageType> WriterType;
	typedef itk::CastImageFilter< InputImageType, MidImageType > CastFilterType;

	
	typename ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName( inputVolume.c_str() );
	reader->Update();


	//Cast input volume to a volume with unsigned char pixel
	typename CastFilterType::Pointer castInputFilter = CastFilterType::New();
	castInputFilter->SetInput(reader->GetOutput());
	castInputFilter->Update();
	MidImageType::Pointer midImage= MidImageType ::New();
	midImage= castInputFilter->GetOutput();

	//Set un-zero pixel to 255, this is a limit of binaryContourImageFilterType
	MidImageType::SizeType midImageSize =
		midImage->GetLargestPossibleRegion().GetSize();
	MidImageType::IndexType indexMid;

	for(unsigned int k=0; k<midImageSize[2]; k++)
	{
		for(unsigned int j=0; j<midImageSize[1]; j++)
		{
			for(unsigned int i=0; i<midImageSize[0]; i++)
			{
				indexMid[0]=i;
				indexMid[1]=j;
				indexMid[2]=k;
				if ( midImage->GetPixel(indexMid)!=0)
				{

					midImage->SetPixel(indexMid, 255);
				}
			}
		}
	}
	
	//Acquire a contour volume
	MidImageType::Pointer  outContour = getVolumeContour(midImage);

	//Get Spacing, size and origin of contour volume
    MidImageType::RegionType outputRegion =
	outContour->GetLargestPossibleRegion();
	MidImageType::SizeType outputSize =
	outputRegion.GetSize();
    MidImageType::SpacingType outputSpacing =
    outContour->GetSpacing();
    MidImageType::PointType outOrigin=
  	outContour->GetOrigin();
 
	//Get points at boundary  	
	MidImageType::IndexType index;

	vtkSmartPointer<vtkPoints> pointsInit = vtkSmartPointer<vtkPoints>::New();
	

	float pointXYZ[3];
	const float matrixVTK2Slicer [3][3]={{-1, 0, 0}, {0,-1, 0}, {0, 0,1}};
	float pointXYZ_S[3];

    MidImageType::PointType contourPoint;
	
	for(unsigned int k=0; k<outputSize[2]; k++)
	{
		for(unsigned int j=0; j<outputSize[1]; j++)
		{
			for(unsigned int i=0; i<outputSize[0]; i++)
			{
				index[0]=i;
				index[1]=j;
				index[2]=k;
				if ( outContour->GetPixel(index)!=0)
				{
                    outContour->TransformIndexToPhysicalPoint(index,contourPoint);
					pointXYZ[2]=contourPoint[2];
					pointXYZ[1]=contourPoint[1];
					pointXYZ[0]=contourPoint[0];

					//transform a point to Slicer space
					vtkMath::Multiply3x3(matrixVTK2Slicer,pointXYZ,pointXYZ_S);

					pointsInit->InsertNextPoint(pointXYZ_S);

				}
				
			}
	
		}
	}

/*	Add the points to a polydata     */    
  vtkSmartPointer<vtkPolyData> polydata = 
  	vtkSmartPointer<vtkPolyData>::New();      
          polydata->SetPoints(pointsInit); 



	

//Adds a normal to each point in polydata. 
   
  // Construct the normal vectors
	vtkSmartPointer<vtkFloatArray> pointNormalsArray = vtkSmartPointer<vtkFloatArray>::New();
	pointNormalsArray->SetNumberOfComponents(3); //3d normals (ie x,y,z)
	pointNormalsArray->SetNumberOfTuples(polydata->GetNumberOfPoints());


	float pT[3]; //Transformed points
	float pN[3]; //Normal of one transformed point
	double pTTrans[3]; //points to be transform
	double bounds[6];

	polydata->GetBounds(bounds);
 
  // Add the data to the normals array
	for (unsigned int i=0; i< polydata->GetNumberOfPoints(); i++)
	{        
	polydata->GetPoint(i, pTTrans);
	pT[0] = pTTrans[0]-(bounds[1]-(bounds[1]-bounds[0])/2);
	pT[1] = pTTrans[1]-(bounds[3]-(bounds[3]-bounds[2])/2);
	pT[2] = pTTrans[2]-(bounds[5]-(bounds[5]-bounds[4])/2);

	pN[0]=pT[0]/sqrt(pT[0]*pT[0]+pT[1]*pT[1]+pT[2]*pT[2]);
	pN[1]=pT[1]/sqrt(pT[0]*pT[0]+pT[1]*pT[1]+pT[2]*pT[2]);
	pN[2]=pT[2]/sqrt(pT[0]*pT[0]+pT[1]*pT[1]+pT[2]*pT[2]);
  
	pointNormalsArray->SetTuple(i, pN);
	}
  
  
	// Add the normals to the points in the polydata
	polydata->GetPointData()->SetNormals(pointNormalsArray); 


	//Make a vtkPolyData with a vertex on each point.
	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputConnection(polydata->GetProducerPort());
	vertexFilter->Update();


	vtkSmartPointer<vtkPolyData> polydataNew = vtkSmartPointer<vtkPolyData>::New();
	polydataNew->ShallowCopy(vertexFilter->GetOutput()); 
  




	//		  
	vtkSmartPointer<vtkXMLPolyDataWriter> writerPoly = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writerPoly->SetInput(polydataNew);
	writerPoly->SetFileName(outputFile.c_str());
	writerPoly->Update();


	vtkSmartPointer<vtkXMLPolyDataReader> readerPoly = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	readerPoly->SetFileName(outputFile.c_str());
	readerPoly->Update();

	//PoissonReconstruction
	vtkSmartPointer<vtkPoissonReconstruction> poissonFilter = vtkSmartPointer<vtkPoissonReconstruction>::New();
	poissonFilter->SetDepth( depth );
	poissonFilter->SetScale(scale);
	poissonFilter->SetSolverDivide(solverDivide);
	poissonFilter->SetIsoDivide(isoDivide);
	poissonFilter->SetSamplesPerNode(samplesPerNode);
	poissonFilter->SetInputConnection(readerPoly->GetOutputPort());
	poissonFilter->Update();

	//Write the file

	vtkSmartPointer<vtkXMLPolyDataWriter> writerSurface = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writerSurface->SetInputConnection(poissonFilter->GetOutputPort());
	writerSurface->SetFileName(outputSurfaceFile.c_str());
	writerSurface->Update(); 



  return EXIT_SUCCESS;
}

} // end of anonymous namespace

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType(inputVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch( componentType )
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0) );
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0) );
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0) );
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0) );
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
