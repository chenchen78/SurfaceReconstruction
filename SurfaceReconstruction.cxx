//#include "itkImageFileWriter.h"

//#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

//#include "CLIModuleTemplateCLP.h"

//add 0814
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"



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

template <class T>
int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  typedef    T InputPixelType;
  typedef    T OutputPixelType;

  typedef itk::Image<InputPixelType,  3> InputImageType;
  typedef itk::Image<OutputPixelType, 3> OutputImageType;

  typedef itk::ImageFileReader<InputImageType>  ReaderType;
  typedef itk::ImageFileWriter<OutputImageType> WriterType;


  typename ReaderType::Pointer reader = ReaderType::New();

  reader->SetFileName( inputVolume.c_str() );
  reader->Update();


  
//Extract points at boundary  
  
//Get Spacing, size and origin of input image
  const  InputImageType::SpacingType& inputSpacing =
    reader->GetOutput()->GetSpacing();
  const  InputImageType::RegionType& inputRegion =
     reader->GetOutput()->GetLargestPossibleRegion();
  const  InputImageType::SizeType& inputSize =
    inputRegion.GetSize();
  const  InputImageType::PointType& inOrigin=
  	 reader->GetOutput()->GetOrigin();


//Get points at boundary  	
   InputImageType::IndexType index;
  
  vtkSmartPointer<vtkPoints> pointsInit = vtkSmartPointer<vtkPoints>::New();

  unsigned int	i,j,k;
  float pointX,pointY,pointZ;

  for(k=0; k<inputSize[2]; k++)
  {
  	for(j=0; j<inputSize[1]; j++)
  	{
  		for(i=0; i<inputSize[0]; i++)
  		{
  			index[0]=i;
  			index[1]=j;
  			index[2]=k;
  			if ( reader->GetOutput()->GetPixel(index)!=0)
				{
					pointZ= inOrigin[2]+k*inputSpacing[2];
					pointY= inOrigin[1]+j*inputSpacing[1];					 	
					pointX= inOrigin[0]+i*inputSpacing[0];
					pointsInit->InsertNextPoint(pointX, pointY, pointZ);
					break;
         }
      }
      for(i=inputSize[0]-1; i>0;i--)
  		{
  			index[0]=i;
  			index[1]=j;
  			index[2]=k;
  			if (reader->GetOutput()->GetPixel(index)!=0)
			  {
					pointZ= inOrigin[2]+k*inputSpacing[2];
					pointY= inOrigin[1]+j*inputSpacing[1];					 	
					pointX= inOrigin[0]+i*inputSpacing[0];
					pointsInit->InsertNextPoint(pointX, pointY, pointZ);
					break;     
			  }

      }
		}

	  for(i=0; i<inputSize[0]; i++)
	  {
			for(j=0; j<inputSize[1]; j++)
  		{
  			index[0]=i;
  			index[1]=j;
  			index[2]=k;
  			if (reader->GetOutput()->GetPixel(index)!=0)
				{
					pointZ= inOrigin[2]+k*inputSpacing[2];
					pointY= inOrigin[1]+j*inputSpacing[1];					 	
					pointX= inOrigin[0]+i*inputSpacing[0];
					pointsInit->InsertNextPoint(pointX, pointY, pointZ);
					break;     
				}
			}
		 	for(j=inputSize[1]-1; j>0;j--)
  		{
  			index[0]=i;
  			index[1]=j;
  			index[2]=k;
  			if(reader->GetOutput()->GetPixel(index)!=0)
				{	
					pointZ= inOrigin[2]+k*inputSpacing[2];
					pointY= inOrigin[1]+j*inputSpacing[1];					 	
					pointX= inOrigin[0]+i*inputSpacing[0];
					pointsInit->InsertNextPoint(pointX, pointY, pointZ);
					break;     
			 	}

      }
		}
  }

  for(j=0; j<inputSize[1]; j++)
  {
	  for(i=0; i<inputSize[0]; i++)
	  {
			for (k=0; k<inputSize[2]; k++)
  		{
  			index[0]=i;
  			index[1]=j;
  			index[2]=k;
  			if (reader->GetOutput()->GetPixel(index)!=0)
				{
					pointZ= inOrigin[2]+k*inputSpacing[2];
					pointY= inOrigin[1]+j*inputSpacing[1];					 	
					pointX= inOrigin[0]+i*inputSpacing[0];
					pointsInit->InsertNextPoint(pointX, pointY, pointZ);
					break;     
				}
			}
		 	for(k=inputSize[2]-1; k>0;k--)
  		{
  			index[0]=i;
  			index[1]=j;
  			index[2]=k;
  			if(reader->GetOutput()->GetPixel(index)!=0)
				{
					pointZ= inOrigin[2]+k*inputSpacing[2];
					pointY= inOrigin[1]+j*inputSpacing[1];					 	
					pointX= inOrigin[0]+i*inputSpacing[0];
					pointsInit->InsertNextPoint(pointX, pointY, pointZ);
					break;     
			 	}

       }
	  }
  }

 
 //Merge points. Removes coincident points from pointsInit
	
	//Set a point point in the pintUniq
  vtkSmartPointer<vtkPoints> pointsSource = vtkSmartPointer<vtkPoints>::New();
  pointsSource->InsertNextPoint(pointsInit->GetPoint(0));
  vtkPolyData* pointsUniq=vtkPolyData::New();
  pointsUniq->SetPoints(pointsSource);

  vtkIdType id;
  vtkSmartPointer<vtkMergePoints> mergePoints = vtkSmartPointer<vtkMergePoints>::New();
  mergePoints->SetDataSet(pointsUniq);
  mergePoints->InitPointInsertion(pointsUniq->GetPoints(), pointsInit->GetBounds());

  for (vtkIdType m = 0; m <pointsInit->GetNumberOfPoints(); m++)
    {
    mergePoints->InsertUniquePoint(pointsInit->GetPoint(m), id);
    }
  
 
//////////////////////////////////////////////////
//Adds a normal to each point in pointsUniq. 
   
  // Construct the normal vectors
  vtkSmartPointer<vtkFloatArray> pointNormalsArray = vtkSmartPointer<vtkFloatArray>::New();
  pointNormalsArray->SetNumberOfComponents(3); //3d normals (ie x,y,z)
  pointNormalsArray->SetNumberOfTuples(pointsUniq->GetNumberOfPoints());




  
  float pT[3]; //Transformed points
  float pN[3]; //Normal of one transformed point
  double pTTrans[3]; //points to be transform
  double bounds[6];
  
  pointsUniq->GetBounds(bounds);
 
  // Add the data to the normals array
  for (i=0; i< pointsUniq->GetNumberOfPoints(); i++)
  {        
		pointsUniq->GetPoint(i, pTTrans);
		pT[0] = pTTrans[0]-(bounds[1]-(bounds[1]-bounds[0])/2);
	 	pT[1] = pTTrans[1]-(bounds[3]-(bounds[3]-bounds[2])/2);
	 	pT[2] = pTTrans[2]-(bounds[5]-(bounds[5]-bounds[4])/2);

		pN[0]=pT[0]/sqrt(pT[0]*pT[0]+pT[1]*pT[1]+pT[2]*pT[2]);
		pN[1]=pT[1]/sqrt(pT[0]*pT[0]+pT[1]*pT[1]+pT[2]*pT[2]);
		pN[2]=pT[2]/sqrt(pT[0]*pT[0]+pT[1]*pT[1]+pT[2]*pT[2]);
	  
		pointNormalsArray->SetTuple(i, pN);
  }
  
  
  // Add the normals to the points in the polydata
  pointsUniq->GetPointData()->SetNormals(pointNormalsArray); 

 
  //Make a vtkPolyData with a vertex on each point.
  vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  vertexFilter->SetInputConnection(pointsUniq->GetProducerPort());
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
  poissonFilter->SetInputConnection(readerPoly->GetOutputPort());
  poissonFilter->Update();
  
  //Write the file
 
  vtkSmartPointer<vtkXMLPolyDataWriter> writerSurface = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writerSurface->SetInputConnection(poissonFilter->GetOutputPort());
  writerSurface->SetFileName(outputSurfaceFile.c_str());
  writerSurface->Update(); 

//  if(pointsUniq)
//	  pointsUniq->Delete();
 // if(poissonFilter)
//	  poissonFilter->Delete();	


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
