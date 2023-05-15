#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkNew.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTrackballActor.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLineSource.h>
#include <vtkPointPicker.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkTransform.h>
#include <vtkPropPicker.h>
#include <vtkCollectionIterator.h>
#include <vtkStringArray.h>
#include <vtkLineSource.h>

#include <QApplication>
#include <QDockWidget>
#include <QGridLayout>
#include <QLabel>
#include <QMainWindow>
#include <QPointer>
#include <QPushButton>
#include <QVBoxLayout>
#include <QInputDialog>
#include <QFileDialog>
#include <QComboBox>
#include <QStandardItem>
#include <QColorDialog>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;
vtkNew<vtkStringArray> nameArray;

namespace {
    class DeleteActorCallback : public vtkCommand
    {
    public:
        static DeleteActorCallback* New()
        {
            return new DeleteActorCallback;
        }

        virtual void Execute(vtkObject* caller, unsigned long eventId, void* callData)
        {
            cout << "hello";
            // Get the picker and the renderer
            vtkPropPicker* picker = vtkPropPicker::SafeDownCast(caller);
            vtkRenderer* renderer = picker->GetRenderer();

            // Get the actor that was picked and remove it from the renderer
            vtkActor* actor = picker->GetActor();
            renderer->RemoveActor(actor);

            // Update the render window
            vtkRenderWindow* renderWindow = renderer->GetRenderWindow();
            renderWindow->Render();
        }
    };
    
    class ScalingInteractorStyle : public vtkInteractorStyleTrackballActor
    {
    public:
        static ScalingInteractorStyle* New();
        vtkTypeMacro(ScalingInteractorStyle, vtkInteractorStyleTrackballActor);
        double scaleFactor = 1.0;
        virtual void OnMouseWheelForward() override
        {
            scaleFactor += 0.1;
            ScaleActor(this->Interactor,scaleFactor);
            vtkInteractorStyleTrackballActor::OnMouseWheelForward();
        }

        virtual void OnMouseWheelBackward() override
        {
            scaleFactor -= 0.1;
            ScaleActor(this->Interactor, scaleFactor);
            vtkInteractorStyleTrackballActor::OnMouseWheelBackward();
        }

    protected:
        vtkSmartPointer<vtkActor> SelectedActor;

        void ScaleActor(vtkRenderWindowInteractor* interactor,double scalefactor)
        {
            vtkRenderer* renderer;
            vtkActor* actor;
            renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
            actor=renderer->GetActors()->GetLastActor();
            if (actor != nullptr) {
                // Create a transform object with the scaling applied.
                vtkNew<vtkTransform> transform;
                transform->PostMultiply();
                transform->Scale(scalefactor, scalefactor, scalefactor);

                // Apply the transform to the actor.
                actor->SetUserTransform(transform.Get());
                actor->Modified();

                // Render the scene to see the scaling effect applied.
                interactor->GetRenderWindow()->Render(); //window instead of renderer?
            }
        }
    };
    vtkStandardNewMacro(ScalingInteractorStyle);

    vtkNew<vtkPoints> drawLine(double px1,double py1,double px2,double py2) {
        vtkNew<vtkPoints> linepoints;
        linepoints->InsertNextPoint(px1, py1, 0.0);
        linepoints->InsertNextPoint(px2, py2, 0.0);
        return linepoints;
    }

    vtkNew<vtkPoints> drawEllipse(double cx, double cy, double rx, double ry ) {
        // Define ellipse parameters
        vtkNew<vtkPoints> ellipsepoints;
        /*double cx = 0.0; // Center X
        double cy = 0.0; // Center Y
        double rx = 0.2; // X-axis radius (W: half - width)
        double ry = 0.1*/; // Y-axis radius (H: half - height)
        // Create points for ellipse vertices
        for (int i = 0; i <= 360; i++) {
            //convert from degree to radian
            //for 0 ≤ t ≤ 2pi.
            double theta = i * vtkMath::Pi() / 180;
            //x(t) = W cos(t) , add cx to shift to the stated center
            double x = cx + rx * cos(theta);
            //y(t) = H sin(t)  , add cy to shift to the stated center
            double y = cy + ry * sin(theta);
            ellipsepoints->InsertNextPoint(x, y, 0.0); // z=0 because its 2D
        }
        return ellipsepoints;
    }

    vtkNew<vtkPoints> drawRegularPolygon(int number_of_sides, double Polygon_raduis, double polygon_center) {
        // Create points for regular polygon vertices
        vtkNew<vtkPoints> regularPolygonPoints;
        //Pointi = ( R cos( 2πi / n ), R sin(2πi / n )),
        int regularPolygonNoOfSides = number_of_sides;
        double regularPolygonRadius = Polygon_raduis;
        double regularPolygonCenter = polygon_center;
        for (int i = 0; i <= regularPolygonNoOfSides; i++) {
            //double theta = i * vtkMath::Pi() / 180;
            //x= R cos( 2πi / n )
            double x = regularPolygonRadius * cos(2 * vtkMath::Pi() * i / regularPolygonNoOfSides);
            //y= R sin(2πi / n )
            double y = regularPolygonRadius * sin(2 * vtkMath::Pi() * i / regularPolygonNoOfSides);
            regularPolygonPoints->InsertNextPoint(x, y, 0.0); // z=0 because its 2D
        }
        return regularPolygonPoints;
    }

    vtkNew<vtkPoints> drawStar(double radius) {
        // Create points for regular polygon vertices
        vtkNew<vtkPoints> starPoints;
        //rosette general equation
        //x(t) = r * cos(k * t * theta)
        //y(t) = r * sin(k * t * theta)
        //r is the radius of the rosette, controlling the size of the curve.
        //k is the number of petals or points on the rosette, controlling the number of "petals" or "points" in the curve.
        //t is the parameter that varies from 0 to 2*pi, controlling the position of the points on the curve.
        //theta is an additional parameter that can be adjusted to control the shape of the rosette. It is typically a constant value.
        double starRadius = radius;
        // This scales the angle by a factor of 14,
        //which determines the number of points or vertices in the star shape. 
        //In this case, it generates a 14-gon shape by using 14 vertices, resulting in a star-like shape.
        double star_K = 14.0;
        int starNumPoints = 10; // number of points to generate (determines smoothness)
        for (int i = 0; i < starNumPoints; ++i) {
            double star_T = 2 * vtkMath::Pi() * i / starNumPoints; // t ranges from 0 to 2*pi , in radians
            //* 14: This scales the angle by a factor of 14,
            //which determines the number of points or vertices in the star shape. 
            //In this case, it generates a 14-gon shape by using 14 vertices, resulting in a star-like shape.
            double x = starRadius * cos(star_T * star_K);
            double y = starRadius * sin(star_T * star_K);
            //flower if numpoints=100 ,phi = 0.2;
            /*double x = starRadius * (cos(t) * cos(starConstantk * t + phi));
            double y = starRadius * (sin(t) * cos(starConstantk * t + phi));*/
            //sparkle or 4 point star if numpoints = 100, phi = 0.2;
            /*double x = starRadius * pow(cos(t) * cos(starConstantk * t + phi), 3);
            double y = starRadius * pow(sin(t) * cos(starConstantk * t + phi), 3);*/
            starPoints->InsertNextPoint(x, y, 0.0); // z=0 because its 2D
        }
        return starPoints;
    }

    vtkNew<vtkPoints>drawPolyline(double PoriginX, double PoriginY,double px0,double py0,double px1,double py1,double px2,double py2, double px3, double py3) {
       // Create five points.
        double origin[3] = { PoriginX, PoriginY, 0.0 };
        double p0[3] = { px0, py0, 0.0 };
        double p1[3] = { px1, py1, 0.0 };
        double p2[3] = { px2, py2, 0.0 };
        double p3[3] = { px3, py3, 0.0 };

        // Create a vtkPoints object and store the points in it
        vtkNew<vtkPoints> points;
        points->InsertNextPoint(origin);
        points->InsertNextPoint(p0);
        points->InsertNextPoint(p1);
        points->InsertNextPoint(p2);
        points->InsertNextPoint(p3);
        return points;
    }

    vtkNew<vtkPoints>drawIrregularPolygon(double px0,double py0,double px1,double py1,double px2,double py2, double px3, double py3, double px4,double py4,double px5,double py5) {
        vtkNew<vtkPoints> irregularPolygonPoints;
        // Generate random coordinates for each vertex

        irregularPolygonPoints->InsertNextPoint(px0, py0, 0);
        irregularPolygonPoints->InsertNextPoint(px1, py1, 0);
        irregularPolygonPoints->InsertNextPoint(px2, py2, 0);
        irregularPolygonPoints->InsertNextPoint(px3, py3, 0);
        irregularPolygonPoints->InsertNextPoint(px4, py4, 0);
        irregularPolygonPoints->InsertNextPoint(px5, py5, 0);
        /* for (int i = 0; i < 6; i++) {
             double x = (double)rand() / RAND_MAX * maxX - maxX / 2;
             double y = (double)rand() / RAND_MAX * maxY - maxY / 2;
             irregularPolygonPoints->InsertNextPoint(x, y, 0.0);
         }*/

        irregularPolygonPoints->InsertNextPoint(px0, py0, 0);
        return irregularPolygonPoints;
    }

    vtkNew<vtkPoints>drawCircle(double radius) {
        vtkNew<vtkPoints> circlepoints;

        double angle_step = 2.0 * vtkMath::Pi() / 100.0;
        for (int i = 0; i <= 100; i++) {
            double x = radius* cos(i * angle_step);
            double y = radius* sin(i * angle_step);
            double z = 0.0;
            circlepoints->InsertNextPoint(x, y, z);
        }
            return circlepoints;

			//vtkNew<vtkIntArray> identifiers;
			//identifiers->SetName("Identifier");
			//identifiers->InsertNextValue(1); // Set the identifier value to 1
   //         circle->GetFieldData()->AddArray(identifiers);
    }

    vtkNew<vtkPoints>drawArc(double centerX, double centerY, double radius, double startAngle, double endAngle ) {
        vtkNew<vtkPoints> arcpoints;
        /*double center[3] = { -0.1 , -0.1 , 0.0 };
        double radius = 0.3;
        double startAngle = 0.0;
        double endAngle = 90.0;*/
        int numSegments = 30;

        // Add the arc points
        for (int i = 0; i <= numSegments; i++)
        {
            double angle = startAngle + (i / static_cast<double>(numSegments)) * (endAngle - startAngle);
            double x = centerX + radius * cos(angle * vtkMath::Pi() / 180.0);
            double y = centerY + radius * sin(angle * vtkMath::Pi() / 180.0);
           
            arcpoints->InsertNextPoint(x, y, 0.0);
        }
        return arcpoints;
    }

    vtkNew<vtkPoints> drawRectangle(double length,double width) {
        vtkNew<vtkPoints> rectanglepoints;
        rectanglepoints->InsertNextPoint(-length/2, -width/2, 0.0);
        rectanglepoints->InsertNextPoint(length / 2, -width / 2, 0.0);
        rectanglepoints->InsertNextPoint(length / 2, width / 2, 0.0);
        rectanglepoints->InsertNextPoint(-length / 2, width / 2, 0.0);
        rectanglepoints->InsertNextPoint(-length/2, -width/2, 0.0);
        return rectanglepoints;
    }

    vtkNew<vtkPoints> drawTriangle(double p1x,double p1y,double p2x,double p2y,double p3x,double p3y) {
        vtkNew<vtkPoints> trianglepoints;
        trianglepoints->InsertNextPoint(-p1x, -p1y, 0.0);
        trianglepoints->InsertNextPoint(p2x, -p2y, 0.0);
        trianglepoints->InsertNextPoint(p3x, p3y, 0.0);
        trianglepoints->InsertNextPoint(-p1x, -p1y, 0.0);
        //trianglepoints->InsertNextPoint(-0.5, 0.25, 0.0);
        //trianglepoints->InsertNextPoint(-0.5, -0.25, 0.0);
        return trianglepoints;
    }

    vtkNew<vtkPoints> drawRhombus(double sidelength) {
        vtkNew<vtkPoints> rhombuspoints;

        double a = 0.1; // side length
        rhombuspoints->InsertNextPoint(-sidelength, 0, 0); // lower left
        rhombuspoints->InsertNextPoint(0, 2* sidelength, 0); // upper left
        rhombuspoints->InsertNextPoint(sidelength, 0, 0); // upper right
        rhombuspoints->InsertNextPoint(0, -2* sidelength, 0); // lower right
        rhombuspoints->InsertNextPoint(-sidelength, 0, 0);
        return rhombuspoints;
    }


    void selectShape(int index, vtkGenericOpenGLRenderWindow* window,vtkRenderer* renderer) {
        vtkNew<vtkPoints> points;
        vtkNew<vtkActor> Actor;
        if (index == 1) {//line
            double x1 = QInputDialog::getDouble(NULL, "Enter first coordinates", "x1 coordinate", 0, -1000, 1000, 2);
            double y1 = QInputDialog::getDouble(NULL, "Enter first coordinates", "y1 coordinate", 0, -1000, 1000, 2);
            double x2 = QInputDialog::getDouble(NULL, "Enter second coordinates", "x2 coordinate", 0, -1000, 1000, 2);
            double y2 = QInputDialog::getDouble(NULL, "Enter second coordinates", "y2 coordinate", 0, -1000, 1000, 2);
            //Actor->SetObjectName("line");
            nameArray->InsertNextValue("line");
            points = drawLine(x1,y1,x2,y2);
        }
        else if (index == 2) {//polyline
            double originX = QInputDialog::getDouble(NULL, "Enter first coordinates", "x1 coordinate", 0, -1000, 1000, 2);
            double originY = QInputDialog::getDouble(NULL, "Enter first coordinates", "y1 coordinate", 0, -1000, 1000, 2);
            double x0 = QInputDialog::getDouble(NULL, "Enter second coordinates", "x2 coordinate", 0, -1000, 1000, 2);
            double y0 = QInputDialog::getDouble(NULL, "Enter second coordinates", "y2 coordinate", 0, -1000, 1000, 2);
            double x1 = QInputDialog::getDouble(NULL, "Enter third coordinates", "x3 coordinate", 0, -1000, 1000, 2);
            double y1 = QInputDialog::getDouble(NULL, "Enter third coordinates", "y3 coordinate", 0, -1000, 1000, 2);
            double x2 = QInputDialog::getDouble(NULL, "Enter fourth coordinates", "x4 coordinate", 0, -1000, 1000, 2);
            double y2 = QInputDialog::getDouble(NULL, "Enter fourth coordinates", "y4 coordinate", 0, -1000, 1000, 2);
            double x3 = QInputDialog::getDouble(NULL, "Enter fifth coordinates", "x5 coordinate", 0, -1000, 1000, 2);
            double y3 = QInputDialog::getDouble(NULL, "Enter fifth coordinates", "y5 coordinate", 0, -1000, 1000, 2);
            Actor->SetObjectName("polyline");
            nameArray->InsertNextValue("polyline");
            points = drawPolyline(originX, originY, x0, y0,x1,y1,x2,y2,x3,y3);
        }
        else if (index == 3) {//irregular polygon
            double x0 = QInputDialog::getDouble(NULL, "Enter first coordinates", "x1 coordinate", 0, -1000, 1000, 2);
            double y0 = QInputDialog::getDouble(NULL, "Enter first coordinates", "y1 coordinate", 0, -1000, 1000, 2);
            double x1 = QInputDialog::getDouble(NULL, "Enter second coordinates", "x2 coordinate", 0, -1000, 1000, 2);
            double y1 = QInputDialog::getDouble(NULL, "Enter second coordinates", "y2 coordinate", 0, -1000, 1000, 2);
            double x2 = QInputDialog::getDouble(NULL, "Enter thrid coordinates", "x3 coordinate", 0, -1000, 1000, 2);
            double y2 = QInputDialog::getDouble(NULL, "Enter third coordinates", "y3 coordinate", 0, -1000, 1000, 2);
            double x3 = QInputDialog::getDouble(NULL, "Enter fourth coordinates", "x4 coordinate", 0, -1000, 1000, 2);
            double y3 = QInputDialog::getDouble(NULL, "Enter fourth coordinates", "y4 coordinate", 0, -1000, 1000, 2);
            double x4 = QInputDialog::getDouble(NULL, "Enter fifth coordinates", "x5 coordinate", 0, -1000, 1000, 2);
            double y4 = QInputDialog::getDouble(NULL, "Enter fifth coordinates", "y5 coordinate", 0, -1000, 1000, 2);
            double x5 = QInputDialog::getDouble(NULL, "Enter sixth coordinates", "x6 coordinate", 0, -1000, 1000, 2);
            double y5 = QInputDialog::getDouble(NULL, "Enter sixth coordinates", "y6 coordinate", 0, -1000, 1000, 2);
            Actor->SetObjectName("irregular polygon");
            nameArray->InsertNextValue("irregular polygon");
            points = drawIrregularPolygon(x0, y0, x1, y1, x2, y2, x3, y3, x4, y4, x5, y5);
        }
        else if (index == 4) {//regular polygon
            int number_of_sides = QInputDialog::getDouble(NULL, "Enter number of sides", "number of sides", 0, -1000, 1000, 2);
            double Polygon_raduis = QInputDialog::getDouble(NULL, "Enter polygon raduis", "polygon raduis", 0, -1000, 1000, 2);
            double polygon_center = QInputDialog::getDouble(NULL, "Enter polygon center", "polygon center", 0, -1000, 1000, 2);
            Actor->SetObjectName("regular polygon");
            nameArray->InsertNextValue("regular polygon");
            points = drawRegularPolygon( number_of_sides, Polygon_raduis, polygon_center);
        }        
        else if (index == 5) {//circle
            double radius = QInputDialog::getDouble(NULL, "Enter radius ", "radius", 0, -1000, 1000, 2);
            Actor->SetObjectName("circle");
            nameArray->InsertNextValue("circle");
            points = drawCircle(radius);
        }        
        else if (index == 6) {//arc
			double centerX = QInputDialog::getDouble(NULL, "Enter point1", "point1", 0, -1000, 1000, 2);
			double centerY = QInputDialog::getDouble(NULL, "Enter point2", "point2", 0, -1000, 1000, 2);
			double radius = QInputDialog::getDouble(NULL, "Enter radius", "radius", 0, -1000, 1000, 2);
            double startangle = QInputDialog::getDouble(NULL, "Enter start angle", "start angle", 0, -1000, 1000, 2);
			double endAngle = QInputDialog::getDouble(NULL, "Enter end angle", "end angle", 0, -1000, 1000, 2);
            Actor->SetObjectName("arc");
            nameArray->InsertNextValue("arc");
            points = drawArc(centerX, centerY, radius, startangle, endAngle );
        }        
        else if (index == 7) {//ellipse
			double centerX = QInputDialog::getDouble(NULL, "Enter point1", "point1", 0, -1000, 1000, 2);
			double centerY = QInputDialog::getDouble(NULL, "Enter point2", "point2", 0, -1000, 1000, 2);
			double radiusA = QInputDialog::getDouble(NULL, "Enter radiusI", "radiusI", 0, -1000, 1000, 2);
			double radiusB = QInputDialog::getDouble(NULL, "Enter radiusII", "radiusII", 0, -1000, 1000, 2);
            Actor->SetObjectName("ellipse");
            nameArray->InsertNextValue("ellipse");
            points = drawEllipse(centerX, centerY, radiusA, radiusB );
        }
        else if (index == 8) {//rectangle
            double length = QInputDialog::getDouble(NULL, "Enter length", "x1 coordinate", 0, -1000, 1000, 2);
            double width = QInputDialog::getDouble(NULL, "Enter width", "y1 coordinate", 0, -1000, 1000, 2);
            Actor->SetObjectName("rectangle");
            nameArray->InsertNextValue("rectangle");
            points = drawRectangle(length,width);
        }
        else if (index == 9) {//triangle
            double p1x = QInputDialog::getDouble(NULL, "Enter point1", "x1 coordinate", 0, -1000, 1000, 2);
            double p1y = QInputDialog::getDouble(NULL, "Enter point1", "y1 coordinate", 0, -1000, 1000, 2);
            double p2x = QInputDialog::getDouble(NULL, "Enter point2", "x2 coordinate", 0, -1000, 1000, 2);
            double p2y = QInputDialog::getDouble(NULL, "Enter point2", "y2 coordinate", 0, -1000, 1000, 2);
            double p3x = QInputDialog::getDouble(NULL, "Enter point3", "x3 coordinate", 0, -1000, 1000, 2);
            double p3y = QInputDialog::getDouble(NULL, "Enter point3", "y3 coordinate", 0, -1000, 1000, 2);
            Actor->SetObjectName("triangle");
            nameArray->InsertNextValue("triangle");
            points = drawTriangle(p1x,p1y,p2x,p2y,p3x,p3y);
        }
        else if (index == 10) {//rhombus
            double sidelength = QInputDialog::getDouble(NULL, "Enter side length", "side length", 0, -1000, 1000, 2);
            Actor->SetObjectName("rhombus");
            nameArray->InsertNextValue("rhombus");
            points = drawRhombus(sidelength);
        }
        else if (index == 11) {//star
            double radius = QInputDialog::getDouble(NULL, "Enter radius", "radius", 0, -1000, 1000, 2);
            Actor->SetObjectName("star");
            nameArray->InsertNextValue("star");
            points = drawStar(radius);
        }
        // Create polyline to connect those points using lines
        //vtkpolyline: type of VTK cell that represents a single polyline in 3D space.
        //vtkNew<vtkPolyLine> polyline;
        vtkNew<vtkLineSource> linesource;
        //ses the number of point IDs in the vtkIdList associated with the polyline object
        //polyline->GetPointIds()->SetNumberOfIds(points->GetNumberOfPoints());
        // iterate through each point in the points object and set the corresponding point ID in the vtkIdList associated with the polyline object. 
       /* for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
            polyline->GetPointIds()->SetId(i, i);
        }*/
        //vtkPolyData:VTK data object that represents a dataset consisting of points, cells, and associated data attributes.
        //vtkNew<vtkPolyData> polydata;
        //set the points object as the points of the polydata object
        //polydata->SetPoints(points);
        linesource->SetPoints(points);
        //allocate memory for the cells in the polydata object. 
        //polydata->Allocate();
        //The polyline(parameter1) will be drawn using lines connecting the points(parameter2) defined by the point IDs in the vtkIdList.
        //polydata->InsertNextCell(polyline->GetCellType(), polyline->GetPointIds());
        vtkNew<vtkPolyDataMapper> mapper;
        //mapper takes data that is going to be rendered
        //mapper->SetInputData(polydata);
         mapper->SetInputConnection(linesource->GetOutputPort());
        //actor is used to change properties
        Actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        Actor->SetVisibility(true);
        Actor->SetMapper(mapper);
        renderer->AddActor(Actor);
        window->Render();
    }

    /*void deleteSelectedShape(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        int numActors = renderer->GetActors()->GetNumberOfItems();
        if (numActors > 0) {
            vtkProp* lastActor = renderer->GetActors()->GetLastProp();
            renderer->RemoveActor(lastActor);
            window->Render();
        }
    }*/

	/*void deleteSelectedShape(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer, int shapeIdentifier) {
		int numActors = renderer->GetActors()->GetNumberOfItems();
		if (numActors > 0) {
			for (int i = 0; i < numActors; i++) {
				vtkProp* actor = renderer->GetActors()->GetItemAsObject(i);
				vtkDataObject* dataObject = actor->GetMapper()->GetInputDataObject(0, 0);
				vtkFieldData* fieldData = dataObject->GetFieldData();
				if (fieldData->HasArray("Identifier")) {
					int identifier = fieldData->GetArray("Identifier")->GetTuple1(0);
					if (identifier == shapeIdentifier) {
						renderer->RemoveActor(actor);
						window->Render();
						return;
					}
				}
			}
			std::cout << "Shape with identifier " << shapeIdentifier << " not found." << std::endl;
		}
		else {
			std::cout << "No shapes to delete." << std::endl;
		}
	}*/

    void applyTranslation(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        double x = QInputDialog::getDouble(NULL, "Enter x translation", "x translation", 0, -1000, 1000, 3);
        double y = QInputDialog::getDouble(NULL, "Enter y translation", "y translation", 0, -1000, 1000, 3);
        vtkProp* prop = renderer->GetActors()->GetLastProp();
        vtkActor* actor = renderer->GetActors()->GetLastActor();
        vtkMapper* mapper = actor->GetMapper();
        vtkPolyData* data = vtkPolyData::SafeDownCast(mapper->GetInput());
        vtkPoints* points = data->GetPoints();
        vtkIdType numPoints = points->GetNumberOfPoints();
        //translate each point Q = P + d
        //Q -> new place
        //P-> old place
        //d-> translation
        for (vtkIdType i = 0; i < numPoints; i++) {
            double* pt = points->GetPoint(i);
            pt[0] += x;
            pt[1] += y;
            points->SetPoint(i, pt);
        }
        points->Modified();
        window->Render();
    }

    void applyRotation(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        double angle = QInputDialog::getDouble(NULL, "Enter angle in degrees", "Rotation angle", 0, -360, 360, 3);
        double angleRad = vtkMath::RadiansFromDegrees(angle);

        vtkProp* prop = renderer->GetActors()->GetLastProp();
        vtkActor* actor = renderer->GetActors()->GetLastActor();
        vtkMapper* mapper = actor->GetMapper();
        vtkPolyData* data = vtkPolyData::SafeDownCast(mapper->GetInput());
        vtkPoints* points = data->GetPoints();
        vtkIdType numPoints = points->GetNumberOfPoints();
        //create an empty matrix
        vtkMatrix4x4* rotation = vtkMatrix4x4::New();
        //identity matrix initially
        rotation->Identity();
        //Qx= Pxcos(θ) − Pysin(θ)
        //Qy = Pxsin(θ) + Pycos(θ)
        //Q -> new place
        //P-> old place
        //set element sets the matrix,its parameters are row,column,value
        rotation->SetElement(0, 0, cos(angleRad));
        rotation->SetElement(0, 1, -sin(angleRad));
        rotation->SetElement(1, 0, sin(angleRad));
        rotation->SetElement(1, 1, cos(angleRad));
        for (vtkIdType i = 0; i < numPoints; i++) {
            double* pt = points->GetPoint(i);
            double p[4] = { pt[0], pt[1], pt[2], 1 };
            double q[4] = { 0, 0, 0, 0 };//store result of array multiplication in q
            rotation->MultiplyPoint(p, q);//multiply matrix of rotation by the matrix of points p and store it in q
            pt[0] = q[0];
            pt[1] = q[1];
            pt[2] = q[2];
            points->SetPoint(i, pt);
        }
        points->Modified();
        window->Render();
    }

    void applyScaling(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        double scalingfactor = QInputDialog::getDouble(NULL, "enter scaling factor", "scaling factor", 0, -1000, 1000, 3);

        vtkProp* prop = renderer->GetActors()->GetLastProp();
        vtkActor* actor = renderer->GetActors()->GetLastActor();
        vtkMapper* mapper = actor->GetMapper();
        vtkPolyData* data = vtkPolyData::SafeDownCast(mapper->GetInput());
        vtkPoints* points = data->GetPoints();
        vtkIdType numPoints = points->GetNumberOfPoints();
        //(Qx, Qy) = (SxPx, SyPy), multiply each old point p by the scaling factor
        for (vtkIdType i = 0; i < numPoints; i++) {
            double* pt = points->GetPoint(i);
            pt[0] *= scalingfactor;
            pt[1] *= scalingfactor;            
            points->SetPoint(i, pt);
        }
        points->Modified();
        window->Render();
    }

    void applyShear(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
            double shearFactorX = QInputDialog::getDouble(NULL, "Enter shear factor X", "shear factor X", 0, -1000, 1000, 3);
            double shearFactorY = QInputDialog::getDouble(NULL, "Enter shear factor Y", "shear factor Y", 0, -1000, 1000, 3);

            vtkProp* prop = renderer->GetActors()->GetLastProp();
            vtkActor* actor = dynamic_cast<vtkActor*>(prop);

            vtkNew<vtkTransform> transform;
            double shearElements[16] = { 1.0, shearFactorX, 0.0, 0.0,
                                        shearFactorY, 1.0, 0.0, 0.0,
                                       0.0, 0.0, 1.0, 0.0,
                                       0.0, 0.0, 0.0, 1.0 };
            // Set the shearing matrix in the transform
            transform->SetMatrix(shearElements);
            if (actor != nullptr) {
                // Set the actor's user transform to be the shearing transform.
                actor->SetUserMatrix(transform->GetMatrix());
                actor->Modified();
                // Render the scene to see the shearing effect applied.
                window->Render();
            }
        }

    void openColorWindow(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        QColorDialog colorDialog;
        if (colorDialog.exec() == QDialog::Accepted) {
            QColor color = colorDialog.currentColor();
            int red = color.red();
            int green = color.green();
            int blue = color.blue();
            vtkActor* lastActor = vtkActor::SafeDownCast(renderer->GetActors()->GetLastProp());
            if (lastActor) {
                lastActor->GetProperty()->SetColor(red / 255.0, green / 255.0, blue / 255.0);
                window->Render();
            }
        }
    }

    void changeLineWidth(int value, vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        vtkActor* lastActor = vtkActor::SafeDownCast(renderer->GetActors()->GetLastProp());
        lastActor->GetProperty()->SetLineWidth(value);
        window->Render();
    }

	void writeInFile(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
		QString fileName = QFileDialog::getSaveFileName(nullptr, "Save File", ".", "Text Files (*.txt)");
		QFile file(fileName);
        vtkSmartPointer<vtkActorCollection> actors = renderer->GetActors();
        //vtkCollectionSimpleIterator* it = actors->NewIterator();
        vtkSmartPointer<vtkCollectionIterator> it;
        int counter = 0;
        for (it = actors->NewIterator(); !it->IsDoneWithTraversal(); it->GoToNextItem())
        {

            vtkSmartPointer<vtkActor> actor = vtkActor::SafeDownCast(it->GetCurrentObject());
            if (actor)
            {
                // Get actor properties
                double* color = actor->GetProperty()->GetColor();
                // Do something with actor properties
                if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                    QTextStream out(&file);
                    vtkSmartPointer<vtkPolyData> polyData = vtkPolyData::SafeDownCast(actor->GetMapper()->GetInput());
                    vtkSmartPointer<vtkPoints> points = polyData->GetPoints();
                    out << nameArray->GetValue(counter) << " ";
                    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
                    {
                        double p[3];
                        points->GetPoint(i, p);
                        out << p[0] << " " << p[1] << " " << p[2] << " ";
                        // do something with the point coordinates of each actor
                    }
                    out << color[0] << " "
                        << color[1] << " "
                        << color[2] << " "
                        << actor->GetProperty()->GetLineWidth() << Qt::endl;
                }
                counter++;
            }
        }
        file.close();
	}

	//void readInputFile(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
	//	QString fileObject = QFileDialog::getOpenFileName(nullptr, "Open File", ".", "Text Files (*.txt)");
	//	if (fileObject.isEmpty()) {
	//		return;  // Dialog was cancelled
	//	}
	//	QFile file(fileObject);
	//	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	//		return;
	//	else {
	//		QTextStream in(&file);
	//		double x1, y1, x2, y2;
	//		if (!in.atEnd()) {
	//			QStringList linepoint1 = in.readLine().split(" ");
	//			linesource->SetPoint1(linepoint1[0].toDouble(), linepoint1[1].toDouble(), 0.0);
	//		}
	//		if (!in.atEnd()) {
	//			QStringList linepoint2 = in.readLine().split(" ");
	//			linesource->SetPoint2(linepoint2[0].toDouble(), linepoint2[1].toDouble(), 0.0);
	//		}
	//		if (!in.atEnd()) {
	//			QStringList rgb = in.readLine().split(" ");
	//			double rgbarr[3];
	//			rgbarr[0] = rgb.at(0).toDouble();
	//			rgbarr[1] = rgb.at(1).toDouble();
	//			rgbarr[2] = rgb.at(2).toDouble();
	//			lineActor->GetProperty()->SetColor(rgbarr);
	//		}
	//		if (!in.atEnd()) {
	//			QStringList lineWidthString = in.readLine().split(" ");;
	//			double lineWidth = lineWidthString.at(0).toDouble();
	//			lineActor->GetProperty()->SetLineWidth(lineWidth);
	//		}
	//		window->Render();
	//		updateTextCoordinates(linesource, TextActor, lineActor);
	//		file.close(); //close the file object.
	//	}
	//}

} // namespace

int main(int argc, char** argv)
{
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());

    QApplication app(argc, argv);

    // main window
    QMainWindow mainWindow;
    mainWindow.resize(1200, 900);

    // control area
    QDockWidget controlDock;
    mainWindow.addDockWidget(Qt::LeftDockWidgetArea, &controlDock);

    QLabel controlDockTitle("Control Dock");
    controlDockTitle.setMargin(20);
    controlDock.setTitleBarWidget(&controlDockTitle);

    QPointer<QVBoxLayout> dockLayout = new QVBoxLayout();
    QWidget layoutContainer;
    layoutContainer.setLayout(dockLayout);
    controlDock.setWidget(&layoutContainer);

    //QComboBox* shapesComboBoxparent ;
    QComboBox shapesComboBox;
    shapesComboBox.addItem(QApplication::tr("Select Shape"));
    shapesComboBox.addItem(QApplication::tr("Line"));
    shapesComboBox.addItem(QApplication::tr("Polyline"));
    shapesComboBox.addItem(QApplication::tr("Irregular Polygon"));
    shapesComboBox.addItem(QApplication::tr("Regular Polygon"));
    shapesComboBox.addItem(QApplication::tr("Circle"));
    shapesComboBox.addItem(QApplication::tr("Arc"));
    shapesComboBox.addItem(QApplication::tr("Ellipse"));
    shapesComboBox.addItem(QApplication::tr("Rectangle"));
    shapesComboBox.addItem(QApplication::tr("Triangle"));
    shapesComboBox.addItem(QApplication::tr("Rhombus"));
    shapesComboBox.addItem(QApplication::tr("Star"));
    shapesComboBox.setItemData(0, QVariant(0), Qt::UserRole - 1);

    dockLayout->addWidget(&shapesComboBox, 1, Qt::AlignTop);

    QPushButton deleteButton("Delete Selected Shape");
    dockLayout->addWidget(&deleteButton, 1, Qt::AlignTop);

    QPushButton translateButton("Apply Translate");
    dockLayout->addWidget(&translateButton, 1, Qt::AlignTop);
    
    QPushButton rotateButton("Apply Rotate");
    dockLayout->addWidget(&rotateButton, 1, Qt::AlignTop);

    QPushButton scaleButton("Apply Scale");
    dockLayout->addWidget(&scaleButton, 1, Qt::AlignTop);

    QPushButton shearButton("Apply shearing");
    dockLayout->addWidget(&shearButton, 1, Qt::AlignTop);

    QPushButton colorButton;
    colorButton.setText("Select color");
    dockLayout->addWidget(&colorButton, 0, Qt::AlignTop);

    QSlider slider;
    slider.setMinimum(0);
    slider.setMaximum(10);
    slider.setValue(0);
    slider.setOrientation(Qt::Horizontal);
    dockLayout->addWidget(&slider, 1, Qt::AlignTop);
    QPushButton readFile;
    readFile.setText("Read Input File");
    dockLayout->addWidget(&readFile, 0, Qt::AlignTop);

    QPushButton writeFile;
    writeFile.setText("Write Input File");
    dockLayout->addWidget(&writeFile, 1, Qt::AlignTop);



     //render area
    QPointer<QVTKOpenGLNativeWidget> vtkRenderWidget = new QVTKOpenGLNativeWidget();
    mainWindow.setCentralWidget(vtkRenderWidget);
    mainWindow.setWindowTitle("VTK Line Example");

    // VTK part
    vtkNew<vtkGenericOpenGLRenderWindow> window;
    vtkRenderWidget->setRenderWindow(window);

    /*---------------------renderers---------------------*/
    // Create a callback for deleting actors

    // Attach the callback to the picker
    vtkNew<vtkPropPicker> picker;
    window->GetInteractor()->SetPicker(picker);
    vtkNew<DeleteActorCallback> callback;
    picker->AddObserver(vtkCommand::EndPickEvent, callback);

    vtkNew<vtkRenderer> renderer;
    window->AddRenderer(renderer);
    
    vtkNew<vtkPointPicker> pointPicker;
    window->SetInteractor(vtkRenderWidget->interactor());

    vtkNew<ScalingInteractorStyle> style;

    //vtkNew<customMouseInteractorStyle> style;
    //style->setLineSource(linesource);
    //style->setVTKActor(lineactor);

    window->GetInteractor()->SetInteractorStyle(style);

       QObject::connect(&shapesComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), [&](int index) {
          ::selectShape(index,window,renderer);
          });

       // connect button to slot
       //QObject::connect(&deleteButton, &QPushButton::clicked, [&]() {
       //    ::deleteSelectedShape(window, renderer);
       //    });
       QObject::connect(&translateButton, &QPushButton::clicked, [&]() {
           ::applyTranslation(window, renderer);
           });
       QObject::connect(&rotateButton, &QPushButton::clicked, [&]() {
           ::applyRotation(window, renderer);
           });
       QObject::connect(&scaleButton, &QPushButton::clicked, [&]() {
           ::applyScaling(window, renderer);
           });
       QObject::connect(&shearButton, &QPushButton::clicked, [&]() {
           ::applyShear(window, renderer);
           });
       QObject::connect(&colorButton, &QPushButton::released,
           [&]() { ::openColorWindow(window, renderer); });
       QObject::connect(&slider, &QSlider::valueChanged, [&](int value) {
           // Do something with the value
           ::changeLineWidth(value,window, renderer);
           });
       /*QObject::connect(&readFile, &QPushButton::released,
           [&]() { ::readInputFile(window, renderer); });*/

       QObject::connect(&writeFile, &QPushButton::released,
           [&]() { ::writeInFile(window, renderer); });

    window->Render();
    mainWindow.show();

    return app.exec();
    //return EXIT_SUCCESS;
}
