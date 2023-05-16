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
#include <vtkMatrix3x3.h>

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
#include <QMessageBox>


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <optional>

using namespace std;
struct Shapes {
    string name;
    optional<double> x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, length, width, radius1, radius2, radius3, startangle, endangle;
    optional<int> numofsides;
};
vector<Shapes> Shapesdata;

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
            ScaleActor(this->Interactor, scaleFactor);
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

        void ScaleActor(vtkRenderWindowInteractor* interactor, double scalefactor)
        {
            vtkRenderer* renderer;
            vtkActor* actor;
            renderer = interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
            actor = renderer->GetActors()->GetLastActor();
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
    /*-------------------------------------------Draw Functions 2D------------------------------------------------------------*/
    vtkNew<vtkPoints> drawLine(double px1, double py1, double px2, double py2) {
        vtkNew<vtkPoints> linepoints;
        linepoints->InsertNextPoint(px1, py1, 0.0);
        linepoints->InsertNextPoint(px2, py2, 0.0);
        return linepoints;
    }

    vtkNew<vtkPoints> drawEllipse(double rx, double ry) {
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
            double x = rx * cos(theta);
            //y(t) = H sin(t)  , add cy to shift to the stated center
            double y = ry * sin(theta);
            ellipsepoints->InsertNextPoint(x, y, 0.0); // z=0 because its 2D
        }

        return ellipsepoints;
    } //i feel like center is useless here

    vtkNew<vtkPoints> drawRegularPolygon(int number_of_sides, double Polygon_raduis) {
        // Create points for regular polygon vertices
        vtkNew<vtkPoints> regularPolygonPoints;
        //Pointi = ( R cos( 2πi / n ), R sin(2πi / n )),
        int regularPolygonNoOfSides = number_of_sides;
        double regularPolygonRadius = Polygon_raduis;
        // double regularPolygonCenter = polygon_center;
        for (int i = 0; i <= regularPolygonNoOfSides; i++) {
            //double theta = i * vtkMath::Pi() / 180;
            //x= R cos( 2πi / n )
            double x = regularPolygonRadius * cos(2 * vtkMath::Pi() * i / regularPolygonNoOfSides);
            //y= R sin(2πi / n )
            double y = regularPolygonRadius * sin(2 * vtkMath::Pi() * i / regularPolygonNoOfSides);
            regularPolygonPoints->InsertNextPoint(x, y, 0.0); // z=0 because its 2D
        }
        return regularPolygonPoints;
    } //i feel like center is useless here

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

    vtkNew<vtkPoints>drawPolyline(double PoriginX, double PoriginY, double px0, double py0, double px1, double py1, double px2, double py2, double px3, double py3) {
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

    vtkNew<vtkPoints>drawIrregularPolygon(double px0, double py0, double px1, double py1, double px2, double py2, double px3, double py3, double px4, double py4, double px5, double py5) {
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
            double x = radius * cos(i * angle_step);
            double y = radius * sin(i * angle_step);
            double z = 0.0;
            circlepoints->InsertNextPoint(x, y, z);
        }
        return circlepoints;

        //vtkNew<vtkIntArray> identifiers;
        //identifiers->SetName("Identifier");
        //identifiers->InsertNextValue(1); // Set the identifier value to 1
//         circle->GetFieldData()->AddArray(identifiers);
    }

    vtkNew<vtkPoints>drawArc(double radius, double startAngle, double endAngle) {
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
            double x = radius * cos(angle * vtkMath::Pi() / 180.0);
            double y = radius * sin(angle * vtkMath::Pi() / 180.0);

            arcpoints->InsertNextPoint(x, y, 0.0);
        }
        return arcpoints;
    }

    vtkNew<vtkPoints> drawRectangle(double length, double width) {
        vtkNew<vtkPoints> rectanglepoints;
        rectanglepoints->InsertNextPoint(-length / 2, -width / 2, 0.0);
        rectanglepoints->InsertNextPoint(length / 2, -width / 2, 0.0);
        rectanglepoints->InsertNextPoint(length / 2, width / 2, 0.0);
        rectanglepoints->InsertNextPoint(-length / 2, width / 2, 0.0);
        rectanglepoints->InsertNextPoint(-length / 2, -width / 2, 0.0);
        return rectanglepoints;
    }

    vtkNew<vtkPoints> drawTriangle(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y) {
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
        rhombuspoints->InsertNextPoint(0, 2 * sidelength, 0); // upper left
        rhombuspoints->InsertNextPoint(sidelength, 0, 0); // upper right
        rhombuspoints->InsertNextPoint(0, -2 * sidelength, 0); // lower right
        rhombuspoints->InsertNextPoint(-sidelength, 0, 0);
        return rhombuspoints;
    }
    /*-------------------------------------------Draw Functions 3D------------------------------------------------------------*/

    vtkNew<vtkPoints> drawSphere(double radius) {
        // Define sphere parameters
        vtkNew<vtkPoints> spherePoints;
        // Create points for sphere vertices
        for (int i = 0; i <= 360; i++) {
            double theta = i * vtkMath::Pi() / 180;
            for (int j = 0; j <= 180; j++) {
                double phi = j * vtkMath::Pi() / 180;
                double x = radius * sin(phi) * cos(theta);
                double y = radius * sin(phi) * sin(theta);
                double z = radius * cos(phi);
                spherePoints->InsertNextPoint(x, y, z);
            }
        }
        return spherePoints;
    }

    vtkNew<vtkPoints> drawEllipsoid(double radiusX, double radiusY, double radiusZ) {

        // Define ellipsoid parameters
        vtkNew<vtkPoints> ellipsoidPoints;
        double cx = 0.0; // Center X
        double cy = 0.0; // Center Y
        double cz = 0.0; // Center Z
        double rx = radiusX; // X-axis radius (W: half - width)  0.2
        double ry = radiusY; // Y-axis radius (H: half - height) 0.2
        double rz = radiusZ; // Z-axis radius (D: half - depth) 0.3
        // Create points for ellipsoid vertices
        const int thetaSteps = 36;
        const int phiSteps = 36;
        const int psiSteps = 72;
        for (int i = 0; i <= thetaSteps; i++) {
            double theta = i * vtkMath::Pi() / thetaSteps;
            for (int j = 0; j <= phiSteps; j++) {
                double phi = j * vtkMath::Pi() / phiSteps;
                for (int k = 0; k <= psiSteps; k++) {
                    double psi = k * vtkMath::Pi() / psiSteps;
                    double x = cx + rx * sin(phi) * cos(theta);
                    double y = cy + ry * sin(phi) * sin(theta);
                    double z = cz + rz * cos(phi);
                    // Rotate the x-y coordinates around the z-axis by psi
                    double x_rot = x * cos(psi) - y * sin(psi);
                    double y_rot = x * sin(psi) + y * cos(psi);
                    ellipsoidPoints->InsertNextPoint(x_rot, y_rot, z);
                }
            }

            return ellipsoidPoints;
        }

    }
    /*vtkNew<vtkPoints> drawcube(double x ,double y , double z ,double length) {*/
    vtkNew<vtkPoints> drawcube(double length) {
        vtkNew<vtkPoints> cubepoints;
        //first face
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        cubepoints->InsertNextPoint(x, y, z);
        cubepoints->InsertNextPoint(x + length, y, z);
        cubepoints->InsertNextPoint(x + length, y + length, z);
        cubepoints->InsertNextPoint(x, y + length, z);
        cubepoints->InsertNextPoint(x, y, z);
        //second face
        cubepoints->InsertNextPoint(x, y, z + length);
        cubepoints->InsertNextPoint(x + length, y, z + length);
        cubepoints->InsertNextPoint(x + length, y, z);
        //starting point
        cubepoints->InsertNextPoint(x, y, z);
        // third face
        cubepoints->InsertNextPoint(x, y + length, z);
        cubepoints->InsertNextPoint(x, y + length, z + length);
        cubepoints->InsertNextPoint(x, y, z + length);
        //fourth face
        cubepoints->InsertNextPoint(x + length, y, z + length);
        cubepoints->InsertNextPoint(x + length, y + length, z + length);
        cubepoints->InsertNextPoint(x, y + length, z + length);
        //fifth
        cubepoints->InsertNextPoint(x + length, y + length, z + length);
        cubepoints->InsertNextPoint(x + length, y + length, z);
        return cubepoints;
    }

    void selectShape(int index, vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        vtkNew<vtkPoints> points;
        vtkNew<vtkActor> Actor;
        /*------------------------------- Select 2D Shapes------------------------------------------------------------------*/
        if (index == 1) {//line
            double x1 = QInputDialog::getDouble(NULL, "Enter first coordinates", "x1 coordinate", 0, -1000, 1000, 2);
            double y1 = QInputDialog::getDouble(NULL, "Enter first coordinates", "y1 coordinate", 0, -1000, 1000, 2);
            double x2 = QInputDialog::getDouble(NULL, "Enter second coordinates", "x2 coordinate", 0, -1000, 1000, 2);
            double y2 = QInputDialog::getDouble(NULL, "Enter second coordinates", "y2 coordinate", 0, -1000, 1000, 2);
            points = drawLine(x1, y1, x2, y2);
            Shapes line;
            line.name = "line";
            line.x1 = x1; line.y1 = y1; line.x2 = x2; line.y2 = y2; line.z1 = 0.0; line.z2 = 0.0;
            Shapesdata.push_back(line);
        }
        else if (index == 2) {//polyline
            double x1 = QInputDialog::getDouble(NULL, "Enter first coordinates", "x1 coordinate", 0, -1000, 1000, 2);
            double y1 = QInputDialog::getDouble(NULL, "Enter first coordinates", "y1 coordinate", 0, -1000, 1000, 2);
            double x2 = QInputDialog::getDouble(NULL, "Enter second coordinates", "x2 coordinate", 0, -1000, 1000, 2);
            double y2 = QInputDialog::getDouble(NULL, "Enter second coordinates", "y2 coordinate", 0, -1000, 1000, 2);
            double x3 = QInputDialog::getDouble(NULL, "Enter third coordinates", "x3 coordinate", 0, -1000, 1000, 2);
            double y3 = QInputDialog::getDouble(NULL, "Enter third coordinates", "y3 coordinate", 0, -1000, 1000, 2);
            double x4 = QInputDialog::getDouble(NULL, "Enter fourth coordinates", "x4 coordinate", 0, -1000, 1000, 2);
            double y4 = QInputDialog::getDouble(NULL, "Enter fourth coordinates", "y4 coordinate", 0, -1000, 1000, 2);
            double x5 = QInputDialog::getDouble(NULL, "Enter fifth coordinates", "x5 coordinate", 0, -1000, 1000, 2);
            double y5 = QInputDialog::getDouble(NULL, "Enter fifth coordinates", "y5 coordinate", 0, -1000, 1000, 2);
            points = drawPolyline(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5);
            Shapes polyline;
            polyline.name = "polyline";
            polyline.x1 = x1; polyline.y1 = y1; polyline.z1 = 0.0; polyline.x2 = x2; polyline.y2 = y2;  polyline.z2 = 0.0; polyline.x3 = x3; polyline.y3 = y3;  polyline.z3 = 0.0;
            polyline.x4 = x4; polyline.y4 = y4;  polyline.z4 = 0.0; polyline.x5 = x5; polyline.y5 = y5;  polyline.z5 = 0.0;
            Shapesdata.push_back(polyline);
        }
        else if (index == 3) {//irregular polygon
            double x1 = QInputDialog::getDouble(NULL, "Enter first coordinates", "x1 coordinate", 0, -1000, 1000, 2);
            double y1 = QInputDialog::getDouble(NULL, "Enter first coordinates", "y1 coordinate", 0, -1000, 1000, 2);
            double x2 = QInputDialog::getDouble(NULL, "Enter second coordinates", "x2 coordinate", 0, -1000, 1000, 2);
            double y2 = QInputDialog::getDouble(NULL, "Enter second coordinates", "y2 coordinate", 0, -1000, 1000, 2);
            double x3 = QInputDialog::getDouble(NULL, "Enter thrid coordinates", "x3 coordinate", 0, -1000, 1000, 2);
            double y3 = QInputDialog::getDouble(NULL, "Enter third coordinates", "y3 coordinate", 0, -1000, 1000, 2);
            double x4 = QInputDialog::getDouble(NULL, "Enter fourth coordinates", "x4 coordinate", 0, -1000, 1000, 2);
            double y4 = QInputDialog::getDouble(NULL, "Enter fourth coordinates", "y4 coordinate", 0, -1000, 1000, 2);
            double x5 = QInputDialog::getDouble(NULL, "Enter fifth coordinates", "x5 coordinate", 0, -1000, 1000, 2);
            double y5 = QInputDialog::getDouble(NULL, "Enter fifth coordinates", "y5 coordinate", 0, -1000, 1000, 2);
            double x6 = QInputDialog::getDouble(NULL, "Enter sixth coordinates", "x6 coordinate", 0, -1000, 1000, 2);
            double y6 = QInputDialog::getDouble(NULL, "Enter sixth coordinates", "y6 coordinate", 0, -1000, 1000, 2);
            Shapes IrregularPolygon;
            IrregularPolygon.name = "irregular polygon";
            IrregularPolygon.x1 = x1; IrregularPolygon.y1 = y1; IrregularPolygon.z1 = 0.0; IrregularPolygon.x2 = x2; IrregularPolygon.y2 = y2;  IrregularPolygon.z2 = 0.0;
            IrregularPolygon.x3 = x3; IrregularPolygon.y3 = y3;  IrregularPolygon.z3 = 0.0; IrregularPolygon.x4 = x4; IrregularPolygon.y4 = y4;  IrregularPolygon.z4 = 0.0;
            IrregularPolygon.x5 = x5; IrregularPolygon.y5 = y5;  IrregularPolygon.z5 = 0.0; IrregularPolygon.x6 = x6; IrregularPolygon.y6 = y6;  IrregularPolygon.z6 = 0.0;
            Shapesdata.push_back(IrregularPolygon);
        }
        else if (index == 4) {//regular polygon
            int number_of_sides = QInputDialog::getDouble(NULL, "Enter number of sides", "number of sides", 0, -1000, 1000, 2);
            double Polygon_raduis = QInputDialog::getDouble(NULL, "Enter polygon raduis", "polygon raduis", 0, -1000, 1000, 2);
            points = drawRegularPolygon(number_of_sides, Polygon_raduis);
            Shapes regularPolygon;
            regularPolygon.name = "regular polygon";
            regularPolygon.numofsides = number_of_sides;
            regularPolygon.radius1 = Polygon_raduis;
            Shapesdata.push_back(regularPolygon);
        }
        else if (index == 5) {//circle
            double radius = QInputDialog::getDouble(NULL, "Enter radius ", "radius", 0, -1000, 1000, 2);
            points = drawCircle(radius);
            Shapes circle;
            circle.name = "circle";
            circle.radius1 = radius;
            Shapesdata.push_back(circle);
        }
        else if (index == 6) {//arc
            //double centerX = QInputDialog::getDouble(NULL, "Enter point1", "point1", 0, -1000, 1000, 2);//useless center
            //double centerY = QInputDialog::getDouble(NULL, "Enter point2", "point2", 0, -1000, 1000, 2);
            double radius = QInputDialog::getDouble(NULL, "Enter radius", "radius", 0, -1000, 1000, 2);
            double startangle = QInputDialog::getDouble(NULL, "Enter start angle", "start angle", 0, -1000, 1000, 2);
            double endAngle = QInputDialog::getDouble(NULL, "Enter end angle", "end angle", 0, -1000, 1000, 2);
            points = drawArc(radius, startangle, endAngle);
            Shapes arc;
            arc.name = "arc";
            arc.radius1 = radius;
            arc.startangle = startangle;
            arc.endangle = endAngle;
            Shapesdata.push_back(arc);
        }
        else if (index == 7) {//ellipse
            //double centerX = QInputDialog::getDouble(NULL, "Enter point1", "point1", 0, -1000, 1000, 2);//remove centers
            //double centerY = QInputDialog::getDouble(NULL, "Enter point2", "point2", 0, -1000, 1000, 2);
            double radiusA = QInputDialog::getDouble(NULL, "Enter radiusI", "radiusI", 0, -1000, 1000, 2);
            double radiusB = QInputDialog::getDouble(NULL, "Enter radiusII", "radiusII", 0, -1000, 1000, 2);
            points = drawEllipse(radiusA, radiusB);
            Shapes ellipse;
            ellipse.name = "ellipse";
            ellipse.radius1 = radiusA;
            ellipse.radius2 = radiusB;
            Shapesdata.push_back(ellipse);
        }
        else if (index == 8) {//rectangle
            double length = QInputDialog::getDouble(NULL, "Enter length", "x1 coordinate", 0, -1000, 1000, 2);
            double width = QInputDialog::getDouble(NULL, "Enter width", "y1 coordinate", 0, -1000, 1000, 2);
            points = drawRectangle(length, width);
            Shapes rectangle;
            rectangle.name = "rectangle";
            rectangle.length = length;
            rectangle.width = width;
            Shapesdata.push_back(rectangle);
        }
        else if (index == 9) {//triangle
            double p1x = QInputDialog::getDouble(NULL, "Enter point1", "x1 coordinate", 0, -1000, 1000, 2);
            double p1y = QInputDialog::getDouble(NULL, "Enter point1", "y1 coordinate", 0, -1000, 1000, 2);
            double p2x = QInputDialog::getDouble(NULL, "Enter point2", "x2 coordinate", 0, -1000, 1000, 2);
            double p2y = QInputDialog::getDouble(NULL, "Enter point2", "y2 coordinate", 0, -1000, 1000, 2);
            double p3x = QInputDialog::getDouble(NULL, "Enter point3", "x3 coordinate", 0, -1000, 1000, 2);
            double p3y = QInputDialog::getDouble(NULL, "Enter point3", "y3 coordinate", 0, -1000, 1000, 2);
            points = drawTriangle(p1x, p1y, p2x, p2y, p3x, p3y);
            Shapes triangle;
            triangle.name = "triangle";
            triangle.x1 = p1x; triangle.y1 = p1y; triangle.z1 = 0.0; triangle.x2 = p2x; triangle.y2 = p2y; triangle.z2 = 0.0; triangle.x3 = p3x; triangle.y3 = p3y; triangle.z3 = 0.0;
            Shapesdata.push_back(triangle);
        }
        else if (index == 10) {//rhombus
            double sidelength = QInputDialog::getDouble(NULL, "Enter side length", "side length", 0, -1000, 1000, 2);
            points = drawRhombus(sidelength);
            Shapes rhombus;
            rhombus.name = "rhombus";
            rhombus.length = sidelength;
            Shapesdata.push_back(rhombus);
        }
        else if (index == 11) {//star
            double radius = QInputDialog::getDouble(NULL, "Enter radius", "radius", 0, -1000, 1000, 2);
            points = drawStar(radius);
            Shapes star;
            star.name = "star";
            star.radius1 = radius;
            Shapesdata.push_back(star);
        }
        /*-------------------------------------select 3D Shapes------------------------------------------------------------------*/
        else if (index == 12) {//sphere
            double radius = QInputDialog::getDouble(NULL, "Enter sphere radius", "radius", 0, -1000, 1000, 2);
            points = drawSphere(radius);
            Shapes sphere;
            sphere.name = "sphere";
            sphere.radius1 = radius;
            Shapesdata.push_back(sphere);
        }
        else if (index == 13) {//cube
            /*double px = QInputDialog::getDouble(NULL, "Enter your x origin ", "x coordinate", 0, -1000, 1000, 2);
            double py = QInputDialog::getDouble(NULL, "Enter your y origin ", "y coordinate", 0, -1000, 1000, 2);
            double pz = QInputDialog::getDouble(NULL, "Enter your z origin", "z coordinate", 0, -1000, 1000, 2);*/
            double cube_length = QInputDialog::getDouble(NULL, "Enter the length ", "length ", 0, -1000, 1000, 2);
            //points = drawcube(px, py, pz, cube_length);
            points = drawcube(cube_length);
            Shapes cube;
            cube.name = "cube";
            cube.length = cube_length;
            Shapesdata.push_back(cube);
        }
        else if (index == 14) {//ellipsoid
            double radiusX = QInputDialog::getDouble(NULL, "Enter Ellipsoid Xradius", "Xradius", 0, -1000, 1000, 2);
            double radiusY = QInputDialog::getDouble(NULL, "Enter Ellipsoid Yradius", "Yradius", 0, -1000, 1000, 2);
            double radiusZ = QInputDialog::getDouble(NULL, "Enter Ellipsoid Zradius", "Zradius", 0, -1000, 1000, 2);
            points = drawEllipsoid(radiusX, radiusY, radiusZ);
            Shapes ellipsoid;
            ellipsoid.name = "ellipsoid";
            ellipsoid.radius1 = radiusX;
            ellipsoid.radius2 = radiusY;
            ellipsoid.radius3 = radiusZ;
            Shapesdata.push_back(ellipsoid);
        }
        vtkNew<vtkLineSource> linesource;
        linesource->SetPoints(points);
        vtkNew<vtkPolyDataMapper> mapper;
        //mapper takes data that is going to be rendered
        mapper->SetInputConnection(linesource->GetOutputPort());
        //actor is used to change properties
        Actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        Actor->SetVisibility(true);
        Actor->SetMapper(mapper);
        renderer->AddActor(Actor);
        window->Render();
    }
    /*-----------------------------delete selected shape function-----------------------------------------*/

    void deleteSelectedShape(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        QMessageBox msgBox;
        msgBox.setText("Delete:");
        msgBox.setStandardButtons(QMessageBox::Close);
        msgBox.setDefaultButton(QMessageBox::Close);
        QAbstractButton* deleteLast = msgBox.addButton("Delete last", QMessageBox::ButtonRole::AcceptRole);
        QAbstractButton* deleteAll = msgBox.addButton("Clear All", QMessageBox::ButtonRole::AcceptRole);
        msgBox.exec();
        int numActor = renderer->GetActors()->GetNumberOfItems();

        if (msgBox.clickedButton() == deleteLast) {
            if (numActor > 0) {
                vtkProp* lastActor = renderer->GetActors()->GetLastProp();
                renderer->RemoveActor(lastActor);
                Shapesdata.pop_back();
            }
        }
        else {
            for (int i = 0; i < numActor; i++)
            {
                vtkPropCollection* actors = renderer->GetActors();
                actors->InitTraversal();
                vtkProp* prop;
                while ((prop = actors->GetNextProp()) != nullptr) {
                    // Check if the prop is an instance of vtkActor
                    if (vtkActor* actor = vtkActor::SafeDownCast(prop)) {
                        renderer->RemoveActor(actor);
                        Shapesdata.pop_back();
                    }
                }
            }
        }
        window->Render();
    }

    /*---------------------------------Transformation functions 2D&3D--------------------------------------------------------*/
    void applyTranslation(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        double x = QInputDialog::getDouble(NULL, "Enter x translation", "x translation", 0, -1000, 1000, 3);
        double y = QInputDialog::getDouble(NULL, "Enter y translation", "y translation", 0, -1000, 1000, 3);
        double z = QInputDialog::getDouble(NULL, "Enter z translation", "z translation", 0, -1000, 1000, 3);
        vtkProp* prop = renderer->GetActors()->GetLastProp();
        vtkActor* actor = renderer->GetActors()->GetLastActor();
        vtkMapper* mapper = actor->GetMapper();
        vtkPolyData* data = vtkPolyData::SafeDownCast(mapper->GetInput());
        vtkPoints* points = data->GetPoints();
        vtkIdType numPoints = points->GetNumberOfPoints();
        vtkMatrix4x4* translation = vtkMatrix4x4::New();
        translation->Identity();
        translation->SetElement(0, 3, x);
        translation->SetElement(1, 3, y);
        translation->SetElement(2, 3, z);
        //translate each point Q = pq
        //Q -> new place
        //p-> old place
        //q-> translation
        for (vtkIdType i = 0; i < numPoints; i++) {
            double* pt = points->GetPoint(i);
            double p[4] = { pt[0], pt[1], pt[2], 1 };
            double q[4] = { 0, 0, 0, 0 };
            translation->MultiplyPoint(p, q);
            pt[0] = q[0];
            pt[1] = q[1];
            pt[2] = q[2];
            points->SetPoint(i, pt);
        }
        points->Modified();
        window->Render();
        translation->Delete();
    }

    void applyRotation(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        Shapes& lastStruct = Shapesdata[Shapesdata.size() - 1];
        string lastElement = lastStruct.name;
        vtkProp* prop = renderer->GetActors()->GetLastProp();
        vtkActor* actor = renderer->GetActors()->GetLastActor();
        vtkMapper* mapper = actor->GetMapper();
        vtkPolyData* data = vtkPolyData::SafeDownCast(mapper->GetInput());
        vtkPoints* points = data->GetPoints();
        vtkIdType numPoints = points->GetNumberOfPoints();
        double angle = QInputDialog::getDouble(NULL, "Enter angle in degrees", "Rotation angle", 0, -360, 360, 3);
        double angleInRadians = vtkMath::RadiansFromDegrees(angle);
        if (lastElement != "sphere" && lastElement != "ellipsoid" && lastElement != "cube") {
            for (vtkIdType i = 0; i < numPoints; i++) {
                double* pt = points->GetPoint(i);
                double x = pt[0];
                double y = pt[1];
                //Qx= Pxcos(θ) − Pysin(θ)
                //Qy = Pxsin(θ) + Pycos(θ)
                //Q -> new place
                //P-> old place
                double transformedX = x * cos(angleInRadians) - y * sin(angleInRadians);
                double transformedY = x * sin(angleInRadians) + y * cos(angleInRadians);

                pt[0] = transformedX;
                pt[1] = transformedY;

                points->SetPoint(i, pt);
            }
        }
        else {
            vtkMatrix4x4* rotationMatrix2D = vtkMatrix4x4::New();
            QMessageBox msgBox;
            msgBox.setText("choose direction of rotation in 3D");
            msgBox.setStandardButtons(QMessageBox::Close);
            msgBox.setDefaultButton(QMessageBox::Close);
            QAbstractButton* xAxisButton = msgBox.addButton("xaxis", QMessageBox::ButtonRole::AcceptRole);
            QAbstractButton* yAxisButton = msgBox.addButton("yaxis", QMessageBox::ButtonRole::AcceptRole);
            QAbstractButton* zAxisButton = msgBox.addButton("zaxis", QMessageBox::ButtonRole::AcceptRole);
            msgBox.exec();
            vtkMatrix4x4* rotationMatrix3D = vtkMatrix4x4::New();
            rotationMatrix3D->Identity();
            if (msgBox.clickedButton() == xAxisButton) {
                rotationMatrix3D->SetElement(1, 1, cos(angleInRadians));
                rotationMatrix3D->SetElement(2, 1, sin(angleInRadians));
                rotationMatrix3D->SetElement(1, 2, -sin(angleInRadians));
                rotationMatrix3D->SetElement(2, 2, cos(angleInRadians));
            }
            else if (msgBox.clickedButton() == yAxisButton) {
                rotationMatrix3D->SetElement(0, 0, cos(angleInRadians));
                rotationMatrix3D->SetElement(0, 2, sin(angleInRadians));
                rotationMatrix3D->SetElement(2, 0, -sin(angleInRadians));
                rotationMatrix3D->SetElement(2, 2, cos(angleInRadians));
            }
            else {
                rotationMatrix3D->SetElement(0, 0, cos(angleInRadians));
                rotationMatrix3D->SetElement(0, 1, -sin(angleInRadians));
                rotationMatrix3D->SetElement(1, 0, sin(angleInRadians));
                rotationMatrix3D->SetElement(1, 1, cos(angleInRadians));
            }
            for (vtkIdType i = 0; i < numPoints; i++) {
                double* pt = points->GetPoint(i);
                double p[4] = { pt[0], pt[1], pt[2],1 };
                double q[4] = { 0, 0, 0, 0 };//store result of array multiplication in q
                rotationMatrix3D->MultiplyPoint(p, q);//multiply matrix of rotation by the matrix of points p and store it in q
                pt[0] = q[0];
                pt[1] = q[1];
                pt[2] = q[2];
                points->SetPoint(i, pt);
            }
        }
        points->Modified();
        window->Render();
    }

    void applyShear(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        Shapes& lastStruct = Shapesdata[Shapesdata.size() - 1];
        string lastElement = lastStruct.name;
        vtkProp* prop = renderer->GetActors()->GetLastProp();
        vtkActor* actor = renderer->GetActors()->GetLastActor();
        vtkMapper* mapper = actor->GetMapper();
        vtkPolyData* data = vtkPolyData::SafeDownCast(mapper->GetInput());
        vtkPoints* points = data->GetPoints();
        vtkIdType numPoints = points->GetNumberOfPoints();
        if (lastElement != "sphere" && lastElement != "ellipsoid" && lastElement != "cube") {
            QMessageBox msgBox;
            msgBox.setText("choose direction of shearing in 2D");
            msgBox.setStandardButtons(QMessageBox::Close);
            msgBox.setDefaultButton(QMessageBox::Close);
            QAbstractButton* xAxisButton = msgBox.addButton("xaxis", QMessageBox::ButtonRole::AcceptRole);
            QAbstractButton* yAxisButton = msgBox.addButton("yaxis", QMessageBox::ButtonRole::AcceptRole);
            msgBox.exec();
            vtkMatrix3x3* shearMatrix2D = vtkMatrix3x3::New();
            shearMatrix2D->Identity();
            double shearFactorX, shearFactorY;
            if (msgBox.clickedButton() == xAxisButton) {
                shearFactorX = QInputDialog::getDouble(NULL, "Enter shear factor X", "shear factor X", 0, -1000, 1000, 3);
            }
            if (msgBox.clickedButton() == yAxisButton) {
                shearFactorY = QInputDialog::getDouble(NULL, "Enter shear factor  Y", "shear factor Y", 0, -1000, 1000, 3);
            }
            for (vtkIdType i = 0; i < numPoints; i++) {
                double* pt = points->GetPoint(i);
                double x = pt[0];
                double y = pt[1];
                double transformedX = x + shearFactorX * y;
                double transformedY = shearFactorY * x + y;
                pt[0] = transformedX;
                pt[1] = transformedY;
                points->SetPoint(i, pt);
            }
        }
        else {
            QMessageBox msgBox;
            msgBox.setText("choose direction of shearing in 3D");
            msgBox.setStandardButtons(QMessageBox::Close);
            msgBox.setDefaultButton(QMessageBox::Close);
            QAbstractButton* xAxisButton = msgBox.addButton("xaxis", QMessageBox::ButtonRole::AcceptRole);
            QAbstractButton* yAxisButton = msgBox.addButton("yaxis", QMessageBox::ButtonRole::AcceptRole);
            QAbstractButton* zAxisButton = msgBox.addButton("zaxis", QMessageBox::ButtonRole::AcceptRole);
            msgBox.exec();
            vtkMatrix4x4* shearMatrix3D = vtkMatrix4x4::New();
            shearMatrix3D->Identity();
            if (msgBox.clickedButton() == xAxisButton) {
                double shearFactorY = QInputDialog::getDouble(NULL, "Enter shear factor Y", "shear factor Y", 0, -1000, 1000, 3);
                double shearFactorZ = QInputDialog::getDouble(NULL, "Enter shear factor Z", "shear factor Z", 0, -1000, 1000, 3);
                shearMatrix3D->SetElement(0, 1, shearFactorY);
                shearMatrix3D->SetElement(0, 2, shearFactorZ);
            }
            else if (msgBox.clickedButton() == yAxisButton) {
                double shearFactorX = QInputDialog::getDouble(NULL, "Enter shear factor X", "shear factor X", 0, -1000, 1000, 3);
                double shearFactorZ = QInputDialog::getDouble(NULL, "Enter shear factor Z", "shear factor Z", 0, -1000, 1000, 3);
                shearMatrix3D->SetElement(1, 0, shearFactorX);
                shearMatrix3D->SetElement(1, 2, shearFactorZ);
            }
            else {
                double shearFactorX = QInputDialog::getDouble(NULL, "Enter shear factor X", "shear factor X", 0, -1000, 1000, 3);
                double shearFactorY = QInputDialog::getDouble(NULL, "Enter shear factor Y", "shear factor Y", 0, -1000, 1000, 3);
                shearMatrix3D->SetElement(2, 0, shearFactorX);
                shearMatrix3D->SetElement(2, 1, shearFactorY);
            }
            for (vtkIdType i = 0; i < numPoints; i++) {
                double* pt = points->GetPoint(i);
                double p[4] = { pt[0], pt[1], pt[2],1 };
                double q[4] = { 0, 0, 0, 0 };//store result of array multiplication in q
                shearMatrix3D->MultiplyPoint(p, q);//multiply matrix of rotation by the matrix of points p and store it in q
                pt[0] = q[0];
                pt[1] = q[1];
                pt[2] = q[2];
                points->SetPoint(i, pt);
            }
        }
        points->Modified();
        window->Render();
    }

    void applyScaling(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        Shapes& lastStruct = Shapesdata[Shapesdata.size() - 1];
        string lastElement = lastStruct.name;
        vtkProp* prop = renderer->GetActors()->GetLastProp();
        vtkActor* actor = renderer->GetActors()->GetLastActor();
        vtkMapper* mapper = actor->GetMapper();
        vtkPolyData* data = vtkPolyData::SafeDownCast(mapper->GetInput());
        vtkPoints* points = data->GetPoints();
        vtkIdType numPoints = points->GetNumberOfPoints();
        if (lastElement != "sphere" && lastElement != "ellipsoid" && lastElement != "cube") {
            double x = QInputDialog::getDouble(NULL, "enter x scaling factor", "x scaling factor", 0, -1000, 1000, 3);
            double y = QInputDialog::getDouble(NULL, "enter y scaling factor", "y scaling factor", 0, -1000, 1000, 3);
            for (vtkIdType i = 0; i < numPoints; i++) {
                double* pt = points->GetPoint(i);
                double xpoint = pt[0];
                double ypoint = pt[1];
                double transformedX = xpoint * x;
                double transformedY = ypoint * y;
                pt[0] = transformedX;
                pt[1] = transformedY;
                points->SetPoint(i, pt);
            }
        }
        else {
            double x = QInputDialog::getDouble(NULL, "enter x scaling factor", "x scaling factor", 0, -1000, 1000, 3);
            double y = QInputDialog::getDouble(NULL, "enter y scaling factor", "y scaling factor", 0, -1000, 1000, 3);
            double z = QInputDialog::getDouble(NULL, "enter z scaling factor", "z scaling factor", 0, -1000, 1000, 3);
            //(Qx, Qy) = (SxPx, SyPy), multiply each old point p by the scaling factor
            vtkMatrix4x4* scaling = vtkMatrix4x4::New();
            scaling->Identity();
            scaling->SetElement(0, 0, x);
            scaling->SetElement(1, 1, y);
            scaling->SetElement(2, 2, z);
            for (vtkIdType i = 0; i < numPoints; i++) {
                double* pt = points->GetPoint(i);
                double p[4] = { pt[0], pt[1], pt[2], 1 };
                double q[4] = { 0, 0, 0, 0 };
                scaling->MultiplyPoint(p, q);
                pt[0] = q[0];
                pt[1] = q[1];
                pt[2] = q[2];
                points->SetPoint(i, pt);
            }
        }
        points->Modified();
        window->Render();
    }
    /*---------------------------------change color function-----------------------------------------------------------*/
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
    /*----------------------change line width function--------------------------------------------------------------*/
    void changeLineWidth(int value, vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        vtkActor* lastActor = vtkActor::SafeDownCast(renderer->GetActors()->GetLastProp());
        lastActor->GetProperty()->SetLineWidth(value);
        window->Render();
    }
    /*------------------------Read & Write Functions---------------------------------------------------------------*/
    void writeInFile(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        QString fileName = QFileDialog::getSaveFileName(nullptr, "Save File", ".", "CSV Files (*.csv)");
        QFile file(fileName);
        vtkSmartPointer<vtkActorCollection> actors = renderer->GetActors();
        vtkSmartPointer<vtkCollectionIterator> it = actors->NewIterator();
        int counter = 0;
        if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            QTextStream out(&file);
            out << "Shape name" << "," << "X1" << "," << "Y1" << "," << "Z1" << "," << "X2" << "," << "Y2" << "," << "Z2" << "," << "X3" << "," << "Y3" << "," << "Z3" << ","
                << "X4" << "," << "Y4" << "," << "Z4" << "," << "X5" << "," << "Y5" << "," << "Z5" << "," << "X6" << "," << "Y6" << "," << "Z6" << "," << "length"
                << "," << "width" << "," << "Radius1" << "," << "Radius2" << "," << "Radius3" << "," << "Start angle" << "," << "End angle" << "," << "No. of sides" << ","
                << "R" << "," << "G" << "," << "B" << "," << "Line width" << "," << Qt::endl;

            for (const auto& shape : Shapesdata) {
                out << QString::fromStdString(shape.name) << ",";
                out << (shape.x1.has_value() ? QString::fromStdString(std::to_string(*shape.x1)) : "none") << ",";
                out << (shape.y1.has_value() ? QString::fromStdString(std::to_string(*shape.y1)) : "none") << ",";
                out << (shape.z1.has_value() ? QString::fromStdString(std::to_string(*shape.z1)) : "none") << ",";
                out << (shape.x2.has_value() ? QString::fromStdString(std::to_string(*shape.x2)) : "none") << ",";
                out << (shape.y2.has_value() ? QString::fromStdString(std::to_string(*shape.y2)) : "none") << ",";
                out << (shape.z2.has_value() ? QString::fromStdString(std::to_string(*shape.z2)) : "none") << ",";
                out << (shape.x3.has_value() ? QString::fromStdString(std::to_string(*shape.x3)) : "none") << ",";
                out << (shape.y3.has_value() ? QString::fromStdString(std::to_string(*shape.y3)) : "none") << ",";
                out << (shape.z3.has_value() ? QString::fromStdString(std::to_string(*shape.z3)) : "none") << ",";
                out << (shape.x4.has_value() ? QString::fromStdString(std::to_string(*shape.x4)) : "none") << ",";
                out << (shape.y4.has_value() ? QString::fromStdString(std::to_string(*shape.y4)) : "none") << ",";
                out << (shape.z4.has_value() ? QString::fromStdString(std::to_string(*shape.z4)) : "none") << ",";
                out << (shape.x5.has_value() ? QString::fromStdString(std::to_string(*shape.x5)) : "none") << ",";
                out << (shape.y5.has_value() ? QString::fromStdString(std::to_string(*shape.y5)) : "none") << ",";
                out << (shape.z5.has_value() ? QString::fromStdString(std::to_string(*shape.z5)) : "none") << ",";
                out << (shape.x6.has_value() ? QString::fromStdString(std::to_string(*shape.x6)) : "none") << ",";
                out << (shape.y6.has_value() ? QString::fromStdString(std::to_string(*shape.y6)) : "none") << ",";
                out << (shape.z6.has_value() ? QString::fromStdString(std::to_string(*shape.z6)) : "none") << ",";
                out << (shape.length.has_value() ? QString::fromStdString(std::to_string(*shape.length)) : "none") << ",";
                out << (shape.width.has_value() ? QString::fromStdString(std::to_string(*shape.width)) : "none") << ",";
                out << (shape.radius1.has_value() ? QString::fromStdString(std::to_string(*shape.radius1)) : "none") << ",";
                out << (shape.radius2.has_value() ? QString::fromStdString(std::to_string(*shape.radius2)) : "none") << ",";
                out << (shape.radius3.has_value() ? QString::fromStdString(std::to_string(*shape.radius3)) : "none") << ",";
                out << (shape.startangle.has_value() ? QString::fromStdString(std::to_string(*shape.startangle)) : "none") << ",";
                out << (shape.endangle.has_value() ? QString::fromStdString(std::to_string(*shape.endangle)) : "none") << ",";
                out << (shape.numofsides.has_value() ? QString::fromStdString(std::to_string(*shape.numofsides)) : "none") << ",";
                vtkActor* actor = vtkActor::SafeDownCast(it->GetCurrentObject());
                if (actor) {
                    double* color = actor->GetProperty()->GetColor();

                    out << color[0] << ","
                        << color[1] << ","
                        << color[2] << ","
                        << actor->GetProperty()->GetLineWidth() << Qt::endl;
                }
                it->GoToNextItem();
            }
        }
        file.close();
    }

    void readInputFile(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        QString fileObject = QFileDialog::getOpenFileName(nullptr, "Open File", ".", "CSV Files (*.csv)");
        if (fileObject.isEmpty()) {
            return;  // Dialog was cancelled
        }
        QFile file(fileObject);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return;
        else {
            QTextStream in(&file);
            QString line = in.readLine();
            while (!in.atEnd()) {
                line = in.readLine();
                QStringList elements = line.split(",");
                vtkNew<vtkPoints> points;
                vtkNew<vtkActor> actor;
                bool ok;
                if (elements[0] == "line") {
                    double x1 = elements[1].toDouble(&ok);
                    double y1 = elements[2].toDouble(&ok);
                    double x2 = elements[4].toDouble(&ok);
                    double y2 = elements[5].toDouble(&ok);
                    points = drawLine(x1, y1, x2, y2);
                }
                else if (elements[0] == "polyline") {
                    double x1 = elements[1].toDouble(&ok);
                    double y1 = elements[2].toDouble(&ok);

                    double x2 = elements[4].toDouble(&ok);
                    double y2 = elements[5].toDouble(&ok);

                    double x3 = elements[7].toDouble(&ok);
                    double y3 = elements[8].toDouble(&ok);

                    double x4 = elements[10].toDouble(&ok);
                    double y4 = elements[11].toDouble(&ok);

                    double x5 = elements[13].toDouble(&ok);
                    double y5 = elements[14].toDouble(&ok);
                    points = drawPolyline(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5);
                }
                else if (elements[0] == "irregular polygon") {
                    double x1 = elements[1].toDouble(&ok);
                    double y1 = elements[2].toDouble(&ok);

                    double x2 = elements[4].toDouble(&ok);
                    double y2 = elements[5].toDouble(&ok);

                    double x3 = elements[7].toDouble(&ok);
                    double y3 = elements[8].toDouble(&ok);

                    double x4 = elements[10].toDouble(&ok);
                    double y4 = elements[11].toDouble(&ok);

                    double x5 = elements[13].toDouble(&ok);
                    double y5 = elements[14].toDouble(&ok);

                    double x6 = elements[16].toDouble(&ok);
                    double y6 = elements[17].toDouble(&ok);
                    points = drawIrregularPolygon(x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6);
                }
                else if (elements[0] == "regular polygon") {
                    double radius = elements[21].toDouble(&ok);
                    int numofsides = elements[25].toDouble(&ok);
                    points = drawRegularPolygon(numofsides, radius);
                }
                else if (elements[0] == "arc") {
                    double radius = elements[21].toDouble(&ok);
                    int startangle = elements[24].toDouble(&ok);
                    int endangle = elements[25].toDouble(&ok);
                    points = drawArc(radius, startangle, endangle);
                }
                else if (elements[0] == "ellipse") {
                    double radius1 = elements[21].toDouble(&ok);
                    double radius2 = elements[22].toDouble(&ok);
                    points = drawEllipse(radius1, radius2);
                }
                else if (elements[0] == "circle") {
                    double radius = elements[21].toDouble(&ok);
                    points = drawCircle(radius);
                }
                else if (elements[0] == "rectangle") {
                    double length = elements[19].toDouble(&ok);
                    double width = elements[20].toDouble(&ok);
                    points = drawRectangle(length, width);
                }
                else if (elements[0] == "star") {
                    double radius = elements[21].toDouble(&ok);
                    points = drawStar(radius);
                }
                else if (elements[0] == "cube") {
                    double length = elements[19].toDouble(&ok);
                    points = drawcube(length);
                }
                else if (elements[0] == "sphere") {
                    double radius = elements[21].toDouble(&ok);
                    points = drawSphere(radius);
                }
                else if (elements[0] == "rhombus") {
                    double length = elements[19].toDouble(&ok);
                    points = drawRhombus(length);
                }
                else if (elements[0] == "triangle") {
                    double x1 = elements[1].toDouble(&ok);
                    double y1 = elements[2].toDouble(&ok);

                    double x2 = elements[4].toDouble(&ok);
                    double y2 = elements[5].toDouble(&ok);

                    double x3 = elements[7].toDouble(&ok);
                    double y3 = elements[8].toDouble(&ok);
                    points = drawTriangle(x1, y1, x2, y2, x3, y3);
                }
                else if (elements[0] == "ellipsoid") {
                    double radius1 = elements[21].toDouble(&ok);
                    double radius2 = elements[22].toDouble(&ok);
                    double radius3 = elements[23].toDouble(&ok);
                    points = drawEllipsoid(radius1,radius2, radius3);
                }
                vtkNew<vtkLineSource> linesource;
                linesource->SetPoints(points);
                vtkNew<vtkPolyDataMapper> mapper;
                //mapper takes data that is going to be rendered
                mapper->SetInputConnection(linesource->GetOutputPort());
                //actor is used to change properties
                actor->GetProperty()->SetColor(elements[27].toDouble(&ok), elements[28].toDouble(&ok), elements[29].toDouble(&ok));
                actor->GetProperty()->SetLineWidth(elements[30].toDouble(&ok));
                actor->SetVisibility(true);
                actor->SetMapper(mapper);
                renderer->AddActor(actor);
                window->Render();
            }
            file.close(); //close the file object.
        }
    }

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
    shapesComboBox.addItem(QApplication::tr("Rectangle/Square"));
    shapesComboBox.addItem(QApplication::tr("Triangle"));
    shapesComboBox.addItem(QApplication::tr("Rhombus"));
    shapesComboBox.addItem(QApplication::tr("Star"));
    shapesComboBox.addItem(QApplication::tr("Sphere"));
    shapesComboBox.addItem(QApplication::tr("Cube"));
    shapesComboBox.addItem(QApplication::tr("Ellipsoid"));
    shapesComboBox.setItemData(0, QVariant(0), Qt::UserRole - 1);

    dockLayout->addWidget(&shapesComboBox, 1, Qt::AlignTop);

    QPushButton deleteButton("Delete Shape");
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

    window->GetInteractor()->SetInteractorStyle(style);

    QObject::connect(&shapesComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), [&](int index) {
        ::selectShape(index, window, renderer);
        });

    // connect button to slot
    QObject::connect(&deleteButton, &QPushButton::clicked, [&]() {
        ::deleteSelectedShape(window, renderer);
        });
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
        ::changeLineWidth(value, window, renderer);
        });
    QObject::connect(&readFile, &QPushButton::released,
        [&]() { ::readInputFile(window, renderer); });

    QObject::connect(&writeFile, &QPushButton::released,
        [&]() { ::writeInFile(window, renderer); });

    window->Render();
    mainWindow.show();

    return app.exec();
    //return EXIT_SUCCESS;
}
