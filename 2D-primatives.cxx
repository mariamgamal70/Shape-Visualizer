﻿#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkNew.h>
//#include <vtkInteractorStyleTrackballCamera.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderWindow.h>
//#include <vtkNamedColors.h>
#include <vtkLineSource.h>
#include <vtkPointPicker.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
//#include <vtkTextActor.h>
//#include <vtkTextProperty.h>
//#include <vtkTextRepresentation.h>
//#include <vtkTextWidget.h>
//#include <vtkParametricEllipsoid.h>
//#include <vtkParametricFunctionSource.h>
//#include <vtkSuperquadricSource.h>
//#include <vtkRegularPolygonSource.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>

#include <QApplication>
#include <QDockWidget>
#include <QGridLayout>
#include <QLabel>
#include <QMainWindow>
#include <QPointer>
#include <QPushButton>
#include <QVBoxLayout>
#include <QInputDialog>
//#include <QFileDialog>
#include <QComboBox>

#include <cmath>
#include <cstdlib>
//#include <random>
#include <iostream>
#include <fstream>
using namespace std;

ofstream myfile;
ifstream readfile;
int filecounter = 0;

namespace {

    void selectShape(int index, vtkGenericOpenGLRenderWindow* window, vtkActor* ellipseActor, vtkActor* regularPolygonActor, vtkActor* starActor, vtk*arcActor, vtk*circleActor, vtkActor* PolylineActor,vtkActor* irregularPolygonActor) {
        if (index == 0) {//line
            cout << "line";
        }
        else if (index == 1) {//polyline
            regularPolygonActor->SetVisibility(false);
            ellipseActor->SetVisibility(false);
            circleActor->SetVisibility(false);
            arcActor->SetVisibility(false);
            PolylineActor->SetVisibility(true);
            irregularPolygonActor->SetVisibility(false);
            starActor->SetVisibility(false);
        }
        else if (index == 2) {//polygon
            regularPolygonActor->SetVisibility(false);
            ellipseActor->SetVisibility(false);
            circleActor->SetVisibility(false);
            arcActor->SetVisibility(false);
            PolylineActor->SetVisibility(false);
            irregularPolygonActor->SetVisibility(true);
            starActor->SetVisibility(false);
        }
        else if (index == 3) {//regular polygon
           regularPolygonActor->SetVisibility(true);
            ellipseActor->SetVisibility(false);
            circleActor->SetVisibility(false);
            arcActor->SetVisibility(false);
            PolylineActor->SetVisibility(false);
            irregularPolygonActor->SetVisibility(false);
            starActor->SetVisibility(false);
        }        
        else if (index == 4) {//circle
            circleActor->SetVisibility(true);
            ellipseActor->SetVisibility(false);
            regularPolygonActor->SetVisibility(false);
            arcActor->SetVisibility(false);
            irregularPolygonActor->SetVisibility(false);
            polylineActor->SetVisibility(false);
            starActor->SetVisibility(false);
        }        
        else if (index == 5) {//arc
            arcActor->SetVisibility(true);
            ellipseActor->SetVisibility(false);
            regularPolygonActor->SetVisibility(false);
            circleActor->SetVisibility(false);
            irregularPolygonActor->SetVisibility(false);
            polylineActor->SetVisibility(false);
            starActor->SetVisibility(false);
        }        
        else if (index == 6) {//ellipse
            ellipseActor->SetVisibility(true);
            regularPolygonActor->SetVisibility(false);
            circleActor->SetVisibility(false);
            arcActor->SetVisibility(false);
           irregularPolygonActor->SetVisibility(false);
            polylineActor->SetVisibility(false);
             starActor->SetVisibility(false);
        }
        else if (index == 7) {//rectangle
            ellipseActor->SetVisibility(true);
            regularPolygonActor->SetVisibility(false);
        }
        else if (index == 8) {//triangle
            ellipseActor->SetVisibility(true);
            regularPolygonActor->SetVisibility(false);
        }
        else if (index == 9) {//rhombus
            ellipseActor->SetVisibility(true);
            regularPolygonActor->SetVisibility(false);
        }
        else if (index == 10) {//star
            starActor->SetVisibility(true);
            ellipseActor->SetVisibility(false);
            regularPolygonActor->SetVisibility(false);
            circleActor->SetVisibility(false);
            arcActor->SetVisibility(false);
           irregularPolygonActor->SetVisibility(false);
            polylineActor->SetVisibility(false);
        }
        window->Render();
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

    QComboBox shapesComboBox ;
    shapesComboBox.addItem(QApplication::tr("Line"));
    shapesComboBox.addItem(QApplication::tr("Polyline"));
    shapesComboBox.addItem(QApplication::tr("Polygon"));
    shapesComboBox.addItem(QApplication::tr("Regular Polygon"));
    shapesComboBox.addItem(QApplication::tr("Circle"));
    shapesComboBox.addItem(QApplication::tr("Arc"));
    shapesComboBox.addItem(QApplication::tr("Ellipse"));
    shapesComboBox.addItem(QApplication::tr("Rectangle"));
    shapesComboBox.addItem(QApplication::tr("Triangle"));
    shapesComboBox.addItem(QApplication::tr("Rhombus"));
    shapesComboBox.addItem(QApplication::tr("Star"));
    dockLayout->addWidget(&shapesComboBox, 1, Qt::AlignTop);

     //render area
    QPointer<QVTKOpenGLNativeWidget> vtkRenderWidget = new QVTKOpenGLNativeWidget();
    mainWindow.setCentralWidget(vtkRenderWidget);
    mainWindow.setWindowTitle("VTK Line Example");

    // VTK part
    vtkNew<vtkGenericOpenGLRenderWindow> window;
    vtkRenderWidget->setRenderWindow(window);

    /*vtkNew<vtkLineSource> linesource;
    vtkNew<vtkPolyDataMapper> linemapper;
    linemapper->SetInputConnection(linesource->GetOutputPort());
    vtkNew<vtkActor> lineactor;
    lineactor->SetMapper(linemapper);*/
    /*----------------------ELLIPSE---------------------*/

    // Define ellipse parameters
    double cx = 0.0; // Center X
    double cy = 0.0; // Center Y
    double rx = 0.3; // X-axis radius (W: half - width)
    double ry = 0.1; // Y-axis radius (H: half - height)
    // Create points for ellipse vertices
    vtkNew<vtkPoints> ellipsepoints;
    for (int i = 0; i <= 360; i++) {
        //convert from degree to radian
        //for 0 ≤ t ≤ 2pi.
        double theta = i* vtkMath::Pi()/180; 
        //x(t) = W cos(t) , add cx to shift to the stated center
        double x = cx + rx * cos(theta);
        //y(t) = H sin(t)  , add cy to shift to the stated center
        double y = cy + ry * sin(theta);
        ellipsepoints->InsertNextPoint(x, y, 0.0); // z=0 because its 2D
    }
    // Create polyline to connect those points using lines
    //vtkpolyline: type of VTK cell that represents a single polyline in 3D space.
    vtkNew<vtkPolyLine> ellipsepolyline;
    //ses the number of point IDs in the vtkIdList associated with the polyline object
    ellipsepolyline->GetPointIds()->SetNumberOfIds(ellipsepoints->GetNumberOfPoints());
    // iterate through each point in the points object and set the corresponding point ID in the vtkIdList associated with the polyline object. 
    for (vtkIdType i = 0; i < ellipsepoints->GetNumberOfPoints(); ++i) {
        ellipsepolyline->GetPointIds()->SetId(i, i);
    }
    //vtkPolyData:VTK data object that represents a dataset consisting of points, cells, and associated data attributes.
    vtkNew<vtkPolyData> ellipsepolydata;
    //set the points object as the points of the polydata object
    ellipsepolydata->SetPoints(ellipsepoints);
    //allocate memory for the cells in the polydata object. 
    ellipsepolydata->Allocate();
    //The polyline(parameter1) will be drawn using lines connecting the points(parameter2) defined by the point IDs in the vtkIdList.
    ellipsepolydata->InsertNextCell(ellipsepolyline->GetCellType(), ellipsepolyline->GetPointIds());
    vtkNew<vtkPolyDataMapper> ellipsemapper;
    //mapper takes data that is going to be rendered
    ellipsemapper->SetInputData(ellipsepolydata);
    //actor is used to change properties
    vtkNew<vtkActor> ellipseActor;
    ellipseActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    ellipseActor->SetVisibility(false);
    ellipseActor->SetMapper(ellipsemapper);
    
    /*---------------------regular polygon----------------------*/

    // Create points for regular polygon vertices
    vtkNew<vtkPoints> regularPolygonPoints;
    //Pointi = ( R cos( 2πi / n ), R sin(2πi / n )),
    double regularPolygonNoOfSides = 6;
    double regularPolygonRadius = 0.2;
    double regularPolygonCenter = 0.0;
    for (int i = 0; i <= regularPolygonNoOfSides; i++) {
        //double theta = i * vtkMath::Pi() / 180;
        //x= R cos( 2πi / n )
        double x = regularPolygonRadius * cos(2 * vtkMath::Pi() * i / regularPolygonNoOfSides);
        //y= R sin(2πi / n )
        double y = regularPolygonRadius * sin(2 * vtkMath::Pi() * i / regularPolygonNoOfSides);
        regularPolygonPoints->InsertNextPoint(x, y, 0.0); // z=0 because its 2D
    }
    // Create polyline to connect those points using lines
    //vtkpolyline: type of VTK cell that represents a single polyline in 3D space.
    vtkNew<vtkPolyLine> regularPolygonPolyline;
    //ses the number of point IDs in the vtkIdList associated with the polyline object
    regularPolygonPolyline->GetPointIds()->SetNumberOfIds(regularPolygonPoints->GetNumberOfPoints());
    // iterate through each point in the points object and set the corresponding point ID in the vtkIdList associated with the polyline object. 
    for (vtkIdType i = 0; i < regularPolygonPoints->GetNumberOfPoints(); ++i) {
        regularPolygonPolyline->GetPointIds()->SetId(i, i);
    }
    //vtkPolyData:VTK data object that represents a dataset consisting of points, cells, and associated data attributes.
    vtkNew<vtkPolyData> regularPolygonPolydata;
    //set the points object as the points of the polydata object
    regularPolygonPolydata->SetPoints(regularPolygonPoints);
    //allocate memory for the cells in the polydata object. 
    regularPolygonPolydata->Allocate();
    //The polyline(parameter1) will be drawn using lines connecting the points(parameter2) defined by the point IDs in the vtkIdList.
    regularPolygonPolydata->InsertNextCell(regularPolygonPolyline->GetCellType(), regularPolygonPolyline->GetPointIds());
    vtkNew<vtkPolyDataMapper> regularPolygonMapper;
    //mapper takes data that is going to be rendered
    regularPolygonMapper->SetInputData(regularPolygonPolydata);
    //actor is used to change properties
    vtkNew<vtkActor> regularPolygonActor;
    regularPolygonActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    regularPolygonActor->SetVisibility(false);
    regularPolygonActor->SetMapper(regularPolygonMapper);

    /*---------------------star----------------------*/
    // Create points for regular polygon vertices
    vtkNew<vtkPoints> starPoints;
    //rosette general equation
    //x(t) = r * cos(k * t * theta)
    //y(t) = r * sin(k * t * theta)
    //r is the radius of the rosette, controlling the size of the curve.
    //k is the number of petals or points on the rosette, controlling the number of "petals" or "points" in the curve.
    //t is the parameter that varies from 0 to 2*pi, controlling the position of the points on the curve.
    //theta is an additional parameter that can be adjusted to control the shape of the rosette. It is typically a constant value.
    double starRadius = 0.2;   
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
    // Create polyline to connect those points using lines
    //vtkpolyline: type of VTK cell that represents a single polyline in 3D space.
    vtkNew<vtkPolyLine> starPolyline;
    //ses the number of point IDs in the vtkIdList associated with the polyline object
    starPolyline->GetPointIds()->SetNumberOfIds(starPoints->GetNumberOfPoints());
    // iterate through each point in the points object and set the corresponding point ID in the vtkIdList associated with the polyline object. 
    for (vtkIdType i = 0; i < starPoints->GetNumberOfPoints(); ++i) {
        starPolyline->GetPointIds()->SetId(i, i);
    }
    //vtkPolyData:VTK data object that represents a dataset consisting of points, cells, and associated data attributes.
    vtkNew<vtkPolyData> starPolydata;
    //set the points object as the points of the polydata object
    starPolydata->SetPoints(starPoints);
    //allocate memory for the cells in the polydata object. 
    starPolydata->Allocate();
    //The polyline(parameter1) will be drawn using lines connecting the points(parameter2) defined by the point IDs in the vtkIdList.
    starPolydata->InsertNextCell(starPolyline->GetCellType(), starPolyline->GetPointIds());
    vtkNew<vtkPolyDataMapper> starMapper;
    //mapper takes data that is going to be rendered
    starMapper->SetInputData(starPolydata);
    //actor is used to change properties
    vtkNew<vtkActor> starActor;
    starActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    starActor->SetVisibility(false);
    starActor->SetMapper(starMapper);

    /*----------------------CIRCLE---------------------*/
    vtkNew<vtkPoints> circlepoints;
    vtkNew<vtkCellArray> circlelines;

    double angle_step = 2.0 * vtkMath::Pi() / 100.0;
    for (int i = 0; i < 100; i++) {
        double x = cos(i * angle_step);
        double y = sin(i * angle_step);
        double z = 0.0;
        circlepoints->InsertNextPoint(x, y, z);
        vtkIdType ids[2] = { i, (i + 1) % 100 };
        circlelines->InsertNextCell(2, ids);
    }

    vtkNew<vtkPolyData> circlePolyData;
    circlePolyData->SetPoints(circlepoints);
    circlePolyData->SetLines(circlelines);

    vtkNew<vtkPolyDataMapper> circleMapper;
    circleMapper->SetInputData(circlePolyData);
   

    vtkNew<vtkActor> circleActor;
    circleActor->GetProperty()->SetColor(0, 0, 1);
    circleActor->SetVisibility(false);
    circleActor->SetMapper(circleMapper);
     /*----------------------ARC---------------------*/

    vtkNew<vtkArcSource> arcSource;
    arcSource->SetCenter(0, 0, 0);
    arcSource->SetPoint1(0.25, 0, 0);
    arcSource->SetPoint2(0, 0.25, 0);
    arcSource->SetAngle(90);
    arcSource->SetResolution(50);
    arcSource->Update();
    vtkNew<vtkPolyDataMapper> arcMapper;
    arcMapper->SetInputConnection(arcSource->GetOutputPort());
    vtkNew<vtkActor> arcActor;
    arcActor->GetProperty()->SetColor(1, 0, 0);
    arcActor->SetVisibility(false);
    arcActor->SetMapper(arcMapper);
    /*---------------------irregular polygone--------------------*/
    vtkNew<vtkPoints> irregularPolygonPoints;

    // Generate random coordinates for each vertex
    srand(time(NULL));
    double maxX = 0.5;
    double maxY = 0.5;
  
    irregularPolygonPoints->InsertNextPoint(0.0, 0.5, 0);
    irregularPolygonPoints->InsertNextPoint(0.5, 0.2, 0);
    irregularPolygonPoints->InsertNextPoint(0.3, -0.2, 0);
    irregularPolygonPoints->InsertNextPoint(0.0, 0.0, 0);
    irregularPolygonPoints->InsertNextPoint(-0.3, -0.3, 0);
    irregularPolygonPoints->InsertNextPoint(-0.5, 0.3, 0);
   /* for (int i = 0; i < 6; i++) {
        double x = (double)rand() / RAND_MAX * maxX - maxX / 2;
        double y = (double)rand() / RAND_MAX * maxY - maxY / 2;
        irregularPolygonPoints->InsertNextPoint(x, y, 0.0);
    }*/
    
    irregularPolygonPoints->InsertNextPoint(0.0, 0.5, 0);

    // Create polyline to connect those points using lines
    vtkNew<vtkPolyLine> irregularPolygonPolyline;
    irregularPolygonPolyline->GetPointIds()->SetNumberOfIds(irregularPolygonPoints->GetNumberOfPoints());
    for (vtkIdType i = 0; i < irregularPolygonPoints->GetNumberOfPoints(); ++i) {
        irregularPolygonPolyline->GetPointIds()->SetId(i, i);
    }

    // Create polydata to represent the irregular polygon
    vtkNew<vtkPolyData> irregularPolygonPolydata;
    irregularPolygonPolydata->SetPoints(irregularPolygonPoints);
    irregularPolygonPolydata->Allocate();
    irregularPolygonPolydata->InsertNextCell(irregularPolygonPolyline->GetCellType(), irregularPolygonPolyline->GetPointIds());

    // Create mapper to map the polydata to graphics primitives
    vtkNew<vtkPolyDataMapper> irregularPolygonMapper;
    irregularPolygonMapper->SetInputData(irregularPolygonPolydata);

    // Create actor to represent the irregular polygon in the scene
    vtkNew<vtkActor> irregularPolygonActor;
    irregularPolygonActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    irregularPolygonActor->SetVisibility(false);
    irregularPolygonActor->SetMapper(irregularPolygonMapper);
 /*----------------------polyline--------------------*/
    

    // Create five points.
    double origin[3] = { -0.2, 0.5, 0.0 };
    double p0[3] = { 0.0, 0.0, 0.0 };
    double p1[3] = { 0.1, 0.0, 0.0 };
    double p2[3] = { 0.2, -0.4, 0.0 };
    double p3[3] = { 0.3, 0.0, 0.0 };

    // Create a vtkPoints object and store the points in it
    vtkNew<vtkPoints> points;
    points->InsertNextPoint(origin);
    points->InsertNextPoint(p0);
    points->InsertNextPoint(p1);
    points->InsertNextPoint(p2);
    points->InsertNextPoint(p3);

    vtkNew<vtkPolyLine> polyLine;
    polyLine->GetPointIds()->SetNumberOfIds(5);
    for (unsigned int i = 0; i < 5; i++)
    {
        polyLine->GetPointIds()->SetId(i, i);
    }

    // Create a cell array to store the lines in and add the lines to it
    vtkNew<vtkCellArray> cells;
    cells->InsertNextCell(polyLine);

    // Create a polydata to store everything in
    vtkNew<vtkPolyData> polyData;

    // Add the points to the dataset
    polyData->SetPoints(points);

    // Add the lines to the dataset
    polyData->SetLines(cells);

    // Setup actor and mapper
    vtkNew<vtkPolyDataMapper> PolylineMapper;
    PolylineMapper->SetInputData(polyData);

    vtkNew<vtkActor> PolylineActor;
    PolylineActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    PolylineActor->SetVisibility(false);
    PolylineActor->SetMapper(PolylineMapper);


    /*---------------------renderers---------------------*/

    vtkNew<vtkRenderer> renderer;
    //renderer->AddActor(lineactor);
    renderer->AddActor(ellipseActor);
    renderer->AddActor(regularPolygonActor);
    renderer->AddActor(starActor);
    renderer->AddActor(circleActor);
    renderer->AddActor(arcActor);
    renderer->AddActor(irregularPolygonActor);
    renderer->AddActor(PolylineActor);

    window->AddRenderer(renderer);
    
    //vtkNew<vtkPointPicker> pointPicker;
    window->SetInteractor(vtkRenderWidget->interactor());
    //window->GetInteractor()->SetPicker(pointPicker);

    //window->SetInteractor(vtkRenderWidget->interactor());

       QObject::connect(&shapesComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), [&](int index) {
          ::selectShape(index,window,ellipseActor,regularPolygonActor,starActor,circleActor, arcActor ,PolylineActor, irregularPolygonActor );
          });
    window->Render();
    mainWindow.show();

    return app.exec();
    //return EXIT_SUCCESS;
}
