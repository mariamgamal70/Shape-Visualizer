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
#include <vtkRegularPolygonSource.h>
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
#include <QFileDialog>
#include <QComboBox>

#include <cmath>
#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
using namespace std;

ofstream myfile;
ifstream readfile;
int filecounter = 0;

namespace {

    void selectShape(int index, vtkGenericOpenGLRenderWindow* window, vtkActor* ellipseActor, vtkActor* regularPolygonActor) {
        if (index == 0) {//line
            cout << "line";
        }
        else if (index == 1) {//polyline

        }
        else if (index == 2) {//polygon

        }
        else if (index == 3) {//regular polygon
            regularPolygonActor->SetVisibility(true);
            ellipseActor->SetVisibility(false);
        }        
        else if (index == 4) {//circle

        }        
        else if (index == 5) {//arc

        }        
        else if (index == 6) {//ellipse
            ellipseActor->SetVisibility(true);
            regularPolygonActor->SetVisibility(false);
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
    double n = 6;
    double r = 0.2;
    double c = 0.0;
    for (int i = 0; i <= n; i++) {
        //double theta = i * vtkMath::Pi() / 180;
        //x= R cos( 2πi / n )
        double x = r * cos(2 * vtkMath::Pi() * i / n);
        //y= R sin(2πi / n )
        double y = r * sin(2 * vtkMath::Pi() * i / n);
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

    /*---------------------renderers---------------------*/

    vtkNew<vtkRenderer> renderer;
    //renderer->AddActor(lineactor);
    renderer->AddActor(ellipseActor);
    renderer->AddActor(regularPolygonActor);
    window->AddRenderer(renderer);
    
    //vtkNew<vtkPointPicker> pointPicker;
    window->SetInteractor(vtkRenderWidget->interactor());
    //window->GetInteractor()->SetPicker(pointPicker);

    //window->SetInteractor(vtkRenderWidget->interactor());

       QObject::connect(&shapesComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), [&](int index) {
          ::selectShape(index,window,ellipseActor,regularPolygonActor);
          });
    window->Render();
    mainWindow.show();

    return app.exec();
    //return EXIT_SUCCESS;
}
