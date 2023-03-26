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

    //void addEllipse(vtkSuperquadricSource* source, vtkActor* actor,vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
    //    /*double width = QInputDialog::getDouble(NULL, "Enter ellipse scaling factors", "width", 0, -1000, 2, 2);
    //    double depth = QInputDialog::getDouble(NULL, "Enter ellipse scaling factors", "depth", 0, -1000, 2, 2);
    //    double thickness = QInputDialog::getDouble(NULL, "Enter ellipse scaling factors", "thickness", 0, -1000, 2, 2);
    //    double centerx = QInputDialog::getDouble(NULL, "Enter ellipse center", "center x", 0, -1000, 2, 2);
    //    double centery = QInputDialog::getDouble(NULL, "Enter ellipse center", "center y", 0, -1000, 2, 2);
    //    double centerz = QInputDialog::getDouble(NULL, "Enter ellipse center", "center z", 0, -1000, 2, 2);
    //    
    //    source->SetPhiRoundness(1);
    //    source->SetThetaRoundness(1);
    //    source->SetScale(width, depth, thickness);
    //    source->Update();

    //    actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    //    actor->SetPosition(centerx, centery, centerz);
    //    renderer->AddActor(actor);
    //    window->Render();*/

    //    double angle = 0;
    //    double r1, r2;
    //    double centerX, centerY;
    //    r1 = 50;
    //    r2 = 30;
    //    centerX = 10.0;
    //    centerY = 5.0;
    //    vtkNew<vtkPoints> points;
    //    int id = 0;
    //    while (angle <= 2.0 * vtkMath::Pi() + (vtkMath::Pi() / 60.0))
    //    {
    //        points->InsertNextPoint(r1 * cos(angle) + centerX,
    //            r2 * sin(angle) + centerY, 0.0);
    //        angle = angle + (vtkMath::Pi() / 60.0);
    //        ++id;
    //    }
    //    vtkNew<vtkPolyLine> line;
    //    line->GetPointIds()->SetNumberOfIds(id);
    //    for (unsigned int i = 0; i < static_cast<unsigned int>(id); ++i)
    //    {
    //        line->GetPointIds()->SetId(i, i);
    //    }

    //    vtkNew<vtkCellArray> lines;
    //    lines->InsertNextCell(line);

    //    vtkNew<vtkPolyData> polyData;
    //    polyData->SetPoints(points);
    //    polyData->SetLines(lines);
    //}
    //void addRegularPolygon(vtkRegularPolygonSource* polygonSource, vtkActor* polyactor,vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
    //    int numsides = QInputDialog::getInt(NULL, "Enter regular polygon info", "number of sides", 0, 0, 10, 2);
    //    double centerx = QInputDialog::getDouble(NULL, "Enter regular polygon center", "center x", 0, -1000, 1000, 2);
    //    double centery = QInputDialog::getDouble(NULL, "Enter regular polygon center", "center y", 0, -1000, 1000, 2);
    //    double centerz = QInputDialog::getDouble(NULL, "Enter regular polygon center", "center z", 0, -1000, 1000, 2);
    //    double radius= QInputDialog::getDouble(NULL, "Enter regular polygon radius", "radius", 0, -1000, 1000, 2);
    //    //vtkNew<vtkRegularPolygonSource> polygonSource;
    //    polygonSource->SetNumberOfSides(numsides);
    //    double center[3] = { centerx, centery,centerz };
    //    polygonSource->SetCenter(center);
    //    polygonSource->SetRadius(radius);          // Horizontal radius
    //    /*vtkNew<vtkPolyDataMapper> polymapper;
    //    polymapper->SetInputConnection(polygonSource->GetOutputPort());
    //    vtkNew <vtkActor> polyactor;
    //    polyactor->SetMapper(polymapper);*/
    //    polyactor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    //    renderer->AddActor(polyactor);
    //    window->Render();
    //}
    void selectShape(int index, vtkGenericOpenGLRenderWindow* window, vtkActor* ellipseActor, vtkActor* regularPolygonActor) {
        if (index == 0) {
            cout << "line";
        }
        else if (index == 1) {//ellipse
            ellipseActor->SetVisibility(true);
            regularPolygonActor->SetVisibility(false);
        }
        else if (index == 2) {//regular polygon
            regularPolygonActor->SetVisibility(true);
            ellipseActor->SetVisibility(false);
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

    /*QPushButton addEllipse;
    addEllipse.setText("Add Ellipse");
    dockLayout->addWidget(&addEllipse, 0, Qt::AlignTop);

    QPushButton addRegularPolygon;
    addRegularPolygon.setText("Add Regular Polygon");
    dockLayout->addWidget(&addRegularPolygon, 0, Qt::AlignTop);*/

    QComboBox shapesComboBox ;
    shapesComboBox.addItem(QApplication::tr("Line"));
    shapesComboBox.addItem(QApplication::tr("Ellipse"));
    shapesComboBox.addItem(QApplication::tr("Regular Polygon"));
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
    double angle = 0;
    double r1, r2;
    double centerX, centerY;
    r1 = 0.25;
    r2 = 0.1;
    centerX = 0;
    centerY = 0;
    vtkNew<vtkPoints> points;
    int id = 0;
    while (angle <= 2.0 * vtkMath::Pi() + (vtkMath::Pi() / 60.0))
    {
        points->InsertNextPoint(r1 * cos(angle) + centerX,
            r2 * sin(angle) + centerY, 0.0);
        angle = angle + (vtkMath::Pi() / 60.0);
        ++id;
    }
    vtkNew<vtkPolyLine> line;
    line->GetPointIds()->SetNumberOfIds(id);
    for (unsigned int i = 0; i < static_cast<unsigned int>(id); ++i)
    {
        line->GetPointIds()->SetId(i, i);
    }

    vtkNew<vtkCellArray> lines;
    lines->InsertNextCell(line);

    vtkNew<vtkPolyData> polyData;
    polyData->SetPoints(points);
    polyData->SetLines(lines);

    vtkNew<vtkPolyDataMapper> ellipsemapper;
    ellipsemapper->SetInputData(polyData);
    vtkNew<vtkActor> ellipseActor;

    ellipseActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    ellipseActor->SetVisibility(false);

    ellipseActor->SetMapper(ellipsemapper);
    
    /*---------------------regular polygon----------------------*/
    vtkNew<vtkRegularPolygonSource> regularPolygonSource;
    regularPolygonSource->SetNumberOfSides(5);
    regularPolygonSource->SetCenter(0,0,0);
    regularPolygonSource->SetRadius(0.2);          // Horizontal radius
    vtkNew<vtkPolyDataMapper> regularPolygonMapper;
    regularPolygonMapper->SetInputConnection(regularPolygonSource->GetOutputPort());
    vtkNew <vtkActor> regularPolygonActor;
    regularPolygonActor->GetProperty()->SetColor(1.0, 0.0, 0.0);
    regularPolygonActor->SetVisibility(false);
    regularPolygonActor->SetMapper(regularPolygonMapper);

    vtkNew<vtkRenderer> renderer;
    //renderer->AddActor(lineactor);
    renderer->AddActor(ellipseActor);
    renderer->AddActor(regularPolygonActor);
    window->AddRenderer(renderer);
    
    vtkNew<vtkPointPicker> pointPicker;
    window->SetInteractor(vtkRenderWidget->interactor());
    window->GetInteractor()->SetPicker(pointPicker);

    window->SetInteractor(vtkRenderWidget->interactor());
    //vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    //vtkNew<vtkPointPicker> pointPicker;
    //window->GetInteractor()->SetPicker(pointPicker);
    //renderWindowInteractor->SetPicker(pointPicker);
    //interactor->SetInteractorStyle(style);

    //// connect the buttons
  /*  QObject::connect(&addEllipse, &QPushButton::released,
        [&]() { ::addEllipse(ellipsesource, ellipseactor,window,renderer); });*/

       QObject::connect(&shapesComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), [&](int index) {
          ::selectShape(index,window,ellipseActor,regularPolygonActor);
          });
    window->Render();
    mainWindow.show();

    return app.exec();
    //return EXIT_SUCCESS;
}
