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
#include <vtkRenderWindowInteractor.h>
#include <vtkLineSource.h>
#include <vtkPointPicker.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkPoints.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkTransform.h>

#include <QApplication>
#include <QDockWidget>
#include <QGridLayout>
#include <QLabel>
#include <QMainWindow>
#include <QPointer>
#include <QPushButton>
#include <QVBoxLayout>
#include <QInputDialog>
#include <QComboBox>

#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

namespace {
    // Define interaction style
    class customMouseInteractorStyle : public vtkInteractorStyleTrackballCamera
    {
    public:
        static customMouseInteractorStyle* New();
        vtkTypeMacro(customMouseInteractorStyle, vtkInteractorStyleTrackballCamera);

        virtual void OnLeftButtonDown() override
        {
            int x = this->Interactor->GetEventPosition()[0];
            int y = this->Interactor->GetEventPosition()[1];
            cout << x<<" " << y<<endl;
            vtkRenderer* renderer = this->Interactor->FindPokedRenderer(x, y);
            this->Interactor->GetPicker()->Pick(x, y, 0, renderer);
            double pickedPoint[3];
            cout << pickedPoint[0] << " " << pickedPoint[1] << " " << pickedPoint[2] << endl;
            this->Interactor->GetPicker()->GetPickPosition(pickedPoint);
            setSelectedActor();
            if (SelectedActor)
            {
                SelectedActor->SetDragable(true);
                SelectedActor->SetPickable(false);
                /*vtkNew <vtkTransform> transform;
                transform->Translate(pickedPoint);
                SelectedActor->SetUserTransform(transform);*/
                LastPosition[0] = pickedPoint[0];
                LastPosition[1] = pickedPoint[1];
                LastPosition[2] = pickedPoint[2];
                cout << LastPosition[0] << " " << LastPosition[1] << " " << LastPosition[2] << " ";
                SelectedActor->SetPosition(pickedPoint);

            }

            // Forward events
            vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
        }
        virtual void OnLeftButtonUp() override
        {
            if (SelectedActor) {
                // Save the new position of the selected actor
                double newPosition[3];
                SelectedActor->GetPosition(newPosition);
                LastPosition[0] = newPosition[0];
                LastPosition[1] = newPosition[1];
                LastPosition[2] = newPosition[2];

                // Reset the actor's properties
                SelectedActor->GetProperty()->SetColor(1, 1, 1);
                SelectedActor->SetDragable(false);
                SelectedActor->SetPickable(true);
                SelectedActor = nullptr;
            }

            // Forward the event to the parent class
            vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
        }
        void setSelectedActor() {
        // safely cast the vtkProp object returned by GetLastProp() to an vtkActor object
            SelectedActor = vtkActor::SafeDownCast(this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastProp());
        }
        void setSelectedPolyData() {
            polyData = vtkPolyData::SafeDownCast(SelectedActor->GetMapper()->GetInputDataObject(0, 0));
        }

    private:
        double LastPosition[3];
        vtkActor* SelectedActor;
        vtkPolyData* polyData;
    };
    vtkStandardNewMacro(customMouseInteractorStyle);

    vtkNew<vtkPoints> drawEllipse() {
        // Define ellipse parameters
        vtkNew<vtkPoints> ellipsepoints;
        double cx = 0.0; // Center X
        double cy = 0.0; // Center Y
        double rx = 0.3; // X-axis radius (W: half - width)
        double ry = 0.1; // Y-axis radius (H: half - height)
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
    vtkNew<vtkPoints> drawRegularPolygon() {
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
        return regularPolygonPoints;
    }
    vtkNew<vtkPoints> drawStar() {
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
        return starPoints;
    }
    vtkNew<vtkPoints>drawPolyline() {
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
        return points;
    }
    vtkNew<vtkPoints>drawIrregularPolygon() {
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
        return irregularPolygonPoints;
    }
    void selectShape(int index, vtkGenericOpenGLRenderWindow* window,vtkRenderer* renderer) {
        vtkNew<vtkPoints> points;
        int numActors = renderer->GetActors()->GetNumberOfItems();
        if (numActors > 0) {
            vtkProp* lastActor = renderer->GetActors()->GetLastProp();
            renderer->RemoveActor(lastActor);
        }
        if (index == 0) {//line

        }
        else if (index == 1) {//polyline
            points = drawPolyline();
        }
        else if (index == 2) {//irregular polygon
            points = drawIrregularPolygon();
        }
        else if (index == 3) {//regular polygon
            points = drawRegularPolygon();
        }        
        else if (index == 4) {//circle

        }        
        else if (index == 5) {//arc

        }        
        else if (index == 6) {//ellipse
            points = drawEllipse();
        }
        else if (index == 7) {//rectangle

        }
        else if (index == 8) {//triangle

        }
        else if (index == 9) {//rhombus

        }
        else if (index == 10) {//star
            points = drawStar();
        }
        // Create polyline to connect those points using lines
        //vtkpolyline: type of VTK cell that represents a single polyline in 3D space.
        vtkNew<vtkPolyLine> polyline;
        //ses the number of point IDs in the vtkIdList associated with the polyline object
        polyline->GetPointIds()->SetNumberOfIds(points->GetNumberOfPoints());
        // iterate through each point in the points object and set the corresponding point ID in the vtkIdList associated with the polyline object. 
        for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
            polyline->GetPointIds()->SetId(i, i);
        }
        //vtkPolyData:VTK data object that represents a dataset consisting of points, cells, and associated data attributes.
        vtkNew<vtkPolyData> polydata;
        //set the points object as the points of the polydata object
        polydata->SetPoints(points);
        //allocate memory for the cells in the polydata object. 
        polydata->Allocate();
        //The polyline(parameter1) will be drawn using lines connecting the points(parameter2) defined by the point IDs in the vtkIdList.
        polydata->InsertNextCell(polyline->GetCellType(), polyline->GetPointIds());
        vtkNew<vtkPolyDataMapper> mapper;
        //mapper takes data that is going to be rendered
        mapper->SetInputData(polydata);
        //actor is used to change properties
        vtkNew<vtkActor> Actor;
        Actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        Actor->SetVisibility(true);
        Actor->SetMapper(mapper);
        renderer->AddActor(Actor);
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
    shapesComboBox.addItem(QApplication::tr("Irregular Polygon"));
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

    /*---------------------renderers---------------------*/

    vtkNew<vtkRenderer> renderer;
    window->AddRenderer(renderer);
    
    vtkNew<vtkPointPicker> pointPicker;
    window->SetInteractor(vtkRenderWidget->interactor());
    window->GetInteractor()->SetPicker(pointPicker);

    vtkNew<customMouseInteractorStyle> style;
    //style->setLineSource(linesource);
    //style->setVTKActor(lineactor);

    window->GetInteractor()->SetInteractorStyle(style);

       QObject::connect(&shapesComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), [&](int index) {
          ::selectShape(index,window,renderer);
          });
    window->Render();
    mainWindow.show();

    return app.exec();
    //return EXIT_SUCCESS;
}
