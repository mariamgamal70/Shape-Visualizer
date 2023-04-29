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
#include <QStandardItem>

#include <cmath>
#include <cstdlib>
#include <iostream>
using namespace std;

namespace {
    //// Define interaction style
    //class customMouseInteractorStyle : public vtkInteractorStyleTrackballCamera
    //{
    //public:
    //    static customMouseInteractorStyle* New();
    //    vtkTypeMacro(customMouseInteractorStyle, vtkInteractorStyleTrackballCamera);

    //    virtual void OnLeftButtonDown() override
    //    {
    //        //get the x and y coordinates of the mouse click.
    //        int x = this->Interactor->GetEventPosition()[0];
    //        int y = this->Interactor->GetEventPosition()[1];
    //        //get the renderer that was clicked on by the mouse.
    //        vtkRenderer* renderer = this->Interactor->FindPokedRenderer(x, y);
    //        //get the actor displayed from the renderer
    //        setSelectedActor(renderer);
    //        //pick the object that was clicked on by the mouse
    //        this->Interactor->GetPicker()->Pick(x, y, 0, renderer);
    //        double pickedPoint[3];
    //        //retrieves the 3D position where the mouse was clicked in the rendering window 
    //        this->Interactor->GetPicker()->GetPickPosition(pickedPoint);
    //        //check if an actor was selected
    //        if (SelectedActor)
    //        {
    //            SelectedActor->SetDragable(true);
    //            SelectedActor->SetPickable(false);
    //            // Create a transform and set the shearing coefficients
    //            vtkNew<vtkTransform> transform;
    //            vtkNew < vtkMatrix4x4> matrix;
    //            matrix->Identity();
    //            matrix->SetElement(0, 1, 0.5);

    //            transform->SetMatrix(matrix);
    //            // Apply the transform to the actor
    //            SelectedActor->SetUserTransform(transform);
    //            // saves the current position of the actor before it is moved.
    //            LastPosition[0] = pickedPoint[0];
    //            LastPosition[1] = pickedPoint[1];
    //            LastPosition[2] = pickedPoint[2];
    //            // sets the new position of the actor to the picked point, which will cause it to move to that location.
    //            SelectedActor->SetPosition(pickedPoint);
    //        }
    //        // Forward events
    //        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    //    }
    //    void setSelectedActor(vtkRenderer* renderer) {
    //    // safely cast the vtkProp object returned by GetLastProp() to an vtkActor object
    //        SelectedActor =renderer->GetActors()->GetLastActor();

    //           //vtkActor::SafeDownCast(this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->GetActors()->GetLastProp());
    //    }
    //    void setSelectedPolyData() {
    //        polyData = vtkPolyData::SafeDownCast(SelectedActor->GetMapper()->GetInputDataObject(0, 0));
    //    }

    //private:
    //    double LastPosition[3];
    //    vtkActor* SelectedActor;
    //    vtkPolyData* polyData;
    //};
    //vtkStandardNewMacro(customMouseInteractorStyle);
    
    
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

    vtkNew<vtkPoints> drawLine() {
        vtkNew<vtkPoints> linepoints;
        linepoints->InsertNextPoint(-0.5, -0.5, 0.0);
        linepoints->InsertNextPoint(0.5, 0.5, 0.0);
        return linepoints;
    }

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
        double regularPolygonRadius = 0.5;
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
        double starRadius = 0.5;
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

    vtkNew<vtkPoints>drawCircle() {
        vtkNew<vtkPoints> circlepoints;

        double angle_step = 2.0 * vtkMath::Pi() / 100.0;
        for (int i = 0; i <= 100; i++) {
            double x = cos(i * angle_step);
            double y = sin(i * angle_step);
            double z = 0.0;
            circlepoints->InsertNextPoint(x, y, z);
        }
            return circlepoints;
    }

    vtkNew<vtkPoints>drawArc() {
        vtkNew<vtkPoints> arcpoints;
        double center[3] = { -0.1 , -0.1 , 0.0 };
        double radius = 0.5;
        double startAngle = 0.0;
        double endAngle = 90.0;
        int numSegments = 20;

        // Add the arc points
        for (int i = 0; i <= numSegments; i++)
        {
            double angle = startAngle + (i / static_cast<double>(numSegments)) * (endAngle - startAngle);
            double x = center[0] + radius * cos(angle * vtkMath::Pi() / 180.0);
            double y = center[1] + radius * sin(angle * vtkMath::Pi() / 180.0);
            double z = center[2];
            arcpoints->InsertNextPoint(x, y, z);
        }
        return arcpoints;
    }

    vtkNew<vtkPoints> drawRectangle() {
        vtkNew<vtkPoints> rectanglepoints;
        rectanglepoints->InsertNextPoint(-0.2, -0.1, 0.0);
        rectanglepoints->InsertNextPoint(0.2, -0.1, 0.0);
        rectanglepoints->InsertNextPoint(0.2, 0.1, 0.0);
        rectanglepoints->InsertNextPoint(-0.2, 0.1, 0.0);
        rectanglepoints->InsertNextPoint(-0.2, -0.1, 0.0);
        return rectanglepoints;
    }

    vtkNew<vtkPoints> drawTriangle() {
        vtkNew<vtkPoints> trianglepoints;
        trianglepoints->InsertNextPoint(-0.2, -0.1, 0.0);
        trianglepoints->InsertNextPoint(0.2, -0.1, 0.0);
        trianglepoints->InsertNextPoint(0.0, 0.1, 0.0);
        trianglepoints->InsertNextPoint(-0.2, -0.1, 0.0);
        //trianglepoints->InsertNextPoint(-0.5, 0.25, 0.0);
        //trianglepoints->InsertNextPoint(-0.5, -0.25, 0.0);
        return trianglepoints;
    }

    vtkNew<vtkPoints> drawRhombus() {
        vtkNew<vtkPoints> rhombuspoints;

        double a = 0.25; // side length
        rhombuspoints->InsertNextPoint(-a, 0, 0); // lower left
        rhombuspoints->InsertNextPoint(0, 2*a, 0); // upper left
        rhombuspoints->InsertNextPoint(a, 0, 0); // upper right
        rhombuspoints->InsertNextPoint(0, -2*a, 0); // lower right
        rhombuspoints->InsertNextPoint(-a, 0, 0);
        return rhombuspoints;
    }

    void selectShape(int index, vtkGenericOpenGLRenderWindow* window,vtkRenderer* renderer) {
        vtkNew<vtkPoints> points;
        int numActors = renderer->GetActors()->GetNumberOfItems();
        if (numActors > 0) {
            vtkProp* lastActor = renderer->GetActors()->GetLastProp();
            renderer->RemoveActor(lastActor);
        }
        if (index == 1) {//line
            points = drawLine();
        }
        else if (index == 2) {//polyline
            points = drawPolyline();
        }
        else if (index == 3) {//irregular polygon
            points = drawIrregularPolygon();
        }
        else if (index == 4) {//regular polygon
            points = drawRegularPolygon();
        }        
        else if (index == 5) {//circle
            points = drawCircle();
        }        
        else if (index == 6) {//arc
            points = drawArc();
        }        
        else if (index == 7) {//ellipse
            points = drawEllipse();
        }
        else if (index == 8) {//rectangle
            points = drawRectangle();
        }
        else if (index == 9) {//triangle
            points = drawTriangle();
        }
        else if (index == 10) {//rhombus
            points = drawRhombus();
        }
        else if (index == 11) {//star
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

    void deleteSelectedShape(vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        int numActors = renderer->GetActors()->GetNumberOfItems();
        if (numActors > 0) {
            vtkProp* lastActor = renderer->GetActors()->GetLastProp();
            renderer->RemoveActor(lastActor);
            window->Render();
        }
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
        for (vtkIdType i = 0; i < numPoints; i++) {
            double* pt = points->GetPoint(i);
            pt[0] *= scalingfactor;
            pt[1] *= scalingfactor;            
            points->SetPoint(i, pt);
        }
        points->Modified();
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

    QPushButton translateButton("Apply Translation");
    dockLayout->addWidget(&translateButton, 1, Qt::AlignTop);
    
    QPushButton rotateButton("Apply rotating");
    dockLayout->addWidget(&rotateButton, 1, Qt::AlignTop);

    QPushButton scaleButton("Apply rotating");
    dockLayout->addWidget(&scaleButton, 1, Qt::AlignTop);

    QPushButton shearButton("Apply shearing");
    dockLayout->addWidget(&shearButton, 1, Qt::AlignTop);

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

    vtkNew<ScalingInteractorStyle> style;

    //vtkNew<customMouseInteractorStyle> style;
    //style->setLineSource(linesource);
    //style->setVTKActor(lineactor);

    window->GetInteractor()->SetInteractorStyle(style);

       QObject::connect(&shapesComboBox, QOverload<int>::of(&QComboBox::currentIndexChanged), [&](int index) {
          ::selectShape(index,window,renderer);
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
    window->Render();
    mainWindow.show();

    return app.exec();
    //return EXIT_SUCCESS;
}
