#include <QVTKOpenGLNativeWidget.h>
#include <vtkActor.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderWindow.h>
#include <vtkLineSource.h>
#include <vtkNamedColors.h>
#include <vtkPointPicker.h>
#include <vtkRendererCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTextRepresentation.h>
#include <vtkTextWidget.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricFunctionSource.h>
#include <vtkRegularPolygonSource.h>
#include <vtkSuperquadricSource.h>



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
    // Define interaction style
    //class customMouseInteractorStyle : public vtkInteractorStyleTrackballCamera
    //{
    //public:
    //    static customMouseInteractorStyle* New();
    //    vtkTypeMacro(customMouseInteractorStyle, vtkInteractorStyleTrackballCamera);

    //    virtual void OnLeftButtonDown() override
    //    {
    //        click++;
    //        vtkRenderWindowInteractor* interactor = this->Interactor;
    //        //int clickpositionone[2];
    //        //this->Interactor->GetEventPosition(clickpositionone); //get mouse coordinates x and y
    //        if (click == 1) {
    //            this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],//pick the first point using mouse x&y
    //                this->Interactor->GetEventPosition()[1],
    //                0, // always zero.
    //                this->Interactor->GetRenderWindow()
    //                ->GetRenderers()
    //                ->GetFirstRenderer());//The renderer in which the picking operation will be performed.
    //            double pickedone[3];
    //            this->Interactor->GetPicker()->GetPickPosition(pickedone);
    //            UpdateFirstPoint(pickedone);
    //        }
    //        if (click == 2) {
    //            double pickedtwo[3];
    //            this->Interactor->GetPicker()->GetPickPosition(pickedtwo);
    //            vtkRenderWindowInteractor* interactor = this->Interactor;
    //            this->Interactor->GetPicker()->Pick(this->Interactor->GetEventPosition()[0],//pick the first point using mouse x&y
    //                this->Interactor->GetEventPosition()[1],
    //                0, // always zero.
    //                this->Interactor->GetRenderWindow()
    //                ->GetRenderers()
    //                ->GetFirstRenderer());//The renderer in which the picking operation will be performed.
    //            UpdateSecondPoint(pickedtwo);
    //            click = 0;
    //        }
    //        double* point1 = LineSource->GetPoint1();
    //        double* point2 = LineSource->GetPoint2();
    //        char text[100];
    //        sprintf(text, "Line coordinates: (%.2f, %.2f) - (%.2f, %.2f)", point1[0], point1[1], point2[0], point2[1]);
    //        TextActor->SetInput(text);
    //        TextActor->Modified();
    //        //writeInFile();
    //        // Forward events
    //        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    //    }
    //    void setLineSource(vtkLineSource* linesource) {
    //        this->LineSource = linesource;
    //    }
    //    void setTextActor(vtkTextActor* actor) {
    //        TextActor = actor;
    //    }
    //    void UpdateFirstPoint(double* pickedone) {
    //        LineSource->SetPoint1(pickedone[0], pickedone[1], pickedone[2]);
    //    }
    //    void UpdateSecondPoint(double* pickedtwo) {
    //        LineSource->SetPoint2(pickedtwo[0], pickedtwo[1], pickedtwo[2]);
    //    }
    //    void setVTKActor(vtkActor* lineActor) {
    //        this->lineActor = lineActor;
    //    }
    //    //void writeInFile() {
    //    //    myfile.open("myfile.txt", ios::out); // Open the file for writing
    //    //    if (myfile.is_open()) { // Check if file opened successfully
    //    //        double* point1 = LineSource->GetPoint1();
    //    //        double* point2 = LineSource->GetPoint2();
    //    //        myfile << point1[0] << " " << point1[1] << endl; // Write data to the file
    //    //        myfile << point2[0] << " " << point2[1] << endl;
    //    //        //myfile << lineActor->GetProperty()->GetColor() << endl;
    //    //        myfile.close(); // Close the file
    //    //    }
    //    //    else {
    //    //        cout << "Unable to create or open the file." << endl;
    //    //    }
    //    //}

    //private:
    //    vtkLineSource* LineSource;
    //    vtkTextActor* TextActor;
    //    vtkActor* lineActor;
    //    int click = 0;
    //};
    //vtkStandardNewMacro(customMouseInteractorStyle);

  //-------------------------------------------------------------------------------------------------------------------------------------------

    //void writeInFile(vtkLineSource* linesource, vtkActor* lineActor) {
    //    QString fileName = QFileDialog::getSaveFileName(nullptr, "Save File", ".", "Text Files (*.txt)");
    //    QFile file(fileName);
    //    if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
    //        QTextStream out(&file);
    //        double* point1 = linesource->GetPoint1();
    //        double* point2 = linesource->GetPoint2();
    //        out << point1[0] << " " << point1[1] << Qt::endl; // write data to the file
    //        out << point2[0] << " " << point2[1] << Qt::endl;
    //        out << lineActor->GetProperty()->GetColor()[0] << " "
    //            << lineActor->GetProperty()->GetColor()[1] << " "
    //            << lineActor->GetProperty()->GetColor()[2] << Qt::endl;
    //        file.close();
    //    }
    //    //filecounter++;
    //    //string filename = "myfile" + to_string(filecounter) + ".txt";
    //    //myfile.open(filename, ios::out); // open the file for writing

    //    //if (myfile.is_open()) { // check if file opened successfully
    //    //    double* point1 = linesource->GetPoint1();
    //    //    double* point2 = linesource->GetPoint2();
    //    //    myfile << point1[0] << " " << point1[1] << endl; // write data to the file
    //    //    myfile << point2[0] << " " << point2[1] << endl;
    //    //    //myfile << lineactor->getproperty()-> getcolor() << endl;
    //    //    myfile.close(); // close the file
    //    //}
    //    //else {
    //    //    cout << "unable to create or open the file." << endl;
    //    //}
    //}

    void updateTextCoordinates(vtkLineSource* linesource, vtkTextActor* TextActor, vtkActor* lineActor) {
        double* point1 = linesource->GetPoint1();
        double* point2 = linesource->GetPoint2();
        char text[100];
        sprintf(text, "Line coordinates: (%.2f, %.2f) - (%.2f, %.2f)", point1[0], point1[1], point2[0], point2[1]);
        TextActor->SetInput(text);
        TextActor->Modified();
        //writeInFile(linesource,lineActor);
    }
    void addEllipse(vtkSuperquadricSource* source, vtkActor* actor,vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        double width = QInputDialog::getDouble(NULL, "Enter ellipse scaling factors", "width", 0, -1000, 2, 2);
        double depth = QInputDialog::getDouble(NULL, "Enter ellipse scaling factors", "depth", 0, -1000, 2, 2);
        double thickness = QInputDialog::getDouble(NULL, "Enter ellipse scaling factors", "thickness", 0, -1000, 2, 2);
        double centerx = QInputDialog::getDouble(NULL, "Enter ellipse center", "center x", 0, -1000, 2, 2);
        double centery = QInputDialog::getDouble(NULL, "Enter ellipse center", "center y", 0, -1000, 2, 2);
        double centerz = QInputDialog::getDouble(NULL, "Enter ellipse center", "center z", 0, -1000, 2, 2);
        
        //vtkNew<vtkSuperquadricSource> source;
        source->SetPhiRoundness(1);
        source->SetThetaRoundness(0.8);
        source->SetScale(width, depth, thickness);
        source->Update();

        //vtkNew<vtkPolyDataMapper> mapper;
        //mapper->SetInputData(source->GetOutput());

        //vtkNew<vtkActor> actor;
        //actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        actor->SetPosition(centerx, centery, centerz);
        renderer->AddActor(actor);
        window->Render();

    }
    void addRegularPolygon(vtkRegularPolygonSource* polygonSource, vtkActor* polyactor,vtkGenericOpenGLRenderWindow* window, vtkRenderer* renderer) {
        int numsides = QInputDialog::getInt(NULL, "Enter regular polygon info", "number of sides", 0, 0, 10, 2);
        double centerx = QInputDialog::getDouble(NULL, "Enter regular polygon center", "center x", 0, -1000, 1000, 2);
        double centery = QInputDialog::getDouble(NULL, "Enter regular polygon center", "center y", 0, -1000, 1000, 2);
        double centerz = QInputDialog::getDouble(NULL, "Enter regular polygon center", "center z", 0, -1000, 1000, 2);
        double radius= QInputDialog::getDouble(NULL, "Enter regular polygon radius", "radius", 0, -1000, 1000, 2);
        //vtkNew<vtkRegularPolygonSource> polygonSource;
        polygonSource->SetNumberOfSides(numsides);
        double center[3] = { centerx, centery,centerz };
        polygonSource->SetCenter(center);
        polygonSource->SetRadius(radius);          // Horizontal radius
        /*vtkNew<vtkPolyDataMapper> polymapper;
        polymapper->SetInputConnection(polygonSource->GetOutputPort());
        vtkNew <vtkActor> polyactor;
        polyactor->SetMapper(polymapper);*/
        polyactor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        renderer->AddActor(polyactor);
        window->Render();
    }

    //void readInputFile(vtkLineSource* linesource, vtkGenericOpenGLRenderWindow* window, vtkTextActor* TextActor, vtkActor* lineActor) {
    //    QString fileObject = QFileDialog::getOpenFileName(nullptr, "Open File", ".", "Text Files (*.txt)");
    //    if (fileObject.isEmpty()) {
    //        return;  // Dialog was cancelled
    //    }
    //    QFile file(fileObject);
    //    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    //        return;
    //    else {
    //        QTextStream in(&file);
    //        double x1, y1, x2, y2;
    //        //QString color, property;
    //        if (!in.atEnd()) {
    //            QStringList linepoint1 = in.readLine().split(" ");
    //            /*x1 = linepoint1[0].toDouble();
    //            y1 = linepoint1[1].toDouble();*/
    //            linesource->SetPoint1(linepoint1[0].toDouble(), linepoint1[1].toDouble(), 0.0);
    //        }
    //        if (!in.atEnd()) {
    //           QStringList linepoint2 = in.readLine().split(" ");
    //           /*x2 = linepoint2[0].toDouble();
    //           y2 = linepoint2[1].toDouble();*/
    //           linesource->SetPoint2(linepoint2[0].toDouble(), linepoint2[1].toDouble(), 0.0);
    //        }
    //        if (!in.atEnd()) {
    //            QStringList rgb = in.readLine().split(" ");
    //            double rgbarr[3];
    //            rgbarr[0]=rgb.at(0).toDouble();
    //            rgbarr[1] = rgb.at(1).toDouble();
    //            rgbarr[2] = rgb.at(2).toDouble();
    //            lineActor->GetProperty()->SetColor(rgbarr);
    //        }
    //        /*if (!in.atEnd()) {
    //            property = in.readLine();
    //        }*/
    //        window->Render();
    //        updateTextCoordinates(linesource, TextActor, lineActor);
    //        file.close(); //close the file object.
    //        }
    //        //readfile.open(fileName, ios::in); //open a file to perform read operation using file object
    //        //if (readfile.is_open()) { //checking whether the file is open
    //        //    double x1, y1, x2, y2;
    //        //    readfile >> x1 >> y1 >> x2 >> y2;
    //        //    linesource->SetPoint1(x1, y1,0.0);
    //        //    linesource->SetPoint2(x2, y2, 0.0);
    //        //    window->Render();
    //        //    updateTextCoordinates(linesource, TextActor, lineActor);
    //        //    readfile.close(); //close the file object.
    //        //}
    //    
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

    QPushButton addEllipse;
    addEllipse.setText("Add Ellipse");
    dockLayout->addWidget(&addEllipse, 0, Qt::AlignTop);

    QPushButton addRegularPolygon;
    addRegularPolygon.setText("Add Regular Polygon");
    dockLayout->addWidget(&addRegularPolygon, 0, Qt::AlignTop);

    /*QPushButton setSecondCoordinate;
    setSecondCoordinate.setText("Set Second Coordinate");
    dockLayout->addWidget(&setSecondCoordinate, 1, Qt::AlignTop);

    QPushButton readFile;
    readFile.setText("Read Input File");
    dockLayout->addWidget(&readFile, 0, Qt::AlignTop);

    QPushButton writeFile;
    writeFile.setText("Write Input File");
    dockLayout->addWidget(&writeFile, 1, Qt::AlignTop);*/

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
    vtkNew<vtkSuperquadricSource> ellipsesource;
    vtkNew<vtkPolyDataMapper> ellipsemapper;
    ellipsemapper->SetInputData(ellipsesource->GetOutput());
    vtkNew<vtkActor> ellipseactor;
    ellipseactor->SetMapper(ellipsemapper);

    vtkNew<vtkRegularPolygonSource> polygonSource;
    vtkNew<vtkPolyDataMapper> polymapper;
    polymapper->SetInputConnection(polygonSource->GetOutputPort());
    vtkNew <vtkActor> polyactor;
    polyactor->SetMapper(polymapper);

    vtkNew<vtkRenderer> renderer;
    //renderer->AddActor(lineactor);
    window->AddRenderer(renderer);

    vtkNew<vtkPointPicker> pointPicker;
    window->SetInteractor(vtkRenderWidget->interactor());
    window->GetInteractor()->SetPicker(pointPicker);

    //vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
    //vtkNew<vtkPointPicker> pointPicker;
    window->SetInteractor(vtkRenderWidget->interactor());
    //window->GetInteractor()->SetPicker(pointPicker);
    //renderWindowInteractor->SetPicker(pointPicker);
    //interactor->SetInteractorStyle(style);

    /*vtkNew<customMouseInteractorStyle> style;
    style->setLineSource(linesource);
    style->setVTKActor(lineactor);
    window->GetInteractor()->SetInteractorStyle(style);*/

    //// connect the buttons
    QObject::connect(&addEllipse, &QPushButton::released,
        [&]() { ::addEllipse(ellipsesource, ellipseactor,window,renderer); });

    QObject::connect(&addRegularPolygon, &QPushButton::released,
        [&]() { ::addRegularPolygon(polygonSource, polyactor,window,renderer); });

   /* QObject::connect(&readFile, &QPushButton::released,
        [&]() { ::readInputFile(linesource, window, textActor, lineactor); });

    QObject::connect(&writeFile, &QPushButton::released,
        [&]() { ::writeInFile(linesource,lineactor); });*/
    
    mainWindow.show();

    return app.exec();
    //return EXIT_SUCCESS;
}
