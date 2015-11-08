#include <QtWidgets>
#include <QtGui>
#include <iostream>

#if defined(__APPLE__) || defined(MACOSX)
  #include <GLUT/glut.h>
#else
  #include "GL/glut.h"
  #include "GL/glu.h"
#endif

#include <QtOpenGL>
#include <math.h>

#include "Visualizer.h"
#include "Model.h"
#include "Sphere.h"
#include "Searcher.h"
#include "Target.h"

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

//! [0]
Visualizer::Visualizer(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    frames = 0;

    setMinimumHeight(500);
    setMinimumWidth(500);

    xRot = 0;
    yRot = 0;
    zRot = 0;
    scale = 0.001f;

    x_min_bound = -1;
    y_min_bound = -1;
    z_min_bound = -1;

    x_max_bound = 1;
    y_max_bound = 1;
    z_max_bound = 1;

    horz_pan = -0.4f;
    vert_pan = -0.4f;

    target_display_radius = 5.0;
    target_detect_radius = 5.0;

}
//! [0]

//! [1]
Visualizer::~Visualizer()
{
    cout << "Visualizer: cleaning up" << endl;
}
//! [1]

//! [2]
QSize Visualizer::minimumSizeHint() const
{
    return QSize(50, 50);
}
//! [2]

//! [3]
QSize Visualizer::sizeHint() const
//! [3] //! [4]
{
    return QSize(1000, 1000);
}
//! [4]

static void qNormalizeAngle(int &angle)
{
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}

//! [5]
void Visualizer::setXRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        emit xRotationChanged(angle);
        updateGL();
    }
}
//! [5]

void Visualizer::setYRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != yRot) {
        yRot = angle;
        emit yRotationChanged(angle);
        updateGL();
    }
}

void Visualizer::setZRotation(int angle)
{
    qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        emit zRotationChanged(angle);
        updateGL();
    }
}

//! [6]
void Visualizer::initializeGL()
{
    //qglClearColor(Qt::white);
    //glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    char *my_argv[] = { "PowerSearch3D", NULL };
    int   my_argc = 1;
    //glutInit(&my_argc, my_argv);
    //glutInitDisplayMode (GLUT_DOUBLE);
    //int window_id = glutCreateWindow("Glut Window");
    //glutSetWindow(window_id);



    //glEnable(GL_DEPTH_TEST);
    //glEnable(GL_CULL_FACE);
    glShadeModel(GL_SMOOTH);
    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
    glEnable(GL_MULTISAMPLE);
   // static GLfloat lightPosition[4] = { 0.5, 5.0, 7.0, 1.0 };
    //glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_DEPTH_CLAMP);



}
//! [6]

void Visualizer::drawRectangularPrism(float x_min, float y_min, float z_min, float x_max, float y_max, float z_max)
{
    //cout << x_min << " " << y_min << " " << z_min << " " << x_max << " " << y_max << " " << z_max << " " << endl;

    glColor4f(1.0, 0.0, 0.0, 0.7);
    glBegin(GL_LINE_LOOP);
    glVertex3f(  x_min,   y_max,   z_max );
    glVertex3f(  x_max,   y_max,   z_max );
    glVertex3f(  x_max,   y_max,  z_min );
    glVertex3f(  x_min,   y_max,  z_min );
    glEnd();

    glColor4f(0.0, 1.0, 0.0, 0.7);
    glBegin(GL_LINE_LOOP);
    glVertex3f( x_min,  y_min,   z_max );
    glVertex3f( x_max,  y_min,   z_max );
    glVertex3f( x_max,  y_min,   z_min );
    glVertex3f( x_min,  y_min,   z_min );
    glEnd();


    glColor4f(1.0, 1.0, 1.0, 0.7);
    glBegin(GL_LINE_LOOP);
    glVertex3f( x_min,  y_max,  z_max );
    glVertex3f( x_min, y_min,  z_max );
    glVertex3f( x_max,  y_min,  z_max );
    glVertex3f( x_max,   y_max,  z_max );
    glEnd();

    glColor4f(0.0, 0.0, 1.0, 0.7);
    glBegin(GL_LINE_LOOP);
    glVertex3f( x_max,  y_max,  z_min );
    glVertex3f( x_min, y_max,  z_min );
    glVertex3f( x_min,  y_min,  z_min );
    glVertex3f( x_max,   y_min,  z_min );
    glEnd();
}

void drawSphere(float x, float y, double r, int lats, int longs)
{
int i, j;
for(i = 0; i <= lats; i++) {
double lat0 = M_PI * (-0.5 + (double) (i - 1) / lats);
double z0 = sin(lat0);
double zr0 = cos(lat0);

double lat1 = M_PI * (-0.5 + (double) i / lats);
double z1 = sin(lat1);
double zr1 = cos(lat1);

glBegin(GL_QUAD_STRIP);
for(j = 0; j <= longs; j++) {
double lng = 2 * M_PI * (double) (j - 1) / longs;
double x = cos(lng);
double y = sin(lng);

glNormal3f(x * zr0, y * zr0, z0);
glColor3f( 1,0,0);
glVertex3f(x * zr0, y * zr0, z0);
glNormal3f(x * zr1, y * zr1, z1);
glColor3f(0,1,0);
glVertex3f(x * zr1, y * zr1, z1);
}
glEnd();
}
}


//! [7]

void Visualizer::paintGL()
{ 
    QPainter painter;
    painter.begin(this);



   // qApp->processEvents();
    //cout << "paintGL called" << endl;
    model->gui_ready = true;

    model->drawing = true;

    if (model->paths_mutex)
    {
        cout << "Visualizer: paths_mutex is set" << endl;
        model->paths_mutex_seen = true;
        return;
    }

    if (model->targets_mutex)
    {
        cout << "Visualizer: targets_mutex is set" << endl;
        model->targets_mutex_seen = true;
        return;
    }

    if (model->searchers_mutex)
    {
        cout << "Visualizer: searchers_mutex is set" << endl;
        model->searchers_mutex_seen = true;
        return;
    }

    x_min_bound = getModel()->getSearchSpace()->getXMinBound();
    y_min_bound = getModel()->getSearchSpace()->getYMinBound();
    z_min_bound = getModel()->getSearchSpace()->getZMinBound();
    x_max_bound = getModel()->getSearchSpace()->getXMaxBound();
    y_max_bound = getModel()->getSearchSpace()->getYMaxBound();
    z_max_bound = getModel()->getSearchSpace()->getZMaxBound();

    saveGLState();

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    //glTranslatef(0.0, 0.0, -10.0);
    glFrustum(-1, 1, -1, 1, 0, 100);
    glViewport(0, 0, 2*width(), 2*height());
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    //glTranslatef(0.0f, 0.0f, -15.0f);
    glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
    glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
    glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);
    glTranslatef(horz_pan, vert_pan, 0.0f);
    //glScalef (1.0, 1.0, 1.0);      /* modeling transformation */
    glScalef(scale, scale, scale);
    glColor4f(1.0, 1.0, 1.0, 0.3);



    drawRectangularPrism(x_min_bound, y_min_bound, z_min_bound, x_max_bound, y_max_bound, z_max_bound);

    QPointF pos;
    QPointF center_pos;
    SearchSpace* space = model->getSearchSpace();
    int n_searchers = space->getNumSearchers();
    Searcher* searchers = space->getSearchers();

    //cout << "Number of searchers: " << n_searchers << endl;

    // draw boundary
    QRect frame_rect = this->contentsRect();
    int frame_height = frame_rect.height();
    int frame_width = frame_rect.width();

    center_pos = frame_rect.center();

    for ( int i = 0; i < n_searchers; i++ )
    {
        float scaled_x = searchers[i].getXPos();//*frame_width/2+center_pos.x();
        float scaled_y = searchers[i].getYPos();//*frame_height/2+center_pos.y();
        float scaled_z = searchers[i].getZPos();//*frame_height/2+center_pos.y();

        glColor3f(searchers[i].red/255.0, searchers[i].blue/255.0, searchers[i].green/255.0);

               GLUquadricObj *quadric=gluNewQuadric();
               gluQuadricNormals(quadric, GLU_SMOOTH);
               glPushMatrix();
               glTranslatef( scaled_x, scaled_y, scaled_z );
               glPushMatrix();
               gluSphere(quadric, target_display_radius, 10,10);
               glPopMatrix();
               gluDeleteQuadric(quadric);

            //cout << "Raw (" << searchers[i].getXPos() << ", " << searchers[i].getYPos() << ", " << searchers[i].getZPos() <<")" << endl;
            //cout << "Scaled (" << scaled_x << ", " << scaled_y << ", " << scaled_z << ")" << endl;


        // Display path
        vector<Coordinate*> path = searchers[i].getPath();

        for(vector<Coordinate*>::iterator it = path.begin(); it != path.end(); ++it)
        {

            float scaled_path_x = (*it)->getX();//*frame_width/2;//+center_pos.x();
            float scaled_path_y = (*it)->getY();//*frame_height/2;//+center_pos.y();
            float scaled_path_z = (*it)->getZ();//*frame_height/2;//+center_pos.y(); // <---- FIX FOR Z

            Coordinate start, end;
            start.setX(scaled_path_x);
            start.setY(scaled_path_y);
            start.setZ(scaled_path_z);


            if (it+1 != path.end())
            {
                scaled_path_x = (*(it+1))->getX();//*frame_width/2;//+center_pos.x();
                scaled_path_y = (*(it+1))->getY();//*frame_height/2;//+center_pos.y();
                scaled_path_z = (*(it+1))->getZ();//*frame_height/2;//+center_pos.y();
            }
            else
            {
                break;
            }

            end.setX(scaled_path_x);
            end.setY(scaled_path_y);
            end.setZ(scaled_path_z);

            //painter.drawLine(path_segment);

            glBegin(GL_LINES);
            glVertex3d(start.getX(), start.getY(), start.getZ());
            glVertex3d(end.getX(), end.getY(), end.getZ());
            glEnd();

        }


    }

    // display targets
    int n_targets = space->getNumTargets();
    Target* targets = space->getTargets();

    for ( int i = 0; i < n_targets; i++ )
    {
        if ( targets[i].isFound() ) glColor4f( 0.0, 1.0, 1.0, 1.0 );
        else glColor4f( 0.0, 1.0, 0.0, 0.2 );
        GLUquadricObj *quadric=gluNewQuadric();
           gluQuadricNormals(quadric, GLU_SMOOTH);
           glPushMatrix();
           glTranslatef( targets[i].getXPos(),targets[i].getYPos() ,targets[i].getZPos() );
           glPushMatrix();
           gluSphere(quadric, target_display_radius, 10,10);
           glPopMatrix();
           gluDeleteQuadric(quadric);
    }

    restoreGLState();
    swapBuffers();


    //emit VisFinishedPaintGL();

    QString framesPerSecond;
       framesPerSecond.setNum(frames /(time.elapsed() / 1000.0), 'f', 2);

       painter.setPen(Qt::white);

       painter.drawText(20, 40, framesPerSecond + " fps");

       painter.end();

       swapBuffers();

       frames++;

       if (!(frames % 100)) {
           time.start();
           frames = 0;
       }


    model->drawing = false;

    //glFinish();
    //qApp->processEvents();

}


void Visualizer::paintGLOld()
{
/*

    saveGLState();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    //glTranslatef(0.0, 0.0, -10.0);
    glFrustum(-1, 1, -1, 1, 0, 100);
    //glViewport(0, 0, 2*width(), 2*height());
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glTranslatef(0.0f+horz_pan, 0.0f+vert_pan, -15.0f);
    glRotatef(xRot / 16.0, 1.0, 0.0, 0.0);
    glRotatef(yRot / 16.0, 0.0, 1.0, 0.0);
    glRotatef(zRot / 16.0, 0.0, 0.0, 1.0);
    //glScalef (1.0, 1.0, 1.0);      // modeling transformation
    glScalef(scale, scale, scale);
    glColor4f(1.0, 1.0, 1.0, 0.3);
    glutWireSphere (1.0, 20, 20);

    QPointF pos;
    QPointF center_pos;
    SearchSphere* space = model->getSearchSpace();
    int n_searchers = space->getNumSearchers();
    Searcher* searchers = space->getSearchers();

    // draw boundary
    QRect frame_rect = this->contentsRect();
    int frame_height = frame_rect.height();
    int frame_width = frame_rect.width();

    center_pos = frame_rect.center();

    for ( int i = 0; i < n_searchers; i++ )
    {
        float scaled_x = searchers[i].getXPos();//*frame_width/2+center_pos.x();
        float scaled_y = searchers[i].getYPos();//*frame_height/2+center_pos.y();
        float scaled_z = searchers[i].getZPos();//*frame_height/2+center_pos.y();

        glColor3f(searchers[i].red/255.0, searchers[i].blue/255.0, searchers[i].green/255.0);


               GLUquadricObj *quadric=gluNewQuadric();
               gluQuadricNormals(quadric, GLU_SMOOTH);
               glPushMatrix();
               glTranslatef( scaled_x, scaled_y, scaled_z );
               glPushMatrix();
               gluSphere(quadric, searchers[i].getRadius(), 10,10);
               glPopMatrix();
               gluDeleteQuadric(quadric);

            //cout << "Raw (" << searchers[i].getXPos() << ", " << searchers[i].getYPos() << ", " << searchers[i].getZPos() <<")" << endl;
            //cout << "Scaled (" << scaled_x << ", " << scaled_y << ", " << scaled_z << ")" << endl;


        // Display path
        vector<Coordinate*> path = searchers[i].getPath();

        for(vector<Coordinate*>::iterator it = path.begin(); it != path.end(); ++it)
        {
            float scaled_path_x = (*it)->getX();//*frame_width/2;//+center_pos.x();
            float scaled_path_y = (*it)->getY();//*frame_height/2;//+center_pos.y();
            float scaled_path_z = (*it)->getZ();//*frame_height/2;//+center_pos.y(); // <---- FIX FOR Z

            Coordinate start, end;
            start.setX(scaled_path_x);
            start.setY(scaled_path_y);
            start.setZ(scaled_path_z);


            if (it+1 != path.end())
            {
                scaled_path_x = (*(it+1))->getX();//*frame_width/2;//+center_pos.x();
                scaled_path_y = (*(it+1))->getY();//*frame_height/2;//+center_pos.y();
                scaled_path_z = (*(it+1))->getZ();//*frame_height/2;//+center_pos.y();
            }
            else
            {
                break;
            }

            end.setX(scaled_path_x);
            end.setY(scaled_path_y);
            end.setZ(scaled_path_z);

            //painter.drawLine(path_segment);

            glBegin(GL_LINES);
            glVertex3d(start.getX(), start.getY(), start.getZ());
            glVertex3d(end.getX(), end.getY(), end.getZ());
            glEnd();

        }


    }

    // display targets
    int n_targets = space->getNumTargets();
    Target* targets = space->getTargets();

    for ( int i = 0; i < n_targets; i++ )
    {
        if ( targets[i].isFound() ) glColor4f( 0.0, 1.0, 0.0, 1.0 );
        else glColor4f( 0.0, 1.0, 0.0, 0.2 );
        GLUquadricObj *quadric=gluNewQuadric();
           gluQuadricNormals(quadric, GLU_SMOOTH);
           glPushMatrix();
           glTranslatef( targets[i].getXPos(),targets[i].getYPos() ,targets[i].getZPos() );
           glPushMatrix();
           gluSphere(quadric, targets[i].getRadius(), 10,10);
           glPopMatrix();
           gluDeleteQuadric(quadric);
    }

    restoreGLState();
    */

}
//! [7]

//! [8]
void Visualizer::resizeGL(int width, int height)
{
    int side = qMin(width, height);
    //glViewport((width - side) / 2, (height - side) / 2, side, side);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
#ifdef QT_OPENGL_ES_1
    glOrthof(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
#else
    glOrtho(-0.5, +0.5, -0.5, +0.5, 4.0, 15.0);
#endif
    glMatrixMode(GL_MODELVIEW);
}
//! [8]

//! [9]
void Visualizer::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
}
//! [9]

//! [10]
void Visualizer::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();

    if (event->buttons() & Qt::LeftButton) {
        setXRotation(xRot + 8 * dy);
        setYRotation(yRot + 8 * dx);
    } else if (event->buttons() & Qt::RightButton) {
        setXRotation(xRot + 8 * dy);
        setZRotation(zRot + 8 * dx);
    }
    lastPos = event->pos();
}

void Visualizer::wheelEvent(QWheelEvent *e)
 {
    //cout << "Wheel event" << endl;
     e->delta() > 0 ? ZoomIn() : ZoomOut();
 }

void Visualizer::ZoomIn()
{
    scale += scale*0.1f;
   // cout << "Zoom: " << scale << endl;
    updateGL();
}

void Visualizer::ZoomOut()
{
    scale -= scale*0.1f;
 //   cout << "Zoom: " << scale << endl;
    updateGL();
}


void Visualizer::saveGLState()
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
}

void Visualizer::restoreGLState()
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glPopAttrib();
}

void Visualizer::setModel(Model* m)
{
    model = m;
}

Model* Visualizer::getModel()
{
    return model;
}

void Visualizer::stop()
{
    emit StopModel();
}

void Visualizer::start()
{
    emit StartModel();
}

QStringList Visualizer::OpenFile()
{
    QString default_path = QString::fromStdString(model->getWorkingPath());

    #if defined(__APPLE__)
    default_path += "/MyFile.txt";
    #endif

    QStringList fileNames = QFileDialog::getOpenFileNames(this, tr("Open File"), default_path, tr("Comma Separated Values File (*.csv)"));

    if (!fileNames.isEmpty())
        getModel()->queueInputFiles( fileNames );

    return fileNames;
}

void Visualizer::displayParseError(QString msg)
{
    QErrorMessage errorMessage;
    errorMessage.showMessage(msg);
    errorMessage.exec();
}

void Visualizer::PanRight()
{
    horz_pan = horz_pan-0.1;
    cout << "Horz Pan: " << horz_pan << endl;
    updateGL();
}

void Visualizer::PanLeft()
{
    horz_pan = horz_pan+0.1;
    cout << "Horz Pan: " << horz_pan << endl;
    updateGL();
}

void Visualizer::PanUp()
{
    vert_pan = vert_pan-0.1;
    cout << "Vert Pan: " << vert_pan << endl;
    updateGL();
}

void Visualizer::PanDown()
{
    vert_pan = vert_pan+0.1;
    cout << "Vert Pan: " << vert_pan << endl;
    updateGL();
}

void Visualizer::GUISetTargetDisplayRadius(const QString &textvalue)
{
    target_display_radius = textvalue.toFloat();
}

void Visualizer::GUISetTargetDetectRadius(const QString &textvalue)
{
    target_detect_radius = textvalue.toFloat();
    getModel()->getSearchSpace()->setTargetDetectRadius(target_detect_radius);
    getModel()->calcTargetsFound();
}

void Visualizer::GUISetClusterRadius(const QString &textvalue)
{
    //cout << "Visualizer::GUISetClusterRadius called" << endl;
    float cluster_radius = textvalue.toFloat();
    emit VisSetClusterRadius(cluster_radius);
}

float Visualizer::getDisplayRadius()
{
    return this->target_display_radius;
}

void Visualizer::GUIgenerate()
{
    cout << "Visualizer::GUIgenerate called" << endl;
    emit generate();
}

void Visualizer::GUIPlaceSearchers()
{
    cout << "Visualizer::GUIPlaceSearchers called" << endl;
    emit VisPlaceSearchers();
}

void Visualizer::GUISetSearchType(int v)
{
    cout << "Visualizer::GUISetSearchType called with " << v << endl;
    emit VisSetSearchType( v );
}

