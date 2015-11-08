#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QString>
#include <QTime>

typedef QVector<float> FloatVector;

class Model;

class Visualizer : public QGLWidget
{
    Q_OBJECT

public:
    Visualizer(QWidget *parent = 0);
    ~Visualizer();

    void setModel(Model* m);
    Model* getModel();

    QSize minimumSizeHint() const;
    QSize sizeHint() const;
    void saveGLState();
    void restoreGLState();
    void paintGLOld();
    void paintGL();
    void drawRectangularPrism(float, float, float, float, float, float);
    float getDisplayRadius();

public slots:
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
    void ZoomIn();
    void ZoomOut();
    void stop();
    void start();

    void displayParseError(QString msg);
    void PanUp();
    void PanDown();
    void PanLeft();
    void PanRight();
    void GUISetTargetDisplayRadius(const QString &);
    void GUISetTargetDetectRadius(const QString &);
    void GUISetClusterRadius(const QString &);
    void GUIgenerate();
    void GUIPlaceSearchers();


    //void updateGL();

signals:
    void xRotationChanged(int angle);
    void yRotationChanged(int angle);
    void zRotationChanged(int angle);
    void StartModel();
    void StopModel();
    void generate();
    void clearPaths();
    void VisPlaceSearchers();
    void VisFinishedPaintGL();
    void VisSetSearchType(int v);
    void VisSetClusterRadius(FloatVector v);

protected:
    void initializeGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void wheelEvent(QWheelEvent *);

private:
    int xRot;
    int yRot;
    int zRot;
    QPoint lastPos;
    QColor qtGreen;
    QColor qtPurple;
    Model* model;
    float scale;
    float x_min_bound;
    float y_min_bound;
    float z_min_bound;
    float x_max_bound;
    float y_max_bound;
    float z_max_bound;
    float horz_pan;
    float vert_pan;
    float target_display_radius;
    float target_detect_radius;
    QTime time;
    int frames;
};

#endif
