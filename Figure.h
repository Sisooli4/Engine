//
// Created by simon on 09.03.22.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H
#include "Line2D.h"
#include "vector3d.h"
#include "easy_image.h"
#include "ini_configuration.h"
#include <cmath>
#include <vector>
#include "l_parser.h"
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

class Face{
public:
    std::vector<int> point_indexes;

    Face(const std::vector<int> &pointIndexes);

    Face();
};

class L3Dsystem {
    std::set<char> alphabet;
    std::map<char, std::string> replacement;
    double angle;
    std::string initiator;
    unsigned int iterations;
    std::map<char, bool> draw;
    std::string DrawString;


public:
    std::vector<Vector3D> pointz;
    std::vector<Face> facez;
    L3Dsystem(std::string inputfile);

    std::string replaceString(std::string, unsigned int iteration);

    std::pair<std::vector<Vector3D>, std::vector<Face>> getP_F();
};

class Figure;
typedef std::vector<Figure> Figures3D;

class Figures{
    std::vector<Figure> figures;
    Mycolor BackGroundcolor;
    Matrix EyePointMatrix;
    Lights3D Lights;
    double dNear;
    double dFar;
    double right;
    double left;
    double top;
    double bottom;

public:
    Figures(const ini::Configuration &configuration);
    img::EasyImage Draw_3Dlines(const ini::Configuration &configuration);
    img::EasyImage Draw_3Dtria(const ini::Configuration &configuration);
    std::vector<Triangle> clipFigures(std::vector<Triangle> triangles);
    double getP(Vector3D p1, Vector3D p2, double val, double dNear, int teller);

    const Lights3D getLights() const;

    const Matrix &getEyePointMatrix() const;

    const std::vector<Figure> &getFigures() const;

    void setFigures(const std::vector<Figure> &figures);

    void generateThickFigures();

    float roundOff(float value, unsigned char prec);
};

class Figure {
public:
    Figure(const std::vector<Vector3D> &points, const std::vector<Face> &faces);

    Figure();

    Figure(ini::Section conf, Figures3D &figures3D, Lights3D lights);

    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Mycolor color;
    Light Reflections;
    double ReflectionCoefficient;
    Matrix TFM;
    int n;
    int m;
    double r;
    bool thick = false;

    Figure createCube();

    Figure createTetrahedron();

    Figure createIcosahedron();

    Figure createOctahedron();

    Figure createDodecahedron();

    Figure createCone(const int n, const double h);

    Figure createCylinder(const int n, const double h, bool thick);

    Figure createSphere(const int n);

    Figure createTorus(const double r,const double R,const int n,const int m);

    Figure createBuckyBall();

    void createMengerSponge(ini::Section conf, Figures3D& fractal, int nrIterations, const Figure& cube, Mycolor Color, Light reflex, double reflexcoef);

    Figure LS3D(std::string inputfile);

    void applyTransformation(Figure &fig, const Matrix &m);

    void generateThickFigure(Figures3D &resultingFigures);
};



void toPolar(const Vector3D &point, double &theta, double &phi, double &r);

Matrix EyePointTrans(const Vector3D &point, bool clipping);

void applyTransformation(Figures3D &figs, const Matrix &m);

Matrix scaleFigure(const double scale);

Matrix rotateX(const double angle);

Matrix rotateY(const double angle);

Matrix rotateZ(const double angle);

Matrix translate(const Vector3D &vector);

Matrix CreateTransformationMatrix(const double angleX, const double angleY, const double angleZ, const double scale, const Vector3D &point);

std::vector<Face> triangulate(const Face& face);

void generateFractal(Figure& fig, Figures3D& fractal, const int nr_iterations, const double scale, ini::Section conf, Mycolor color, Light l, double reflexcoef);

void applyTransformation(Figure &fig, const Matrix &m);




#endif //ENGINE_FIGURE_H
