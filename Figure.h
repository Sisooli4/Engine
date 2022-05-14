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
};



class Figure {
public:
    Figure(const std::vector<Vector3D> &points, const std::vector<Face> &faces);

    Figure();

    Figure(ini::Section conf, Figures3D &figures3D);

    std::vector<Vector3D> points;
    std::vector<Face> faces;
    Mycolor color;
    Matrix TFM;

    Figure createCube();

    Figure createTetrahedron();

    Figure createIcosahedron();

    Figure createOctahedron();

    Figure createDodecahedron();

    Figure createCone(const int n, const double h);

    Figure createCylinder(const int n, const double h);

    Figure createSphere(const int n);

    Figure createTorus(const double r,const double R,const int n,const int m);

    Figure createBuckyBall();

    Figure LS3D(std::string inputfile);


};



void toPolar(const Vector3D &point, double &theta, double &phi, double &r);

Matrix EyePointTrans(const Vector3D &point);

void applyTransformation(Figures3D &figs, const Matrix &m);

Matrix scaleFigure(const double scale);

Matrix rotateX(const double angle);

Matrix rotateY(const double angle);

Matrix rotateZ(const double angle);

Matrix translate(const Vector3D &vector);

Matrix CreateTransformationMatrix(const double angleX, const double angleY, const double angleZ, const double scale, const Vector3D &point);

std::vector<Face> triangulate(const Face& face);

void generateFractal(Figure& fig, Figures3D& fractal, const int nr_iterations, const double scale, ini::Section conf);

void applyTransformation(Figure &fig, const Matrix &m);


#endif //ENGINE_FIGURE_H
