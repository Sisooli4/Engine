//
// Created by simon on 02.03.22.
//

#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H
#include "easy_image.h"
#include <cmath>
#include <cstdio>
#include <vector>
#include "ini_configuration.h"
#include "vector3d.h"



class Point2D {

public:
    double x;
    double y;

    Point2D(double x, double y);

    Point2D();

    virtual ~Point2D();

    double getX() const;

    void setX(double x);

    double getY() const;

    void setY(double y);


};

class Mycolor{
    double red;
    double green;
    double blue;

public:
    Mycolor();

    Mycolor(double red, double green, double blue);

    virtual ~Mycolor();

    double getRed() const;

    void setRed(double red);

    double getGreen() const;

    void setGreen(double green);

    double getBlue() const;

    void setBlue(double blue);
};


class Light {
public:
    //de ambiente licht component
    Mycolor ambientLight;
    //de diffuse licht component
    Mycolor diffuseLight;
    //de diffuse licht component
    Mycolor specularLight;

    Light(const Mycolor &ambientLight);

    Light(const Mycolor &ambientLight, const Mycolor &diffuseLight);

    Light(const Mycolor &ambientLight, const Mycolor &diffuseLight, const Mycolor &specularLight);
};

typedef std::vector<Light> Lights3D;

class InfLight: public Light
{
public:
    //de richting waarin het
    //licht schijnt
    Vector3D ldVector;

    InfLight(const Mycolor &ambientLight, const Mycolor &diffuseLight, const Vector3D &ldVector);
};

class PointLight: public Light
{
public:
    //de locatie van de puntbron
    Vector3D location;
    //de hoek van een spotlicht
    double spotAngle;
};

class Line2D {
    Point2D p1;
    Point2D p2;
    Mycolor color;



public:
    double z1;
    double z2;

    Line2D(const Point2D &p1, const Point2D &p2, const Mycolor &color);

    virtual ~Line2D();

    const Point2D &getP1() const;

    void setP1(const Point2D &p1);

    const Point2D &getP2() const;

    void setP2(const Point2D &p2);

    const Mycolor &getColor() const;

    void setColor(const Mycolor &color);

    void setP1(double x, double y);

    void setP2(double x, double y);

    double get_P1X();

    double get_P1Y();

    double get_P2X();

    double get_P2Y();

};
Point2D doProjection(const Vector3D &point, const double d=1);

class ZBuffer: public std::vector<std::vector<double> >
{
public:
    std::vector<std::vector<double> > zbuffer;
    //Constructor: maakt een Z-Buffer van de correcte
    //grootte aan en initialiseert alle velden op +inf
    ZBuffer(const int width, const int height);
};

class Triangle{
    Vector3D A;
    Vector3D B;
    Vector3D C;
    Point2D Aa;
    Point2D Ba;
    Point2D Ca;
    Mycolor color;
    Lights3D lights;

public:
    Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, const Mycolor color);

    Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, const Point2D &aa, const Point2D &ba,
             const Point2D &ca, const Mycolor &color);

    const Vector3D &getA() const;

    const Vector3D &getB() const;

    const Vector3D &getC() const;

    const Point2D &getAa() const;

    const Point2D &getBa() const;

    const Point2D &getCa() const;

    const Mycolor &getColor() const;

    void setA(const Vector3D &a);

    void setB(const Vector3D &b);

    void setC(const Vector3D &c);

    void setAa(const Point2D &aa);

    void setBa(const Point2D &ba);

    void setCa(const Point2D &ca);

    void setC1(const Mycolor &c);

};

using Lines2D = std::vector<Line2D>;

std::vector<double> Get_Coordinates(std::vector<Line2D*>& lines);

void draw_zbuf_line(ZBuffer & Z, img::EasyImage & Im, unsigned int x0, unsigned int y0, double z0,
                    unsigned int x1, unsigned int y1, double z1, const Mycolor &color);

img::EasyImage Draw_lines(const ini::Configuration& configuration, std::vector<Line2D*>& lines);

std::vector<double> Get_Coordinates(std::vector<Triangle>& triangles);

img::EasyImage Draw_tria(const ini::Configuration& configuration, std::vector<Triangle>);

void draw_zbuf_triag(ZBuffer&, img::EasyImage&, Vector3D const& A, Vector3D const& B,Vector3D const& C,Vector3D const& Aa, Vector3D const& Ba,Vector3D const& Ca,
                     double d, double dx, double dy, const Mycolor& color);






#endif //ENGINE_LINE2D_H
