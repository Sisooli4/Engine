//
// Created by simon on 02.03.22.
//

#include "Line2D.h"
#include <limits>
#include <assert.h>
#include <algorithm>
#include <math.h>
#include <iostream>

Point2D::Point2D(double x, double y) : x(x), y(y) {}

double Point2D::getX() const {
    return x;
}

void Point2D::setX(double x) {
    Point2D::x = x;
}

double Point2D::getY() const {
    return y;
}

void Point2D::setY(double y) {
    Point2D::y = y;
}

Point2D::~Point2D() {

}

Point2D::Point2D() {}

Light::Light(const Mycolor &ambientLight) : ambientLight(ambientLight) {}

Light::Light(const Mycolor &ambientLight, const Mycolor &diffuseLight) : ambientLight(ambientLight),
                                                                         diffuseLight(diffuseLight) {}

Light::Light(const Mycolor &ambientLight, const Mycolor &diffuseLight, const Mycolor &specularLight) : ambientLight(
        ambientLight), diffuseLight(diffuseLight), specularLight(specularLight) {}

InfLight::InfLight(const Mycolor &ambientLight, const Mycolor &diffuseLight, const Vector3D &ldVector) : Light(
        ambientLight, diffuseLight), ldVector(ldVector) {}


const Point2D &Line2D::getP1() const {
    return p1;
}

void Line2D::setP1(const Point2D &p1) {
    Line2D::p1 = p1;
}

const Point2D &Line2D::getP2() const {
    return p2;
}

void Line2D::setP2(const Point2D &p2) {
    Line2D::p2 = p2;
}



Line2D::~Line2D() {

}

double Line2D::get_P1X() {
    return p1.x;
}

double Line2D::get_P1Y() {
    return p1.y;
}

double Line2D::get_P2X() {
    return p2.x;
}

double Line2D::get_P2Y() {
    return p2.y;
}

void Line2D::setP1(double x, double y) {
    p1.x = x;
    p1.y = y;
}

void Line2D::setP2(double x, double y) {
    p2.x = x;
    p2.y = y;
}


const Mycolor &Line2D::getColor() const {
    return color;
}

void Line2D::setColor(const Mycolor &color) {
    Line2D::color = color;
}

Line2D::Line2D(const Point2D &p1, const Point2D &p2, const Mycolor &color) : p1(p1), p2(p2), color(color) {}


Mycolor::Mycolor(double red, double green, double blue) : red(red), green(green), blue(blue) {

}
Mycolor::Mycolor():red(0), green(0), blue(0) {}

Mycolor::~Mycolor() {}

double Mycolor::getRed() const {
    return red;
}

void Mycolor::setRed(double red) {
    Mycolor::red = red;
}

double Mycolor::getGreen() const {
    return green;
}

void Mycolor::setGreen(double green) {
    Mycolor::green = green;
}

double Mycolor::getBlue() const {
    return blue;
}

void Mycolor::setBlue(double blue) {
    Mycolor::blue = blue;
}
Point2D doProjection(const Vector3D &point, const double d){

    double x = (d*point.x)/(-point.z);
    double y = (d*point.y)/(-point.z);

    return Point2D(x,y);
}

std::vector<double> Get_Coordinates(std::vector<Line2D*> &lines) {
    double BigX = lines[0]->getP1().getX();
    double BigY = lines[0]->getP1().getY();
    double SmallX = lines[0]->getP1().getX();
    double SmallY = lines[0]->getP1().getY();


    for (auto i:lines){
        if (i->getP1().getX() > BigX){
            BigX = i->getP1().getX();


        }
        else if(i->getP1().getX() < SmallX){
            SmallX = i->getP1().getX();
        }
        if (i->getP2().getX() > BigX){
            BigX = i->getP2().getX();
        }
        else if(i->getP2().getX() < SmallX){
            SmallX = i->getP2().getX();
        }
        if (i->getP1().getY() > BigY){
            BigY = i->getP1().getY();
        }
        else if(i->getP1().getY() < SmallY) {
            SmallY = i->getP1().getY();
        }
        if (i->getP2().getY() > BigY){
            BigY = i->getP2().getY();
        }
        else if(i->getP2().getY() < SmallY) {
            SmallY = i->getP2().getY();
        }
    }
    return std::vector<double> {BigX, BigY, SmallX, SmallY};

}

ZBuffer::ZBuffer(const int width, const int height) {
    double posInf = std::numeric_limits<double>::infinity();
    std::vector<std::vector<double>> v1;
    for (int i=0; i<width-1; i++){
        std::vector<double> v2;
        for(int j=0; j<height-1; j++){
            v2.push_back(posInf);
        }
        v1.push_back(v2);
    }
    zbuffer = v1;
}

void draw_zbuf_line(ZBuffer & Z, img::EasyImage & Im,  unsigned int x0,  unsigned int y0,  double z0,
                    unsigned int x1,  unsigned int y1,  double z1, const Mycolor &mycolor){
    assert(x0 < Im.get_width() && y0 < Im.get_height());
    assert(x1 < Im.get_width() && y1 < Im.get_height());
    img::Color color = img::Color(lround(mycolor.getRed()*255), lround(mycolor.getGreen()*255), lround(mycolor.getBlue()*255));
    double z;


    if(x0==x1 and y0==y1){
        z = 1/ std::max(z0,z1);
        if (z<=Z.zbuffer[x0][y0]){
            Z.zbuffer[x0][y0] = z;
            (Im)(x0, y0) = color;
        }
    }
    else if (x0 == x1){
        //special case for x0 == x1
        if(y1<y0){
            std::swap(y0,y1);
            std::swap(z0,z1);
        }
        double a = y1-y0;
        unsigned int k = 0;
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++){
            double cur_z_value = Z.zbuffer[x0][i];
            double iOverA = (a - k) / a;
            double new_z_value = (iOverA / z0) + ((1 - iOverA) / z1);
            k++;
            if (new_z_value < cur_z_value) {
                (Im)(x0, i) = color;
                Z.zbuffer[x0][i] = new_z_value;
            }
        }
    }
    //special case for y0 == y1
    else if (y0 == y1){
        if(x1<x0){
            std::swap(x0,x1);
            std::swap(z0,z1);
        }
        double a = x1 - x0;
        unsigned int k = 0;
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
            double cur_z_value = Z.zbuffer[i][y0];
            double iOverA = (a - k) / a;
            double new_z_value = (iOverA / z0) + ((1 - iOverA) / z1);
            k++;
            if (new_z_value < cur_z_value) {
                (Im)(i, y0) = color;
                Z.zbuffer[i][y0] = new_z_value;
            }
        }
        return;
    }
    else
    {
        if (x0 > x1)
        {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0)
        {
            double a = x1 - x0;
            unsigned int k = 0;
            for (unsigned int i = 0; i <= (x1 - x0); i++)
            {
                double cur_z_value = Z.zbuffer[x0 + i][(unsigned int) lround(y0 + m * i)];
                double iOverA = (a - k) / a;
                double new_z_value = (iOverA / z0) + ((1 - iOverA) / z1);
                k++;
                if (new_z_value < cur_z_value) {
                    (Im)(x0 + i, (unsigned int) lround(y0 + m * i)) = color;
                    Z.zbuffer[x0 + i][(unsigned int) lround(y0 + m * i)] = new_z_value;
                }
            }
            return;
        }
        else if (m > 1.0)
        {
            double a = y1 - y0;
            unsigned int k = 0;
            for (unsigned int i = 0; i <= (y1 - y0); i++) {
                double cur_z_value = Z.zbuffer[(unsigned int) lround(x0 + (i / m))][y0 + i];
                double iOverA = (a - k) / a;
                double new_z_value = (iOverA / z0) + ((1 - iOverA) / z1);
                k++;
                if (new_z_value < cur_z_value) {
                    (Im)((unsigned int) lround(x0 + (i / m)), y0 + i) = color;
                    Z.zbuffer[(unsigned int) lround(x0 + (i / m))][y0 + i] = new_z_value;
                }
            }
            return;
        }
        else if (m < -1.0)
        {
            double a = y0 - y1;
            unsigned int k = 0;
            for (unsigned int i = 0; i <= (y0 - y1); i++)
            {
                double cur_z_value = Z.zbuffer[(unsigned int) lround(x0 - (i / m))][y0 - i];
                double iOverA = (a - k) / a;
                double new_z_value = (iOverA / z0) + ((1 - iOverA) / z1);
                k++;

                if (new_z_value < cur_z_value) {
                    (Im)((unsigned int) lround(x0 - (i / m)), y0 - i) = color;
                    Z.zbuffer[(unsigned int) lround(x0 - (i / m))][y0 - i] = new_z_value;
                }
            }
            return;
        }
    }
}
img::EasyImage Draw_lines(const ini::Configuration& configuration, std::vector<Line2D*> &lines) {
    std::vector<double> Coordinates = Get_Coordinates(lines);
    double rangeX = Coordinates[0]-Coordinates[2];
    double rangeY = Coordinates[1]-Coordinates[3];
    double imageX = (configuration["General"]["size"].as_double_or_die()*rangeX)/std::max(rangeX,rangeY);
    double imageY = (configuration["General"]["size"].as_double_or_die()*rangeY)/std::max(rangeX,rangeY);
    double d = 0.95*(imageX/rangeX);

    for (auto i:lines) {
        i->setP1(i->get_P1X()*d, i->get_P1Y()*d);
        i->setP2(i->get_P2X()*d, i->get_P2Y()*d);
    }

    double DCX = d*((Coordinates[0]+Coordinates[2])/2);
    double DCY = d*((Coordinates[1]+Coordinates[3])/2);
    double dx = imageX/2 - DCX;
    double dy = imageY/2 - DCY;

    for (auto i:lines) {
        i->setP1(i->get_P1X()+dx, i->get_P1Y()+dy);
        i->setP2(i->get_P2X()+dx, i->get_P2Y()+dy);
    }
    std::vector<double> color = configuration["General"]["backgroundcolor"];
    img::EasyImage image = img::EasyImage (int(lround(imageX)), int(lround(imageY)), img::Color(int(lround(color[0]*255)),
                                                                                                int(lround(color[1]*255)),int(lround(color[2]*255))));
    if (configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe"){
        ZBuffer ZB = ZBuffer(int(lround(imageX)), int(lround(imageY)));
        for (auto i:lines){
            draw_zbuf_line(ZB, image, int(lround(i->get_P1X())), int(lround(i->get_P1Y())), i->z1, int(lround(i->get_P2X())), int(lround(i->get_P2Y())), i->z2, Mycolor(i->getColor().getRed(), i->getColor().getGreen(), i->getColor().getBlue()));
        }
        return image;
    }
    else{
        for(auto i:lines){
            image.draw_line(int(lround(i->get_P1X())), int(lround(i->get_P1Y())), int(lround(i->get_P2X())), int(lround(i->get_P2Y())), img::Color(lround(i->getColor().getRed()*255), lround(i->getColor().getGreen()*255), lround(i->getColor().getBlue()*255)));
        }
        return image;
    }
}

void draw_zbuf_triag(ZBuffer& Z, img::EasyImage& Im, Vector3D const&A, Vector3D const&B, Vector3D const&C,Vector3D const& Aa, Vector3D const& Ba,Vector3D const& Ca,
                     double d, double dx, double dy, const Mycolor& color){
    img::Color c = img::Color(lround(color.getRed()*255), lround(color.getGreen()*255), lround(color.getBlue()*255));
    double posInf = std::numeric_limits<double>::infinity();
    double negInf = -std::numeric_limits<double>::infinity();
    Point2D G = Point2D((Aa.x+Ba.x+Ca.x)/3, (Aa.y+Ba.y+Ca.y)/3);
    double eenOzg = ((1.0/(3*Aa.z)) + (1.0/(3*Ba.z)) + (1.0/(3*Ca.z)));
    Vector3D u = B-A;
    Vector3D v = C-A;
    Vector3D w = Vector3D::point(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
    Vector3D n = w;
    n.normalise();

    double k = w.x*A.x+w.y*A.y+w.z*A.z;
    double dzdx = w.x/(-d*k);
    double dzdy = w.y/(-d*k);

    int Ymin = std::min(Aa.y, std::min(Ba.y, Ca.y));
    int Ymax = std::max(Aa.y, std::max(Ba.y, Ca.y));

    for(int i=Ymin; i<=Ymax; i++){
        double xlab = posInf;
        double xlbc = posInf;
        double xlac = posInf;
        double xrab = negInf;
        double xrbc = negInf;
        double xrac = negInf;

        if((i-Aa.y)*(i-Ba.y) <= 0 and Aa.y != Ba.y){
            double xi = Ba.x + (Aa.x-Ba.x)*(i-Ba.y)/(Aa.y-Ba.y);
            xlab = xi;
            xrab = xi;
        }

        if((i-Ba.y)*(i-Ca.y) <= 0 and Ba.y != Ca.y){
            double xi = Ca.x + (Ba.x-Ca.x)*(i-Ca.y)/(Ba.y-Ca.y);
            xlbc = xi;
            xrbc = xi;
        }

        if((i-Aa.y)*(i-Ca.y) <= 0 and Aa.y != Ca.y){
            double xi = Ca.x + (Aa.x-Ca.x)*(i-Ca.y)/(Aa.y-Ca.y);
            xlac = xi;
            xrac = xi;
        }

        if(xlab == xlac and xlac == xlbc and xrab == xrac and xrac == xrbc){
            continue;
        }
        else{
            int xl = lround(std::min(xlab, std::min(xlbc, xlac))+0.5);
            int xr = lround(std::max(xrab, std::max(xrbc, xrac))-0.5);
            for(int j=xl; j<=xr; j++){
                double eenOverZ = 1.0001*eenOzg+(j-G.x)*dzdx+(i-G.y)*dzdy;
                if(eenOverZ<Z.zbuffer[j][i]){
                    Z.zbuffer[j][i] = eenOverZ;
                    (Im)(j, i) = c;
                }
            }
        }
    }
}
std::vector<double> Get_Coordinates(std::vector<Triangle> &triangles) {
    double BigX = triangles[0].getAa().getX();
    double BigY = triangles[0].getAa().getY();
    double SmallX = triangles[0].getAa().getX();
    double SmallY = triangles[0].getAa().getY();


    for (auto i:triangles){
        if (i.getAa().getX() > BigX){
            BigX = i.getAa().getX();
        }
        else if(i.getAa().getX() < SmallX){
            SmallX = i.getAa().getX();
        }
        if (i.getBa().getX() > BigX){
            BigX = i.getBa().getX();
        }
        else if(i.getBa().getX() < SmallX){
            SmallX = i.getBa().getX();
        }
        if (i.getCa().getX() > BigX){
            BigX = i.getCa().getX();
        }
        else if(i.getCa().getX() < SmallX){
            SmallX = i.getBa().getX();
        }
        if (i.getAa().getY() > BigY){
            BigY = i.getAa().getY();
        }
        else if(i.getAa().getY() < SmallY) {
            SmallY = i.getAa().getY();
        }
        if (i.getBa().getY() > BigY){
            BigY = i.getBa().getY();
        }
        else if(i.getBa().getY() < SmallY) {
            SmallY = i.getBa().getY();
        }
        if (i.getCa().getY() > BigY){
            BigY = i.getCa().getY();
        }
        else if(i.getCa().getY() < SmallY) {
            SmallY = i.getCa().getY();
        }
    }
    return std::vector<double> {BigX, BigY, SmallX, SmallY};

}

img::EasyImage Draw_tria(const ini::Configuration& configuration, std::vector<Triangle> triangles){
    std::vector<double> Coordinates = Get_Coordinates(triangles);
    double rangeX = Coordinates[0]-Coordinates[2];
    double rangeY = Coordinates[1]-Coordinates[3];
    double imageX = (configuration["General"]["size"].as_double_or_die()*rangeX)/std::max(rangeX,rangeY);
    double imageY = (configuration["General"]["size"].as_double_or_die()*rangeY)/std::max(rangeX,rangeY);
    std::cout<<"imageX:"<<imageX<<std::endl;
    std::cout<<"imageY:"<<imageY<<std::endl;
    double d = 0.95*(imageX/rangeX);
    std::cout<<"d:"<<d<<std::endl;
    for (auto& i:triangles) {
        i.setAa(doProjection(i.getA(), d));
        i.setBa(doProjection(i.getB(), d));
        i.setCa(doProjection(i.getC(), d));
    }

    double DCX = d*((Coordinates[0]+Coordinates[2])/2);
    double DCY = d*((Coordinates[1]+Coordinates[3])/2);
    double dx = imageX/2 - DCX;
    double dy = imageY/2 - DCY;
    std::cout<<"dx:"<<dx<<std::endl;
    std::cout<<"dy:"<<dy<<std::endl;
    for (auto& i:triangles) {
        i.setAa(Point2D(i.getAa().x+dx, i.getAa().y+dy));
        i.setBa(Point2D(i.getBa().x+dx, i.getBa().y+dy));
        i.setCa(Point2D(i.getCa().x+dx, i.getCa().y+dy));
    }
    std::vector<double> color = configuration["General"]["backgroundcolor"];
    img::EasyImage image = img::EasyImage (int(lround(imageX)), int(lround(imageY)), img::Color(int(lround(color[0]*255)),
                                                                                                int(lround(color[1]*255)),int(lround(color[2]*255))));
    ZBuffer ZB = ZBuffer(int(lround(imageX)), int(lround(imageY)));
    for (auto i:triangles){
        Vector3D Aa = Vector3D::point(i.getAa().x, i.getAa().y, i.getA().z);
        Vector3D Ba = Vector3D::point(i.getBa().x, i.getBa().y, i.getB().z);
        Vector3D Ca = Vector3D::point(i.getCa().x, i.getCa().y, i.getC().z);
        draw_zbuf_triag(ZB, image, i.getA(),i.getB(),i.getC(), Aa, Ba, Ca, d, dx, dy, i.getColor());
    }
    return image;
}

Triangle::Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, const Point2D &aa, const Point2D &ba,
                   const Point2D &ca, const Mycolor &c1) : A(a), B(b), C(c), Aa(aa), Ba(ba), Ca(ca), color(c1), lights({}) {}

Triangle::Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, Mycolor Color) : A(a), B(b), C(c), color(Color), lights({}){}

const Vector3D &Triangle::getA() const {
    return A;
}

const Vector3D &Triangle::getB() const {
    return B;
}

const Vector3D &Triangle::getC() const {
    return C;
}

const Point2D &Triangle::getAa() const {
    return Aa;
}

const Point2D &Triangle::getBa() const {
    return Ba;
}

const Point2D &Triangle::getCa() const {
    return Ca;
}

const Mycolor &Triangle::getColor() const{
    return color;
}

void Triangle::setA(const Vector3D &a) {
    A = a;
}

void Triangle::setB(const Vector3D &b) {
    B = b;
}

void Triangle::setC(const Vector3D &c) {
    C = c;
}

void Triangle::setAa(const Point2D &aa) {
    Aa = aa;
}

void Triangle::setBa(const Point2D &ba) {
    Ba = ba;
}

void Triangle::setCa(const Point2D &ca) {
    Ca = ca;
}

void Triangle::setC1(const Mycolor &c) {
    Triangle::color = c;
}





