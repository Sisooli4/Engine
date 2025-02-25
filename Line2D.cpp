//
// Created by simon on 02.03.22.
//

#include "Line2D.h"
#include <limits>
#include <cassert>
#include <algorithm>
#include <cmath>
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

Light::Light(const Mycolor &ambientLight) : ambientLight(ambientLight) {
    diffuseLight = Mycolor(-1,-1,-1);
    specularLight = Mycolor(-1,-1,-1);
}

Light::Light(const Mycolor &ambientLight, const Mycolor &diffuseLight) : ambientLight(ambientLight),
                                                                         diffuseLight(diffuseLight) {
    specularLight = Mycolor(-1,-1,-1);
}

Light::Light(const Mycolor &ambientLight, const Mycolor &diffuseLight, const Mycolor &specularLight) : ambientLight(
        ambientLight), diffuseLight(diffuseLight), specularLight(specularLight) {}

Light::Light() {}

Light::~Light() {

}

InfLight::InfLight(const Mycolor &ambientLight, const Mycolor &diffuseLight, const Vector3D &ldVector) : Light(
        ambientLight, diffuseLight), ldVector(ldVector) {
    specularLight = Mycolor(-1,-1,-1);
}

InfLight::InfLight() {
    ambientLight = Mycolor(-1,-1,-1);
    diffuseLight = Mycolor(-1,-1,-1);
    specularLight = Mycolor(-1,-1,-1);
}

const Vector3D &InfLight::getLdVector() const {
    return ldVector;
}


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


Mycolor::Mycolor(double red, double green, double blue) {
    if(red>1){
        Mycolor::red = 1;
    }
    else{
        Mycolor::red=red;
    }
    if(green>1){
        Mycolor::green = 1;
    }
    else{
        Mycolor::green = green;
    }
    if(blue>1){
        Mycolor::blue = 1;
    }
    else{
        Mycolor::blue = blue;
    }
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

void draw_zbuf_triag(ZBuffer& Z, img::EasyImage& Im, Triangle triangle, double d, double dx, double dy){
    double posInf = std::numeric_limits<double>::infinity();
    double negInf = -std::numeric_limits<double>::infinity();
    bool back = false;
    Point2D G = Point2D((triangle.getAa().x+triangle.getBa().x+triangle.getCa().x)/3, (triangle.getAa().y+triangle.getBa().y+triangle.getCa().y)/3);
    double eenOzg = ((1.0/(3*triangle.getA().z)) + (1.0/(3*triangle.getB().z)) + (1.0/(3*triangle.getC().z)));
    Vector3D u = triangle.getB()-triangle.getA();
    Vector3D v = triangle.getC()-triangle.getA();
    Vector3D w = Vector3D::point(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
    double k = w.x*triangle.getA().x+w.y*triangle.getA().y+w.z*triangle.getA().z;
    if(triangle.getClipping() and k>0){
        back = true;
    }

    Vector3D n = Vector3D::normalise(Vector3D::cross(u,v));
    if(back){
        n = -n;
    }

    double cos;
    Vector3D l;
    if(triangle.getReflections().diffuseLight.getRed() != -1){
        for (auto* light:triangle.getLights()){
            auto* infLight = dynamic_cast<InfLight*>(light);
            if (infLight != nullptr){
                l = infLight->ldVector*triangle.getEyePointMatrix();
                l = -l;
                l.normalise();
                cos = n.dot(l);
                if (cos > 0){
                    triangle.setC1(Mycolor((triangle.getColor().getRed()+ (light->diffuseLight.getRed() *triangle.getReflections().diffuseLight.getRed()*cos)),
                                           (triangle.getColor().getGreen()+ (light->diffuseLight.getGreen()*triangle.getReflections().diffuseLight.getGreen()*cos)),
                                           (triangle.getColor().getBlue()+(light->diffuseLight.getBlue()*triangle.getReflections().diffuseLight.getBlue()*cos))));

                }
                }
        }
    }
    img::Color c = img::Color(lround(triangle.getColor().getRed()*255), lround(triangle.getColor().getGreen()*255), lround(triangle.getColor().getBlue()*255));



    double dzdx = w.x/(-d*k);
    double dzdy = w.y/(-d*k);

    int Ymin = lround(std::min(triangle.getAa().y, std::min(triangle.getBa().y, triangle.getCa().y))+0.5);
    int Ymax = lround(std::max(triangle.getAa().y, std::max(triangle.getBa().y, triangle.getCa().y))-0.5);

    for(int i=Ymin; i<=Ymax; i++) {
        double xlab = posInf;
        double xlbc = posInf;
        double xlac = posInf;
        double xrab = negInf;
        double xrbc = negInf;
        double xrac = negInf;
        Mycolor co = Mycolor(-1,-1,-1);
        Mycolor col = Mycolor(-1,-1,-1);

        if ((i - triangle.getAa().y) * (i - triangle.getBa().y) <= 0 and triangle.getAa().y != triangle.getBa().y) {
            double xi = triangle.getBa().x + (triangle.getAa().x - triangle.getBa().x) * (i - triangle.getBa().y) /
                                             (triangle.getAa().y - triangle.getBa().y);
            xlab = xi;
            xrab = xi;
        }

        if ((i - triangle.getBa().y) * (i - triangle.getCa().y) <= 0 and triangle.getBa().y != triangle.getCa().y) {
            double xi = triangle.getCa().x + (triangle.getBa().x - triangle.getCa().x) * (i - triangle.getCa().y) /
                                             (triangle.getBa().y - triangle.getCa().y);
            xlbc = xi;
            xrbc = xi;
        }

        if ((i - triangle.getAa().y) * (i - triangle.getCa().y) <= 0 and triangle.getAa().y != triangle.getCa().y) {
            double xi = triangle.getCa().x + (triangle.getAa().x - triangle.getCa().x) * (i - triangle.getCa().y) /
                                             (triangle.getAa().y - triangle.getCa().y);
            xlac = xi;
            xrac = xi;
        }

        if (xlab == xlac and xlac == xlbc and xrab == xrac and xrac == xrbc) {
            continue;
        } else {
            int xl = lround(std::min(xlab, std::min(xlbc, xlac)) + 0.5);
            int xr = lround(std::max(xrab, std::max(xrbc, xrac)) - 0.5);
            for (int j = xl; j <= xr; j++) {
                double eenOverZ = 1.0001 * eenOzg + (j - G.x) * dzdx + (i - G.y) * dzdy;
                if (eenOverZ < Z.zbuffer[j][i]) {
                    double z = 1.0 / eenOverZ;
                    double x = -z * (j - dx) / d;
                    double y = -z * (i - dy) / d;
                    int teller = 0;
                    if (triangle.getReflections().diffuseLight.getRed() != -1) {
                        for (auto *light: triangle.getLights()) {
                            auto *pointLight = dynamic_cast<PointLight *>(light);
                            if (pointLight != nullptr) {
                                l = pointLight->location * triangle.getEyePointMatrix() -
                                    Vector3D::point(x, y, z);
                                l.normalise();
                                cos = n.dot(l);
                                if (cos > 0) {
                                    if (pointLight->spotAngle == -1) {
                                        if (teller == 0){
                                            teller++;
                                            co = Mycolor((triangle.getColor().getRed() +
                                                          (light->diffuseLight.getRed() *
                                                           triangle.getReflections().diffuseLight.getRed() *
                                                           cos)),
                                                         (triangle.getColor().getGreen() +
                                                          (light->diffuseLight.getGreen() *
                                                           triangle.getReflections().diffuseLight.getGreen() *
                                                           cos)),
                                                         (triangle.getColor().getBlue() +
                                                          (light->diffuseLight.getBlue() *
                                                           triangle.getReflections().diffuseLight.getBlue() *
                                                           cos)));

                                        }
                                        else{
                                            co = Mycolor((co.getRed() +
                                                          (light->diffuseLight.getRed() *
                                                           triangle.getReflections().diffuseLight.getRed() *
                                                           cos)),
                                                         (co.getGreen() +
                                                          (light->diffuseLight.getGreen() *
                                                           triangle.getReflections().diffuseLight.getGreen() *
                                                           cos)),
                                                         (co.getBlue() +
                                                          (light->diffuseLight.getBlue() *
                                                           triangle.getReflections().diffuseLight.getBlue() *
                                                           cos)));
                                        }

                                    } else {
                                        double p = pointLight->spotAngle;
                                        if (cos > std::cos(p)) {
                                            if(teller == 0){
                                                teller++;
                                                co = Mycolor((triangle.getColor().getRed() +
                                                              (light->diffuseLight.getRed() *
                                                               triangle.getReflections().diffuseLight.getRed() *
                                                               (1 - (1 - cos) / (1 - std::cos(p))))),
                                                             (triangle.getColor().getGreen() +
                                                              (light->diffuseLight.getGreen() *
                                                               triangle.getReflections().diffuseLight.getGreen() *
                                                               (1 - (1 - cos) / (1 - std::cos(p))))),
                                                             (triangle.getColor().getBlue() +
                                                              (light->diffuseLight.getBlue() *
                                                               triangle.getReflections().diffuseLight.getBlue() *
                                                               (1 - (1 - cos) / (1 - std::cos(p))))));

                                            }
                                            else{
                                                co = Mycolor((co.getRed() +
                                                              (light->diffuseLight.getRed() *
                                                               triangle.getReflections().diffuseLight.getRed() *
                                                               (1 - (1 - cos) / (1 - std::cos(p))))),
                                                             (co.getGreen() +
                                                              (light->diffuseLight.getGreen() *
                                                               triangle.getReflections().diffuseLight.getGreen() *
                                                               (1 - (1 - cos) / (1 - std::cos(p))))),
                                                             (co.getBlue() +
                                                              (light->diffuseLight.getBlue() *
                                                               triangle.getReflections().diffuseLight.getBlue() *
                                                               (1 - (1 - cos) / (1 - std::cos(p))))));
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (triangle.getReflections().specularLight.getRed() != -1) {
                        if(cos>0) {
                            for (auto *light: triangle.getLights()) {
                                auto *infLight = dynamic_cast<InfLight *>(light);
                                if (infLight != nullptr) {
                                    Vector3D r = 2 * cos * n - l;
                                    r.normalise();
                                    Vector3D cam = -Vector3D::point(x, y, z);
                                    cam.normalise();
                                    double cosb = r.dot(cam);
                                    if (cosb > 0) {
                                        if (teller == 0) {
                                            teller += 2;
                                            col = Mycolor((triangle.getColor().getRed() +
                                                           (light->specularLight.getRed() *
                                                            triangle.getReflections().specularLight.getRed() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (triangle.getColor().getGreen() +
                                                           (light->specularLight.getGreen() *
                                                            triangle.getReflections().specularLight.getGreen() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (triangle.getColor().getBlue() +
                                                           (light->specularLight.getBlue() *
                                                            triangle.getReflections().specularLight.getBlue() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))));

                                        } else if (teller == 1) {
                                            teller++;
                                            col = Mycolor((co.getRed() +
                                                           (light->specularLight.getRed() *
                                                            triangle.getReflections().specularLight.getRed() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (co.getGreen() +
                                                           (light->specularLight.getGreen() *
                                                            triangle.getReflections().specularLight.getGreen() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (co.getBlue() +
                                                           (light->specularLight.getBlue() *
                                                            triangle.getReflections().specularLight.getBlue() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))));

                                        } else {
                                            col = Mycolor((col.getRed() +
                                                           (light->specularLight.getRed() *
                                                            triangle.getReflections().specularLight.getRed() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (col.getGreen() +
                                                           (light->specularLight.getGreen() *
                                                            triangle.getReflections().specularLight.getGreen() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (col.getBlue() +
                                                           (light->specularLight.getBlue() *
                                                            triangle.getReflections().specularLight.getBlue() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))));
                                        }
                                    }
                                }
                            }
                            for (auto *light: triangle.getLights()) {
                                auto *pointLight = dynamic_cast<PointLight *>(light);
                                if (pointLight != nullptr) {
                                    Vector3D r = 2 * cos * n - l;
                                    r.normalise();
                                    Vector3D cam = Vector3D::point(0, 0, 0) - Vector3D::point(x, y, z);
                                    cam.normalise();
                                    double cosb = r.dot(cam);
                                    if (cosb > 0) {
                                        if (teller == 0) {
                                            teller += 2;
                                            col = Mycolor((triangle.getColor().getRed() +
                                                           (light->specularLight.getRed() *
                                                            triangle.getReflections().specularLight.getRed() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (triangle.getColor().getGreen() +
                                                           (light->specularLight.getGreen() *
                                                            triangle.getReflections().specularLight.getGreen() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (triangle.getColor().getBlue() +
                                                           (light->specularLight.getBlue() *
                                                            triangle.getReflections().specularLight.getBlue() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))));

                                        } else if (teller == 1) {
                                            teller++;
                                            col = Mycolor((co.getRed() +
                                                           (light->specularLight.getRed() *
                                                            triangle.getReflections().specularLight.getRed() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (co.getGreen() +
                                                           (light->specularLight.getGreen() *
                                                            triangle.getReflections().specularLight.getGreen() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (co.getBlue() +
                                                           (light->specularLight.getBlue() *
                                                            triangle.getReflections().specularLight.getBlue() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))));

                                        } else {
                                            col = Mycolor((col.getRed() +
                                                           (light->specularLight.getRed() *
                                                            triangle.getReflections().specularLight.getRed() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (col.getGreen() +
                                                           (light->specularLight.getGreen() *
                                                            triangle.getReflections().specularLight.getGreen() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))),
                                                          (col.getBlue() +
                                                           (light->specularLight.getBlue() *
                                                            triangle.getReflections().specularLight.getBlue() *
                                                            pow(cosb, triangle.getReflectionCoefficient()))));
                                        }
                                    }
                                }
                            }
                        }
                    }
                    Z.zbuffer[j][i] = eenOverZ;
                        if(col.getRed() != -1 and col.getBlue() != -1 and col.getGreen() != -1){
                            (Im)(j, i) = img::Color(col.getRed()*255, col.getGreen()*255, col.getBlue()*255);
                        }
                        else if(co.getRed() != -1 and co.getBlue() != -1 and co.getGreen() != -1){
                            (Im)(j, i) = img::Color(co.getRed()*255, co.getGreen()*255, co.getBlue()*255);
                        }
                        else{
                            (Im)(j, i) = c;
                        }

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


    for (auto& i:triangles){
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
            SmallX = i.getCa().getX();
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
    double d = 0.95*(imageX/rangeX);
    for (auto& i:triangles) {
        i.setAa(doProjection(i.getA(), d));
        i.setBa(doProjection(i.getB(), d));
        i.setCa(doProjection(i.getC(), d));
    }

    double DCX = d*((Coordinates[0]+Coordinates[2])/2);
    double DCY = d*((Coordinates[1]+Coordinates[3])/2);
    double dx = imageX/2 - DCX;
    double dy = imageY/2 - DCY;
    for (auto& i:triangles) {
        i.setAa(Point2D(i.getAa().x+dx, i.getAa().y+dy));
        i.setBa(Point2D(i.getBa().x+dx, i.getBa().y+dy));
        i.setCa(Point2D(i.getCa().x+dx, i.getCa().y+dy));
    }
    std::vector<double> color = configuration["General"]["backgroundcolor"];
    img::EasyImage image = img::EasyImage (int(lround(imageX)), int(lround(imageY)), img::Color(int(lround(color[0]*255)),
                                                                                                int(lround(color[1]*255)),int(lround(color[2]*255))));
    ZBuffer ZB = ZBuffer(int(lround(imageX)), int(lround(imageY)));
    for (auto& i:triangles){
        draw_zbuf_triag(ZB, image, i, d, dx, dy);
    }
    return image;
}

Triangle::Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, const Point2D &aa, const Point2D &ba,
                   const Point2D &ca, const Mycolor &c1) : A(a), B(b), C(c), Aa(aa), Ba(ba), Ca(ca), color(c1) {}

Triangle::Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, Mycolor Color) : A(a), B(b), C(c), color(Color){}

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

Triangle::Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, const Point2D &aa, const Point2D &ba,
                   const Point2D &ca, const Lights3D& lights, const Light &reflections) : A(a), B(b), C(c), Aa(aa),
                                                                                          Ba(ba), Ca(ca),
                                                                                          lights(lights),
                                                                                          reflections(reflections) {}

Triangle::Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, const Point2D &aa, const Point2D &ba,
                   const Point2D &ca, const Lights3D& lights, const Light &reflections, double reflectionCoefficient)
        : A(a), B(b), C(c), Aa(aa), Ba(ba), Ca(ca), lights(lights), reflections(reflections),
          reflectionCoefficient(reflectionCoefficient) {}

const Lights3D &Triangle::getLights() const {
    return lights;
}

const Light &Triangle::getReflections() const {
    return reflections;
}

double Triangle::getReflectionCoefficient() const {
    return reflectionCoefficient;
}

void Triangle::setEyePointMatrix(const Matrix &eyePointMatrix) {
    EyePointMatrix = eyePointMatrix;
}

const Matrix &Triangle::getEyePointMatrix() const {
    return EyePointMatrix;
}

void Triangle::setClipping(bool clipping) {
    Triangle::clipping = clipping;
}

bool Triangle::getClipping() const {
    return clipping;
}

Triangle::Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, const Lights3D &lights,
                   const Light &reflections, double reflectionCoefficient) : A(a), B(b), C(c), lights(lights),
                                                                             reflections(reflections),
                                                                             reflectionCoefficient(
                                                                                     reflectionCoefficient) {}

Triangle::Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, const Mycolor &color,
                   const Lights3D &lights, const Light &reflections, double reflectionCoefficient) : A(a), B(b), C(c),
                                                                                                     color(color),
                                                                                                     lights(lights),
                                                                                                     reflections(
                                                                                                             reflections),
                                                                                                     reflectionCoefficient(
                                                                                                             reflectionCoefficient) {}


PointLight::PointLight(const Mycolor &ambientLight, const Mycolor &diffuseLight, const Vector3D &location) : Light(
        ambientLight, diffuseLight), location(location) {
    specularLight = Mycolor(-1,-1,-1);
    spotAngle = -1;
}

PointLight::PointLight(const Mycolor &ambientLight, const Mycolor &diffuseLight, const Mycolor &specularLight,
                       const Vector3D &location) : Light(ambientLight, diffuseLight, specularLight),
                                                   location(location) {
    spotAngle = -1;
}

PointLight::PointLight() {
    spotAngle = -1;
    specularLight = Mycolor(-1,-1,-1);
    ambientLight = Mycolor(-1,-1,-1);
    diffuseLight = Mycolor(-1,-1,-1);
}
