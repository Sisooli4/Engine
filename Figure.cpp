//
// Created by simon on 09.03.22.
//

#include "Figure.h"
#include <bits/stdc++.h>
#include <stack>
#include <cmath>

void toPolar(const Vector3D &point, double &theta, double &phi, double &r) {
    r = sqrt(pow(point.x, 2) + pow(point.y, 2) + pow(point.z, 2));
    theta = std::atan2(point.y, point.x);
    phi = std::acos(point.z / r);
}

Matrix EyePointTrans(const Vector3D &point, bool clipping = false){
double theta;
double phi;
double r;
if(clipping){
    phi = point.x;
    theta = point.y;
    r = point.z;
}
else{
    toPolar(point, theta, phi, r);
}
Matrix ETM = Matrix();
ETM(1, 1) = -sin(theta);
ETM(1, 2) = -cos(theta) * cos(phi);
ETM(1, 3) = cos(theta) * sin(phi);
ETM(2, 1) = cos(theta);
ETM(2, 2) = -sin(theta) * cos(phi);
ETM(2, 3) = sin(theta) * sin(phi);
ETM(3, 2) = sin(phi);
ETM(3, 3) = cos(phi);
ETM(4, 3) = -r;
return ETM;

}

void applyTransformation(Figures3D &figs, const Matrix &m) {
    for (auto& figure: figs) {
        for (auto &point: figure.points) {
            point = point * m;
        }
    }
}


Matrix scaleFigure(const double scale){
    Matrix SF = Matrix();
    SF(1,1) = scale;
    SF(2,2) = scale;
    SF(3,3) = scale;
    return SF;
}

Matrix rotateX(const double angle){
    Matrix RX = Matrix();
    RX(2,2) = cos(angle);
    RX(2,3) = sin(angle);
    RX(3,2) = -sin(angle);
    RX(3,3) = cos(angle);
    return RX;
}

Matrix rotateY(const double angle){
    Matrix RY = Matrix();
    RY(1,1) = cos(angle);
    RY(1,3) = -sin(angle);
    RY(3,1) = sin(angle);
    RY(3,3) = cos(angle);
    return RY;
}

Matrix rotateZ(const double angle){
    Matrix RZ = Matrix();
    RZ(1,1) = cos(angle);
    RZ(1,2) = sin(angle);
    RZ(2,1) = -sin(angle);
    RZ(2,2) = cos(angle);
    return RZ;
}

Matrix translate(const Vector3D &vector){
    Matrix T = Matrix();
    T(4,1) = vector.x;
    T(4,2) = vector.y;
    T(4,3) = vector.z;
    return T;
}

Matrix CreateTransformationMatrix(const double angleX, const double angleY, const double angleZ, const double scale, const Vector3D &point){
    Matrix m1 = scaleFigure(scale);
    Matrix m2 = rotateX(angleX);
    Matrix m3 = rotateY(angleY);
    Matrix m4 = rotateZ(angleZ);
    Matrix m5 = translate(point);

    return m1*m2*m3*m4*m5;
}
Face::Face(const std::vector<int> &pointIndexes) : point_indexes(pointIndexes) {}

Face::Face() {}

Figures::Figures(const ini::Configuration &configuration) {
    std::vector<double> color = configuration["General"]["backgroundcolor"];
    BackGroundcolor = Mycolor(color[0], color[1], color[2]);
    std::vector<double> Eye = configuration["General"]["eye"];
    EyePointMatrix = EyePointTrans(Vector3D::point(Eye[0], Eye[1], Eye[2]));
    if(configuration["General"]["clipping"].exists()) {
        if (configuration["General"]["clipping"].as_bool_or_die()) {
            double theta1;
            double phi1;
            double r1;
            double theta2;
            double phi2;
            double r2;
            std::vector<double> view = configuration["General"]["viewDirection"];
            Vector3D viewDirection = Vector3D::vector(-view[0], -view[1], -view[2]);
            toPolar(viewDirection, theta1, phi1, r1);
            std::vector<double> Eye1 = configuration["General"]["eye"];
            toPolar(Vector3D::point(Eye1[0], Eye1[1], Eye1[2]), theta2, phi2, r2);
            EyePointMatrix = EyePointTrans(Vector3D::point(phi1, theta1, r2), true);
            dNear = configuration["General"]["dNear"].as_double_or_die();
            dFar = configuration["General"]["dFar"].as_double_or_die();
            right = dNear * tan(configuration["General"]["hfov"].as_double_or_die()*M_PI/360);
            left = -right;
            top = right / configuration["General"]["aspectRatio"].as_double_or_die();
            bottom = -top;
        }
    }
    if(!configuration["General"]["nrLights"].exists()){
        for (auto i=0 ; i<int(configuration["General"]["nrFigures"]); i++){
            Figures3D figures3D = {};
            Figure(configuration["Figure" + std::to_string(i)], figures3D, {});
            Figures3D fig = figures3D;
            figures.insert(figures.end(), fig.begin(), fig.end());
    }
    }
    else {
        Lights3D lights = {};
        for (auto i = 0; i < int(configuration["General"]["nrLights"]); i++) {
            if (!configuration["Light" + std::to_string(i)]["infinity"].exists()) {
                auto * l = new Light();
                std::vector<double> light = configuration["Light" + std::to_string(i)]["ambientLight"];
                l->ambientLight = Mycolor(light[0], light[1], light[2]);
                lights.push_back(l);
            }
            else if(!configuration["Light" + std::to_string(i)]["specularLight"].exists()){
                if(configuration["Light" + std::to_string(i)]["infinity"].as_bool_or_die()){
                    auto * l = new InfLight();
                    std::vector<double> alight = configuration["Light" + std::to_string(i)]["ambientLight"];
                    std::vector<double> dlight = configuration["Light" + std::to_string(i)]["diffuseLight"];
                    std::vector<double> direction = configuration["Light" + std::to_string(i)]["direction"];
                    l->ambientLight = Mycolor(alight[0], alight[1], alight[2]);
                    l->diffuseLight = Mycolor(dlight[0], dlight[1], dlight[2]);
                    l->ldVector = Vector3D::vector(direction[0], direction[1], direction[2]);
                    lights.push_back(l);
                }
                else{
                    auto * l = new PointLight();
                    std::vector<double> alight = configuration["Light" + std::to_string(i)]["ambientLight"];
                    std::vector<double> dlight = configuration["Light" + std::to_string(i)]["diffuseLight"];
                    std::vector<double> location = configuration["Light" + std::to_string(i)]["location"];
                    l->ambientLight = Mycolor(alight[0], alight[1], alight[2]);
                    l->diffuseLight = Mycolor(dlight[0], dlight[1], dlight[2]);
                    l->location = Vector3D::point(location[0], location[1], location[2]);
                    if(configuration["Light" + std::to_string(i)]["spotAngle"].exists()){
                        double spotAngle = configuration["Light" + std::to_string(i)]["spotAngle"].as_double_or_die()*M_PI/180;
                        l->spotAngle = spotAngle;
                    }
                    lights.push_back(l);
                }

            }
            else{
                if(configuration["Light" + std::to_string(i)]["infinity"].as_bool_or_die()){
                    auto * l = new InfLight();
                    std::vector<double> alight = configuration["Light" + std::to_string(i)]["ambientLight"];
                    std::vector<double> dlight = configuration["Light" + std::to_string(i)]["diffuseLight"];
                    std::vector<double> slight = configuration["Light" + std::to_string(i)]["specularLight"];
                    std::vector<double> direction = configuration["Light" + std::to_string(i)]["direction"];
                    l->ambientLight = Mycolor(alight[0], alight[1], alight[2]);
                    l->diffuseLight = Mycolor(dlight[0], dlight[1], dlight[2]);
                    l->ldVector = Vector3D::vector(direction[0], direction[1], direction[2]);
                    l->specularLight = Mycolor(slight[0], slight[1], slight[2]);
                    lights.push_back(l);
                }
                else{
                    auto * l = new PointLight();
                    std::vector<double> alight = configuration["Light" + std::to_string(i)]["ambientLight"];
                    std::vector<double> dlight = configuration["Light" + std::to_string(i)]["diffuseLight"];
                    std::vector<double> slight = configuration["Light" + std::to_string(i)]["specularLight"];
                    std::vector<double> location = configuration["Light" + std::to_string(i)]["location"];
                    l->ambientLight = Mycolor(alight[0], alight[1], alight[2]);
                    l->diffuseLight = Mycolor(dlight[0], dlight[1], dlight[2]);
                    l->location = Vector3D::point(location[0], location[1], location[2]);
                    l->specularLight = Mycolor(slight[0], slight[1], slight[2]);
                    if(configuration["Light" + std::to_string(i)]["spotAngle"].exists()){
                        double spotAngle = configuration["Light" + std::to_string(i)]["spotAngle"];
                        l->spotAngle = spotAngle;
                    }
                    lights.push_back(l);
                }

            }

        }
        Lights = lights;
        for (auto i = 0; i < int(configuration["General"]["nrFigures"]); i++) {
            Figures3D figures3D = {};
            Figure(configuration["Figure" + std::to_string(i)], figures3D, lights);
            Figures3D fig = figures3D;
            figures.insert(figures.end(), fig.begin(), fig.end());
        }
    }
}

img::EasyImage Figures::Draw_3Dlines(const ini::Configuration &configuration) {
    std::vector<Line2D*> Lines;
    for (auto& i: figures){
        Matrix CompMatrix = i.TFM*EyePointMatrix;
        for (auto &point: i.points) {
            point = point * CompMatrix;
        }
        for(auto &vlak:i.faces){
            for (int j=0; j<vlak.point_indexes.size()-1; j++){
                double z1 = i.points[vlak.point_indexes[j]].z;
                double z2 = i.points[vlak.point_indexes[j+1]].z;
                Line2D* l = new Line2D(doProjection(i.points[vlak.point_indexes[j]]), doProjection(i.points[vlak.point_indexes[j+1]]), i.color);
                if (std::string(configuration["General"]["type"]) == "2DLSystem" or std::string(configuration["General"]["type"]) == "Wireframe" or std::string(configuration["General"]["type"]) == "ZBufferedWireframe"){
                l->z1 = z1;
                l->z2 = z2;}
                Lines.push_back(l);

            }
            double z1 = i.points[vlak.point_indexes[vlak.point_indexes.size()-1]].z;
            double z2 = i.points[vlak.point_indexes[0]].z;
            Line2D* l = new Line2D(doProjection(i.points[vlak.point_indexes[vlak.point_indexes.size()-1]]), doProjection(i.points[vlak.point_indexes[0]]), i.color);
            if (std::string(configuration["General"]["type"]) == "2DLSystem" or std::string(configuration["General"]["type"]) == "Wireframe" or std::string(configuration["General"]["type"]) == "ZBufferedWireframe"){
            l->z1 = z1;
            l->z2 = z2;
            Lines.push_back(l);}
        }
    }
    if(Lines.empty()){
        return img::EasyImage();
    }
    else{
        return Draw_lines(configuration, Lines);
    }
}

img::EasyImage Figures::Draw_3Dtria(const ini::Configuration &configuration) {
    std::vector<Triangle> Triangles;
    for (auto& i: figures){
        Mycolor c = Mycolor(i.color.getRed(), i.color.getGreen(), i.color.getBlue());
        Matrix CompMatrix = i.TFM*EyePointMatrix;
        for (auto &point: i.points) {
            point = point * CompMatrix;
        }
        for(auto &vlak:i.faces){
            std::vector<Face> tri = triangulate(vlak);
            for (auto& tr:tri){
                Vector3D p1 = Vector3D::point(i.points[tr.point_indexes[0]].x, i.points[tr.point_indexes[0]].y,i.points[tr.point_indexes[0]].z);
                Vector3D p2 = Vector3D::point(i.points[tr.point_indexes[1]].x, i.points[tr.point_indexes[1]].y,i.points[tr.point_indexes[1]].z);
                Vector3D p3 = Vector3D::point(i.points[tr.point_indexes[2]].x, i.points[tr.point_indexes[2]].y,i.points[tr.point_indexes[2]].z);
                Point2D p1a = doProjection(p1);
                Point2D p2a = doProjection(p2);
                Point2D p3a = doProjection(p3);
                if(configuration["General"]["type"].as_string_or_die() != "LightedZBuffering"){
                    Triangle triangle = Triangle(p1, p2, p3, p1a, p2a, p3a, c);
                    Triangles.push_back(triangle);
                }
                else if(!configuration["Figure0"]["diffuseReflection"].exists()){
                    Triangle triangle = Triangle(p1, p2, p3, p1a, p2a, p3a, Lights, i.Reflections);
                    triangle.setC1(c);
                    triangle.setEyePointMatrix(EyePointMatrix);
                    Triangles.push_back(triangle);
                }
                else if(!configuration["Figure0"]["specularReflection"].exists()){
                    Triangle triangle = Triangle(p1, p2, p3, p1a, p2a, p3a, Lights, i.Reflections, i.ReflectionCoefficient);
                    triangle.setC1(c);
                    triangle.setEyePointMatrix(EyePointMatrix);
                    Triangles.push_back(triangle);
                }
                else{
                    Triangle triangle = Triangle(p1, p2, p3, p1a, p2a, p3a, Lights, i.Reflections, i.ReflectionCoefficient);
                    triangle.setC1(c);
                    triangle.setEyePointMatrix(EyePointMatrix);
                    Triangles.push_back(triangle);
                }

            }
        }
        }
    if(configuration["General"]["clipping"].exists()){
        if(configuration["General"]["clipping"]){
            Triangles = clipFigures(Triangles);
        }
    }
    return Draw_tria(configuration, Triangles);
}

std::vector<Triangle> Figures::clipFigures(std::vector<Triangle> triangles) {

    /*
    for(auto& i:triangles){

        std::cout<<"A.x = "<<doProjection(i.getA()).x<<"  A.y = "<<doProjection(i.getA()).y<<"  A.z = "<<i.getA().z<<std::endl;
        std::cout<<"B.x = "<<doProjection(i.getB()).x<<"  B.y = "<<doProjection(i.getB()).y<<"  B.z = "<<i.getB().z<<std::endl;
        std::cout<<"C.x = "<<doProjection(i.getC()).x<<"  C.y = "<<doProjection(i.getC()).y<<"  C.z = "<<i.getC().z<<std::endl;
        std::cout<<std::endl;
    }*/
    Matrix eye = triangles[0].getEyePointMatrix();
    std::vector<double> edges = {-dNear, -dFar, right, left, top, bottom};
    int teller = 0;
    double A;
    double B;
    double C;
    double p;
    for (auto value: edges) {
        std::vector<Triangle> Triangles = {};
        for (auto& i: triangles) {
            if (teller == 0 or teller == 1){
                A = i.getA().z;
                B = i.getB().z;
                C = i.getC().z;
            }
            else if(teller == 2 or teller == 3) {
                A = doProjection(i.getA(), dNear).x;
                B = doProjection(i.getB(), dNear).x;
                C = doProjection(i.getC(), dNear).x;
            }
            else{
                A = doProjection(i.getA(), dNear).y;
                B = doProjection(i.getB(), dNear).y;
                C = doProjection(i.getC(), dNear).y;
            }
            if(teller == 0 or teller == 2 or teller == 4) {
                if (A <= value and B <= value and C <= value) {
                    Triangles.push_back(i);
                } else if (A <= value and B <= value and C > value) {
                    Vector3D A1 = i.getC();
                    Vector3D B1 = i.getA();
                    Vector3D C1 = i.getB();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor()));
                    } else{
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                    } else if (A <= value and B > value and C <= value) {
                    Vector3D A1 = i.getB();
                    Vector3D B1 = i.getC();
                    Vector3D C1 = i.getA();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor()));
                    } else{
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                } else if (A > value and B <= value and C <= value) {
                    Vector3D A1 = i.getA();
                    Vector3D B1 = i.getB();
                    Vector3D C1 = i.getC();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor()));
                    } else{
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                } else if (A <= value and B > value and C > value) {
                    Vector3D A1 = i.getA();
                    Vector3D B1 = i.getB();
                    Vector3D C1 = i.getC();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor()));
                    }
                    else{
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                } else if (A > value and B <= value and C > value) {
                    Vector3D A1 = i.getB();
                    Vector3D B1 = i.getC();
                    Vector3D C1 = i.getA();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor()));
                    }
                    else{
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                } else if (A > value and B > value and C <= value) {
                    Vector3D A1 = i.getC();
                    Vector3D B1 = i.getA();
                    Vector3D C1 = i.getB();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor()));
                    }
                    else{
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                }
            }
            else{
                if (A >= value and B >= value and C >= value) {
                    Triangles.push_back(i);
                } else if (A >= value and B >= value and C < value) {
                    Vector3D A1 = i.getC();
                    Vector3D B1 = i.getA();
                    Vector3D C1 = i.getB();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor()));
                    } else{
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }

                } else if (A >= value and B < value and C >= value) {
                    Vector3D A1 = i.getB();
                    Vector3D B1 = i.getC();
                    Vector3D C1 = i.getA();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor()));
                    } else{
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                } else if (A < value and B >= value and C >= value) {
                    Vector3D A1 = i.getA();
                    Vector3D B1 = i.getB();
                    Vector3D C1 = i.getC();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor()));
                    } else{
                        Triangles.push_back(Triangle(p1, B1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                        Triangles.push_back(Triangle(B1, C1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                } else if (A >= value and B < value and C < value) {
                    Vector3D A1 = i.getA();
                    Vector3D B1 = i.getB();
                    Vector3D C1 = i.getC();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor()));
                    }
                    else{
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                } else if (A < value and B >= value and C < value) {
                    Vector3D A1 = i.getB();
                    Vector3D B1 = i.getC();
                    Vector3D C1 = i.getA();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor()));
                    }
                    else{
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                } else if (A < value and B < value and C >= value) {
                    Vector3D A1 = i.getC();
                    Vector3D B1 = i.getA();
                    Vector3D C1 = i.getB();
                    p = getP(A1, B1, value, dNear, teller);
                    Vector3D p1 = p*A1+(1-p)*B1;
                    p = getP(A1, C1, value, dNear, teller);
                    Vector3D p2 = p*A1+(1-p)*C1;
                    if(i.getLights().empty()){
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor()));
                    }
                    else{
                        Triangles.push_back(Triangle(A1, p1, p2, i.getColor(), i.getLights(), i.getReflections(), i.getReflectionCoefficient()));
                    }
                }
        }
        }
        triangles = Triangles;
        teller += 1;
    }
    for(auto &i:triangles){
        i.setAa(doProjection(i.getA()));
        i.setBa(doProjection(i.getB()));
        i.setCa(doProjection(i.getC()));
        i.setEyePointMatrix(eye);
        i.setClipping(true);
    }
/*
    for(auto& i:triangles){
        std::cout<<"A.x = "<<doProjection(i.getA()).x<<"  A.y = "<<doProjection(i.getA()).y<<"  A.z = "<<i.getA().z<<std::endl;
        std::cout<<"B.x = "<<doProjection(i.getB()).x<<"  B.y = "<<doProjection(i.getB()).y<<"  B.z = "<<i.getB().z<<std::endl;
        std::cout<<"C.x = "<<doProjection(i.getC()).x<<"  C.y = "<<doProjection(i.getC()).y<<"  C.z = "<<i.getC().z<<std::endl;
        std::cout<<std::endl;
    }
*/
    return triangles;
}

double Figures::getP(Vector3D p1, Vector3D p2, double val, double dNear, int teller) {
    double p;
    if (teller == 0 or teller == 1){
        p = (val-p2.z)/(p1.z-p2.z);
    }
    else if(teller == 2 or teller == 3){
        p = (p2.x*dNear+p2.z*val)/((p2.x-p1.x)*dNear+(p2.z-p1.z)*val);
    }
    else{
        p = (p2.y*dNear+p2.z*val)/((p2.y-p1.y)*dNear+(p2.z-p1.z)*val);
    }
    return p;
}

const Lights3D Figures::getLights() const {
    return Lights;
}

const Matrix &Figures::getEyePointMatrix() const {
    return EyePointMatrix;
}

const std::vector<Figure> &Figures::getFigures() const {
    return figures;
}

void Figures::setFigures(const std::vector<Figure> &figures) {
    Figures::figures = figures;
}



Figure::Figure(ini::Section conf, Figures3D &figures3D, Lights3D lights = {}) {
    std::string str = conf["type"].as_string_or_die();
    if(str.substr(0,5) == "Thick"){
        m = conf["m"];
        n = conf["n"];
        r = conf["radius"];
        str = str.substr(5,str.size());
        thick = true;
    }

    if (str == "3DLSystem") {
        L3Dsystem l3D = L3Dsystem(conf["inputfile"]);
        this->points = l3D.pointz;
        this->faces = l3D.facez;

    }
    else if (str == "LineDrawing"){
    for (int i=0 ; i<conf["nrPoints"].as_int_or_die(); i++){
        std::vector<double> p = conf["point" + std::to_string(i)];
        Vector3D point = Vector3D::point(p[0], p[1], p[2]);
        points.push_back(point);
    }
    for (int i=0 ; i<conf["nrLines"].as_int_or_die(); i++){
        std::vector<int> l = conf["line" + std::to_string(i)];
        Face f = Face(l);
        faces.push_back(f);
    }}
    else{
        if (str == "Cube"){
            Figure f = createCube();
            points = f.points;
            faces = f.faces;
        }
        else if(str == "Tetrahedron"){
            Figure f = createTetrahedron();
            points = f.points;
            faces = f.faces;
        }

        else if(str == "Icosahedron"){

            Figure f = createIcosahedron();
            points = f.points;
            faces = f.faces;
        }

        else if(str == "Octahedron"){
            Figure f = createOctahedron();
            points = f.points;
            faces = f.faces;
        }

        else if(str == "Dodecahedron"){
            Figure f = createDodecahedron();
            points = f.points;
            faces = f.faces;
        }

        else if(str == "Cone"){
            Figure f = createCone(conf["n"], conf["height"]);
            points = f.points;
            faces = f.faces;
        }

        else if(str == "Cylinder"){
            Figure f = createCylinder(conf["n"], conf["height"], false);
            points = f.points;
            faces = f.faces;
        }

        else if(str == "Sphere"){
            Figure f = createSphere(conf["n"]);
            points = f.points;
            faces = f.faces;
        }

        else if(str == "Torus"){
            Figure f = createTorus(conf["r"], conf["R"], conf["n"], conf["m"]);
            points = f.points;
            faces = f.faces;
        }

        else if(str == "BuckyBall"){
            Figure f = createBuckyBall();
            points = f.points;
            faces = f.faces;
        }}

        Reflections = Light();
        ReflectionCoefficient = -1.0;
        if (lights.empty()){
                std::vector<double> C = conf["color"];
                color = Mycolor(C[0], C[1], C[2]);
            }
        else{
            Mycolor Color = Mycolor();
            std::vector<double> ambientref = conf["ambientReflection"];
            for(auto &light:lights){
                Color.setRed(Color.getRed()+(light->ambientLight.getRed()*ambientref[0]));
                Color.setGreen(Color.getGreen()+(light->ambientLight.getGreen()*ambientref[1]));
                Color.setBlue(Color.getBlue()+(light->ambientLight.getBlue()*ambientref[2]));
                std::vector<double> ar = conf["ambientReflection"];
                Reflections = Light(Mycolor(ar[0], ar[1], ar[2]));
            }
            color = Color;
            if(conf["diffuseReflection"].exists()){
            if(!conf["specularReflection"].exists()){
                std::vector<double> ar = conf["ambientReflection"];
                std::vector<double> dr = conf["diffuseReflection"];
                Reflections = Light(Mycolor(ar[0], ar[1], ar[2]), Mycolor(dr[0], dr[1], dr[2]));

            }
            else{
                std::vector<double> ar = conf["ambientReflection"];
                std::vector<double> dr = conf["diffuseReflection"];
                std::vector<double> sr = conf["specularReflection"];
                Reflections = Light(Mycolor(ar[0], ar[1], ar[2]), Mycolor(dr[0], dr[1], dr[2]), Mycolor(sr[0], sr[1], sr[2]));
                ReflectionCoefficient = conf["reflectionCoefficient"];
            }
            }
        }
        if(str != "FractalTetrahedron" and str != "FractalCube" and str != "FractalIcosahedron" and str != "FractalOctahedron" and str != "FractalDodecahedron" and str != "FractalBuckyBall" and str != "MengerSponge") {

            std::vector<double> center = conf["center"];
            TFM = CreateTransformationMatrix(conf["rotateX"].as_double_or_die() * M_PI / 180,
                                             conf["rotateY"].as_double_or_die() * M_PI / 180,
                                             conf["rotateZ"].as_double_or_die() * M_PI / 180, conf["scale"],
                                             Vector3D::point(center[0], center[1], center[2]));
            figures3D.push_back(*this);
        }

        else if(str == "MengerSponge"){
            std::vector<Figure> figures;
            createMengerSponge(conf, figures, conf["nrIterations"], createCube(), color, Reflections, ReflectionCoefficient);
            figures3D = figures;
        }

        else if(str == "FractalTetrahedron"){
            Figure f = createTetrahedron();
            std::vector<Figure> figures = {};

            generateFractal(f, figures, conf["nrIterations"].as_int_or_die(), conf["fractalScale"].as_double_or_die(), conf, color, Reflections, ReflectionCoefficient);
            figures3D = figures;
        }
        else if(str == "FractalCube"){
            Figure f = createCube();
            std::vector<Figure> figures = {};
            generateFractal(f, figures, conf["nrIterations"].as_int_or_die(), conf["fractalScale"].as_double_or_die(), conf, color, Reflections, ReflectionCoefficient);
            figures3D = figures;
        }

        else if(str == "FractalIcosahedron"){
            Figure f = createIcosahedron();
            std::vector<Figure> figures = {};
            generateFractal(f, figures, conf["nrIterations"].as_int_or_die(), conf["fractalScale"].as_double_or_die(), conf, color,Reflections, ReflectionCoefficient);
            figures3D = figures;
        }

        else if(str == "FractalOctahedron"){
            Figure f = createOctahedron();
            std::vector<Figure> figures = {};
            generateFractal(f, figures, conf["nrIterations"].as_int_or_die(), conf["fractalScale"].as_double_or_die(), conf, color,Reflections, ReflectionCoefficient);
            figures3D = figures;
        }

        else if(str == "FractalDodecahedron") {
            Figure f = createDodecahedron();
            std::vector<Figure> figures = {};
            generateFractal(f, figures, conf["nrIterations"].as_int_or_die(), conf["fractalScale"].as_double_or_die(), conf, color,Reflections, ReflectionCoefficient);
            figures3D = figures;
        }
        else if(str == "FractalBuckyBall") {
            Figure f = createBuckyBall();
            std::vector<Figure> figures = {};
            generateFractal(f, figures, conf["nrIterations"].as_int_or_die(), conf["fractalScale"].as_double_or_die(), conf, color,Reflections, ReflectionCoefficient);
            figures3D = figures;
        }
        else{
            Figure f = Figure();
            figures3D.push_back(f);
        }
    }

Figure::Figure(const std::vector<Vector3D> &points, const std::vector<Face> &faces) : points(points), faces(faces) {}

Figure::Figure() {}

Figure Figure::createCube() {

    points = std::vector<Vector3D>{Vector3D::point(1,-1,-1),Vector3D::point(-1,1,-1),Vector3D::point(1,1,1),
                                 Vector3D::point(-1,-1,1),Vector3D::point(1,1,-1),Vector3D::point(-1,-1,-1),
                                 Vector3D::point(1,-1,1),Vector3D::point(-1,1,1)};
    faces = std::vector<Face> {Face({0,4,2,6}),Face({4,1,7,2}),Face({1,5,3,7}),
                               Face({5,0,6,3}),Face({6,2,7,3}),Face({0,5,1,4})

    };
    return Figure(points,faces);
}



Figure Figure:: createTetrahedron(){
    points = std::vector<Vector3D>{Vector3D::point(1,-1,-1),Vector3D::point(-1,1,-1),Vector3D::point(1,1,1),
                                   Vector3D::point(-1,-1,1)};
    faces = std::vector<Face> {Face({0,1,2}),Face({1,3,2}),Face({0,3,1}),
                               Face({0,2,3})};
    return Figure(points,faces);
}

Figure Figure:: createIcosahedron(){
    points.push_back(Vector3D::point(0,0, sqrt(5)/2));
    for (int i=2;i<7;i++){
        points.push_back(Vector3D::point(cos((i-2)*2*M_PI/5),sin((i-2)*2*M_PI/5),0.5));
    }
    for (int i = 7; i<12; i++){
        points.push_back(Vector3D::point(cos(M_PI/5+(i-7)*2*M_PI/5),sin(M_PI/5+(i-7)*2*M_PI/5),-0.5));
    }
    points.push_back(Vector3D::point(0,0, -sqrt(5)/2));

    faces = std::vector<Face> {Face({0,1,2}),Face({0,2,3}),Face({0,3,4}),
                               Face({0,4,5}),Face({0,5,1}),Face({1,6,2}),Face({2,6,7}),
            Face({2,7,3}),Face({3,7,8}),Face({3,8,4}),Face({4,8,9}),
            Face({4,9,5}),Face({5,9,10}),Face({5,10,1}),Face({1,10,6}),
            Face({11,7,6}),Face({11,8,7}),Face({11,9,8}),Face({11,10,9}),
                               Face({11,6,10})};
    return Figure(points,faces);
}

Figure Figure:: createOctahedron(){
    points = std::vector<Vector3D>{Vector3D::point(1,0,0),Vector3D::point(0,1,0),Vector3D::point(-1,0,0),
                                   Vector3D::point(0,-1,0),Vector3D::point(0,0,-1),Vector3D::point(0,0,1),};
    faces = std::vector<Face> {Face({0,1,5}),Face({1,2,5}),Face({2,3,5}),
                               Face({3,0,5}),Face({1,0,4}),Face({2,1,4}), Face({3,2,4}),Face({0,3,4})

    };
    return Figure(points,faces);
}

Figure Figure:: createDodecahedron(){
    Figure f = createIcosahedron();
    std::vector<Vector3D> p;
    for (auto vlak:f.faces){
        p.push_back(Vector3D::point((f.points[vlak.point_indexes[0]].x+f.points[vlak.point_indexes[1]].x+f.points[vlak.point_indexes[2]].x)/3,
                                         (f.points[vlak.point_indexes[0]].y+f.points[vlak.point_indexes[1]].y+f.points[vlak.point_indexes[2]].y)/3,
                                         (f.points[vlak.point_indexes[0]].z+f.points[vlak.point_indexes[1]].z+f.points[vlak.point_indexes[2]].z)/3));
    }
    faces = std::vector<Face> {Face({0,1,2,3,4}),Face({0,5,6,7,1}),Face({1,7,8,9,2}),
                               Face({2,9,10,11,3}),Face({3,11,12,13,4}),Face({4,13,14,5,0}),
                               Face({19,18,17,16,15}), Face({19,14,13,12,18}),Face({18,12,11,10,17}),
                               Face({17,10,9,8,16}),Face({16,8,7,6,15}),
                               Face({15,6,5,14,19})};
    return Figure(p,faces);
}

Figure Figure:: createCone(const int n, const double h){
    for (int i=0; i<n; i++){
        points.push_back(Vector3D::point(cos(2*i*M_PI/n), sin(2*i*M_PI/n), 0));
    }
    points.push_back(Vector3D::point(0,0,h));
    for (int i=0;i<n;i++){
        faces.push_back(Face({i, (i+1)%n, n}));
    }
    Face f = Face();
    for(int i=n-1;i>-1;i--){
        f.point_indexes.push_back(i);
    }
    faces.push_back(f);
    return Figure(points,faces);
}

Figure Figure:: createCylinder(const int n, const double h, bool thick){
    for (int i=0; i<n; i++){
        points.push_back(Vector3D::point(cos(2*i*M_PI/n), sin(2*i*M_PI/n), 0));
    }
    for (int i=0; i<n; i++){
        points.push_back(Vector3D::point(cos(2*i*M_PI/n), sin(2*i*M_PI/n), h));
    }
    for (int i=0;i<n-1;i++){
        faces.push_back(Face({i, (i+1)%n, (i+n+1)%(2*n),i+n }));
    }
    faces.push_back(Face({n-1,0,n, 2*n-1}));
    if(!thick){
        for(int i=0; i<2; i++){
            Face f = Face();
            for(int j=i*n;j<(i+1)*n;j++){
                f.point_indexes.push_back(j);
            }
            faces.push_back(f);
        }
    }
    return Figure(points,faces);
}


Figure Figure::createSphere(const int n){
    Figure f = createIcosahedron();
    faces = f.faces;
    points = f.points;
    for (int i = 0; i<n; i++){
        std::vector<Face> Faces;
        for (auto vlak:faces){
            points.push_back(Vector3D::point((points[vlak.point_indexes[0]].x+points[vlak.point_indexes[1]].x)/2,
                                        (points[vlak.point_indexes[0]].y+points[vlak.point_indexes[1]].y)/2,
                                        (points[vlak.point_indexes[0]].z+points[vlak.point_indexes[1]].z)/2));

            points.push_back(Vector3D::point((points[vlak.point_indexes[1]].x+points[vlak.point_indexes[2]].x)/2,
                                          (points[vlak.point_indexes[1]].y+points[vlak.point_indexes[2]].y)/2,
                                          (points[vlak.point_indexes[1]].z+points[vlak.point_indexes[2]].z)/2));

            points.push_back(Vector3D::point((points[vlak.point_indexes[2]].x+points[vlak.point_indexes[0]].x)/2,
                                          (points[vlak.point_indexes[2]].y+points[vlak.point_indexes[0]].y)/2,
                                          (points[vlak.point_indexes[2]].z+points[vlak.point_indexes[0]].z)/2));

            int size = points.size();
            Faces.push_back(Face({vlak.point_indexes[0],size-3,size-1}));
            Faces.push_back( Face({vlak.point_indexes[1], size-2, size-3}));
            Faces.push_back(Face({vlak.point_indexes[2], size-1, size-2}));
            Faces.push_back(Face({size-3, size-2, size-1}));
        }
        faces = Faces;
    }
    std::vector<Vector3D> P;
    for(auto i:points){
        double r = sqrt(pow(i.x,2)+ pow(i.y,2)+ pow(i.z,2));
        P.push_back(Vector3D::point(i.x/r, i.y/r, i.z/r));
    }
    points = P;
    return Figure(points,faces);
}

Figure Figure:: createTorus(const double r, const double R, const int n,const int m){
    for (int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            double u = 2*i*M_PI/n;
            double v = 2*j*M_PI/n;
            points.push_back(Vector3D::point((R+r*cos(v))*cos(u), (R+r*cos(v))*sin(u), r*sin(v)));
        }
    }

    for (int i=0; i<n; i++){
        for(int j=0; j<m; j++) {
            faces.push_back(Face({i*m+j, (i+1)%n*m+j, ((i+1)%(n))*m+(j+1)%m, i*m+(j+1)%m}));
        }
        }
    return Figure(points,faces);
}

Figure Figure::LS3D(std::string inputfile) {
    return Figure();
}

Figure Figure::createBuckyBall() {
    Figure f = createIcosahedron();
    std::vector<Face> Faces;
    std::vector<Vector3D> Points;
        for(auto tria:f.faces) {
            Vector3D p1 = Vector3D::point(points[tria.point_indexes[0]].x, points[tria.point_indexes[0]].y,
                                          points[tria.point_indexes[0]].z);
            Vector3D p2 = Vector3D::point(points[tria.point_indexes[1]].x, points[tria.point_indexes[1]].y,
                                          points[tria.point_indexes[1]].z);
            Vector3D p3 = Vector3D::point(points[tria.point_indexes[2]].x, points[tria.point_indexes[2]].y,
                                          points[tria.point_indexes[2]].z);


            Points.push_back(Vector3D::point(p1.x * 2 / 3 + p2.x / 3,
                                             p1.y * 2 / 3 + p2.y / 3,
                                             p1.z * 2 / 3 + p2.z / 3));

            Points.push_back(Vector3D::point(p1.x * 2 / 3 + p2.x / 3 + p1.x * 2 / 3 + p2.x / 3 - p1.x,
                                             p1.y * 2 / 3 + p2.y / 3 + p1.y * 2 / 3 + p2.y / 3 - p1.y,
                                             p1.z * 2 / 3 + p2.z / 3 + p1.z * 2 / 3 + p2.z / 3 - p1.z));

            Points.push_back(Vector3D::point(p2.x * 2 / 3 + p3.x / 3,
                                             p2.y * 2 / 3 + p3.y / 3,
                                             p2.z * 2 / 3 + p3.z / 3));

            Points.push_back(Vector3D::point(p2.x * 2 / 3 + p3.x / 3 + p2.x * 2 / 3 + p3.x / 3 - p2.x,
                                             p2.y * 2 / 3 + p3.y / 3 + p2.y * 2 / 3 + p3.y / 3 - p2.y,
                                             p2.z * 2 / 3 + p3.z / 3 + p2.z * 2 / 3 + p3.z / 3 - p2.z));

            Points.push_back(Vector3D::point(p3.x * 2 / 3 + p1.x / 3,
                                             p3.y * 2 / 3 + p1.y / 3,
                                             p3.z * 2 / 3 + p1.z / 3));

            Points.push_back(Vector3D::point(p3.x * 2 / 3 + p1.x / 3 + p3.x * 2 / 3 + p1.x / 3 - p3.x,
                                             p3.y * 2 / 3 + p1.y / 3 + p3.y * 2 / 3 + p1.y / 3 - p3.y,
                                             p3.z * 2 / 3 + p1.z / 3 + p3.z * 2 / 3 + p1.z / 3 - p3.z));

            int size = Points.size();

            /*faces.push_back(Face({tria.point_indexes[0],size-6,size-1}));
            faces.push_back(Face({tria.point_indexes[1],size-4,size-5}));
            faces.push_back(Face({tria.point_indexes[2],size-2,size-3}));*/
            Faces.push_back(Face({size - 6, size - 5, size - 4, size - 3, size - 2, size - 1}));
        }
    Faces.push_back(Face({ 0, 6, 12, 18, 24}));
    Faces.push_back(Face({41,8,4,34,33}));
    Faces.push_back(Face({10,9,45,53,14}));
    Faces.push_back(Face({16,15,57,65,20}));
    Faces.push_back(Face({22,21,69,77,26}));
    Faces.push_back(Face({28,27,81,89,2}));
    Faces.push_back(Face({44,43,92,91,99}));
    Faces.push_back(Face({56,55,98,97,105}));
    Faces.push_back(Face({68,67,104,103,111}));
    Faces.push_back(Face({80,79,110,109,117}));
    Faces.push_back(Face({32,31,116,115,93}));
    Faces.push_back(Face({90,96,102,108,114}));
    return Figure(Points, Faces);
}

void applyTransformation(Figure &fig, const Matrix &m) {
    for (auto &point: fig.points) {
        point = point * m;
    }
}

void Figure::applyTransformation(Figure &fig, const Matrix &m) {
    for (auto &point: fig.points) {
        point = point * m;
    }
}

void Figure::createMengerSponge(ini::Section conf,Figures3D& fractal, int nr_iterations, const Figure& fig, Mycolor Color, Light reflex, double reflexcoef) {
    Figure F = fig;
    if (nr_iterations!=0){
        Matrix Ms = scaleFigure(1.0/3.0);
        applyTransformation(F, Ms);
        int teller = 0;
        for(auto i = 0; i<(double)F.points.size(); i ++) {
            if (teller == 0){
                for (int j=0; j<4;j++){
                    for (int k=0; k<3; k++){
                        Figure Fig = F;
                        Matrix Mt = translate((fig.points[faces[j].point_indexes[k]]+fig.points[faces[j].point_indexes[k+1]])/2.0 - (Fig.points[faces[j].point_indexes[k]]+Fig.points[faces[j].point_indexes[k+1]])/2.0);
                        applyTransformation(Fig, Mt);
                        teller = 1;
                        createMengerSponge(conf,fractal,nr_iterations-1,Fig, color, reflex, reflexcoef);
                    }
                }
            }
            Figure Fi = F;
            Matrix Mt = translate(fig.points[i] - Fi.points[i]);
            applyTransformation(Fi, Mt);
            createMengerSponge(conf, fractal, nr_iterations - 1, Fi, color, reflex, reflexcoef);
        }
    }
    else{
        F.color = color;

        std::vector<double> center = conf["center"];
        F.TFM = CreateTransformationMatrix(conf["rotateX"].as_double_or_die()*M_PI/180, conf["rotateY"].as_double_or_die()*M_PI/180, conf["rotateZ"].as_double_or_die()*M_PI/180, conf["scale"], Vector3D::point(center[0], center[1], center[2]));

        F.Reflections = reflex;
        F.ReflectionCoefficient = reflexcoef;
        fractal.push_back(F);
        return;
    }
}



void Figures::generateThickFigures() {
    std::vector<Figure> figs;
    Figures3D resultingFigures;
    for (auto &i: figures){
        if(i.thick){
            i.generateThickFigure(resultingFigures);
            for (auto& fa:resultingFigures){
                fa.Reflections = i.Reflections;
                fa.ReflectionCoefficient = i.ReflectionCoefficient;
                figs.push_back(fa);
        }
        }
        else{
            figs.push_back(i);
        }
    }
    figures = figs;
}

float Figures::roundOff(float value, unsigned char prec) {
    float pow_10 = pow(10.0f, (float)prec);
    return round(value * pow_10) / pow_10;
}

void Figure::generateThickFigure(Figures3D &resultingFigures) {
    for(auto& point:points){
        Figure f;
        f.createSphere(m);
        Matrix Ms = scaleFigure(r);
        applyTransformation(f,Ms);
        Matrix Mt = translate(Vector3D::vector(point.x, point.y, point.z));
        applyTransformation(f,Mt);
        f.color = color;
        resultingFigures.push_back(f);
    }
    for(auto& i:faces){
        for (int j=0; j<i.point_indexes.size()-1; j++){
            Figure fig;
            fig.createCylinder(n, (points[i.point_indexes[j+1]] - points[i.point_indexes[j]]).length()/r, true);
            Matrix Ms = scaleFigure(r);
            applyTransformation(fig,Ms);
            Vector3D Pr = points[i.point_indexes[j+1]]-points[i.point_indexes[j]];
            double phi;
            double theta;
            double r;
            toPolar(Pr, theta, phi, r);
            Matrix Y = rotateY(phi);
            Matrix Z = rotateZ(theta);
            Matrix Mt = translate(points[i.point_indexes[j]]);
            applyTransformation(fig, Y*Z*Mt);
            fig.color = color;
            resultingFigures.push_back(fig);
        }
        Figure fig;
        fig.createCylinder(n, (points[i.point_indexes[0]] - points[i.point_indexes[i.point_indexes.size()-1]]).length()/r, true);
        Matrix Ms = scaleFigure(r);
        applyTransformation(fig,Ms);
        Vector3D Pr = points[i.point_indexes[0]]-points[i.point_indexes[i.point_indexes.size()-1]];
        double phi;
        double theta;
        double r;
        toPolar(Pr, theta, phi, r);
        Matrix Y = rotateY(phi);
        Matrix Z = rotateZ(theta);
        Matrix Mt = translate(points[i.point_indexes[i.point_indexes.size()-1]]);
        applyTransformation(fig, Y*Z*Mt);
        fig.color = color;
        resultingFigures.push_back(fig);
        }

}


L3Dsystem::L3Dsystem(std::string inputfile) {
    LParser::LSystem3D l_system ;
    std::ifstream input_stream(inputfile);
    input_stream >> l_system;
    alphabet = l_system.get_alphabet();
    for (auto i:alphabet){
        replacement[i] =  l_system.get_replacement(i);
    }
    angle = l_system.get_angle()*M_PI/180;
    initiator = l_system.get_initiator();
    iterations = l_system.get_nr_iterations();
    for (auto i:alphabet){
        draw[i] = l_system.draw(i);
    }
    DrawString = replaceString(initiator, iterations);
    std::pair<std::vector<Vector3D>, std::vector<Face>> p_f = getP_F();
    pointz = p_f.first;
    facez = p_f.second;
    input_stream.close();
}

std::string L3Dsystem::replaceString(std::string string, unsigned int iteration) {
    std::string str;
    if (iteration >= 0){
        iteration -= 1;
        for(auto i:string){
            if (alphabet.find(i) != alphabet.end()){
                str += replacement.at(i);
            }
            else{
                str += i;
            }
        }
        if(iteration > 0){
            str = replaceString(str,iteration);
        }
    }
    return str;
}

std::pair<std::vector<Vector3D>, std::vector<Face>> L3Dsystem::getP_F() {

    std::vector<Vector3D> p;
    std::vector<Face> f;
Vector3D H = Vector3D::vector(1,0,0);
Vector3D L = Vector3D::vector(0,1,0);
Vector3D U = Vector3D::vector(0,0,1);

Vector3D Point = Vector3D::point(0,0,0);
std::stack<std::pair<std::vector<Vector3D>, int>> stack;
p.push_back(Point);
int last_index = 0;
int teller = 0;
    for (auto i:DrawString){
        if (alphabet.find(i) != alphabet.end()){
            if(draw[i]){
                if(alphabet.find(DrawString[teller+1]) != alphabet.end()){
                    if(draw[DrawString[teller+1]]){
                        Point = Point + H;
                        teller++;
                        continue;
                    }
                }
            }
            Point = Point + H;
            p.push_back(Point);

            if (draw[i]){
                f.push_back(Face({last_index, int(p.size()-1)}));
            }
            last_index =  int(p.size()-1);
            teller++;
            }
        else{
            if (i == '+'){
                Vector3D H2 = H;
                H = H*cos(angle)+L* sin(angle);
                L = -H2* sin(angle) + L* cos(angle);
            }
            else if(i == '-'){
                Vector3D H2 = H;
                H = H*cos(-angle)+L* sin(-angle);
                L = -H2* sin(-angle) + L* cos(-angle);
            }
            else if(i == '^'){
                Vector3D H2 = H;
                H = H*cos(angle)+U* sin(angle);
                U = -H2* sin(angle) + U* cos(angle);
            }
            else if(i == '&'){
                Vector3D H2 = H;
                H = H*cos(-angle)+U* sin(-angle);
                U = -H2* sin(-angle) + U* cos(-angle);
            }
            else if(i == '\\'){
                Vector3D L2 = L;
                L = L*cos(angle)-U* sin(angle);
                U = L2* sin(angle) + U* cos(angle);
            }
            else if(i == '/'){
                Vector3D L2 = L;
                L = L*cos(-angle)-U* sin(-angle);
                U = L2* sin(-angle) + U* cos(-angle);
            }
            else if(i == '|'){
                H = -H;
                L = -L;
            }
            else if (i == '('){
                stack.push(std::pair<std::vector<Vector3D>, int >({H,L,U,Point}, last_index));
            }

            else if (i == ')'){
                std::pair<std::vector<Vector3D>,int> top;
                top = stack.top();
                stack.pop();
                H = top.first[0];
                L = top.first[1];
                U = top.first[2];
                Point = top.first[3];
                last_index = top.second;
            }
            teller++;
        }
    }
    return std::pair<std::vector<Vector3D>, std::vector<Face>> (p,f);
}

std::vector<Face> triangulate(const Face& face){
    std::vector<Face> F;
    if(face.point_indexes.size()>3){
        for(int i=1; i<=face.point_indexes.size()-2;i++){
            F.push_back(Face({face.point_indexes[0], face.point_indexes[i], face.point_indexes[i+1]}));
        }
    }
    else{
        F.push_back(face);
    }
    return F;
}

void generateFractal(Figure& fig, Figures3D& fractal, const int nr_iterations, const double scale, ini::Section conf, Mycolor color, Light reflex, double reflexcoef){
    Figure F = fig;
    if (nr_iterations!=0){
        Matrix Ms = scaleFigure(1/scale);
        applyTransformation(F, Ms);
    for(int i = 0; i<fig.points.size(); i++) {
        Figure Fi = F;
        Matrix Mt = translate(fig.points[i] - Fi.points[i]);
        applyTransformation(Fi, Mt);
        generateFractal(Fi, fractal, nr_iterations - 1, scale, conf, color, reflex,reflexcoef);
    }
    }
    else{
        F.color = color;

        std::vector<double> center = conf["center"];
        F.TFM = CreateTransformationMatrix(conf["rotateX"].as_double_or_die()*M_PI/180, conf["rotateY"].as_double_or_die()*M_PI/180, conf["rotateZ"].as_double_or_die()*M_PI/180, conf["scale"], Vector3D::point(center[0], center[1], center[2]));

        F.Reflections = reflex;
        F.ReflectionCoefficient = reflexcoef;
        fractal.push_back(F);
    }
}



