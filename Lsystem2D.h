//
// Created by simon on 06.03.22.
//

#ifndef ENGINE_LSYSTEM2D_H
#define ENGINE_LSYSTEM2D_H
#include "easy_image.h"
#include "ini_configuration.h"
#include <cmath>
#include "Line2D.h"
#include <vector>
#include "l_parser.h"
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

class Line2D;
class Lsystem2D {
    std::set<char> alphabet;
    std::map<char, std::string> replacement;
    double angle;
    std::string initiator;
    unsigned int iterations;
    double starting_angle;
    std::map<char, bool> draw;
    std::string DrawString;
    Mycolor BackGroundcolor;
    Mycolor LineColor;
    std::vector<Line2D*> lines;

public:
    Lsystem2D(const ini::Configuration &configuration);

    img::EasyImage Draw_Lsystem(const ini::Configuration &configuration);

    std::string replaceString(std::string string, int iterations);

    virtual ~Lsystem2D();

};


#endif //ENGINE_LSYSTEM2D_H
