//
// Created by simon on 06.03.22.
//

#include "Lsystem2D.h"
#include <stack>

Lsystem2D::Lsystem2D(const ini::Configuration &configuration) {
    std::string FLSystem = configuration["2DLSystem"]["inputfile"];
    std::vector<double> color = configuration["General"]["backgroundcolor"];
    BackGroundcolor = Mycolor(color[0], color[1], color[2]);
    color = configuration["2DLSystem"]["color"];
    LineColor = Mycolor(color[0], color[1], color[2]);
    LParser::LSystem2D l_system ;

    std::ifstream input_stream(FLSystem);
    input_stream >> l_system;
    alphabet = l_system.get_alphabet();
    for (auto i:alphabet){
        replacement[i] =  l_system.get_replacement(i);
    }
    angle = l_system.get_angle()*M_PI/180;
    initiator = l_system.get_initiator();
    iterations = l_system.get_nr_iterations();
    starting_angle = l_system.get_starting_angle()*M_PI/180;
    for (auto i:alphabet){
        draw[i] = l_system.draw(i);
    }
    input_stream.close();
}

img::EasyImage Lsystem2D::Draw_Lsystem(const ini::Configuration &configuration) {
    double CurrentAngle = starting_angle;
    std::stack<std::pair<Point2D, double>> stack;
    Point2D p1 = Point2D(0.0,0.0);
    DrawString = replaceString(initiator, iterations);
    for (auto i:DrawString){
        if (alphabet.find(i) != alphabet.end()){
            Point2D p2 = Point2D(p1.x+std::cos(CurrentAngle), p1.y+std::sin(CurrentAngle));
            if(draw[i]){
                lines.push_back(new Line2D(p1,p2,LineColor));
                p1 = p2;
            }
            else{
                p1 = p2;
            }
        }
        else if (i == '+'){
            CurrentAngle += angle;
        }
        else if(i == '-'){
            CurrentAngle -= angle;
        }
        else if (i == '('){
            stack.push(std::pair<Point2D,double>(p1, CurrentAngle));
        }
        else if (i == ')'){
            std::pair<Point2D, double> update = stack.top();
            stack.pop();
            p1 = update.first;
            CurrentAngle = update.second;
        }
    }
    return Draw_lines(configuration, lines);;
}

Lsystem2D::~Lsystem2D() {
    for (auto i:lines){
        delete i;
    }
}


std::string Lsystem2D::replaceString(std::string string, int iterations) {
    std::string str;
    if (iterations > -1){
        iterations -= 1;
    for(auto i:string){
        if (alphabet.find(i) != alphabet.end()){
                str += replacement.at(i);
        }
        else{
            str += i;
        }
    }
    if(iterations > 0){
        str = replaceString(str,iterations);
    }
    }
    return str;
}
