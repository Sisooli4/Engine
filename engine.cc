#include "easy_image.h"
#include "ini_configuration.h"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include "Lsystem2D.h"
#include "vector3d.h"
#include "Figure.h"
#include <chrono>

img::EasyImage generate_image(const ini::Configuration &configuration)
{
    std::string str = configuration["General"]["type"];
    if (str == "2DLSystem"){
        Lsystem2D l = Lsystem2D(configuration);

        return l.Draw_Lsystem(configuration);
    }
    else if(str == "Wireframe"){
            Figures f = Figures(configuration);

            f.generateThickFigures();

            return f.Draw_3Dlines(configuration);
        }
    else if(str == "ZBufferedWireframe"){
        Figures f = Figures(configuration);

        f.generateThickFigures();

        return f.Draw_3Dlines(configuration);
    }
    else if(str == "ZBuffering"){
        Figures f = Figures(configuration);

        f.generateThickFigures();

        return f.Draw_3Dtria(configuration);
    }

    else if(str == "LightedZBuffering"){
        Figures f = Figures(configuration);

        f.generateThickFigures();

        return f.Draw_3Dtria(configuration);
    }

    else{
        return img::EasyImage();
    }

}

int main(int argc, char const* argv[])
{
    int retVal = 0;
    auto start = std::chrono::high_resolution_clock::now();
    auto timeAll = start - start;
    try
    {
        std::vector<std::string> args =std::vector<std::string>(argv+1, argv+argc);
        if (args.empty()) {
            std::ifstream fileIn("filelist");
            std::string filelistName;
            while (std::getline(fileIn, filelistName)) {
                args.push_back(filelistName);
            }
        }
        for(std::string fileName : args)
        {
            ini::Configuration conf;
            try
            {
                std::ifstream fin(fileName);
                fin >> conf;
                fin.close();
            }
            catch(ini::ParseException& ex)
            {
                std::cerr << "Error parsing file: " << fileName << ": " << ex.what() << std::endl;
                retVal = 1;
                continue;
            }
            auto start = std::chrono::high_resolution_clock::now();
            img::EasyImage image = generate_image(conf);
            auto end = std::chrono::high_resolution_clock::now();
            auto time = end - start;
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(time);
            timeAll += ms;
            std::cout<<fileName<<"  "<<ms.count()/1000.0<<std::endl;
            if(image.get_height() > 0 && image.get_width() > 0)
            {
                std::string::size_type pos = fileName.rfind('.');
                if(pos == std::string::npos)
                {
                    //filename does not contain a '.' --> append a '.bmp' suffix
                    fileName += ".bmp";
                }
                else
                {
                    fileName = fileName.substr(0,pos) + ".bmp";
                }
                try
                {
                    std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                    f_out << image;

                }
                catch(std::exception& ex)
                {
                    std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                    retVal = 1;
                }
            }
            else
            {
                std::cout << "Could not generate image for " << fileName << std::endl;
            }
        }
    }
    catch(const std::bad_alloc &exception)
    {
        //When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
        //Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
        //(Unless of course you are already consuming the maximum allowed amount of memory)
        //If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
        //mark the test as failed while in reality it just needed a bit more memory
        std::cerr << "Error: insufficient memory" << std::endl;
        retVal = 100;
    }
    auto msAll = std::chrono::duration_cast<std::chrono::milliseconds>(timeAll);
    std::cout<<msAll.count()/1000.0<<std::endl;
    return retVal;
}
