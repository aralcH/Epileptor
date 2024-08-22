#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <random>
#include <ctime>

using namespace std;

int main()
{
    int N=90, i;
    double x0, dx0, x1, x2, y1, y2, z;

    //x0=-2.0;
    //dx0=-1.5;
    x1=-1.0;
    y1=-10.0;
    x2=-0.5;
    y2=0.0;
    z=3.0;

    ofstream fileX1("x1_ini.txt");
    ofstream fileX2("x2_ini.txt");
    ofstream fileY1("y1_ini.txt");
    ofstream fileY2("y2_ini.txt");
    ofstream fileZ ("z_ini.txt");
    ofstream fileX0 ("x0.txt");

    default_random_engine gen(time(0));
    normal_distribution<double> gauss(0.0, 0.5);
    uniform_real_distribution<double> distribution(-2.4, -2.0);

    for (i=0; i<N; i++)
    {
        fileX1 << gauss(gen)+x1 << endl;
        fileY1 << gauss(gen)+y1 << endl;
        fileX2 << gauss(gen)+x2 << endl;
        fileY2 << gauss(gen)+y2 << endl;
        fileZ << gauss(gen)+z << endl;
        fileX0 << distribution(gen) << endl;
    }


    return 0;
}
