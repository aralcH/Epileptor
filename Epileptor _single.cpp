#include <iostream>
#include <stdlib.h>
#include <chrono>
#include <algorithm>
#include <random>
#include <fstream>

using namespace std;


int main()
{
    /** Variable declaration **/
    int i, k, Nsteps=1000000, n, N=10000;
    double x1[2], x2[2], y1[2], y2[2], z[2], xn[N];
    double  h, i1, i2, x0, raizh, g, t0, desv;

    /** Random number generation **/
    desv=0.5;
    mt19937 e(chrono::steady_clock::now().time_since_epoch().count());
    normal_distribution<double> gauss(0.0, desv);

    /** Creating/Opening files **/
    ofstream Psi("Psi.txt");
    ofstream Psi1("Psi1.txt");
    ofstream Psi2("Psi2.txt");
    ofstream y_1("y_1.txt");
    ofstream y_2("y_2.txt");
    ofstream z_z("z_z.txt");
    ofstream X1("x1.txt");
    ofstream X2("x2.txt");
    ofstream Z("z.txt");
    ofstream Y1("y1.txt");
    ofstream Y2("y2.txt");
    ofstream tri("tri_05.txt");

    /************* Value initialization **************/

    x1[0]=-1.8;
    y1[0]=-15.0;
    x2[0]=-1.0;
    y2[0]=0.0;
    z[0]=3.0;
    i1=3.1;
    i2=0.45;
    x0=-1.6;
    for (i=0; i<N; i++)
    {
        xn[i]=x1[0];
    }
    t0=0;
    h=1.0/256.0;
    raizh=sqrt(h);

    /*************************************************/

    /******************** NUMERICAL INTEGRATION **********************/
    for (k=0; k<Nsteps; k++)
    {
        // g function integration
        g=0;
        for (n=0; n<N; n++) g=g+h*exp(-0.01*h*(N-n))*xn[n];

        // Integration of the equations

        x1[1]=x1[0]+h*(y1[0]-z[0]+i1)+raizh*gauss(e);
  
        if (x1[0]<0) x1[1]=x1[1]-h*(x1[0]*x1[0]*x1[0]-3*x1[0]*x1[0]);
        else x1[1]=x1[1]-h*(x2[0]-0.6*(z[0]-4)*(z[0]-4))*x1[0];


        y1[1]=y1[0]+h*(1-5*x1[0]*x1[0]-y1[0])+raizh*gauss(e);

        x2[1]=x2[0]+h*(-y2[0]+x2[0]-x2[0]*x2[0]*x2[0]+i2+0.002*g-0.3*(z[0]-3.5))+raizh*gauss(e);

        y2[1] = y2[0]-h*(0.1*y2[0])+raizh*gauss(e);
        if (x2[0]>=-0.25) y2[1]=y2[1]+h*0.6*(x2[0]+0.25);

        z[1]=z[0]+h*0.00035*(4*(x1[0]-x0)-z[0])+raizh*gauss(e);
        if (z[0]<0) z[1]=z[1]-h*(0.00035*4*0.1*pow(z[0],7));

        // xn vector actualization
        for (i=0; i<N-1; i++) xn[i+1]=xn[i];
        xn[0]=x1[1];

        /** Writting current equation value in files **/
        Psi << k*h << "    " << -x1[1]+x2[1] << endl;
        Psi1 << k*h << "    " << x1[1] << endl;
        Psi2 << k*h << "    " << x2[1] << endl;
        y_1 << k * h << "    " << y1[1] << endl;
        y_2 << k * h << "    " << y2[1] << endl;
        z_z << k * h << "    " << z[1] << endl;

        X1 << k * h << "    " << x1[1] << endl;
        X2 << k * h << "    " << x2[1] << endl;
        Y1 << k*h << "    " << y1[1] << endl;
        Y2 << k*h << "    " << y2[1] << endl;
        Z << k * h << "    " << z[1] << endl;

        tri << -x1[1]+x2[1] << "   " << z[1] << "   " << y1[1] << endl;


        /** Initial value actualization **/
        x1[0]=x1[1];
        y1[0]=y1[1];
        x2[0]=x2[1];
        y2[0]=y2[1];
        z[0]=z[1];

    }

    /*****************************************************************/

    return 0;
}




