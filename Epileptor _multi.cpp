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
    int i, k, Nsteps= 500000, n, N=5000, M=90, m, Fsamp;
    double x1[M][2], x2[M][2], y1[M][2], y2[M][2], z[M][2];
    double** xn = new double*[M];
    for (i=0; i<M; i++) xn[i] = new double[N];
    double x0[M], K[M][M];
    double  h, i1, i2, raizh, g, t0, desv, S;

    /** Random number generation **/
    desv=0.05;
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
    ofstream Y1("y1.txt");
    ofstream Y2 ("y2.txt");
    ofstream Z("z.txt");

    ifstream x1_i("x1_ini.txt");
    ifstream x2_i("x2_ini.txt");
    ifstream z_i("z_ini.txt");
    ifstream y1_i("y1_ini.txt");
    ifstream y2_i("y2_ini.txt");
    ifstream x0_i("x0.txt");
    ifstream conectoma("conectoma.txt");

 /************* Value initialization **************/
    for (m=0; m<M; m++)
    {
        x1_i >> x1[m][0];
        y1_i >> y1[m][0];
        x2_i >> x2[m][0];
        y2_i >> y2[m][0];
        z_i >> z[m][0];
        x0_i >> x0[m];
        for (i=0; i<N; i++) xn[m][i]=x1[m][0];
        for (i=0; i<M; i++) conectoma >> K[m][i];
    }

    i1=3.1;
    i2=0.45;
    t0=0;

    Fsamp=100;
    h=0.01;
    raizh=sqrt(h);

 /*************************************************/


 /******************** NUMERICAL INTEGRATION **********************/
    for (k=0; k<Nsteps; k++)
    {

        // Writting time in files
        if (k%Fsamp==0)
        {
            Psi << k*h << "    " ;
            Psi1 << k*h << "    " ;
            Psi2 << k*h << "    " ;
            y_1 << k*h << "    " ;
            y_2 << k*h << "    " ;
            z_z << k*h << "    " ;
        }

     /******************** Epileptor Loops ********************/
        for (m=0; m<M; m++)
        {
            // g function integration
            g=0;
            for (n=0; n<N; n++) g=g+h*exp(-h*(N-n)/100)*xn[m][n];

            // Integration of the equations

            x1[m][1]=x1[m][0]+h*(y1[m][0]-z[m][0]+i1);
            if (x1[m][0]<0) x1[m][1]=x1[m][1]-h*(x1[m][0]*x1[m][0]*x1[m][0]-3*x1[m][0]*x1[m][0]);
            else x1[m][1]=x1[m][1]-h*(x2[m][0]-0.6*(z[m][0]-4)*(z[m][0]-4))*x1[m][0];


            y1[m][1]=y1[m][0]+h*(1-5*x1[m][0]*x1[m][0]-y1[m][0]);

            x2[m][1]=x2[m][0]+h*(-y2[m][0]+x2[m][0]-x2[m][0]*x2[m][0]*x2[m][0]+i2+g/500-(z[m][0]-3.5)/3)+raizh*gauss(e);

            y2[m][1] = y2[m][0]-h*(y2[m][0])/10+raizh*gauss(e);
            if (x2[m][0]>=-0.25) y2[m][1]=y2[m][1]+h*2*(x2[m][0]+0.25)/3;

            S=0;
            for (i=0; i<M; i++) S+=K[m][i]*(x1[i][0]-x1[m][0]);

            z[m][1]=z[m][0]+h*(4*(x1[m][0]-x0[m])-z[m][0]-S)/2857;
            z[m][1]+=raizh*gauss(e);

            // xn vector actualization
            for (i=0; i<N-1; i++) xn[m][i+1]=xn[m][i];
            xn[m][0]=x1[m][1];

            /** Writting current equation value in files **/
            if (k%Fsamp==0)
            {
                // Sequential files
                Psi << -x1[m][1]+x2[m][1] << "    " ;
                Psi1 << x1[m][1] << "    " ;
                Psi2 << x2[m][1] << "    " ;
                y_1 << y1[m][1] << "    " ;
                y_2 << y2[m][1] << "    " ;
                z_z << z[m][1] << "    " ;

                // Heatmaps
                X1 << k*h <<"    "<< m+1 << "    " << x1[m][1] << endl;
                X2 << k*h <<"    "<< m+1 << "    " << x2[m][1] << endl;
                Y1 << k*h <<"    "<< m+1 << "    " << y1[m][1] << endl;
                Y2 << k*h <<"    "<< m+1 << "    " << y2[m][1] << endl;
                Z << k*h <<"    "<< m+1 << "    " << z[m][1] << endl;

            }

        }
     /*********************************************************/

        //Adding a line break after the measurement
        if (k%Fsamp==0)
        {
            Psi << endl;
            Psi1 << endl;
            Psi2 << endl;
            y_1 << endl;
            y_2 << endl;
            z_z << endl;

            X1 << endl;
            X2 << endl;
            Y1 << endl;
            Y2 << endl;
            Z << endl;
        }

        /** Initial value actualization **/
        for (m=0; m<M; m++)
        {   x1[m][0]=x1[m][1];
            y1[m][0]=y1[m][1];
            x2[m][0]=x2[m][1];
            y2[m][0]=y2[m][1];
            z[m][0]=z[m][1];
        }
    }
 /*****************************************************************/


 /***************** DELETING THE POINTER VALUES *******************/
    for (i=0; i<M; i++) delete[] xn[i];
    delete[] xn;
 /*****************************************************************/

    return 0;
}



