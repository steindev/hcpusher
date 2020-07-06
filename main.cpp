/**
 * A (single) particle simulation to test different pushers.
 **/

#include <iostream>
#include <cmath>
#include "vector_ops.hpp"
#include "pusher.hpp"

using namespace std;
using Vec3D = vector<double>;
int main ()
{
    double m = 9.1e-31; // Mass in kg
    double q = -1.6022e-19; // Charge in As
    double c = 299792458; // Speed of light in m/s

    cout << "m*c/q = " << (m*c)/q << endl << endl;

    // Starting conditions
    Vec3D x_0 = {0., 0., 0.};
    Vec3D x_i = x_0;

    double beta = 0.9999;
    Vec3D u_vector_0 = { 0. , beta * c / sqrt(1. - beta*beta) , 0. }; // u = p/m = gamma*v, u bei 0
    Vec3D u_vector_i = u_vector_0; // p/m
    Vec3D u_vector_i_by_c = scalmultip(u_vector_i,1./c);
    double gamma_i = sqrt(1.+ innerprod(u_vector_i_by_c,u_vector_i_by_c));


    // Define fields
    // B-Field
    Vec3D B = {0., 0., 1.};
    cout << "B = " << vec_print(B) << "T\n";
    // E-Field
    Vec3D E = { -beta * c * B[2], 0., 0.};
    cout << "E = " << vec_print(E) << "V/m\n";

    cout << "Force / q = " << vec_print( vec_add( E , crossprod( scalmultip( u_vector_0 , 1. / gamma_i ) , B ) ) ) << "V/m\n";

    // Simulation defines
    double Dt = 1.; // Time step in s 
    //double Dx = c * Dt * sqrt(3.) / 0.995; // Spatial step in m
    int n = 10; // Number of time steps to simulate
    int output_step = 1;

    // Choose pusher
    // Pusher< HC >: Higuera-Cary pusher
    // Pusher< Boris >: Boris pusher
    Pusher< HC > push(Dt, q, m, c);

    // Initial diagnostic print
    cout << "Dt = " << Dt << "s\n";
    //cout << "Dx = " << Dx << "m\n";
    cout << "u(t=0) = " << vec_print(u_vector_0) << "            ";
    cout << "x(t=0) = " << vec_print(x_0) << "            ";
    cout << "gamma(t=0) = " << gamma_i << endl << endl;

    // Calculate particle trajectory
    for (int i=0; i < n; i++ )
    {
        // First half of position update
        Vec3D x_first_half = vec_add( x_i, (scalmultip(u_vector_i, Dt/(2.*gamma_i))));

            //cout << "x(i+1/2) = " << vec_print(x_first_half) << endl;

        // Solution for new momentum u = gamma * v = p/m
        u_vector_i = push(u_vector_i, E, B);

        // Second half of position update
        u_vector_i_by_c = scalmultip(u_vector_i,1./c);

        gamma_i = sqrt(1.+ innerprod(u_vector_i_by_c,u_vector_i_by_c));

            //cout << "gamma_i = " << gamma_i << endl;

        x_i = vec_add(x_first_half , scalmultip(u_vector_i, Dt/(2.*gamma_i)));

        // Output
        if ( (i+1) % output_step == 0 )
        {
            cout << "u(t=" << (i+1)*Dt << ") = " << vec_print(u_vector_i) << "            ";
        
            cout << "x(t=" << (i+1)*Dt << ") = " << vec_print(x_i) << "            ";
            
            cout << "gamma = " << gamma_i << "            ";
            
            cout << "gamma*beta_x = " <<  u_vector_i_by_c[0] << "\n";
        }
    }
    
}

