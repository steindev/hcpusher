/**
* Header file implementing different pusher methods for a particle simulation.
**/

#pragma once

#include <iostream>
#include <cmath>
#include "vector_ops.hpp"

// Generic functor for different pusher implementations via template specializations
template < class PusherType >
struct Pusher
{
    Pusher(const float& Dt, const float& q, const float& m, const float& c) : 
        Dt(Dt), q(q), m(m), c(c)
    {}

    Vec3D operator ()(const Vec3D& u_vector_i, const Vec3D& E, const Vec3D& B);

  private:
    const float Dt;
    const float q;
    const float m;
    const float c;
};

/**
 * Implementation of the Higuera-Cary pusher as presented in doi:10.1063/1.4979989.
 * A correction is applied to the given formulas as documented by the WarpX team:
 * (https://github.com/ECP-WarpX/WarpX/issues/320).
 * 
 * Further references:
 * [Higuera's article on arxiv](https://arxiv.org/abs/1701.05605)
 * [Riperda's comparison of relativistic particle integrators](https://doi.org/10.3847/1538-4365/aab114)
 */
struct HC {};

template < >
Vec3D Pusher< HC >::operator ()(const Vec3D& u_vector_i, const Vec3D& E, const Vec3D& B)
{
    // First half electric field acceleration
    Vec3D u_minus = vec_add(u_vector_i, scalmultip(E,(q*Dt)/(2.*m)));

        //using namespace std;
        //cout << "u_minus = " << vec_print(u_minus) << endl;
        //cout << "u_minus / c = " << vec_print(scalmultip( u_minus , 1./c )) << endl;

    // Auxiliary quantities
    double gamma_minus = sqrt(1. + (innerprod(u_minus, u_minus))/(pow(c,2.)));

        //cout << "gamma_minus = " << gamma_minus << endl;

    Vec3D tau = scalmultip( B, (q * Dt ) / ( 2. * m ));

        //cout << "tau = " << vec_print(tau) << endl;

    double sigma = pow(gamma_minus,2.) - innerprod(tau,tau);

        //cout << "sigma = " << sigma << endl;

    double u_star = innerprod( u_minus, scalmultip(tau,1./c) );

        //cout << "u_star = " << u_star << endl;

        //cout << "sigma^2 = " << pow(sigma,2.) << endl;

        //cout << "4tau^2 = " << 4. * (innerprod(tau,tau)) << endl;

        //cout << "u*^2 = " << pow(u_star,2.) << endl;

        //cout << "(sigma^2 + 4tau^2 + u*^2)^(1/2) = " << sqrt( pow(sigma,2.) + 4. * (innerprod(tau,tau)+ pow(u_star,2.))) << endl;

        //cout << "sigma + (sigma^2 + 4tau^2 + u*^2)^(1/2) = " << (sigma + sqrt( pow(sigma,2.) + 4. * (innerprod(tau,tau)+ pow(u_star,2.)))) << endl;

    double gamma_plus = sqrt( .5 * ( sigma + sqrt( pow(sigma,2.) + 4. * ( innerprod(tau,tau) + pow(u_star,2.) ) ) ) );

        //cout << "gamma_plus = " << gamma_plus << endl;

    Vec3D t_vector = scalmultip(tau,1./gamma_plus);

        //cout << "t_vector = " << vec_print(t_vector) << endl;

    double s = 1./(1.+ innerprod(t_vector,t_vector));

        //cout << "s = " << s << endl;

    // Rotation step
    Vec3D u_plus = scalmultip(vec_add( vec_add( u_minus , scalmultip( t_vector , innerprod( u_minus , t_vector))), crossprod(u_minus,t_vector)), s);

        //cout << "u_plus = " << vec_print(u_plus) << endl;
        //cout << "u_plus / c = " << vec_print(scalmultip( u_plus , 1./c )) << endl;

    // Second half electric field acceleration
    Vec3D u_prime1 = scalmultip( E, ( q*Dt) / (2.*m) );
    Vec3D u_prime2 = crossprod(u_plus,t_vector);
    Vec3D u_prime = vec_add( u_prime1 , u_prime2 );
    
        //cout << "u_primeE / c = " << vec_print(scalmultip( u_prime1 , 1./c )) << endl;
        //cout << "u_primet / c = " << vec_print(scalmultip( u_prime2 , 1./c )) << endl;
        //cout << "u_prime / c = " << vec_print(scalmultip( u_prime , 1./c )) << endl;
    
    Vec3D u_ip1 = vec_add( u_plus , u_prime );
    
        //cout << "u_new / c = " << vec_print(scalmultip( u_ip1 , 1./c )) << endl;
        //cout << "==============================================================="<< endl;
    
    return u_ip1;
}



/**
 * Implementation of the Boris pusher as presented in doi:10.3847/1538-4365/aab114.
 */
struct Boris {};

template < >
Vec3D Pusher< Boris >::operator ()(const Vec3D& u_i, const Vec3D& E, const Vec3D& B)
{
        //using namespace std;
        //cout << "u_i = " << vec_print(u_i) << endl;
        //cout << "E = " << vec_print(E) << endl;
        //cout << "B = " << vec_print(B) << endl;
        //cout << "epsilon = " << vec_print( scalmultip( E , q * Dt / (2.*m) ) ) << endl;
    
    // First half electric field acceleration
    Vec3D u_minus = vec_add( u_i, scalmultip( E , q * Dt / (2.*m) ) );
        //cout << "u_minus = " << vec_print(u_minus) << endl;

    // Auxiliary quantities
    double gamma_minus = sqrt(1. + (innerprod(u_minus, u_minus))/(pow(c,2.)));
        //cout << "gamma_minus = " << gamma_minus << endl;

    Vec3D t_vector = scalmultip( B, (q * Dt ) / ( 2. * m * gamma_minus ) );
        //cout << "t_vector = " << vec_print(t_vector) << endl;

    Vec3D s_vector = scalmultip( t_vector , 2./(1.+ innerprod(t_vector,t_vector)) ) ;
        //cout << "s_vector = " << vec_print(s_vector) << endl;

    // Rotation step
    Vec3D u_prime = vec_add( u_minus , crossprod( u_minus , t_vector ) );
    Vec3D u_plus = vec_add( u_minus , crossprod( u_prime , s_vector ) );
        //cout << "u_plus = " << vec_print(u_plus) << endl;

    // Second half electric field acceleration
    Vec3D u_ip1 = vec_add(u_plus , scalmultip( E, ( q*Dt) / (2.*m) ) );
        //cout << "u_i+1 = " << vec_print(u_ip1) << endl << endl;
    
    return u_ip1;
}
