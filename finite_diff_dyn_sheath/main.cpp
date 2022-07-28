/* Two-species Vlasov-Poisson solver on [-1,1] \times \RR using
 
 with :
 - non-emitting boundary conditions on the distribution functions
 - Source term for the ions
 - Dynamic electric field at the boundary
 
 
 Simulation of the symmetric sheath problem with ionization term.
 
 Discretization :
 - Upwind for the transport
 - Integration for the Poisson equation -\lambda^2 dx E = rho with trapeze formula for rho.  The integration  preserves the symmetry of rho if it is even.
 
 
 Initial data :
 
 - Initial data are in the file initial_data.hpp
 
 Author : Mehdi BASDI
 Date : 27/07/2022
 
 
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include "Matrix.h"
#include "initial_data.hpp"

using namespace std;

/*Physical parameter*/
double nu    = 1.0;
double lambda = 0.1;
double  T     = 2.0;

/* Discretization parameter */
double Nx    = 200;
double Nv    = 400;
double xmin  =-1.0;
double xmax  = 1.0;
// The grid in velocity corresponds to the electronic one which contains the support of boths ions and electrond density function
double vmin_e  = -10./sqrt(mu);
double  vmax_e = +10./sqrt(mu);
double  dx     = (xmax - xmin)/Nx;
double  dv     = (vmax_e  -vmin_e)/Nv;
double  CFL_x  = 0.5* dx/vmax_e;
double  CFL_v  = mu*dv/10;
int Nt         = floor(T/min(CFL_x,CFL_v))+1;


/* The discretized quantities*/
Matrix fi(Nx+1,Nv+1);
Matrix fe(Nx+1,Nv+1);
Vector rho(Nx+1);
Vector E(Nx+1);
Vector tE(Nx+1);
double J_xmax, J_xmin;
double Me;
double dt = T/Nt;


int main(int argc, char ** argv)
{
    
    double tn(0);
    
    /*Initialisation */
    J_xmin = 0; J_xmax = 0;
    for(int i = 0 ; i < Nx+1 ; i++)
    {
        
        E[i] = E0(xmin + i * dx);
        for(int j = 0 ; j < Nv +1; j++)
        {
            fi[i][j] = fi_0(xmin + i * dx, vmin_e + j*dv);
            fe[i][j] = fe_0(xmin + i * dx, vmin_e + j*dv);
        }
    }
    
    /* Prototypes of the functions defined here under*/
    Vector compute_rho(Matrix &, Matrix &,double&);
    Vector compute_eletric_field(Vector &,Vector &,double&,double &,double &);
    void compute_current_at_boundaries(Matrix &, Matrix &,double &,double&);
    Matrix update_fi(Matrix &,Matrix &, Vector&); Matrix update_fe(Matrix&,Vector&);
    void plot_sol(Matrix&,Matrix& , Vector &, Vector&, int &);
    
    /* Main temporal loop*/
    
    
    int iter(0);
    cout << " CFL_x = " << CFL_x << "    CFL_v = " << CFL_v << endl;
    while( tn < T )
    {
        rho = compute_rho(fi,fe,Me);
        compute_current_at_boundaries(fi,fi,J_xmin,J_xmax);
        E   = compute_eletric_field(rho,E,J_xmin,J_xmax,Me);tE = -(1/mu)*E;
        plot_sol(fi,fe,rho,E,iter);
        fi  = update_fi(fi,fe,E);
        fe =  update_fe(fe,tE);
        tn+=dt;
        cout << "Compute sol at time tn = " << tn << "  iter = " << iter <<  endl;
        //cout << "|| rho(tn) || = " << rho.linfty() << endl;
        iter++;
        plot_sol(fi,fe,rho,E,iter);
    }
    
    
    return 0;
}

void plot_sol(Matrix & fii,Matrix & fee, Vector & rho, Vector & E, int & iter)
{
    
    ofstream file("sol_phase_space"+to_string(iter)+".dat");
    for(int i = 0 ; i < Nx +1; i++)
    {
        for(int j = 0 ; j < Nv +1 ; j++)
        {
            file  << xmin + i * dx << "  " << vmin_e + j*dv <<"  " << fii[i][j] << "  "<< fee[i][j] <<endl;
        }
        file << endl;
    }
    file.close();
    file.open("sol_space"+to_string(iter)+".dat");
    for(int i = 0 ; i < Nx +1; i++)
    {
        file << xmin + i*dx << "  " << rho[i] << "  " << E[i] << endl;
    }
    file.close();
    
}

/* Compute the charge density the total eletronic mass :
 rho(x) = \int_{\RR} (fi - fe)(x,v)dv;
 Me     = \int_{-1}^{1} \int_{\RR} fe(x,v)dvdx
 */
Vector compute_rho(Matrix &fi, Matrix &fe, double & Me)
{
    Vector ro(Nx+1);
    Me = 0.;
    for(int i = 0 ; i < Nx+1 ; i++)
    {
        for(int j = 0 ; j < Nv +1; j++)
        {
            ro[i] += fi[i][j] * dv - fe[i][j] *dv;
        }
        Me += rho[i] * dx;
    }
    return ro;
}

/* Compute the current at x = xmax and x = xmin
 
 J(t,x) = \int_{\RR} (f_i - f_e)(t,x,v)vdv
 
 */

void compute_current_at_boundaries(Matrix &fi, Matrix &fe, double & J_l, double & J_r)
{
    J_r = 0; J_l = 0;
    for(int j = 0 ; j < Nv +1 ; j++)
        {
            J_r += dv* (fi[Nx][j] * (vmin_e + j * dv) - fe[Nx][j] * (vmin_e + j *dv));
            J_l += dv* (fi[0][j] * (vmin_e + j * dv) - fe[0][j] * (vmin_e + j *dv));
        }
    
}

/* Compute the eletric field with the symmetric formula:
 E(x) = 1/(2 Lambda^2) * (int_{-1}^{x} rho ds - int_{x}^{1} rho ds) + 0.5 (E(1) + E(-1))
 
 where
 
- E^n+1(1) = E^n(1)   + (\nu dt/2 Lambda^2) * M_e^n   -(dt/Lambda^2) J^n(1),
- E^n+1(-1) = E^n(-1) - (\nu dt/2 Lambda^2) * M_e^n   -(dt/Lambda^2) J^n(-1)

 */
Vector compute_eletric_field(Vector& ro,Vector & E,double & J_left, double & J_right, double & Mass_e)
{
    double E_xmax,E_xmin;
    
    // Compute the new boundary conditions on the electric_field;
    E_xmax = E[Nx] + 0.5 * (nu * dt/(lambda*lambda)) * Mass_e  - (dt/(lambda*lambda)) * J_right;
    E_xmin = E[0]  - 0.5 * (nu * dt/(lambda*lambda)) * Mass_e  - (dt/(lambda*lambda)) * J_left;
    
    double mass_rho_minus; //Quadrature for int_{-1}^{x_i} rho(x)dx.
    double mass_rho_plus;  //Quadrature for int_{x_i}^{1} rho(x)dx
    for(int i = 0 ; i < Nx+1 ; i++)
    {
        mass_rho_plus= 0; mass_rho_minus= 0;
        for(int k = 0; k <= i-1 ; k++)
        {
            mass_rho_minus += 0.5*dx * (ro[k+1] + ro[k]);
        }
        for(int  k = i; k<= Nx-1;k++)
        {
            mass_rho_plus += 0.5 *dx * (ro[k+1]+ro[k]);
        }
        E[i] = (1/(2*lambda*lambda)) * (mass_rho_minus - mass_rho_plus) + 0.5 * (E_xmax + E_xmin);
    }
    return E;
    
}

/* Update the distribution function */

Matrix update_fi(Matrix & fi_0,Matrix & fe_0, Vector & E)
{
    Matrix new_fi(Nx+1,Nv+1);
    double v_plus, v_minus;
    double E_plus, E_minus;
    /* Evole the discrete density in the interior domain*/
    for(int i = 1 ; i < Nx ; i++)
    {
        E_plus =  0.5  * (E[i] + abs(E[i]));
        E_minus = 0.5  * (E[i] - abs(E[i]));
        for(int j =1 ; j < Nv; j++)
        {
            v_plus =  0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
            v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
            new_fi[i][j] = fi_0[i][j]
            - (dt/dx) * ( v_plus * (fi_0[i][j] - fi_0[i-1][j])  + v_minus *(fi_0[i+1][j] - fi_0[i][j]) )
            - (dt/dv) * ( E_plus * (fi_0[i][j] - fi_0[i][j-1]) + E_minus * (fi_0[i][j+1] - fi_0[i][j]) )
            + nu * dt * fe_0[i][j];
                                                
        }
    }
    /* Apply the boundary condition at x = 0*/
    for(int j = 0; j < Nv  ; j++)
    {
        v_plus =  0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
        v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
        if( vmin_e + j*dv > 0)
        {
            new_fi[0][j] = 0.;
        }
        else
        {
            E_plus =  0.5  * (E[0] + abs(E[0]));
            E_minus = 0.5  * (E[0] - abs(E[0]));
            new_fi[0][j] = fi_0[0][j]
            - (dt/dx) * ( v_minus *(fi_0[1][j] - fi_0[0][j]) )
            - (dt/dv) * ( E_plus * (fi_0[0][j] - fi_0[0][j-1]) + E_minus * (fi_0[0][j+1] - fi_0[0][j]) )
            + nu * dt * fe_0[0][j];
        }
    }
    /* Apply the boundary condition at  x= 1*/
    for(int j = 0; j < Nv ; j++)
    {
        v_plus = 0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
        v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
        if( vmin_e + j*dv < 0)
        {
            new_fi[Nx][j] = 0.;
        }
        else
        {
            E_plus =  0.5  * (E[Nx] + abs(E[Nx]));
            E_minus = 0.5  * (E[Nx] - abs(E[Nx]));
            new_fi[Nx][j] = fi_0[Nx][j]
            - (dt/dx) * ( v_plus * (fi_0[Nx][j] - fi_0[Nx-1][j])  )
            - (dt/dv) * ( E_plus * (fi_0[Nx][j] - fi_0[Nx][j-1]) + E_minus * (fi_0[Nx][j+1] - fi_0[Nx][j]) )
            + nu * dt * fe_0[Nx][j];
        }
        
        /* Apply the boundary condition at v = v_max and v = vmin*/
        for(int i = 1 ; i < Nx ; i++)
        {
            
            
            new_fi[i][Nv] = 0;
            new_fi[i][0]  = 0;
         }
    }
    

    return new_fi;
}

// Udpate the distribution function

Matrix update_fe(Matrix & fe_0, Vector & E)
{
    Matrix new_fe(Nx+1,Nv+1);
    double v_plus, v_minus;
    double E_plus, E_minus;
    /* Evole the discrete density in the interior domain*/
    for(int i = 1 ; i < Nx ; i++)
    {
        E_plus =  0.5  * (E[i] + abs(E[i]));
        E_minus = 0.5  * (E[i] - abs(E[i]));
        for(int j =1 ; j < Nv; j++)
        {
            v_plus = 0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
            v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
            new_fe[i][j] = fe_0[i][j]
            - (dt/dx) * ( v_plus * (fe_0[i][j] - fe_0[i-1][j])  +v_minus *(fe_0[i+1][j] - fe_0[i][j]) )
            - (dt/dv) * ( E_plus * (fe_0[i][j] - fe_0[i][j-1]) + E_minus * (fe_0[i][j+1] - fe_0[i][j]) );
            
                                                
        }
    }
    
    /* Apply the boundary condition at x = 0*/
    for(int j = 0; j < Nv  ; j++)
    {
        v_plus =  0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
        v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
        if( vmin_e + j*dv > 0)
        {
            new_fe[0][j] = 0.;
        }
        else
        {
            E_plus =  0.5  * (E[0] + abs(E[0]));
            E_minus = 0.5  * (E[0] - abs(E[0]));
            new_fe[0][j] = fe_0[0][j]
            - (dt/dx) * ( v_minus *(fe_0[1][j] - fe_0[0][j]) )
            - (dt/dv) * ( E_plus * (fe_0[0][j] - fe_0[0][j-1]) + E_minus * (fe_0[0][j+1] - fe_0[0][j]) );
            
        }
    }
    /* Apply the boundary condition at  x= 1*/
    for(int j = 0; j < Nv ; j++)
    {
        v_plus = 0.5 * (vmin_e + j * dv + abs(vmin_e + j *dv));
        v_minus = 0.5 * (vmin_e + j * dv - abs(vmin_e + j *dv));
        if( vmin_e + j*dv < 0)
        {
            new_fe[Nx][j] = 0.;
        }
        else
        {
            E_plus =  0.5  * (E[Nx] + abs(E[Nx]));
            E_minus = 0.5  * (E[Nx] - abs(E[Nx]));
            new_fe[Nx][j] = fe_0[Nx][j]
            - (dt/dx) * ( v_plus * (fe_0[Nx][j] - fe_0[Nx-1][j])  )
            - (dt/dv) * ( E_plus * (fe_0[Nx][j] - fe_0[Nx][j-1]) + E_minus * (fe_0[Nx][j+1] - fe_0[Nx][j]) );
            
        }
        
        /* Apply the boundary condition at v = v_max and v = vmin*/
        for(int i = 1 ; i < Nx ; i++)
        {
            
            
            new_fe[i][Nv] = 0;
            new_fe[i][0]  = 0;
         }
    }
    

    return new_fe;
}


