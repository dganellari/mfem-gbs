#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include "mfem.hpp"
#include <thread>

double pi = 3.14159265358979323846;
mfem::real_t u_0(const mfem::Vector &x);
mfem::real_t p_0(const mfem::Vector &x);
void bfield(const mfem::Vector &x, mfem::Vector &v);

struct Parameters {
    int ref_lvls = 4;
    int Nt       = 64;
    double tmax  = pi*2.*std::sqrt(2.);
    double dt    = tmax/Nt;
    int order    = 1;
    double tol   = 1e-14;
    const char* mesh_file = "./ref-cube.mesh";
    double dt_over_two = dt/2.0;
    int iter     = 1000;
    const char* path_save = "./out/classic/";
};

int main(int argc, char *argv[]) {

    // timer
    auto start = std::chrono::high_resolution_clock::now();

    // simulation parameters
    Parameters param;
    int ref_lvls = param.ref_lvls;
    double dt    = param.dt;
    double tmax  = param.tmax;
    int order    = param.order;
    double tol   = param.tol;
    double dt_over_two = param.dt_over_two;
    int iter     = param.iter;
    const char* path_save = param.path_save;

    // glvis
    // char vishost[] = "localhost";
    // int  visport1   = 19916;
    // int  visport2   = 19917; // call with ./glvis -p 19917
    // mfem::socketstream u_sock(vishost, visport1);
    // mfem::socketstream p_sock(vishost, visport1);
    // u_sock.precision(8);
    // p_sock.precision(8);
    // mfem::socketstream u_sock_init(vishost, visport1);
    // mfem::socketstream p_sock_init(vishost, visport1);
    // u_sock_init.precision(8);
    // p_sock_init.precision(8);

    // mesh
    const char *mesh_file = param.mesh_file;
    mfem::Mesh mesh(mesh_file, 1, 1); 
    int dim = mesh.Dimension();
    for (int i =0; i<ref_lvls; i++){
        mesh.UniformRefinement();
    }

    // FE spaces
    mfem::FiniteElementCollection *fec_CG = new mfem::H1_FECollection(order,dim);
    mfem::FiniteElementCollection *fec_DG = new mfem::L2_FECollection(order,dim);
    mfem::FiniteElementSpace CG_u(&mesh, fec_CG);
    // mfem::FiniteElementSpace CG_u(&mesh, fec_DG);
    mfem::FiniteElementSpace CG_p(&mesh, fec_CG);

    // essential true dofs
    mfem::Array<int> ess_tdof_p;
    CG_p.GetBoundaryTrueDofs(ess_tdof_p);
    for (int i=0; i<ess_tdof_p.Size(); i++) {
            ess_tdof_p[i] += CG_u.GetNDofs();
    }
    
    // unkown functions
    mfem::GridFunction u(&CG_u);
    mfem::GridFunction p(&CG_p);

    // initial condition
    mfem::FunctionCoefficient u_0_coeff(u_0);
    mfem::FunctionCoefficient p_0_coeff(p_0);
    u.ProjectCoefficient(u_0_coeff);
    p.ProjectCoefficient(p_0_coeff);

    // glvis: init cond
    // u_sock_init << "solution\n" << mesh << u << "window_title 'init cond u'" << std::endl;
    // p_sock_init << "solution\n" << mesh << p << "window_title 'init cond p'" << std::endl;

    // exporting tools to paraview
    mfem::ParaViewDataCollection *pd = new mfem::ParaViewDataCollection(path_save, &mesh);
    pd->RegisterField("u" , &u);
    pd->RegisterField("p" , &p);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(mfem::VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    int Nit = 0;
    pd->SetCycle(Nit);
    pd->SetTime(0.0);
    pd->Save();

    // old time step values
    mfem::Vector u_old(u.Size()); u_old = 0.;
    mfem::Vector p_old(p.Size()); p_old = 0.;
    
    // system size
    int ssize = u.Size() + p.Size();
    std::cout << "size: " << ssize << std::endl;

    // vector x: the "full" one (not the tdof one)
    // here we implicitely also set the BC, 
    // as they are given by the u_0_coeff,p_0_coeff
    mfem::Vector x(ssize); // maybe use a blockvector instead
    x.SetVector(u,0);
    x.SetVector(p,u.Size());

    // identify dofs of u and phi
    mfem::Array<int> u_dofs (u.Size());
    mfem::Array<int> p_dofs (p.Size());
    std::iota(&u_dofs[0], &u_dofs[u.Size()], 0);
    std::iota(&p_dofs[0], &p_dofs[p.Size()], u.Size());

    // matrices
    // [ M -E] [u] _ [Mu + Ep + f] 
    // [-F -N] [p] ¯ [Fu + Np + g] 

    // Matrix M
    mfem::BilinearForm blf_M(&CG_u);
    blf_M.AddDomainIntegrator(new mfem::MassIntegrator()); //=(u,v)
    blf_M.Assemble();
    blf_M.Finalize();
    mfem::SparseMatrix M(blf_M.SpMat());
    M.Finalize();

    // Matrix N: parallel diffusion
    mfem::VectorFunctionCoefficient b_gfcoeff(dim, bfield);
    mfem::OuterProductCoefficient K(b_gfcoeff, b_gfcoeff);
    mfem::BilinearForm blf_N_par(&CG_p);
    blf_N_par.AddDomainIntegrator(new mfem::DiffusionIntegrator(K)); // (b·∇p,b·∇q)
    blf_N_par.Assemble();
    blf_N_par.Finalize();
    mfem::SparseMatrix N_par(blf_N_par.SpMat());
    N_par.Finalize();

    // Matrix N: full diffusion
    mfem::BilinearForm blf_N_full(&CG_p);
    blf_N_full.AddDomainIntegrator(new mfem::DiffusionIntegrator);  // (∇p,∇q)
    blf_N_full.Assemble();
    blf_N_full.Finalize();
    mfem::SparseMatrix N_full(blf_N_full.SpMat());

    // Matrix N: perp diffusion
    mfem::SparseMatrix N_n(N_full);
    N_n.Add(-1.0, N_par); // (∇p,∇q) - (b·∇p, b·∇q)
    N_n *= -1.;
    N_n.Finalize();

    // Matrix E and F
    mfem::MixedBilinearForm blf_E(&CG_p, &CG_u);
    blf_E.AddDomainIntegrator(new mfem::MixedDirectionalDerivativeIntegrator(b_gfcoeff)); //=(b·∇p,u)
    blf_E.Assemble();
    blf_E.Finalize();
    mfem::SparseMatrix E(blf_E.SpMat());
    E *= -1.*dt_over_two;

    mfem::SparseMatrix *F;
    F = Transpose(E);
    E.Finalize();
    F->Finalize();
    
    // initialize system matrices
    mfem::Array<int> offsets (3);
    offsets[0] = 0;
    offsets[1] = u.Size();
    offsets[2] = p.Size();
    offsets.PartialSum(); // exclusive scan
    mfem::BlockOperator A(offsets);
    A.SetBlock(0,0, &M);
    A.SetBlock(0,1, &E);
    A.SetBlock(1,0, F);
    A.SetBlock(1,1, &N_n);

    // rhs
    mfem::Vector b(ssize);
    mfem::Vector bsub1(u.Size()); // subvectors to update rhs in time loop
    mfem::Vector bsub2(p.Size()); // subvectors to update rhs in time loop

    // enforce BC
    mfem::Vector X(ssize); 
    mfem::Vector B(ssize); 
    mfem::Operator *A_constrained = nullptr;

    // TODO: we should probably call 
    // A.FormLinearSystem(ess_tdof_p, x, b, A_constrained, X, B);
    // outside the time loop only once

    // time loop
    double t;
    for (t = dt ; t < tmax+dt ; t+=dt) {
    
        // update old values before computing new ones
        u_old = u;
        p_old = p;

        // update rhs
        b = 0.0;
        bsub1 = 0.0;
        bsub2 = 0.0;
        M   .AddMult(u_old,bsub1);
        E   .AddMult(p_old,bsub1,-1);
        N_n .AddMult(p_old,bsub2,1);
        F  ->AddMult(u_old,bsub2,-1);
        b.AddSubVector(bsub1,0);
        b.AddSubVector(bsub2,u.Size());
        
        // enforce BC
        A.FormLinearSystem(ess_tdof_p, x, b, A_constrained, X, B);
        
        // solve 
        mfem::MINRES(*A_constrained, B, X, 0, iter, tol, tol);  
        A_constrained->RecoverFEMSolution(X, b, x);
    
        // extract sol from x and store in u,p
        x.GetSubVector(u_dofs, u);
        x.GetSubVector(p_dofs, p);

        // energy
        double energy = blf_M.InnerProduct(u_old,u_old)
                       + blf_N_full.InnerProduct(p_old,p_old)
                       - blf_N_par.InnerProduct(p_old,p_old);
        std::cout << "t = " << t << ", energy = " << energy << std::endl;

        // stream to glvis
        // u_sock << "solution\n" << mesh << u
        //     << "window_title 'u'" << std::endl;
        // p_sock << "solution\n" << mesh << p
        //     << "window_title 'p'" << std::endl;
        // std::this_thread::sleep_for(std::chrono::milliseconds(500));
        
        // Export data to paraview
        Nit ++;
        pd->SetCycle(Nit);
        pd->SetTime(t);
        pd->Save();
        
    } // time loop
    
    // error
    mfem::real_t up_err = u.ComputeL2Error(u_0_coeff);
    std::cout << "u   L2 error " << up_err << std::endl;
    mfem::real_t p_err = p.ComputeL2Error(p_0_coeff);
    std::cout << "phi L2 error " << p_err << std::endl;
    
    // glvis
    // u_sock << "solution\n" << mesh << u << "window_title 'SOLu'" << std::endl;
    // p_sock << "solution\n" << mesh << p << "window_title 'SOLp'" << std::endl;

    // free memory
    delete fec_CG;

    // timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";

} // main

mfem::real_t u_0(const mfem::Vector &x) {
    
    double X = x(0);
    double Y = x(1);
    double Z = x(2);
    // return std::sin(pi*X) * std::sin(pi*Y) * std::sin(pi*Z);
    return std::sin(pi*X) * std::sin(pi*Y) * (1.- std::cos(pi*Z));
}

mfem::real_t p_0(const mfem::Vector &x) {
    return 0.;
}

void bfield(const mfem::Vector &x, mfem::Vector &returnvalue) { 
   
    returnvalue(0) = 0.0; //+ 0.0*x(0);
    returnvalue(1) = 0.0; //+ 0.0*x(0);
    returnvalue(2) = 1.0; //+ 0.0*x(0);
}
