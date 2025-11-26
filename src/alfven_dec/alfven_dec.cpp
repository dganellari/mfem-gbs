

// decoupled alfven wave equation 
// 
// ∂t^2 div ∇⊥ϕ + b·∇(b·∇ϕ) = 0
//
// with CG elems and Newmark in time


#include <fstream>
#include <iostream>
#include <algorithm>
#include <chrono>
#include "mfem.hpp"
#include <thread>


double pi = 3.14159265358979323846;
double eps = 0.0;
double sq1meps2 = std::sqrt(1.-eps*eps);
int test_case = 0; // test 0 and 7 are implemented (terminology from fenicsx code)


struct Parameters {
    int ref_lvls = 4;
    int Nt       = 10;
    double tmax  = 2*pi*std::sqrt(2.)*(3./4.)*0.1; // = 6.6 = 3/4 of a period
    double dt    = tmax/Nt;
    int order    = 1;
    double tol   = 1e-14;
    const char* mesh_file = "./ref-cube.mesh";
    double dt_over_two = dt/2.0;
    int iter     = 1000;
    const char* path_save = "./out/dec_par"; // NO TRAILING SLASH!, otherwise pvd file doesnt get generated
    int kdim     = 100;
};


class GeneralResidualMonitor : public mfem::IterativeSolverMonitor {
public:
    GeneralResidualMonitor(MPI_Comm comm, const std::string &prefix_,
                           int print_lvl, std::vector<mfem::real_t> &resvec_)
        : prefix(prefix_), residuals(resvec_)

    {
#ifndef MFEM_USE_MPI
        print_level = print_lvl;
#else
        int rank;
        MPI_Comm_rank(comm, &rank);
        if (rank == 0) {
            print_level = print_lvl;
        }
        else {
            print_level = -1;
        }
#endif
    }

    void MonitorResidual(int it, mfem::real_t norm,
                         const mfem::Vector &r, bool final) override;

private:
    const std::string prefix;
    int print_level;
    mutable mfem::real_t norm0;
    std::vector<mfem::real_t> &residuals;
};

void GeneralResidualMonitor::MonitorResidual(int it, mfem::real_t norm,
                                             const mfem::Vector &r, bool final) {
    if (print_level == 1 || (print_level == 3 && (final || it == 0))) {
        std::cout << prefix << " iteration " << std::setw(2) << it
                  << " : ||r|| = " << norm;
        if (it > 0) {
            std::cout << ",  ||r||/||r_0|| = " << norm / norm0;
        }
        else {
            norm0 = norm;
        }
        std::cout << '\n';
    }
    residuals.push_back(norm); // store residual for every iteration
}


mfem::real_t p_0(const mfem::Vector &x);
void bfield(const mfem::Vector &x, mfem::Vector &v);
void bfield_lhs(const mfem::Vector &x, mfem::Vector &v);


int main(int argc, char *argv[]) {

    // timer
    auto start = std::chrono::high_resolution_clock::now();

    // Initialize MPI and HYPRE.
    mfem::Mpi::Init();
    mfem::Hypre::Init();
    int myid = mfem::Mpi::WorldRank();

    // simulation parameters
    Parameters param;
    int ref_lvls = param.ref_lvls;
    double dt    = param.dt;
    int Nt       = param.Nt;
    double tmax  = param.tmax;
    int order    = param.order;
    double tol   = param.tol;
    double dt_over_two = param.dt_over_two;
    int iter     = param.iter;
    int kdim     = param.kdim;
    const char* path_save = param.path_save;

    // mesh
    const char *mesh_file = param.mesh_file;
    mfem::Mesh mesh(mesh_file, 1, 1); 
    int dim = mesh.Dimension();
    for (int i =0; i<ref_lvls; i++){
        mesh.UniformRefinement();
    }

    // parallel mesh: partitioning of the original one
    mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    // FE spaces
    mfem::FiniteElementCollection *fec_CG = new mfem::H1_FECollection(order,dim);
    mfem::ParFiniteElementSpace CG_p(&pmesh, fec_CG);

    // essential true dofs
    mfem::Array<int> ess_tdof_p;
    CG_p.GetBoundaryTrueDofs(ess_tdof_p);
    
    // unknown
    mfem::ParGridFunction p(&CG_p);

    // initial condition
    mfem::FunctionCoefficient p_0_coeff(p_0);
    p.ProjectCoefficient(p_0_coeff);

    // old time step values
    mfem::Vector p_old(p.Size());     p_old = 0.;
    mfem::Vector p_old_old(p.Size()); p_old_old = 0.;
    mfem::Vector p_diff(p.Size());    p_diff = 0.;
    mfem::Vector p_mid(p.Size());     p_mid = 0.;
    p_old_old = p; // as initial velocity is zero
    p_old_old *= std::pow(1/2.,0.5) * dt;
    p = 0.;

    // p_old = p;
    // p_old_old = 0.;

    // exporting tools to paraview
    mfem::ParaViewDataCollection *pd = new mfem::ParaViewDataCollection(path_save, &pmesh);
    pd->RegisterField("p" , &p);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(mfem::VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    int Nit = 0;
    pd->SetCycle(Nit);
    pd->SetTime(0.0);
    pd->Save();
    
    // system size
    int ssize = p.Size();

    // vector x: the "full" one (not the tdof one), enforce BC
    mfem::Vector x(ssize); 
    x.SetVector(p,0);

    // identify dofs of u and phi
    mfem::Array<int> p_dofs (p.Size());
    std::iota(&p_dofs[0], &p_dofs[p.Size()], 0);

    // parallel and perpendicular projections
    mfem::VectorFunctionCoefficient b_gfcoeff(dim, bfield);
    mfem::OuterProductCoefficient K(b_gfcoeff, b_gfcoeff); // Matrix (b. )b, parallel projection
    mfem::CrossCrossCoefficient Q(1., b_gfcoeff);          // Matrix bxbx, perpendicular projection
    double dt2 = std::pow(dt,2);
    
    // LHS bilinearform
    mfem::MatrixSumCoefficient Q_dt(Q, Q, 1./dt2, 0.); // scaling 1/dt^2
    mfem::MatrixSumCoefficient K_dt(K, K, 1./4, 0.);   // scaling 1/4
    mfem::ParBilinearForm N_lhs(&CG_p);
    N_lhs.AddDomainIntegrator(new mfem::DiffusionIntegrator(Q_dt)); // 1/Δt^2 Δ_perp
    N_lhs.AddDomainIntegrator(new mfem::DiffusionIntegrator(K_dt)); // 1/4 Δ_par
    N_lhs.Assemble();
    N_lhs.Finalize();

    // RHS par
    mfem::ParBilinearForm blf_N_par(&CG_p);
    blf_N_par.AddDomainIntegrator(new mfem::DiffusionIntegrator(K)); // (b·∇p,b·∇q)
    blf_N_par.Assemble();
    blf_N_par.Finalize();

    // RHS perp
    mfem::ParBilinearForm blf_N_perp(&CG_p);
    blf_N_perp.AddDomainIntegrator(new mfem::DiffusionIntegrator(Q)); // (∇⊥p, ∇⊥q) := (∇p,∇q) - (b·∇p, b·∇q)
    blf_N_perp.Assemble();
    blf_N_perp.Finalize();

    // matrices for energy computation
    mfem::SparseMatrix par(blf_N_par.SpMat());
    mfem::SparseMatrix perp(blf_N_perp.SpMat());
    par.Finalize();
    perp.Finalize();
    
    // RHS vector
    mfem::Vector b(ssize);

    // BC-constrained ops and vectors
    mfem::OperatorPtr A_constrained;
    mfem::Vector X(ssize); X=0.;
    mfem::Vector B(ssize); B=0.;

    // solver prep
    mfem::CGSolver pcg(MPI_COMM_WORLD);
    pcg.SetPrintLevel(0);
    pcg.SetMaxIter(200);
    pcg.SetRelTol(1e-10);
    pcg.SetAbsTol(1e-10);

    

    

    // residual monitor
    std::vector<mfem::real_t> resvec;
    GeneralResidualMonitor mon(MPI_COMM_WORLD, "CG", 1, resvec);








    // nr of iterations array
    int numiter[Nt + 1];
    int index = 0;
    
    // AMG     
    mfem::HypreBoomerAMG amg;
    amg.SetPrintLevel(0);      // turn off Hypre output
    // amg.SetTol(0.0);           // no iteration, pure preconditioning
    // amg.SetOperator(A_constrained);
    // Optionally tweak AMG parameters
    // amg.SetCoarsening(6);       // e.g., HMIS coarsening
    // amg.SetInterpolation(6);    // e.g., extended+i interpolation
    // amg.SetRelaxType(6);        // e.g., symmetric Gauss-Seidel
    // amg.SetMaxLevels(25);

    // time loop
    double t = 0.0;
    for (t = dt ; t < tmax+dt ; t+=dt) {
        
        // update rhs
        b = 0.0;
        blf_N_perp.AddMult(p_old, b, +2/dt2);     // + 2 phi^n 
        blf_N_perp.AddMult(p_old_old, b, -1/dt2); // - phi^n-1
        blf_N_par.AddMult(p_old, b, -2/4.);       // - 2 phi^n / 4
        blf_N_par.AddMult(p_old_old, b, -1/4.);   // - phi^n-1 / 4
        
        // enforce BC
        N_lhs.FormLinearSystem(ess_tdof_p, x, b, A_constrained, X, B);
        
        // CG
        pcg.SetPreconditioner(amg);
        if (index > -1) {pcg.SetMonitor(mon);}
        pcg.SetOperator(*A_constrained);
        pcg.Mult(B, X);
        numiter[index] = pcg.GetNumIterations();

        // FEM sol into x
        N_lhs.RecoverFEMSolution(X, b, x);

        // sol from x to p
        x.GetSubVector(p_dofs, p);

        // (energy) TODO
        p_diff = p;
        p_diff -= p_old;
        p_diff *= 1/dt;
        p_mid = p;
        p_mid += p_old;
        p_mid *= 0.5;
        // double energy = 0.;
        double energy = perp.InnerProduct(p_diff,p_diff)
                       + par.InnerProduct(p_mid,p_mid);
        if (myid == 0) {
            std::cout << "t = " << t << ", energy = " << energy << ", iter = " << numiter[index] << std::endl;
        }
        index ++;

        // paraview TODO
        Nit ++;
        pd->SetCycle(Nit);
        pd->SetTime(t);
        pd->Save();

        // update old values
        p_old_old = p_old;
        p_old = p;

    } // time loop
    
    // error
    mfem::real_t p_err = p.ComputeL2Error(p_0_coeff);

    // print iter
    double avg, max, min;
    max=numiter[0];
    min=numiter[0];
    for (int i=0; i<Nt; i++) {
        avg+=numiter[i];
        if (max < numiter[i]) {
            max = numiter[i];
        }
        if (min > numiter[i]) {
            min = numiter[i];
        }
    }
    avg = avg/(Nt+1);




    // save in csv file
    std::ofstream fout("residuals.csv");
    if (!fout.is_open()) { std::cerr << "Error opening file\n"; }
    for (size_t i = 0; i < resvec.size(); ++i) {
        fout << i << "," << resvec[i] << "\n";
    }
    fout.close();






    // free memory
    delete fec_CG, pd;

    // timer and prints
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    if (myid == 0) {
        std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
        std::cout << "matrix size: " << ssize << std::endl;
        std::cout << "Ni: " << std::pow(pmesh.GetNE(),1./3.) << std::endl;
        std::cout << "Nt: " << Nt << std::endl;
        
        // print info for tex tables
        std::cout << "avg/max/min iter, runtime, error" << std::endl;
        std::cout << std::round(avg*10.0)/10.0             << " & " 
        << max                                   << " & "
        << min                                   << " & "
        << std::round(elapsed.count()*10.0)/10.0 << " & "
        << p_err                                 << " \\\\"
        << std::endl;
    }
        
} // main

mfem::real_t p_0(const mfem::Vector &x) {
    
    
    if (test_case==0) {
        double X = x(0);
        double Y = x(1);
        double Z = x(2);
        // return 0.;
        // 
        return std::sin(pi*X) * std::sin(pi*Y) * std::sin(pi*Z);
    }
    else if (test_case==7) {
        double X = x(0);
        double Y = x(1);
        double Z = x(2);
        return std::sin(pi*X) * std::sin(pi*Y) * std::sin(pi*Z);
    }
    else {
        std::cout << "Error: test case not implemented!" << std::endl;
        return -1.;
    }

}

void bfield(const mfem::Vector &x, mfem::Vector &returnvalue) { 
   
    if (test_case==0) {
        returnvalue(0) = 0.0; //+ 0.0*x(0);
        returnvalue(1) = 0.0; //+ 0.0*x(0);
        returnvalue(2) = 1.0; //+ 0.0*x(0);
    }
    else if (test_case==7) {
        returnvalue(0) = eps;
        returnvalue(1) = 0.;
        returnvalue(2) = sq1meps2;
    }
    else {
        std::cout << "Error: test case not implemented!" << std::endl;
    }
}

void bfield_lhs(const mfem::Vector &x, mfem::Vector &returnvalue) { 
   
    if (test_case==0) {
        Parameters param;
        double dt    = param.dt;
        double factor = -std::pow(dt,-2) + 1/4.;
        returnvalue(0) = +0.0 * factor;
        returnvalue(1) = +0.0 * factor;
        returnvalue(2) = +1.0 * factor;
    }
    else {
        std::cout << "Error: test case not implemented!" << std::endl;
    }
}

