// decoupled alfven wave equation with CUDA support
// 
// ∂t^2 div ∇⊥ϕ + b·∇(b·∇ϕ) = 0
//
// with CG elems and Newmark in time
// CUDA-enabled for GPU acceleration

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
        : prefix(prefix_), residuals(resvec_), comm_(comm)
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
    MPI_Comm comm_;
};

void GeneralResidualMonitor::MonitorResidual(int it, mfem::real_t norm,
                                             const mfem::Vector &r, bool final) {
    int rank;
    MPI_Comm_rank(comm_, &rank);
    if (rank == 0) {
        std::cout << "DEBUG Monitor: it=" << it << ", norm=" << norm 
                  << ", final=" << final << ", r.Size()=" << r.Size() << std::endl;
    }
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
mfem::real_t rhs_source(const mfem::Vector &x);
void bfield(const mfem::Vector &x, mfem::Vector &v);
void bfield_lhs(const mfem::Vector &x, mfem::Vector &v);

int main(int argc, char *argv[]) {
    // timer
    auto start = std::chrono::high_resolution_clock::now();

    // Initialize MPI and HYPRE.
    mfem::Mpi::Init();
    mfem::Hypre::Init();
    int myid = mfem::Mpi::WorldRank();

    // Parse command-line options
    mfem::OptionsParser args(argc, argv);
    std::string mesh_file = "./ref-cube.mesh";
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    std::string precon_type = "amg";
    args.AddOption(&precon_type, "-p", "--preconditioner", "Preconditioner type: none, amg, smoother");
    args.Parse();
    if (!args.Good()) {
        if (myid == 0) args.PrintUsage(std::cout);
        return 1;
    }

    // NEW: Enable CUDA device for GPU acceleration
    mfem::Device device("cuda");
    if (myid == 0) {
        device.Print();  // Print device info (e.g., CUDA version, GPU name)
    }

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
    const char *mesh_file_c = mesh_file.c_str();
    mfem::Mesh mesh(mesh_file_c, 1, 1); 
    int dim = mesh.Dimension();
    // 5. Refine the serial mesh on all processors to increase the resolution. In
    //    this example we do 'ref_levels' of uniform refinement. We choose
    //    'ref_levels' to be the largest number that gives a final mesh with no
    //    more than 10,000 elements.
    {
       int ref_levels =
          (int)floor(log(10000./mesh.GetNE())/log(2.)/dim);
       for (int l = 0; l < ref_levels; l++)
       {
          mesh.UniformRefinement();
       }
    }

    // 6. Define a parallel mesh by a partitioning of the serial mesh. Refine
    //    this mesh further in parallel to increase the resolution. Once the
    //    parallel mesh is defined, the serial mesh can be deleted.
    mfem::ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    {
       int par_ref_levels = 3;  // Increased from 2 for larger problem (~2M DOF)
       for (int l = 0; l < par_ref_levels; l++)
       {
          pmesh.UniformRefinement();
       }
    }

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
    p_old_old *= std::pow(1.0/2.0, 0.5) * dt;  // FIXED: Avoid integer division
    p = 0.;

    // exporting tools to paraview (transfer to CPU for I/O)
    mfem::Vector p_cpu = p;  // Copy to CPU for ParaView
    mfem::ParFiniteElementSpace CG_p_host(&pmesh, fec_CG);
    mfem::ParGridFunction p_gf(&CG_p_host);  // Create host GridFunction for ParaView
    p_gf = p_cpu;  // Set data from CPU vector
    mfem::ParaViewDataCollection *pd = new mfem::ParaViewDataCollection(path_save, &pmesh);
    pd->RegisterField("p" , &p_gf);  // Use host GridFunction
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
    
    // Variable coefficient to make operator ill-conditioned (100:1 contrast)
    mfem::FunctionCoefficient var_coeff([](const mfem::Vector &x) -> double {
        double X = x(0), Y = x(1), Z = x(2);
        // High contrast localized around center (like material interface)
        return 1.0 + 99.0 * std::exp(-50*((X-0.5)*(X-0.5) + (Y-0.5)*(Y-0.5) + (Z-0.5)*(Z-0.5)));
    });
    
    // LHS bilinearform
    mfem::MatrixSumCoefficient Q_dt(Q, Q, 1./dt2, 0.); // scaling 1/dt^2
    mfem::ScalarMatrixProductCoefficient Q_dt_var(var_coeff, Q_dt); // Apply variable coefficient to perpendicular term
    mfem::MatrixSumCoefficient K_dt(K, K, 1./4, 0.);   // scaling 1/4
    mfem::ParBilinearForm N_lhs(&CG_p);
    N_lhs.AddDomainIntegrator(new mfem::DiffusionIntegrator(Q_dt_var)); // 1/Δt^2 Δ_perp with variable coefficient
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
    
    // Optional: Add a source term to make RHS more challenging
    mfem::FunctionCoefficient source_coeff(rhs_source);
    mfem::ParLinearForm source_lf(&CG_p);
    source_lf.AddDomainIntegrator(new mfem::DomainLFIntegrator(source_coeff));
    source_lf.Assemble();

    // BC-constrained ops and vectors
    mfem::OperatorPtr A_constrained;
    mfem::Vector X(ssize); X=0.;
    mfem::Vector B(ssize); B=0.;

    // solver prep
    mfem::CGSolver pcg(MPI_COMM_WORLD);
    pcg.SetPrintLevel(0);
    pcg.SetMaxIter(500);   // Increased from 200 to allow convergence
    pcg.SetRelTol(1e-10);  // Relaxed from 1e-12
    pcg.SetAbsTol(1e-10);  // Relaxed from 1e-12

    // residual monitor
    std::vector<mfem::real_t> resvec;
    GeneralResidualMonitor mon(MPI_COMM_WORLD, "CG", 1, resvec);

    // nr of iterations array
    int numiter[Nt + 1];
    int index = 0;
    
    // AMG     
    mfem::Solver *precon = nullptr;
    if (precon_type == "amg") {
        precon = new mfem::HypreBoomerAMG();
        static_cast<mfem::HypreBoomerAMG*>(precon)->SetPrintLevel(0);
    } else if (precon_type == "smoother") {
        precon = new mfem::HypreSmoother();
        static_cast<mfem::HypreSmoother*>(precon)->SetType(mfem::HypreSmoother::l1Jacobi);
    } else if (precon_type == "none") {
        precon = nullptr;
    } else {
        if (myid == 0) std::cout << "Invalid preconditioner type: " << precon_type << std::endl;
        return 1;
    }

    // time loop
    double t = 0.0;
    for (t = dt ; t < tmax+dt ; t+=dt) {
        
        // update rhs
        b = 0.0;
        blf_N_perp.AddMult(p_old, b, +2/dt2);     // + 2 phi^n 
        blf_N_perp.AddMult(p_old_old, b, -1/dt2); // - phi^n-1
        blf_N_par.AddMult(p_old, b, -2/4.);       // - 2 phi^n / 4
        blf_N_par.AddMult(p_old_old, b, -1/4.);   // - phi^n-1 / 4
        // Add challenging source term
        b.Add(1.0, source_lf);
        
        // enforce BC
        N_lhs.FormLinearSystem(ess_tdof_p, x, b, A_constrained, X, B);
        
        // CG
        if (precon) pcg.SetPreconditioner(*precon);
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
        double energy = perp.InnerProduct(p_diff,p_diff)
                       + par.InnerProduct(p_mid,p_mid);
        if (myid == 0) {
            std::cout << "t = " << t << ", energy = " << energy << ", iter = " << numiter[index] << std::endl;
        }
        index ++;

        // paraview (transfer to CPU)
        p_cpu = p;  // Copy GPU data to CPU for output
        p_gf = p_cpu;  // Update GridFunction with new data
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

    // free memory (FIXED: Separate deletes)
    delete fec_CG;
    delete pd;
    if (precon) delete precon;

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

//AI generated source term to challenge CG solver
mfem::real_t rhs_source(const mfem::Vector &x) {
    // MUCH stronger high-frequency source to force deep Krylov subspace exploration
    double X = x(0);
    double Y = x(1);
    double Z = x(2);

    // Aggressive multi-scale forcing with large amplitude
    // This should force CG to explore many Krylov directions
    double source = 0.0;
    
    // High-frequency components (wavenumbers 20π, 30π, 40π)
    source += 5.0 * std::sin(20*pi*X) * std::cos(20*pi*Y) * std::sin(20*pi*Z);
    source += 3.0 * std::cos(30*pi*X) * std::sin(30*pi*Y) * std::cos(30*pi*Z);
    source += 2.0 * std::sin(40*pi*X) * std::sin(40*pi*Y) * std::sin(40*pi*Z);
    
    // Add some mid-frequency chaos
    source += 4.0 * std::sin(15*pi*X) * std::cos(15*pi*Y) * std::sin(15*pi*Z);
    
    // Mix with spatial gradients (creates anisotropy)
    source *= (1.0 + 0.5*X) * (1.0 + 0.3*Y);
    
    return source;
}

void bfield(const mfem::Vector &x, mfem::Vector &returnvalue) { 
    if (test_case==0) {
        returnvalue(0) = 0.0;
        returnvalue(1) = 0.0;
        returnvalue(2) = 1.0;
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