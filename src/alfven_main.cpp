#include "mfem.hpp"
#include "core/alfven_solver.hpp"
#include "utils/functions.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <iomanip>

int main(int argc, char *argv[])
{
    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    // Parameters (default values)
    int ref_lvls = 2;
    int Nt = 16;
    double tmax = M_PI * 2.0 * sqrt(2.0);
    double dt = tmax / Nt;
    int order = 1;
    double tol = 1e-14;
    std::string mesh_file = std::string(DATA_DIR) + "/ref-cube.mesh";
//  double dt_over_two = dt / 2.0;
    int iter = 1000;
    const char* path_save = "./out/classic/";

    // Parse command-line arguments
    mfem::OptionsParser args(argc, argv);
    args.AddOption(&ref_lvls, "-r", "--refine", "Number of refinements.");
    args.AddOption(&Nt, "-nt", "--ntsteps", "Number of timesteps.");
    args.Parse();

    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout); // Print all options used

    dt = tmax / Nt;

    // Load mesh
    const char *mesh_file_cstr = mesh_file.c_str();
    mfem::Mesh mesh(mesh_file_cstr, 1, 1);
    int dim = mesh.Dimension();
    std::cout << "Mesh Dimension : " << dim << std::endl;

    // Refine mesh
    for (int i = 0; i < ref_lvls; i++) {
        std::cout << "Refinement: " << i + 1 << std::endl;
        mesh.UniformRefinement();
    }

    // Print mesh info
    std::cout << "\nNumber of elements : " << mesh.GetNE() << std::endl;
    std::cout <<   "Number of vertices : " << mesh.GetNV() << std::endl;

    // FE spaces
    mfem::H1_FECollection fec_CG(order, dim);
    mfem::FiniteElementSpace CG_u(&mesh, &fec_CG);
    mfem::FiniteElementSpace CG_p(&mesh, &fec_CG);

    std::cout << "Number of u unknowns: " << CG_u.GetTrueVSize() << std::endl;
    std::cout << "Number of p unknowns: " << CG_p.GetTrueVSize() << std::endl;

    // unknown grid functions
    mfem::GridFunction u(&CG_u);
    mfem::GridFunction p(&CG_p);

    // Initial conditions
    mfem::FunctionCoefficient u_0_coeff(mfem::u_0);
    mfem::FunctionCoefficient p_0_coeff(mfem::p_0);
    u.ProjectCoefficient(u_0_coeff);
    p.ProjectCoefficient(p_0_coeff);

    // Export to ParaView
    mfem::ParaViewDataCollection *pd = new mfem::ParaViewDataCollection(path_save, &mesh);
    pd->RegisterField("u", &u);
    pd->RegisterField("p", &p);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(mfem::VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    int Nit = 0;
    pd->SetCycle(Nit);
    pd->SetTime(0.0);
    pd->Save();

    // Old time step values
    mfem::Vector u_old(u.Size());
    mfem::Vector p_old(p.Size());
    u_old = 0.0;
    p_old = 0.0;

    // System size (u + p)
    int ssize = u.Size() + p.Size();
    std::cout << "size: " << ssize << std::endl << std::endl;

    // Vector x: the "full" one (not the tdof one)
    mfem::Vector x(ssize);
    x.SetVector(u, 0);
    x.SetVector(p, u.Size());

    // Identify dofs of u and p
    mfem::Array<int> u_dofs(u.Size());
    mfem::Array<int> p_dofs(p.Size());
    std::iota(&u_dofs[0], &u_dofs[u.Size()], 0);
    std::iota(&p_dofs[0], &p_dofs[p.Size()], u.Size());

    // Create Alfven operator
    mfem::AlfvenOperator oper(CG_u, CG_p, dt);

    // RHS vector
    mfem::Vector b(ssize);

    // Vectors for constrained system (with BC)
    mfem::Vector X(ssize);
    mfem::Vector B(ssize);
    mfem::Operator *A_constrained = nullptr;

    // Time loop
    double t = 0.0;
    mfem::real_t energy;
    int jj = 0;

    for (t = dt; t < tmax + dt; t += dt) {
        jj += 1;

        // Update old values before computing new ones
        u_old = u;
        p_old = p;

        // Form RHS
        oper.FormRHS(u_old, p_old, b);

        // Enforce BC
        oper.GetSystemOperator().FormLinearSystem(oper.GetEssentialTDofs(), x, b, 
                                                   A_constrained, X, B);

        // Solve
        MINRES(*A_constrained, B, X, 0, iter, tol, tol);
        A_constrained->RecoverFEMSolution(X, b, x);

        // Extract solution from x and store in u,p
        x.GetSubVector(u_dofs, u);
        x.GetSubVector(p_dofs, p);

        // Compute energy
        energy = oper.ComputeEnergy(u_old, p_old);
        std::cout << "step:\t" << jj << "\tt = " << t << "\tenergy = " << energy << std::endl;

        // Export data to ParaView
        Nit++;
        pd->SetCycle(Nit);
        pd->SetTime(t);
        pd->Save();
    }

    // Compute errors
    // u:         FE solution 
    // u_0_coeff: reference/exact solution 
    mfem::real_t up_err = u.ComputeL2Error(u_0_coeff);
    std::cout << "u   L2 error " << std::scientific << std::setprecision(16) << up_err << std::endl;

    mfem::real_t p_err = p.ComputeL2Error(p_0_coeff);
    std::cout << "phi L2 error " << std::scientific << std::setprecision(16) << p_err << std::endl;

    // Free memory
    delete pd;

    // Timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";

    return 0;
}
