#include "mfem.hpp"
#include "core/alfven_solver.hpp"
#include "utils/functions.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <algorithm>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    // Start timer
    auto start = chrono::high_resolution_clock::now();

    // Parameters (default values)
    int ref_lvls = 2;
    int Nt = 16;
    double tmax = M_PI * 2.0 * sqrt(2.0);
    double dt = tmax / Nt;
    int order = 1;
    double tol = 1e-14;
    string mesh_file = string(DATA_DIR) + "/ref-cube.mesh";
//  double dt_over_two = dt / 2.0;
    int iter = 1000;
    const char* path_save = "./out/classic/";

    // Parse command-line arguments
    OptionsParser args(argc, argv);
    args.AddOption(&ref_lvls, "-r", "--refine", "Number of refinements.");
    args.AddOption(&Nt, "-nt", "--ntsteps", "Number of timesteps.");
    args.Parse();

    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout); // Print all options used

    dt = tmax / Nt;

    // Load mesh
    const char *mesh_file_cstr = mesh_file.c_str();
    Mesh mesh(mesh_file_cstr, 1, 1);
    int dim = mesh.Dimension();
    cout << "Mesh Dimension : " << dim << endl;

    // Refine mesh
    for (int i = 0; i < ref_lvls; i++) {
        cout << "Refinement: " << i + 1 << endl;
        mesh.UniformRefinement();
    }

    // Print mesh info
    cout << "\nNumber of elements : " << mesh.GetNE() << endl;
    cout <<   "Number of vertices : " << mesh.GetNV() << endl;

    // FE spaces
    H1_FECollection fec_CG(order, dim);
    FiniteElementSpace CG_u(&mesh, &fec_CG);
    FiniteElementSpace CG_p(&mesh, &fec_CG);

    cout << "Number of u unknowns: " << CG_u.GetTrueVSize() << endl;
    cout << "Number of p unknowns: " << CG_p.GetTrueVSize() << endl;

    // unknown grid functions
    GridFunction u(&CG_u);
    GridFunction p(&CG_p);

    // Initial conditions
    FunctionCoefficient u_0_coeff(u_0);
    FunctionCoefficient p_0_coeff(p_0);
    u.ProjectCoefficient(u_0_coeff);
    p.ProjectCoefficient(p_0_coeff);

    // Export to ParaView
    ParaViewDataCollection *pd = new ParaViewDataCollection(path_save, &mesh);
    pd->RegisterField("u", &u);
    pd->RegisterField("p", &p);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    int Nit = 0;
    pd->SetCycle(Nit);
    pd->SetTime(0.0);
    pd->Save();

    // Old time step values
    Vector u_old(u.Size());
    Vector p_old(p.Size());
    u_old = 0.0;
    p_old = 0.0;

    // System size (u + p)
    int ssize = u.Size() + p.Size();
    cout << "size: " << ssize << endl << endl;

    // Vector x: the "full" one (not the tdof one)
    Vector x(ssize);
    x.SetVector(u, 0);
    x.SetVector(p, u.Size());

    // Identify dofs of u and p
    Array<int> u_dofs(u.Size());
    Array<int> p_dofs(p.Size());
    iota(&u_dofs[0], &u_dofs[u.Size()], 0);
    iota(&p_dofs[0], &p_dofs[p.Size()], u.Size());

    // Create Alfven operator
    AlfvenOperator oper(CG_u, CG_p, dt);

    // RHS vector
    Vector b(ssize);

    // Vectors for constrained system (with BC)
    Vector X(ssize);
    Vector B(ssize);
    Operator *A_constrained = nullptr;

    // Time loop
    double t = 0.0;
    real_t energy;
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
        cout << "step:\t" << jj << "\tt = " << t << "\tenergy = " << energy << endl;

        // Export data to ParaView
        Nit++;
        pd->SetCycle(Nit);
        pd->SetTime(t);
        pd->Save();
    }

    // Compute errors
    // u:         FE solution 
    // u_0_coeff: reference/exact solution 
    real_t up_err = u.ComputeL2Error(u_0_coeff);
    cout << "u   L2 error " << up_err << endl;

    real_t p_err = p.ComputeL2Error(p_0_coeff);
    cout << "phi L2 error " << p_err << endl;

    // Free memory
    delete pd;

    // Timer
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time: " << elapsed.count() << " seconds\n";

    return 0;
}
