#include "mfem.hpp"
#include "core/wave_solver.hpp"
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;

int main(int argc, char *argv[])
{
    // Parameters (defaults)
    std::string mesh_file = std::string(DATA_DIR) + "/star.mesh";
    int ref_levels = 1;
    int order = 2;
    int ode_solver_type = 10;
    real_t t_final = 0.1;
    real_t dt = 1.0e-2;
    real_t speed = 1.0;
    bool dirichlet = true;
    bool visualization = false;  // Add visualization option
    int vis_steps = 1;           // Visualization frequency

    // Parse command-line arguments
    OptionsParser args(argc, argv);
    args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use.");
    args.AddOption(&ref_levels, "-r", "--refine", "Number of refinements.");
    args.AddOption(&order, "-o", "--order", "Finite element order.");
    args.AddOption(&ode_solver_type, "-s", "--solver", "ODE solver type.");
    args.AddOption(&t_final, "-tf", "--t-final", "Final time.");
    args.AddOption(&dt, "-dt", "--time-step", "Time step size.");
    args.AddOption(&speed, "-c", "--speed", "Wave speed.");
    args.AddOption(&dirichlet, "-d", "--dirichlet", "-n", "--neumann", 
                   "Dirichlet boundary conditions.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization", "Enable/disable GLVis visualization.");
    args.AddOption(&vis_steps, "-vs", "--visualization-steps",
                   "Visualize every n-th timestep.");
    args.Parse();
    
    if (!args.Good())
    {
        args.PrintUsage(std::cout);
        return 1;
    }
    args.PrintOptions(std::cout);  // Print all options used

    cout << "\n=== MFEM Wave Mini-App ===" << endl;

    // Validate mesh file exists
    ifstream mesh_check(mesh_file);
    if (!mesh_check.good())
    {
        cerr << "Error: Mesh file '" << mesh_file << "' not found." << endl;
        return 1;
    }
    mesh_check.close();

    // Read mesh
    Mesh mesh(mesh_file.c_str(), 1, 1);
    int dim = mesh.Dimension();
    cout << "Mesh dimension: " << dim << endl;

    // Refine mesh
    for (int lev = 0; lev < ref_levels; lev++)
    {
        mesh.UniformRefinement();
    }
    cout << "Number of elements: " << mesh.GetNE() << endl;

    // FE space
    H1_FECollection fe_coll(order, dim);
    FiniteElementSpace fespace(&mesh, &fe_coll);
    cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;

    // Initial conditions
    GridFunction u_gf(&fespace);
    GridFunction dudt_gf(&fespace);
    
    FunctionCoefficient u_0(InitialSolution);
    u_gf.ProjectCoefficient(u_0);
    Vector u;
    u_gf.GetTrueDofs(u);

    FunctionCoefficient dudt_0(InitialRate);
    dudt_gf.ProjectCoefficient(dudt_0);
    Vector dudt;
    dudt_gf.GetTrueDofs(dudt);

    // Wave operator
    Array<int> ess_bdr(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0);
    if (dirichlet) ess_bdr = 1; else ess_bdr = 0;
    WaveOperator oper(fespace, ess_bdr, speed);

    // ODE solver
    SecondOrderODESolver *ode_solver = SecondOrderODESolver::Select(ode_solver_type);
    if (!ode_solver)
    {
        cerr << "Error: Invalid ODE solver type " << ode_solver_type << endl;
        return 1;
    }
    ode_solver->Init(oper);

    // Time integration
    real_t t = 0.0;
    int num_steps = static_cast<int>(ceil(t_final / dt));
    
    cout << "\n=== Starting time integration ===" << endl;
    cout << "Time steps: " << num_steps << ", dt: " << dt << endl;
    
    for (int ti = 1; ti <= num_steps; ti++)
    {
        ode_solver->Step(u, dudt, t, dt);
        
        if (ti % 10 == 0 || ti == num_steps)
        {
            cout << "Step " << ti << "/" << num_steps << ", t = " << t << endl;
        }
        
        oper.SetParameters(u);
        
        // TODO: Add visualization support here when needed
    }

    // Save final solution
    u_gf.SetFromTrueDofs(u);
    ofstream osol("final_solution.gf");
    osol.precision(8);
    u_gf.Save(osol);

    cout << "\n=== Simulation complete ===" << endl;
    cout << "Solution saved to final_solution.gf" << endl;

    delete ode_solver;

    return 0;
}