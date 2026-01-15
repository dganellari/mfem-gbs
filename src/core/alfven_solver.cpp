#include "utils/functions.hpp"
#include "alfven_solver.hpp"

namespace mfem {

AlfvenOperator::AlfvenOperator(FiniteElementSpace &fespace_u_,
                               FiniteElementSpace &fespace_p_,
                               real_t dt)
    : fespace_u(fespace_u_), fespace_p(fespace_p_),
      M(NULL), N_par(NULL), N_full(NULL), blf_E(NULL), 
      F_mat(NULL), A(NULL)
{
    // AssembleSystem needs the updated dt_over_two
    dt_over_two = dt / 2.0;

    // Get essential boundary true dofs for pressure
    fespace_p.GetBoundaryTrueDofs(ess_tdof_p);
    // Offset by number of u dofs (p is the second block)
    for (int i = 0; i < ess_tdof_p.Size(); i++) {
        ess_tdof_p[i] += fespace_u.GetNDofs();
    }

    // Setup block offsets
    offsets.SetSize(3);
    offsets[0] = 0;
    offsets[1] = fespace_u.GetNDofs();
    offsets[2] = fespace_p.GetNDofs();
    offsets.PartialSum();
    
    AssembleSystem();
}

void AlfvenOperator::AssembleSystem()
{
    int dim = fespace_u.GetMesh()->Dimension();
    
    // Matrix M: mass matrix for u
    M = new BilinearForm(&fespace_u);
    M->AddDomainIntegrator(new MassIntegrator());
    M->Assemble();
    M->Finalize();
    M_mat = M->SpMat();
    M_mat.Finalize();
    
    // Magnetic field coefficient
    VectorFunctionCoefficient b_gfcoeff(dim, bfield);
    OuterProductCoefficient K(b_gfcoeff, b_gfcoeff);
    
    // Matrix N_par: parallel diffusion (b·∇p, b·∇q)
    N_par = new BilinearForm(&fespace_p);
    N_par->AddDomainIntegrator(new DiffusionIntegrator(K));
    N_par->Assemble();
    N_par->Finalize();
    SparseMatrix N_par_mat(N_par->SpMat());
    N_par_mat.Finalize();
    
    // Matrix N_full: full diffusion (∇p, ∇q)
    N_full = new BilinearForm(&fespace_p);
    N_full->AddDomainIntegrator(new DiffusionIntegrator());
    N_full->Assemble();
    N_full->Finalize();
    SparseMatrix N_full_mat(N_full->SpMat());
    
    // Matrix N_n: perpendicular diffusion = -(N_full - N_par)
    N_n_mat = N_full_mat;
    N_n_mat.Add(-1.0, N_par_mat); // (∇p,∇q) - (b·∇p, b·∇q)
    N_n_mat *= -1.0;
    N_n_mat.Finalize();
    
    // Matrix E: mixed form (b·∇p, u)
    blf_E = new MixedBilinearForm(&fespace_p, &fespace_u);
    blf_E->AddDomainIntegrator(new MixedDirectionalDerivativeIntegrator(b_gfcoeff));
    blf_E->Assemble();
    blf_E->Finalize();
    E_mat = blf_E->SpMat();
    E_mat *= -1.0 * dt_over_two;
    E_mat.Finalize();
    
    // Matrix F: transpose of E
    F_mat = Transpose(E_mat);
    F_mat->Finalize();
    
    // Build block operator
    // [ M   E ]
    // [ F  N_n]
    A = new BlockOperator(offsets);
    A->SetBlock(0, 0, &M_mat);
    A->SetBlock(0, 1, &E_mat);
    A->SetBlock(1, 0,  F_mat);
    A->SetBlock(1, 1, &N_n_mat);
}

void AlfvenOperator::FormRHS(const Vector &u_old, const Vector &p_old, Vector &b) const
{
    int u_size = fespace_u.GetNDofs();
    int p_size = fespace_p.GetNDofs();
    
    b = 0.0;
    Vector bsub1(u_size);
    Vector bsub2(p_size);
    
    bsub1 = 0.0;
    bsub2 = 0.0;
    
    // First block: M*u_old - E*p_old
    M_mat.AddMult(u_old, bsub1);
    E_mat.AddMult(p_old, bsub1, -1.0);
    
    // Second block: -F*u_old + N_n*p_old
    F_mat ->AddMult(u_old, bsub2, -1.0);
    N_n_mat.AddMult(p_old, bsub2,  1.0);
    
    // Add subvectors to b (which is zero, so this sets them)
    b.AddSubVector(bsub1, 0);
    b.AddSubVector(bsub2, u_size);
}

real_t AlfvenOperator::ComputeEnergy(const Vector &u, const Vector &p)
{
    real_t energy =      M->InnerProduct(u, u) +
                    N_full->InnerProduct(p, p) -
                     N_par->InnerProduct(p, p);
    return energy;
}

AlfvenOperator::~AlfvenOperator()
{
    delete A;
    delete F_mat;
    delete blf_E;
    delete N_full;
    delete N_par;
    delete M;
}

} // namespace mfem