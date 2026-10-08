#include "diffusionSystem.h"
#include <assert.h>

DiffusionSystem::DiffusionSystem()
    : yInf(0)
    , wallConst(0)
{
}

void DiffusionSystem::get_A(dvec& a, dvec& b, dvec& c)
{
    assert(mathUtils::notnan(D));
    assert(mathUtils::notnan(B));

    for (size_t j=1; j<=N-2; j++) {
        c1[j] = 0.5*B[j]/(dlj[j]*r[j]);
        c2[j] = rphalf[j]*(D[j]+D[j+1])/hh[j];
    }
    // Conductance of the first face (x[0]..x[1]), used by a[1] below. The
    // j = 0 boundary rows already exchange flux across this face with j = 1;
    // without c2[0], j = 1 never sees the other side of that exchange
    // (a[1] = 0), so the flux is created/destroyed there instead of
    // conserved. Harmless for a zero-gradient far field, but at a burner
    // face it drained ~1/3 of the inlet H and cut off the wall heat loss.
    c2[0] = rphalf[0]*(D[0]+D[1])/hh[0];

    assert(mathUtils::notnan(c1));
    assert(mathUtils::notnan(c2));

    // Left boundary value
    size_t jStart;
    if (grid.leftBC == BoundaryCondition::FixedValue) {
        jStart = 1;
        // Clear coefficients left over if this system used a flux BC before
        b[0] = 0;
        c[0] = 0;
    } else if (grid.leftBC == BoundaryCondition::ControlVolume ||
               grid.leftBC == BoundaryCondition::InletFlux) {
        jStart =  1;
        double c0 = B[0] * (grid.alpha + 1) * (D[0]+D[1]) / (2 * hh[0] * hh[0]);
        b[0] = -c0;
        c[0] = c0;
    } else if (grid.leftBC == BoundaryCondition::WallFlux) {
        jStart = 1;
        double c0 = B[0] * (grid.alpha + 1) / hh[0];
        double d = 0.5 * (D[0]+D[1]);
        b[0] = - c0 * (d / hh[0] + wallConst);
        c[0] = d * c0 / hh[0];
    } else  { // (leftBC == BoundaryCondition::ZeroGradient)
        // In the case of a zero gradient boundary condition, the boundary value
        // is not computed, and the value one point in is computed by substituting
        // y[0] = y[1] in the finite difference formula.
        jStart = 2;
        b[1] = -c1[1]*c2[1];
        c[1] = c1[1]*c2[1];
    }

    // Right boundary value
    size_t jStop;
    if (grid.rightBC == BoundaryCondition::FixedValue) {
        jStop = N-1;
        // Clear coefficients left over if this system used a flux BC before
        a[N-1] = 0;
        b[N-1] = 0;
    } else if (grid.rightBC == BoundaryCondition::Wall) {
        // Zero-flux wall: half control volume from x[N-3/2] to x[N-1], whose
        // only flux is across the face shared (via c2[N-2]) with row N-2
        jStop = N-1;
        double c0 = B[N-1] * c2[N-2] / (r[N-1] * hh[N-2]);
        a[N-1] = c0;
        b[N-1] = -c0;
    } else { // (rightBC == BoundaryCondition::ZeroGradient)
        // In the case of a zero gradient boundary condition, the boundary value
        // is not computed, and the value one point in is computed by substituting
        // y[N-1] = y[N-2] in the finite difference formula.
        jStop = N-2;
        a[N-2] = c1[N-2]*c2[N-3];
        b[N-2] = -c1[N-2]*c2[N-3];
    }

    // Intermediate points
    for (size_t j=jStart; j<jStop; j++) {
        a[j] = c1[j]*c2[j-1];
        b[j] = -c1[j]*(c2[j-1] + c2[j]);
        c[j] = c1[j]*c2[j];
    }
}

void DiffusionSystem::get_k(dvec& k)
{
    assert(mathUtils::notnan(splitConst));
    k = splitConst;
    if (grid.leftBC == BoundaryCondition::WallFlux) {
        k[0] += B[0] * (grid.alpha + 1) / hh[0] * wallConst * yInf;
    }
    assert(mathUtils::notnan(k));
}

void DiffusionSystem::resize(size_t N_)
{
    N = N_;
    B.setConstant(N, NaN);
    D.setConstant(N, NaN);
    splitConst.setConstant(N, NaN);
    c1.setConstant(N, 0);
    c2.setConstant(N, 0);
}

void DiffusionSystem::resetSplitConstants()
{
    splitConst.setZero(N);
}
