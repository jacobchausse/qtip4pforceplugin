/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ReferenceqTIP4PKernels.h"
#include "qTIP4PForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/Vec3.h"
#include "openmm/reference/ReferencePlatform.h"
#include <cmath>
#include <vector>
#include <omp.h>

using namespace qTIP4PPlugin;
using namespace OpenMM;
using namespace std;


static vector<Vec3>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->positions);
}

static vector<Vec3>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<Vec3>*) data->forces);
}

void ReferenceCalcqTIP4PForceKernel::initialize(const System& system, const qTIP4PForce& force) {
    force.getParticles(particles_O, particles_H1, particles_H2, particles_M);
}

double ReferenceCalcqTIP4PForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<Vec3>& pos = extractPositions(context);
    vector<Vec3>& force = extractForces(context);

    int numWaters = particles_O.size();

    double energy = 0, V;

    Vec3 F_O, F_H1, F_H2, F_HM, F_MM;
    
    for (int index1 = 0; index1 < numWaters; index1++) {
        // intramolecular forces
        int water1_O = particles_O[index1];
        int water1_H1 = particles_H1[index1];
        int water1_H2 = particles_H2[index1];
        int water1_M = particles_M[index1];

        Vec3 pos1_O = pos[water1_O];
        Vec3 pos1_H1 = pos[water1_H1];
        Vec3 pos1_H2 = pos[water1_H2];
        Vec3 pos1_M = pos[water1_M];

        Vec3 r1 = pos1_H1 - pos1_O;
        Vec3 r2 = pos1_H2 - pos1_O;

        double d1 = sqrt(r1.dot(r1));
        double d2 = sqrt(r2.dot(r2));
        
        OHBondEnergyForces(r1, r2, d1, d2, F_O, F_H1, F_H2, V);

        force[water1_O] += F_O;
        force[water1_H1] += F_H1;
        force[water1_H2] += F_H2;

        energy += V;
        
        HOHAngleEnergyForces(r1, r2, d1, d2, F_O, F_H1, F_H2, V);

        force[water1_O] += F_O;
        force[water1_H1] += F_H1;
        force[water1_H2] += F_H2;

        energy += V;
        
        for (int index2 = index1 + 1; index2 < numWaters; index2++) {
            
            // intermolecular forces
            
            int water2_O = particles_O[index2];
            int water2_H1 = particles_H1[index2];
            int water2_H2 = particles_H2[index2];
            int water2_M = particles_M[index2];

            Vec3 pos2_O = pos[water2_O];
            Vec3 pos2_H1 = pos[water2_H1];
            Vec3 pos2_H2 = pos[water2_H2];
            Vec3 pos2_M = pos[water2_M];

            // LJ OO
            Vec3 rOO = pos2_O - pos1_O;
            LennardJonesEnergyForces(water2_O, water1_O, rOO, COO, DOO, force, energy);
            
            // Coulombic HH
            Vec3 rH1H1 = pos2_H1 - pos1_H1;
            CoulombicEnergyForces(water2_H1, water1_H1, rH1H1, keqHqH, force, energy);

            Vec3 rH1H2 = pos2_H1 - pos1_H2;
            CoulombicEnergyForces(water2_H1, water1_H2, rH1H2, keqHqH, force, energy);

            Vec3 rH2H1 = pos2_H2 - pos1_H1;
            CoulombicEnergyForces(water2_H2, water1_H1, rH2H1, keqHqH, force, energy);

            Vec3 rH2H2 = pos2_H2 - pos1_H2;
            CoulombicEnergyForces(water2_H2, water1_H2, rH2H2, keqHqH, force, energy);

            // Coulombic HM
            Vec3 rH1M = pos2_H1 - pos1_M;
            CoulombicEnergyForces(water2_H1, water1_M, rH1M, keqHqM, force, energy);

            Vec3 rH2M = pos2_H2 - pos1_M;
            CoulombicEnergyForces(water2_H2, water1_M, rH2M, keqHqM, force, energy);

            Vec3 rMH1 = pos2_M - pos1_H1;
            CoulombicEnergyForces(water2_M, water1_H1, rMH1, keqHqM, force, energy);

            Vec3 rMH2 = pos2_M - pos1_H2;
            CoulombicEnergyForces(water2_M, water1_H2, rMH2, keqHqM, force, energy);

            // Coulombic MM
            Vec3 rMM = pos2_M - pos1_M;
            CoulombicEnergyForces(water2_M, water1_M, rMM, keqMqM, force, energy);
        }
        
    }
    
    return energy;
}


void ReferenceCalcqTIP4PForceKernel::OHBondEnergyForces(Vec3& r1, Vec3& r2, double& d1, double& d2, Vec3& F_O, Vec3& F_H1, Vec3& F_H2, double& Vxyz) {

    double d1d = d1-r0;
    double d1d2 = d1d*d1d;
    double d1d3 = d1d*d1d2;
    double d1d4 = d1d*d1d3;

    double d2d = d2-r0;
    double d2d2 = d2d*d2d;
    double d2d3 = d2d*d2d2;
    double d2d4 = d2d*d2d3;

    Vxyz = Astr*(d1d2+d2d2) + Bstr*(d1d3+d2d3) + Cstr*(d1d4+d2d4);

    F_H1 = -(2*Astr*d1d + 3*Bstr*d1d2 + 4*Cstr*d1d3)*r1/d1;
    F_H2 = -(2*Astr*d2d + 3*Bstr*d2d2 + 4*Cstr*d2d3)*r2/d2;
    F_O = -(F_H1 + F_H2);

}

void ReferenceCalcqTIP4PForceKernel::HOHAngleEnergyForces(Vec3& r1, Vec3& r2, double& d1, double& d2, Vec3& F_O, Vec3& F_H1, Vec3& F_H2, double& Vxyz) {
    
    double cosangle = r1.dot(r2)/(d1*d2);
    double angle = acos(cosangle);

    Vxyz = 0.5*k*pow(angle-angle0, 2);

    double factor = k*(angle-angle0)/sin(angle);

    F_H1 = factor*(r2/d2-r1/d1*cosangle)/d1;
    F_H2 = factor*(r1/d1-r2/d2*cosangle)/d2;
    F_O = -(F_H1 + F_H2);
}

void ReferenceCalcqTIP4PForceKernel::LennardJonesEnergyForces(int& atomA, int& atomB, Vec3& r, const double& C, const double& D, vector<Vec3>& forces, double& energy) {
    double d2inv = 1./r.dot(r);
    double beta = C*d2inv;

    double beta3 = beta*beta*beta;
    double beta6 = beta3*beta3;
    
    energy += D*(beta3 - 1.0)*beta3;;

    double dVdr = D*(12.0*beta3 - 6.0)*beta3*d2inv;

    for (int k=0; k<3; k++){
        double force = dVdr*r[k];
        forces[atomA][k] += force;
        forces[atomB][k] -= force;
    }
}

void ReferenceCalcqTIP4PForceKernel::CoulombicEnergyForces(int& atomA, int& atomB, Vec3& r, const double& keqq, vector<Vec3>& forces, double& energy) {
    double d2inv = 1./r.dot(r);
    double dinv = sqrt(d2inv);
    double V = keqq*dinv;

    energy += V;

    double dVdr = V*d2inv;

    for (int k=0; k<3; k++){
        double force = dVdr*r[k];
        forces[atomA][k] += force;
        forces[atomB][k] -= force;
    }
}