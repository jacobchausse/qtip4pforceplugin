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

    double energy = 0;

    int particle_O, particle_H1, particle_H2, particle_M;

    Vec3 F_O, F_H1, F_H2, r1, r2, pos1_O, pos1_H1, pos1_H2, pos1_M, pos2_O, pos2_H1, pos2_H2, pos2_M;
    double d1, d2, V;
    
    int water1_O, water1_H1, water1_H2, water1_M, water2_O, water2_H1, water2_H2, water2_M;
    Vec3 rOO, rHH, rHM, rMM, F_HM, F_MM;
    
    for (int index1 = 0; index1 < numWaters; index1++) {
        // intramolecular forces
        water1_O = particles_O[index1];
        water1_H1 = particles_H1[index1];
        water1_H2 = particles_H2[index1];
        water1_M = particles_M[index1];

        pos1_O = pos[water1_O];
        pos1_H1 = pos[water1_H1];
        pos1_H2 = pos[water1_H2];
        pos1_M = pos[water1_M];

        r1 = pos1_H1 - pos1_O;
        r2 = pos1_H2 - pos1_O;

        d1 = sqrt(r1.dot(r1));
        d2 = sqrt(r2.dot(r2));
        
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
            
            water2_O = particles_O[index2];
            water2_H1 = particles_H1[index2];
            water2_H2 = particles_H2[index2];
            water2_M = particles_M[index2];

            pos2_O = pos[water2_O];
            pos2_H1 = pos[water2_H1];
            pos2_H2 = pos[water2_H2];
            pos2_M = pos[water2_M];

            // LJ OO
            rOO = pos2_O - pos1_O;
            OOLennardJonesEnergyForces(rOO, F_O, V);
            force[water2_O] += F_O;
            force[water1_O] -= F_O;
            energy += V;
            
            // Coulombic HH
            rHH = pos2_H1 - pos1_H1;
            HHCoulombicEnergyForces(rHH, F_H1, V);
            force[water2_H1] += F_H1;
            force[water1_H1] -= F_H1;
            energy += V;

            rHH = pos2_H1 - pos1_H2;
            HHCoulombicEnergyForces(rHH, F_H1, V);
            force[water2_H1] += F_H1;
            force[water1_H2] -= F_H1;
            energy += V;

            rHH = pos2_H2 - pos1_H1;
            HHCoulombicEnergyForces(rHH, F_H2, V);
            force[water2_H2] += F_H2;
            force[water1_H1] -= F_H2;
            energy += V;

            rHH = pos2_H2 - pos1_H2;
            HHCoulombicEnergyForces(rHH, F_H2, V);
            force[water2_H2] += F_H2;
            force[water1_H2] -= F_H2;
            energy += V;

            // Coulombic HM
            rHM = pos2_H1 - pos1_M;
            HMCoulombicEnergyForces(rHM, F_HM, V);
            force[water2_H1] += F_HM;
            force[water1_M] -= F_HM;
            energy += V;

            rHM = pos2_H2 - pos1_M;
            HMCoulombicEnergyForces(rHM, F_HM, V);
            force[water2_H2] += F_HM;
            force[water1_M] -= F_HM;
            energy += V;

            rHM = pos2_M - pos1_H1;
            HMCoulombicEnergyForces(rHM, F_HM, V);
            force[water2_M] += F_HM;
            force[water1_H1] -= F_HM;
            energy += V;

            rHM = pos2_M - pos1_H2;
            HMCoulombicEnergyForces(rHM, F_HM, V);
            force[water2_M] += F_HM;
            force[water1_H2] -= F_HM;
            energy += V;

            // Coulombic MM
            rMM = pos2_M - pos1_M;
            MMCoulombicEnergyForces(rMM, F_MM, V);
            force[water2_M] += F_MM;
            force[water1_M] -= F_MM;
            energy += V;     
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

void ReferenceCalcqTIP4PForceKernel::OOLennardJonesEnergyForces(Vec3& rOO, Vec3& F_OO, double& Vxyz) {
    double dOO2 = rOO.dot(rOO);

    double dOO6inv = pow(dOO2, -3);
    double dOO12inv = dOO6inv*dOO6inv;
    
    Vxyz = AOO*dOO12inv - BOO*dOO6inv;

    F_OO = (12*dOO12inv*AOO - 6*dOO6inv*BOO)/dOO2*rOO;
}

void ReferenceCalcqTIP4PForceKernel::HHCoulombicEnergyForces(Vec3& rHH, Vec3& F_HH, double& Vxyz) {
    double dHH2inv = 1./rHH.dot(rHH);

    Vxyz = keqHqH*sqrt(dHH2inv);

    F_HH = Vxyz*dHH2inv*rHH;
}

void ReferenceCalcqTIP4PForceKernel::HMCoulombicEnergyForces(Vec3& rHM, Vec3& F_HM, double& Vxyz) {
    double dHM2inv = 1./rHM.dot(rHM);

    Vxyz = keqHqM*sqrt(dHM2inv);

    F_HM = Vxyz*dHM2inv*rHM;
}

void ReferenceCalcqTIP4PForceKernel::MMCoulombicEnergyForces(Vec3& rMM, Vec3& F_MM, double& Vxyz) {
    double dMM2inv = 1./rMM.dot(rMM);

    Vxyz = keqMqM*sqrt(dMM2inv);

    F_MM = Vxyz*dMM2inv*rMM;
}