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

#include "qTIP4PForce.h"
#include "internal/qTIP4PForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"

using namespace qTIP4PPlugin;
using namespace OpenMM;
using namespace std;

qTIP4PForce::qTIP4PForce() {
}

int qTIP4PForce::getNumWaters() const {
    return particles_O.size();
}

int qTIP4PForce::addWater(int particle_O, int particle_H1, int particle_H2, int particle_M) {
    particles_O.push_back(particle_O);
    particles_H1.push_back(particle_H1);
    particles_H2.push_back(particle_H2);
    particles_M.push_back(particle_M);
    return particles_O.size()-1;
}

void qTIP4PForce::getWater(int index, int& particle_O, int& particle_H1, int& particle_H2, int& particle_M) const {
    ASSERT_VALID_INDEX(index, particles_O);
    particle_O = particles_O[index];
    particle_H1 = particles_H1[index];
    particle_H2 = particles_H2[index];
    particle_M = particles_M[index];
}

void qTIP4PForce::getParticles(std::vector<int>& particles_O, std::vector<int>& particles_H1, std::vector<int>& particles_H2, std::vector<int>& particles_M) const {
    particles_O = this->particles_O;
    particles_H1 = this->particles_H1;
    particles_H2 = this->particles_H2;
    particles_M = this->particles_M;
}

ForceImpl* qTIP4PForce::createImpl() const {
    return new qTIP4PForceImpl(*this);
}
