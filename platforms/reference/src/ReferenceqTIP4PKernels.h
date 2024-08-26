#ifndef REFERENCE_EXAMPLE_KERNELS_H_
#define REFERENCE_EXAMPLE_KERNELS_H_

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

#include "qTIP4PKernels.h"
#include "openmm/Platform.h"
#include "openmm/reference/ReferencePlatform.h"
#include <vector>
#include <cmath>

using namespace qTIP4PPlugin;
using namespace OpenMM;
using namespace std;

/**
 * This kernel is invoked by qTIP4PForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcqTIP4PForceKernel : public CalcqTIP4PForceKernel {
public:
    ReferenceCalcqTIP4PForceKernel(string name, const Platform& platform) : CalcqTIP4PForceKernel(name, platform) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the qTIP4PForce this kernel will be used for
     */
    void initialize(const System& system, const qTIP4PForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the qTIP4PForce to copy the parameters from
     */
private:

    double executeAnalytical(ContextImpl& context);

    double executeTabulated(ContextImpl& context);

    void OHBondEnergyForces(Vec3& r1, Vec3& r2, double& d1, double& d2, Vec3& F_O, Vec3& F_H1, Vec3& F_H2, double& Vxyz);

    void HOHAngleEnergyForces(Vec3& r1, Vec3& r2, double& d1, double& d2, Vec3& F_O, Vec3& F_H1, Vec3& F_H2, double& Vxyz);
    
    void LennardJonesEnergyForces(int& atomA, int& atomB, Vec3& r, const double& C, const double& D, vector<Vec3>& forces, double& energy);

    void CoulombicEnergyForces(int& atomA, int& atomB, Vec3& r, const double& keqq, vector<Vec3>& forces, double& energy);

    void tabulateHH_MH();

    void tabulatedCoulombicEnergyForcesHH(int& atomA, int& atomB, Vec3& r, vector<Vec3>& forces, double& energy);

    void tabulatedCoulombicEnergyForcesMH(int& atomA, int& atomB, Vec3& r, vector<Vec3>& forces, double& energy);

    double d_to_q(double d) {return d*d;}
    double q_to_d(double q) {return sqrt(q);}

    double d_to_k(double d) {return d;}
    double k_to_d(double k) {return k;}

    double q_to_k(double q) {return d_to_k(q_to_d(q));}
    double k_to_q(double k) {return d_to_q(k_to_d(k));}

    int index_4d_to_1d(int i, int j, int k, int l, int N) {return l+N*(k+N*(j+N*i));}    
    
    bool tabulated;

    std::vector<double> tabulated_V_HH;
    std::vector<double> tabulated_F_HH;

    std::vector<double> tabulated_V_MH;
    std::vector<double> tabulated_F_MH;
    
    std::vector<int> tabulated_rescaling_HH_MH;

    double one_by_dq_HH_MH, qmin;
    double dmin, dmax;
    int prec_N, rescaling_N;
    
    
    std::vector<int> particles_O, particles_H1, particles_H2, particles_M;

    const double angle0=1.87448361664;
    const double k=367.6;

    const double r0=0.09419;
    const double D=485.7;
    const double a=22.87;

    const double Astr=D*a*a;
    const double Bstr=-D*a*a*a;
    const double Cstr=7./12.*D*a*a*a*a;

    //double charge_O=0;
    const double sigma_O=0.31589;
    const double epsilon_O=0.7749;

    const double charge_H=0.5564;
    //double sigma_H=1;
    //double epsilon_H=0;

    const double charge_M=-1.1128;
    //double sigma_M=1;
    //double epsilon_M=0;

    const double AOO = 4*epsilon_O*pow(sigma_O, 12);
    const double BOO = 4*epsilon_O*pow(sigma_O, 6);
    const double COO = sigma_O*sigma_O;
    const double DOO = 4*epsilon_O;
    const double ke = 1/(4*M_PI*EPSILON0);
    const double keqHqH = ke*charge_H*charge_H;
    const double keqHqM = ke*charge_H*charge_M;
    const double keqMqM = ke*charge_M*charge_M;
};

#endif /*REFERENCE_EXAMPLE_KERNELS_H_*/
