%module qtip4pforceplugin

%import(module="openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
};

%{
#include "qTIP4PForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import openmm as mm
import openmm.unit as unit
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


namespace qTIP4PPlugin {

class qTIP4PForce : public OpenMM::Force {
public:
    qTIP4PForce();

    int qTIP4PForce::getNumWaters() const;

    int qTIP4PForce::addWater(int particle_O, int particle_H1, int particle_H2, int particle_M);

    /*
     * The reference parameters to this function are output values.
     * Marking them as such will cause swig to return a tuple.
    */
    %apply int& OUTPUT {int& particle_O};
    %apply int& OUTPUT {int& particle_H1};
    %apply int& OUTPUT {int& particle_H2};
    %apply int& OUTPUT {int& particle_M};
    void qTIP4PForce::getWater(int index, int& particle_O, int& particle_H1, int& particle_H2, int& particle_M) const;
    %clear int& particle_O;
    %clear int& particle_H1;
    %clear int& particle_H2;
    %clear int& particle_M;

    /*
     * Add methods for casting a Force to an qTIP4PForce.
    */
    %extend {
        static qTIP4PPlugin::qTIP4PForce& cast(OpenMM::Force& force) {
            return dynamic_cast<qTIP4PPlugin::qTIP4PForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<qTIP4PPlugin::qTIP4PForce*>(&force) != NULL);
        }
    }
};

}

%pythoncode %{

import warnings

def replaceWaters(system: mm.System, topology: mm.app.Topology):
    """Replaces all the water-water and intra-water forces in the system with the qTIP4P force plugin. Does not affect other water-system interactions

    Parameters
    ----------
    system : System
        System where the water forces should be replaced.
    topology : Topology
        Topology used to create the system.
    """

    warnings.warn('\'replaceWaters\' is not the most robust function, double-checking that accurate forces are being calculated is encouraged.')

    qtip4pforce = qTIP4PForce()
    system.addForce(qtip4pforce)

    energy_expression_OH_bond = 'D*(a^2*(r-r0)^2-a^3*(r-r0)^3+(7/12)*a^4*(r-r0)^4)'

    ind_list_O = []
    ind_list_H1 = []
    ind_list_H2 = []
    ind_list_M = []

    for res in topology.residues():
        if res.name != 'HOH':
            continue
        
        ind_O = None
        ind_H1 = None
        ind_H2 = None
        ind_M = None

        for atom in res.atoms():
            if atom.name == 'O':
                ind_O = atom.index
            elif atom.name == 'H1':
                ind_H1 = atom.index
            elif atom.name == 'H2':
                ind_H2 = atom.index
            elif atom.name == 'M':
                ind_M = atom.index
            else:
                raise Exception(f'HOH residue contains unknown particle: {atom.name}')
        
        if ind_M is None:
            raise Exception('One or more HOH residues does not contain the M virtual site required by qTIP4P')
        
        qtip4pforce.addWater(ind_O, ind_H1, ind_H2, ind_M)
        
        ind_list_O.append(ind_O)
        ind_list_H1.append(ind_H1)
        ind_list_H2.append(ind_H2)
        ind_list_M.append(ind_M)

    ind_list_water = ind_list_O + ind_list_H1 + ind_list_H2 + ind_list_M

    # remove custom bond force interaction (O-H bond)
    for i, force in enumerate(system.getForces()):
        if type(force) == mm.CustomBondForce and force.getEnergyFunction() == energy_expression_OH_bond:
            system.removeForce(i)
    
    # remove nonbonded force interactions (exclude all water-water non bonded interactions)
    for i, force in enumerate(system.getForces()):
        if type(force) != mm.NonbondedForce:
            continue
        
        nparticles_nonbonded = force.getNumParticles()

        # if there is the same number of particles in the NonBondedForce as there is atoms in the waters, just remove the force
        if len(ind_list_water) == nparticles_nonbonded:
            system.removeForce(i)
            break
        
        # otherwise create exceptions for all the water-water nonbonded forces
        for index1 in range(len(ind_list_water)):
            for index2 in range(index2, len(ind_list_water)):
                particle1 = ind_list_water[index1]
                particle2 = ind_list_water[index2]
                force.addException(particle1, particle2, 0, 0, 0, True)
        
    # remove hamonic angle force (H-O-H angle)
    for i, force in enumerate(system.getForces()):
        if type(force) != mm.HarmonicAngleForce:
            continue
        
        nangles_harmonic = force.getNumAngles()
        
        # otherwise reconstruct the harmonic angle force without the water angles
        harmonic_angle_force = mm.HarmonicAngleForce()

        for j in range(nangles_harmonic):
            particle1, particle2, particle3, angle, k = force.getAngleParameters(j)

            # we only need to check 1 particle
            if (particle1 in ind_list_water):
                continue

            harmonic_angle_force.addAngle(particle1, particle2, particle3, angle, k)
        
        # if this force did not contain any water harmonic angle forces the go to the next force and don't touch this one
        if harmonic_angle_force.getNumAngles() == nangles_harmonic:
            break

        # remove the old harmonic angle force
        system.removeForce(i)

        # add the new force to the system only if there was any angles added
        if harmonic_angle_force.getNumAngles() > 0:

            # setting some stuff to be equal
            harmonic_angle_force.setUsesPeriodicBoundaryConditions(force.usesPeriodicBoundaryConditions())
            harmonic_angle_force.setForceGroup(force.getForceGroup())
            harmonic_angle_force.setName(force.getName())

            system.addForce(harmonic_angle_force)

%}
