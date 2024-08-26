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
    qTIP4PForce(bool tabulate=false, double dmin=0.0, double dmax=0.0, int prec_N=0, int rescaling_N=0);

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

    bool qTIP4PForce::isTabulated() const;

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
