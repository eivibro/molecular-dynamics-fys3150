#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System()
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    for(Atom *atom : m_atoms){
        if(this->m_systemSize[0]<atom->position[0]){atom->position[0]-=this->m_systemSize[0];}
        if(this->m_systemSize[1]<atom->position[1]){atom->position[1]-=this->m_systemSize[1];}
        if(this->m_systemSize[2]<atom->position[2]){atom->position[2]-=this->m_systemSize[2];}
        if(0>atom->position[0]){atom->position[0]+=this->m_systemSize[0];}
        if(0>atom->position[1]){atom->position[1]+=this->m_systemSize[1];}
        if(0>atom->position[2]){atom->position[2]+=this->m_systemSize[2];}

    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    vector<double> momentum(3,0.0);
    for(Atom *atom : m_atoms){
        momentum[0] += atom->mass()*atom->velocity[0];
        momentum[1] += atom->mass()*atom->velocity[1];
        momentum[2] += atom->mass()*atom->velocity[2];
    }
    int numberOfAtoms = m_atoms.size();
    for(Atom *atom : m_atoms){
        atom->velocity.setX(atom->velocity.x()-momentum[0]/atom->mass()/numberOfAtoms);
        atom->velocity.setY(atom->velocity.y()-momentum[1]/atom->mass()/numberOfAtoms);
        atom->velocity.setZ(atom->velocity.z()-momentum[2]/atom->mass()/numberOfAtoms);
    }
    for(int i = 0; i < 3; i++){
        momentum[i] = 0;
    }
    for(Atom *atom : m_atoms){
        momentum[0] += atom->mass()*atom->velocity[0];
        momentum[1] += atom->mass()*atom->velocity[1];
        momentum[2] += atom->mass()*atom->velocity[2];
    }

}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).
    for(int i=0; i<15; i++) {
        Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
        double x = Random::nextDouble()*10; // random number in the interval [0,10]
        double y = Random::nextDouble()*10;
        double z = Random::nextDouble()*10;
        atom->position.set(x,y,z);
        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);
    }
    setSystemSize(vec3(10, 10, 10));
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_time += dt;
}
