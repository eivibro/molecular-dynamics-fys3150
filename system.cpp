#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include <cmath>

//void System::setRi2(double ri2)
//{
//    m_ri2 = ri2;
//}

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
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_
    //conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    for(Atom *atom : m_atoms){
        if(m_systemSize[0]<atom->position[0]){
            atom->position[0]-=m_systemSize[0];
            atom->initialPosition[0] -= m_systemSize[0];
        }
        if(m_systemSize[1]<atom->position[1]){
            atom->position[1]-=m_systemSize[1];
            atom->initialPosition[1] -= m_systemSize[1];
        }
        if(m_systemSize[2]<atom->position[2]){
            atom->position[2]-=m_systemSize[2];
            atom->initialPosition[2] -= m_systemSize[2];
        }
        if(0>atom->position[0]){
            atom->position[0]+=m_systemSize[0];
            atom->initialPosition[0] += m_systemSize[0];
        }
        if(0>atom->position[1]){
            atom->position[1]+=m_systemSize[1];
            atom->initialPosition[1] += m_systemSize[1];
        }
        if(0>atom->position[2]){
            atom->position[2]+=m_systemSize[2];
            atom->initialPosition[2] += m_systemSize[2];
        }
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so
    //the total momentum becomes zero.
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
}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant,
                              double temperature) {
    std::cout << "Temp: " << temperature << std::endl;
    for(int i=0; i<numberOfUnitCellsEachDimension; i++){
        for(int j=0; j<numberOfUnitCellsEachDimension; j++){
            for(int k=0; k<numberOfUnitCellsEachDimension; k++){
                Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                double x = i*latticeConstant;
                double y = j*latticeConstant;
                double z = k*latticeConstant;
                atom->position.set(x,y,z);
                atom->initialPosition.set(x,y,z);
                atom->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom);
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = i*latticeConstant+latticeConstant/2;
                y = j*latticeConstant+latticeConstant/2;
                z = k*latticeConstant;
                atom1->position.set(x,y,z);
                atom1->initialPosition.set(x,y,z);
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);
                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = i*latticeConstant;
                y = j*latticeConstant+latticeConstant/2;
                z = k*latticeConstant+latticeConstant/2;
                atom2->position.set(x,y,z);
                atom2->initialPosition.set(x,y,z);
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);
                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                x = i*latticeConstant+latticeConstant/2;
                y = j*latticeConstant;
                z = k*latticeConstant+latticeConstant/2;
                atom3->position.set(x,y,z);
                atom3->initialPosition.set(x,y,z);
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);
            }
        }
    }
    double length = numberOfUnitCellsEachDimension*latticeConstant;//+latticeConstant/2;
    setSystemSize(vec3(length, length, length));
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
