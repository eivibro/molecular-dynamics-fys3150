#include "eulercromer.h"
#include "../system.h"

void EulerCromer::integrate(System *system, double dt)
{
    system->calculateForces();
    for(Atom *atom : system->atoms()) {
        atom->velocity += atom->force*dt / atom->mass();
        //std::cout << "Force: " << atom->force << std::endl;
        atom->position += atom->velocity*dt;
    }
    system->applyPeriodicBoundaryConditions();
}
