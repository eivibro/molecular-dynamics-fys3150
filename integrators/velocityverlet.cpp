#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"

void VelocityVerlet::integrate(System *system, double dt)
{
    vector<vec3> vsAtHalfDt;
    int i = 0;
    for(Atom *atom : system->atoms()){
        double dtOver2m = dt/(2*atom->mass());
        vec3 vAtHalfDt = atom->velocity+atom->force*dtOver2m;
        vsAtHalfDt.push_back(vAtHalfDt);
        atom->position += vAtHalfDt*dt;
        atom->velocity = vAtHalfDt+atom->force*dtOver2m;
    }
    system->calculateForces();
    for(Atom *atom : system->atoms()){
        double dtOver2m = dt/(2*atom->mass());
        atom->velocity = vsAtHalfDt[i]+atom->force*dtOver2m;
        i++;
    }
    system->applyPeriodicBoundaryConditions();
}
