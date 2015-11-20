#include "lennardjones.h"
#include <iostream>
LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{
    m_potentialEnergy = 0; // Remember to compute this in the loop
    vector<Atom *> atoms = system->atoms(); //Fetching vector of atoms
    double length = system->systemSize()[0]; //Cubic system

    //Declaring constants to reduce number of computations in the loops
    double halfLength = length/2;
    double sigmaPower6 = m_sigma*m_sigma*m_sigma*m_sigma*m_sigma*m_sigma;
    double sigmaPower12 = sigmaPower6*sigmaPower6;
    double sigmaPower12Times2 = 2.0*sigmaPower12;
    double epsilonTimes4 = 4*m_epsilon;
    double epsilonTimes24 = 24*m_epsilon;

    double rCut = 2.5*m_sigma;
    double potentialEnergyAtRCut = epsilonTimes4*(sigmaPower12/pow(rCut,12) - sigmaPower6/pow(rCut,6));

    for(int i = 0; i < atoms.size(); i++){
        for(int j = i+1; j < atoms.size(); j++){
            //The coordinate difference between the considered atoms
            double xij = atoms[i]->position[0]-atoms[j]->position[0];
            double yij = atoms[i]->position[1]-atoms[j]->position[1];
            double zij = atoms[i]->position[2]-atoms[j]->position[2];
            //Accounting for boundary conditions
            if(xij > halfLength){xij=xij-length;}
            if(xij < -halfLength){xij=xij+length;}
            if(yij > halfLength){yij=yij-length;}
            if(yij < -halfLength){yij=yij+length;}
            if(zij > halfLength){zij=zij-length;}
            if(zij < -halfLength){zij=zij+length;}

            //Constants to reduce number of computations
            double rij2 = xij*xij+yij*yij+zij*zij;
            if(rij2>rCut*rCut){continue;}

            double oneOverRij2 = 1.0/rij2;
            double oneOverRij6 = oneOverRij2*oneOverRij2*oneOverRij2;
            double oneOverRij12 = oneOverRij6*oneOverRij6;

            //Calculating forces
            double F = epsilonTimes24*(sigmaPower12Times2*oneOverRij12-sigmaPower6*oneOverRij6)*oneOverRij2;
            atoms[i]->force[0] += F*xij;
            atoms[i]->force[1] += F*yij;
            atoms[i]->force[2] += F*zij;
            atoms[j]->force[0] -= F*xij;
            atoms[j]->force[1] -= F*yij;
            atoms[j]->force[2] -= F*zij;

            //Calculating the potential energy
            m_potentialEnergy += epsilonTimes4*(sigmaPower12*oneOverRij12 - sigmaPower6*oneOverRij6) - potentialEnergyAtRCut;
        }
    }
}
