#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"
#include "unitconverter.h"

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::open(const char *filename) {
    if(file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file "
                  << filename << ", but some file is already open." << std::endl;
        exit(1);
    }

    file.open(filename);
}

void StatisticsSampler::close() {
    if(file.is_open()) {
        file.close();
    }
}

void StatisticsSampler::saveToFile(System &system)
{
    if(file.is_open()) {
        file << system.time() << "    " << m_kineticEnergy
                  << "    " << m_potentialEnergy << "    "
                  << m_temperature<< "    "<< m_diffusionConstant << std::endl;
    }
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusion(system);
    saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential()->potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    m_temperature = m_temperatureCoefficient*UnitConverter::energyToSI(m_kineticEnergy);
}

void StatisticsSampler::sampleDensity(System &system)
{

}

void StatisticsSampler::sampleDiffusion(System &system)
{
    m_ri2;
    m_diffusionConstant = 0;
    for(Atom *atom:system.atoms()){
        m_ri2 += (atom->position-atom->initialPosition).lengthSquared();
        m_diffusionConstant = m_ri2/(6*system.atoms().size()*system.time());
    }
}

