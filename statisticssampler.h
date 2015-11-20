#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
#include <fstream>
class System;
class StatisticsSampler
{
private:
    double m_kineticEnergy = 0;
    double m_potentialEnergy = 0;
    double m_temperature = 0;
    double m_density = 0;
    double const m_temperatureCoefficient = 3.0/2.0/500.0/(1.381e-23);
    double m_diffusionConstant = 0;
    double m_ri2;
    std::ofstream file;


public:
    StatisticsSampler();
    void saveToFile(System &system);
    void open(const char *filename);
    void close();
    void sample(System &system);
    void sampleKineticEnergy(System &system);
    void samplePotentialEnergy(System &system);
    void sampleTemperature(System &system);
    void sampleDensity(System &system);
    void sampleDiffusion(System &system);
    double kineticEnergy() { return m_kineticEnergy; }
    double potentialEnergy() { return m_potentialEnergy; }
    double totalEnergy() { return m_kineticEnergy+m_potentialEnergy; }
    double temperature() { return m_temperature; }
    double density() { return m_density; }
    double diffusionConstant() { return m_diffusionConstant;}
};
#endif
