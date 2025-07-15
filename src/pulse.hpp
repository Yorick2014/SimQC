#pragma once

#include <vector>
#include "simulation_params.hpp"
#include "constants.hpp"

class Pulse {
    public:
    // spectrum
    std::vector<double> frequency;
    std::vector<double> intensity;

    // time
    std::vector<double> time;
    std::vector<double> t_intensity;
    
    void get_spectrum(const Laser &laser);
    void spectrum_to_csv(const std::string& path);
    void spectrum_to_time_domain(const Laser &laser, const QuantumChannel &channel);
        
};
double gaussian_spectrum(double nu, double nu0, double sigma_nu);