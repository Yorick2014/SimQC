#pragma once

#include <vector>
#include "simulation_params.hpp"
#include "constants.hpp"

class Pulse {
    public:
    std::vector<double> frequency;
    std::vector<double> intensity;
    Laser laser;
    
    void get_spectrum();
    // void print_frequency();
    // void print_intensity();
        
};
double gaussian_spectrum(double nu, double nu0, double sigma_nu);