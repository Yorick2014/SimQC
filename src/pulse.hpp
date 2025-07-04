#pragma once

#include <vector>

class Pulse {

    public:
        struct Spectrum {
        std::vector<double> frequency;
        std::vector<double> intensity;
        };

        Spectrum get_spectrum();
};