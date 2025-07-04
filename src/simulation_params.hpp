#pragma once
#include "../external/json/json.hpp"

struct SimulationParams {
    double central_wavelength;
    double pulse_duration;
    int number_points;
    double dispersion;
    double attenuation;
    double quantum_efficiency;
    double repetition_rate;
    int detector_dark_count;
    double channel_length;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(SimulationParams,
    central_wavelength,
    pulse_duration,
    number_points,
    dispersion,
    attenuation,
    quantum_efficiency,
    repetition_rate,
    detector_dark_count,
    channel_length
)
