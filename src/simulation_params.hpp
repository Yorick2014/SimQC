#pragma once

#include "nlohmann/json.hpp"
#include <fstream>

template<typename T>
T load_config(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Ошибка: не удалось открыть файл " + path);
    }

    nlohmann::json json;
    file >> json;

    return json.get<T>();
}

struct Laser{
    double central_wavelength;
    double pulse_duration;
    double avg_count_photons;
    unsigned int number_points; //число точек для спектра
    double repeat_rate; // частота повторения импульсов
};
struct QuantumChannel{
    double channel_length;
    double chromatic_dispersion;
    double channel_attenuation;
    bool isAtt;
    bool isCromDisp; //хроматическая дисперсия
};
struct Photodetector{
    double quantum_efficiency;
    double dead_time;
    double time_slot;
};

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Laser,
    central_wavelength,
    pulse_duration,
    avg_count_photons,
    number_points,
    repeat_rate
)
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(QuantumChannel,
    channel_length,
    chromatic_dispersion,
    channel_attenuation,
    isAtt,
    isCromDisp
)

NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Photodetector,
    quantum_efficiency,
    dead_time,
    time_slot
)
