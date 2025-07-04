#include <iostream>
#include <fstream>

// #include "pulse.hpp"
#include "simulation_params.hpp"

SimulationParams load_config(const std::string& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Ошибка: не удалось открыть файл " + path);
    }

    nlohmann::json j;
    file >> j;

    return j.get<SimulationParams>();
}

int main() {
    try {
        auto params = load_config("cfg/params.json");
        std::cout << "Центральная длина волны: " << params.central_wavelength << " м\n";
        std::cout << "Test " << params.number_points + 1.5 << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Ошибка конфигурации: " << e.what() << '\n';
    }

    return 0;
}