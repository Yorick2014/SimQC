#include <iostream>
#include <string>
#include "simulation_params.hpp"
#include "pulse.hpp"


int main() {
    std::cout << "Enter path to cfg laser:" << std::endl;
    
    std::string path;
    std::cin >> path;

    try {
        Laser laser = load_config<Laser>(path); // cfg/laser_params.json

        Pulse pulse;
        pulse.get_spectrum();

        for (double f : pulse.frequency) std::cout << f << "\n";
        for (double i : pulse.intensity) std::cout << i << "\n";
        // std::cout << data.frequency << std::endl;
        // write_spectrum_to_csv(data, "spectrum.csv");
    } catch (const std::exception& e) {
        std::cerr << "Ошибка конфигурации: " << e.what() << '\n';
    }

    return 0;
}