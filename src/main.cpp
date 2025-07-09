#include <iostream>
#include <string>
#include "simulation_params.hpp"
#include "pulse.hpp"


int main() {
    std::cout << "Enter path to cfg laser:" << std::endl;
    
    // std::string path;
    // std::cin >> path;
    Laser laser;

    try {
        laser = load_config<Laser>("cfg/laser_params.json"); // cfg/laser_params.json
        
    } catch (const std::exception& e) {
        std::cerr << "Ошибка конфигурации: " << e.what() << '\n';
    }
    Pulse pulse;
    pulse.laser = laser;
    pulse.get_spectrum();
    

    for (double f : pulse.frequency) std::cout << f << "\n";
    for (double i : pulse.intensity) std::cout << i << "\n";

    std::cout << "End" << std::endl;
    return 0;
}