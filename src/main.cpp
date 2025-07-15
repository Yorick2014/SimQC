#include <iostream>
#include <string>
#include "simulation_params.hpp"
#include "pulse.hpp"


int main() {
    // std::cout << "Enter path to cfg laser:" << std::endl;
    std::cout << "Start" << std::endl;
    
    // std::string path;
    // std::cin >> path;
    Laser laser;
    QuantumChannel q_channel;

    try {
        laser = load_config<Laser>("cfg/laser_params.json"); // cfg/laser_params.json
        q_channel = load_config<QuantumChannel>("cfg/channel_params.json");
        
    } catch (const std::exception& e) {
        std::cerr << "Ошибка конфигурации: " << e.what() << '\n';
    }
    Pulse pulse;
    pulse.get_spectrum(laser);
    // pulse.spectrum_to_csv("spectrum.csv");
    pulse.spectrum_to_time_domain(laser, q_channel);

    double total_time = 0;
    for (int i = 0; i < pulse.time.size(); i++)
    {
        total_time = total_time + pulse.time[i];
    }
    std::cout << "Total time (ns): " << total_time * 1e9 << std::endl;

    // for (double f : pulse.frequency) std::cout << f << "\n";
    // for (double i : pulse.intensity) std::cout << i << "\n";

    std::cout << "End" << std::endl;
    return 0;
}