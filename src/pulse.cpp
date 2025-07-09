#include <iostream>
#include "pulse.hpp"

void Pulse::get_spectrum(){
    std::cout << "get spectrum" << std::endl;
    // Центральная частота (Гц)
    double nu0 = SPEED_LIGHT / laser.central_wavelength;
    std::cout << "laser.central_wavelength = " << laser.central_wavelength << std::endl;
    std::cout << "get spectrum. nu0 = " << nu0 << std::endl;

    // Стандартное отклонение во времени
    const double sigma_t = laser.pulse_duration / (2 * sqrt(2 * log(2)));
    // Стандартное отклонение в частотной области
    const double sigma_nu = 1 / (2 * M_PI * sigma_t);

    // Ширина спектра (FWHM) в Гц
    const double delta_nu = GAUS_K * (1 / laser.pulse_duration);

    // Ширина спектра (FWHM) в нм
//    const double delta_lamda = (pow(laser.centralWavelength, 2) / SPEED_LIGHT) * delta_nu;

    // Определяем диапазон частот
    double nu_min = nu0 - 5 * delta_nu;
    double nu_max = nu0 + 5 * delta_nu;
    int N = laser.number_points;  // Количество точек
    double step = (nu_max - nu_min) / N;

    frequency.clear();
    intensity.clear();
    // Рассчитываем спектр и заполняем структуру
    for (int i = 0; i < N; i++) {
        double nu = nu_min + i * step;
        double intensity = gaussian_spectrum(nu, nu0, sigma_nu);
        
        this->frequency.push_back(nu);
        this->intensity.push_back(intensity);
    }
};

double gaussian_spectrum(double nu, double nu0, double sigma_nu)
{
    return exp(-pow((nu - nu0), 2) / (2 * pow(sigma_nu, 2)));
}