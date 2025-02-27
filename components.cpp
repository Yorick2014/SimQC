#include "components.h"
#include "subsystemparam.h"
#include <algorithm>

const double SPEED_LIGHT = 3e8; // Скорость света, м/с
const double GAUS_K = 0.441;

Components::Components()
{

}

Components::~Components()
{

}

SpectrumData Components::get_spectrum(const Laser &laser) {
    SpectrumData data;

    // Центральная частота (Гц)
    const double nu0 = SPEED_LIGHT / laser.centralWavelength;
    // std::cout << "Центральная частота: " << nu0 / 1e12 << " ТГц\n";

    // Стандартное отклонение во времени
    const double sigma_t = laser.pulseDuration / (2 * sqrt(2 * log(2)));
    // Стандартное отклонение в частотной области
    const double sigma_nu = 1 / (2 * M_PI * sigma_t);

    // Ширина спектра (FWHM) в Гц
    const double delta_nu = GAUS_K * (1 / laser.pulseDuration);

    // Ширина спектра (FWHM) в нм
//    const double delta_lamda = (pow(laser.centralWavelength, 2) / SPEED_LIGHT) * delta_nu;

    // Определяем диапазон частот
    double nu_min = nu0 - 5 * delta_nu;
    double nu_max = nu0 + 5 * delta_nu;
    int N = 30;  // Количество точек
    double step = (nu_max - nu_min) / N;

    // Рассчитываем спектр и заполняем структуру
    for (int i = 0; i < N; i++) {
        double nu = nu_min + i * step;
        double intensity = gaussian_spectrum(nu, nu0, sigma_nu);
        data.frequency.append(nu);
        data.intensity.append(intensity);
    }

    return data;
}

double Components::gaussian_spectrum(double nu, double nu0, double sigma_nu)
{
    return exp(-pow((nu - nu0), 2) / (2 * pow(sigma_nu, 2)));
}
