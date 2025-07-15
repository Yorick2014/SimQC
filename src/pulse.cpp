#include <iostream>
#include <cmath>
#include <complex>
#include <random>
#include <vector>
#include "pulse.hpp"

void Pulse::get_spectrum(const Laser &laser){
    // Центральная частота (Гц)
    double nu0 = SPEED_LIGHT / laser.central_wavelength;

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

void Pulse::spectrum_to_csv(const std::string& path) {
    std::ofstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл: " + path);
    }

    file << std::fixed << std::setprecision(3); // округление до тысячных

    size_t N = frequency.size();
    for (size_t i = 0; i < N; ++i) {
        file << frequency[i] << "," << intensity[i] << "\n";
    }
}

void Pulse::spectrum_to_time_domain(const Laser &laser, const QuantumChannel &channel)
{
    const int N_time = laser.number_points;

    // Диапазон частот спектра
    double& nu_Min = frequency.front();
    double& nu_Max = frequency.back();
    double nu0 = 0.5 * (nu_Min + nu_Max);

    int N = frequency.size(); // или laser.numberPoints
    double dnu = (nu_Min - nu_Max) / (N - 1); // шаг по частоте
    double dt = 1.0 / (N * dnu);

    double t_min = -0.5 * N * dt;

    double lambda0 = laser.central_wavelength;   // м

    // Массив комплексной амплитуды во временной области
    std::vector<std::complex<double>> E_time(N_time, std::complex<double>(0.0, 0.0));

    for (int i = 0; i < N_time; ++i) {
        double t = t_min + i * dt;
        std::complex<double> sum(0.0, 0.0);

        for (int j = 0; j < N; ++j) {
            double nu = frequency[j];
            double I_nu = intensity[j];
            double amp = std::sqrt(I_nu);

            double lambda = SPEED_LIGHT / nu;
            double deltaLambda = lambda - lambda0;
            double deltaLambda_nm = deltaLambda * 1e9;

            // Хроматическая дисперсия как задержка
            double delay = 0.0;
            if (channel.is_crom_disp) {
                // β [ps/(nm·km)] * Δλ [нм] * L [km] = delay [ps]
                delay = channel.chromatic_dispersion * deltaLambda_nm * channel.channel_length; // в пс
                delay *= 1e-12; // перевод в секунды
            }

            double t_shifted = t + delay;
            std::complex<double> phase = std::exp(std::complex<double>(0.0, 2.0 * M_PI * (nu - nu0) * t_shifted));

            sum += amp * phase * dnu;
        }
        E_time[i] = sum;
    }

    time.reserve(N_time);
    t_intensity.reserve(N_time);
    for (int i = 0; i < N_time; ++i) {
        double t = t_min + i * dt;
        double I_t = std::norm(E_time[i]);  // norm(x) = |x|^2
        time.push_back(t);
        t_intensity.push_back(I_t);
    }

    double sum_energy = 0.0;
    for (int i = 0; i < N_time - 1; ++i) {
        double I_mid = 0.5 * (t_intensity[i] + t_intensity[i+1]);
        sum_energy += I_mid * dt;
    }
    double pulse_energy = laser.avg_count_photons * PLANCK_CONSTANT * (SPEED_LIGHT / lambda0);
    if (sum_energy > 0.0) {
        double scale_factor = pulse_energy / sum_energy;
        for (int i = 0; i < N_time; ++i) {
            t_intensity[i] *= scale_factor;
        }
    }
}