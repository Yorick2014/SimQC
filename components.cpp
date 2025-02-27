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
    int N = laser.numberPoints;  // Количество точек
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

TimeDomainData Components::get_time_domain(const SpectrumData &spectrum, const Laser &laser) {
    TimeDomainData timeData;

    // Количество точек по времени для численного интегрирования
    const int N_time = laser.numberPoints;

    // Вычисляем центральную частоту (nu0) как середину диапазона спектра
    double nu0 = (spectrum.frequency.first() + spectrum.frequency.last()) / 2.0;

    // Оценим "ширину" спектра. Зная, что в get_spectrum диапазон частот:
    // nu_min = nu0 - 5*delta_nu, nu_max = nu0 + 5*delta_nu,
    // полная ширина: 10*delta_nu.
    // По определению Фурье для Гауссового сигнала обратное преобразование даёт огибающую:
    // E_envelope(t) ~ ∫ exp(-x^2/(2σ_ν^2)) cos(2π x t) dx = √(2πσ_ν^2) exp(-2π²σ_ν²t²)
    // где x = ν - ν0.
    // Если FWHM спектра delta_nu = GAUS_K/T, то можно оценить σ_ν через:
    // FWHM = 2√(2 ln2) σ_ν   ⇒   σ_ν = (nu_max - nu_min) / (2√(2 ln2)*5)
    double fullWidth = spectrum.frequency.last() - spectrum.frequency.first(); // = 10*delta_nu
    double sigma_nu = fullWidth / (10.0 * 2.0 * sqrt(2 * log(2))); // оценка σ_ν

    // Оценим характерную длительность огибающей.
    // Для Гауссового поля обратное преобразование даёт E_envelope(t) ~ exp(-2π²σ_ν²t²).
    // Тогда FWHM во времени: t_FWHM = sqrt(ln2)/(πσ_ν)
    double t_FWHM = sqrt(log(2)) / (M_PI * sigma_nu);

    // Зададим диапазон времени, например, от -5*t_FWHM до 5*t_FWHM
    double t_min = -5 * t_FWHM;
    double t_max =  5 * t_FWHM;
    double dt = (t_max - t_min) / (N_time - 1);

    // Определяем шаг по частоте для численного интегрирования
    double dnu = (spectrum.frequency.last() - spectrum.frequency.first()) / (spectrum.frequency.size() - 1);

    // Для каждого t вычисляем интеграл:
    // E(t) = ∫ S(ν) exp(2π i ν t)dν = exp(2π i ν0 t) ∫ S(ν) exp(2π i (ν-ν0)t)dν
    // Принимая, что спектральная фаза нулевая, считаем только реальную часть (огибающую):
    // E_envelope(t) = ∫ S(ν) cos[2π (ν - ν0)t] dν
    for (int i = 0; i < N_time; ++i) {
        double t = t_min + i * dt;
        double sum = 0.0;
        for (int j = 0; j < spectrum.frequency.size(); ++j) {
            double nu = spectrum.frequency[j];
            double value = spectrum.intensity[j];
            double x = nu - nu0; // сдвиг для удаления несущей
            sum += value * cos(2 * M_PI * x * t) * dnu;
        }
        timeData.time.append(t);
        timeData.intensity.append(sum);
    }

    return timeData;
}
