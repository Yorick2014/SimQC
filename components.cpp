#include "components.h"
#include "subsystemparam.h"
#include <algorithm>

const double SPEED_LIGHT = 3e8; // Скорость света, м/с
const double GAUS_K = 0.441;
static const double PLANCK_CONSTANT = 6.62607015e-34; // Дж·с

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

TimeDomainData Components::get_time_domain(const SpectrumData &spectrum, const Laser &laser)
{
    TimeDomainData timeData;

    // Количество точек по времени для численного интегрирования
    const int N_time = laser.numberPoints;

    // Определяем центральную частоту (nu0) как середину диапазона
    double nu_Min = spectrum.frequency.first();
    double nu_Max = spectrum.frequency.last();
    double nu0   = 0.5 * (nu_Min + nu_Max);

    // Оценка ширины спектра (fullWidth = 10 * delta_nu при вашем get_spectrum)
    double full_Width = nu_Max - nu_Min;
    // Оцениваем сигма_ν (sigma_nu) для гауссового спектра
    // FWHM = 2 * sqrt(2 ln2) * sigma_nu => sigma_nu = fullWidth / (10 * 2 sqrt(2 ln2))
    double sigma_nu = full_Width / (10.0 * 2.0 * std::sqrt(2.0 * std::log(2.0)));

    // Оценка FWHM во времени (для гауссового импульса при нулевой фазе)
    // t_FWHM ~ sqrt(ln2) / (π * sigma_nu)
    double t_FWHM = std::sqrt(std::log(2.0)) / (M_PI * sigma_nu);

    // Задаём диапазон времени
    double t_min = -5.0 * t_FWHM;
    double t_max =  5.0 * t_FWHM;
    double dt    = (t_max - t_min) / (N_time - 1);

    // Шаг по частоте для численного интегрирования
    int N_freq = spectrum.frequency.size();
    double dnu = (nu_Max - nu_Min) / (N_freq - 1);

    std::vector<double> E_time(N_time, 0.0); // временная амплитуда (до возведения в квадрат)

     for (int i = 0; i < N_time; ++i) {
         double t = t_min + i * dt;
         double sum = 0.0;
         for (int j = 0; j < N_freq; ++j) {
             double nu    = spectrum.frequency[j];
             double I_nu  = spectrum.intensity[j];  // интенсивность в частотной области
             double ampNu = std::sqrt(I_nu);        // амплитуда в частотной области
             double x     = nu - nu0;               // сдвиг к центральной частоте

             sum += ampNu * std::cos(2.0 * M_PI * x * t) * dnu;
         }
         E_time[i] = sum; // амплитуда E(t)
     }

     // 3) Превращаем амплитуду в интенсивность:
     //    I(t) = E(t)^2
     //    (пока без нормировки на число фотонов)
     timeData.time.reserve(N_time);
     timeData.intensity.reserve(N_time);
     for (int i = 0; i < N_time; ++i) {
         double t = t_min + i * dt;
         double I_t = E_time[i]*E_time[i]; // интенсивность во времени
         timeData.time.push_back(t);
         timeData.intensity.push_back(I_t);
     }

     // --- 4) Нормировка на N фотонов ---
     // Число фотонов = laser.averageCountPhotons
     // Энергия импульса должна быть E_pulse = N * h * nu0
     // Посчитаем интеграл I(t) dt (численная интеграция методом трапеций или прямоугольников).
     // Это даст суммарную энергию (в Дж), если I(t) в Вт.

     double sum_energy = 0.0;
     for (int i = 0; i < N_time - 1; ++i) {
         double I1 = timeData.intensity[i];
         double I2 = timeData.intensity[i+1];
         // метод трапеций
         double I_mid = 0.5 * (I1 + I2);
         sum_energy += I_mid * dt;
     }

     double pulse_energy = laser.averageCountPhotons * PLANCK_CONSTANT * nu0; // (Дж)

     // Чтобы итоговая энергия стала равна pulse_energy, масштабируем интенсивность
     // scaleFactor = desiredEnergy / sumEnergy
     // Тогда I_scaled(t) = scaleFactor * I(t)
     // и интеграл по времени даст desiredEnergy.
     if (sum_energy > 0.0) {
         double scale_factor = pulse_energy / sum_energy;
         for (int i = 0; i < N_time; ++i) {
             timeData.intensity[i] *= scale_factor;
         }
     }

    return timeData;
}

