#include "components.h"
#include <algorithm>
#include <QDebug>
#include <cmath>
#include <random>
#include <iomanip>
#include <complex>

const double SPEED_LIGHT = 3e8; // Скорость света, м/с
const double GAUS_K = 0.441; // Cвязь спектра ширины импульса с его продолжительностью во времени
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

double calculate_fwhm_time(const std::vector<std::complex<double>> &E_time, double dt)
{
    int N = E_time.size();
    std::vector<double> intensity(N);

    for (int i = 0; i < N; ++i) {
        intensity[i] = std::norm(E_time[i]);
    }

    double maxI = *std::max_element(intensity.begin(), intensity.end());
    double halfMax = maxI * 0.5;

    int left = -1, right = -1;
    for (int i = 1; i < N; ++i) {
        if (intensity[i - 1] < halfMax && intensity[i] >= halfMax && left == -1) {
            left = i;
        }
        if (intensity[i - 1] >= halfMax && intensity[i] < halfMax) {
            right = i;
        }
    }

    if (left != -1 && right != -1 && right > left) {
        double t_fwhm_sec = (right - left) * dt;
        return t_fwhm_sec * 1e9; // в наносекундах
    }

    return 0.0;
}


TimeDomainData Components::get_time_domain(const SpectrumData &spectrum, const Laser &laser, const QuantumChannel &quantumChannel)
{
    TimeDomainData timeData;
    const int N_time = laser.numberPoints;  // Число временных точек

    // Диапазон частот спектра
    double nu_Min = spectrum.frequency.first();
    double nu_Max = spectrum.frequency.last();
    double nu0 = 0.5 * (nu_Min + nu_Max);

    // Точная ширина по уровню FWHM
    double full_Width = nu_Max - nu_Min;

    double sigma_nu = full_Width / (10.0 * 2.0 * std::sqrt(2.0 * std::log(2.0)));
    double t_FWHM = std::sqrt(std::log(2.0)) / (M_PI * sigma_nu);

    int N = spectrum.frequency.size(); // или laser.numberPoints
    double dnu = (spectrum.frequency.last() - spectrum.frequency.first()) / (N - 1);
    double dt = 1.0 / (N * dnu);

    double t_min = -0.5 * N * dt;
//    double t_max =  0.5 * N * dt;

    double lambda0 = laser.centralWavelength;   // м

    // Массив комплексной амплитуды во временной области
    std::vector<std::complex<double>> E_time(N_time, std::complex<double>(0.0, 0.0));

    for (int i = 0; i < N_time; ++i) {
        double t = t_min + i * dt;
        std::complex<double> sum(0.0, 0.0);

        for (int j = 0; j < N; ++j) {
            double nu = spectrum.frequency[j];
            double I_nu = spectrum.intensity[j];
            double amp = std::sqrt(I_nu);

            double lambda = SPEED_LIGHT / nu;
            double deltaLambda = lambda - lambda0;
            double deltaLambda_nm = deltaLambda * 1e9;
            // double omega = 2.0 * M_PI * nu;

            // Хроматическая дисперсия как задержка
            double delay = 0.0;
            if (quantumChannel.isCromDisp) {
                // β [ps/(nm·km)] * Δλ [нм] * L [km] = delay [ps]
                delay = quantumChannel.chromaticDispersion * deltaLambda_nm * quantumChannel.channelLength; // в пс
                delay *= 1e-12; // перевод в секунды
            }

            double t_shifted = t + delay;
            std::complex<double> phase = std::exp(std::complex<double>(0.0, 2.0 * M_PI * (nu - nu0) * t_shifted));

            sum += amp * phase * dnu;
        }
        E_time[i] = sum;

        if (quantumChannel.isCromDisp) {
            double fwhm_time_ns = calculate_fwhm_time(E_time, dt);
            qDebug() << "Pulse temporal FWHM (after dispersion):" << fwhm_time_ns << "ns";
        }
    }

    // Вычисляем интенсивность как квадрат модуля комплексного поля.
    timeData.time.reserve(N_time);
    timeData.intensity.reserve(N_time);
    for (int i = 0; i < N_time; ++i) {
        double t = t_min + i * dt;
        double I_t = std::norm(E_time[i]);  // norm(x) = |x|^2
        timeData.time.push_back(t);
        timeData.intensity.push_back(I_t);
    }

    double sum_energy = 0.0;
    for (int i = 0; i < N_time - 1; ++i) {
        double I_mid = 0.5 * (timeData.intensity[i] + timeData.intensity[i+1]);
        sum_energy += I_mid * dt;
    }
    double pulse_energy = laser.averageCountPhotons * PLANCK_CONSTANT * (SPEED_LIGHT / lambda0);
    if (sum_energy > 0.0) {
        double scale_factor = pulse_energy / sum_energy;
        for (int i = 0; i < N_time; ++i) {
            timeData.intensity[i] *= scale_factor;
        }
    }

    return timeData;
}

double generate_random_0_to_1() {
    // Инициализация генератора
    static std::random_device rd;
    static std::mt19937 gen(rd());

    // Распределение для целых чисел от 0 до 9999
    static std::uniform_int_distribution<int> dist(0, 9999);

    // Генерация числа и преобразование в диапазон [0.0, 1.0)
    return dist(gen) / 10000.0;
}

double get_att (double att, double lengthChannel) {
    return 1 - pow(10, -0.1 * att * lengthChannel);
}

bool is_photon_loss (double num1, double num2){

    if (num1 > num2){
        return false;
    }
    else
        return true;
}

// Функция для генерации композиции импульсов с заданным периодом между ними
TimeDomainData Components::generateCompositePulse(const TimeDomainData &singlePulse, const Laser &laser, int numPulses, double dt, const QuantumChannel &quantumChannel)
{
    int N_single = laser.numberPoints;

    // Перевод частоты генерации из МГц в Гц и вычисление периода
    double repRateHz = laser.repRate * 1e6;
    double T_sep = (repRateHz > 0.0) ? (1.0 / repRateHz) : 0.0;
    int shiftSamples = static_cast<int>(std::round(T_sep / dt));

    // Общая длина результирующего сигнала
    int N_composite = (numPulses - 1) * shiftSamples + N_single;

    TimeDomainData composite;
    composite.time.resize(N_composite);
    composite.intensity.resize(N_composite, 0.0);

    // Инициализация генератора случайных чисел
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::poisson_distribution<int> dist(laser.averageCountPhotons);

    int global_count_p = 0;;
    for (int p = 0; p < numPulses; ++p) {
        int n_p = dist(gen); // Число фотонов в одном импульсе
        int count_p = n_p;
//        qDebug() << "В импульсе " << p + 1 << "было фотонов: " << n_p;
        if (n_p > 0 && quantumChannel.isAtt == true){
            for (int i = 0; i < n_p; i++)
            {
                // потеря фотонов (затухание)
                if (is_photon_loss(generate_random_0_to_1(),
                                   get_att(quantumChannel.channelAttenuation, quantumChannel.channelLength)) == true) count_p--;
            }
        }

//        qDebug() << "В импульсе " << p + 1 << "осталось фотонов: " << count_p;
        if (count_p > 0)
            global_count_p++;
        double scale_factor = (laser.averageCountPhotons != 0.0) ? (static_cast<double>(count_p) / laser.averageCountPhotons) : 0.0;

        int offset = p * shiftSamples;
        for (int j = 0; j < N_single; ++j) {
            if (offset + j < N_composite) {
                composite.intensity[offset + j] += singlePulse.intensity[j] * scale_factor;
            }
        }
    }

    qDebug() << "Кол-во отправленных импульсов: " << numPulses;
    qDebug() << "Кол-во дошедших импульсов: " << global_count_p; // Дошедших, но не факт, что принятых
    double pulseRelation { static_cast<double>(global_count_p) / static_cast<double>(numPulses)}; // Отношение импульсов
    qDebug() << "Отношение отправленных импульсов к дошедшим " << pulseRelation;

    // массив временных точек
    for (int i = 0; i < N_composite; ++i) {
        composite.time[i] = i * dt;
    }
    return composite;
}
