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

//        if (quantumChannel.isCromDisp) {
//            double fwhm_time_ns = calculate_fwhm_time(E_time, dt);
//            qDebug() << "Pulse temporal FWHM (after dispersion):" << fwhm_time_ns << "ns";
//        }
    }

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
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<int> dist(0, 9999);

    // Диапазон [0.0, 1.0)
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

int Components::get_photons (const Laser &laser, const QuantumChannel &quantumChannel){

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::poisson_distribution<int> dist(laser.averageCountPhotons);

    int n_p = dist(gen); // Число фотонов в одном импульсе
    int count_p = n_p;
    if (n_p > 0 && quantumChannel.isAtt == true){
        for (int i = 0; i < n_p; i++)
        {
            // потеря фотонов (затухание)
            if (is_photon_loss(generate_random_0_to_1(),
                               get_att(quantumChannel.channelAttenuation, quantumChannel.channelLength)) == true) count_p--;
        }
    }

    return count_p;
}

// Функция для генерации композиции импульсов с заданным периодом между ними
TimeDomainData Components::generateCompositePulse(const TimeDomainData &singlePulse, const Laser &laser, int num_pulses, double dt, const QuantumChannel &quantumChannel)
{
    int N_single = laser.numberPoints;

    // Перевод частоты генерации из МГц в Гц и вычисление периода
    double repRateHz = laser.repRate * 1e6;
    double T_sep = (repRateHz > 0.0) ? (1.0 / repRateHz) : 0.0;
    int shiftSamples = static_cast<int>(std::round(T_sep / dt));

    // Общая длина результирующего сигнала
    int N_composite = (num_pulses - 1) * shiftSamples + N_single;

    TimeDomainData composite;
    composite.time.resize(N_composite);
    composite.intensity.resize(N_composite, 0.0);

    // Инициализация генератора случайных чисел
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::poisson_distribution<int> dist(laser.averageCountPhotons);

    int res_num_pulse = 0;
    for (int p = 0; p < num_pulses; ++p) {

        int res_num_ph = get_photons(laser, quantumChannel); // получение числа фотонов
        if (res_num_ph > 0)
            res_num_pulse++;
        double scale_factor = (laser.averageCountPhotons != 0.0) ? (static_cast<double>(res_num_ph) / laser.averageCountPhotons) : 0.0;

        int offset = p * shiftSamples;
        for (int j = 0; j < N_single; ++j) {
            if (offset + j < N_composite) {
                composite.intensity[offset + j] += singlePulse.intensity[j] * scale_factor;
            }
        }
    }

    qDebug() << "Кол-во отправленных импульсов: " << num_pulses;
    qDebug() << "Кол-во дошедших импульсов: " << res_num_pulse; // Дошедших, но не факт, что принятых
//    double pulseRelation { static_cast<double>(num_photons) / static_cast<double>(numPulses)}; // Отношение импульсов
//    qDebug() << "Отношение отправленных импульсов к дошедшим " << pulseRelation;

    // массив временных точек
    for (int i = 0; i < N_composite; ++i) {
        composite.time[i] = i * dt;
    }
    return composite;
}

void Components::gen_ph_timelabel(unsigned int numPulses, const std::vector<unsigned int>& numPhotons,
                                  std::vector<std::vector<double>>& ph_time, const Photodetector &detector,
                                  std::vector<double> &time_slots, double &time){
    ph_time.resize(numPulses);
    double threshold = detector.quantum_efficiency / 100.0; // перевод QE из процентов

    for (unsigned int i = 0; i < numPulses; ++i) {
        ph_time[i].resize(numPhotons[i]);

        for (unsigned int j = 0; j < ph_time[i].size(); ++j) {
            double rand_num = generate_random_0_to_1();

            // проверка QE и присваивание временной метки
            if (rand_num <= threshold && rand_num > 0){
                ph_time[i][j] = time_slots[i] + rand_num * time;
            }
//            qDebug() << "Pulse" << i + 1 << "Photon" << j + 1 << "time label:" << ph_time[i][j];
//            qDebug() << time;
        }

        std::sort(ph_time[i].begin(), ph_time[i].end());
    }

//    qDebug() << "Out ph_time:" << ph_time;
}

void Components::get_time_slot(unsigned int num_pulses, const Laser &laser, std::vector<double>& time_slots, const Photodetector &detector){

    double step = 1.0 / (laser.repRate * 1e6);

    time_slots.clear();

    for (unsigned int i = 0; i < num_pulses; ++i) {
        time_slots.push_back(i * step);
    }
}

int Components::reg_pulses(std::vector<std::vector<double>>& ph_time,
                            const Photodetector &detector, std::vector<double> &time_slots, std::vector<double> &time_reg){

    int num_reg_pulses = 0;
    double detector_dead = detector.dead_time;
    double slot_duration = detector.time_slot;
    double last_registered_time = -1e9;

    for (unsigned int i = 0; i <  ph_time.size(); ++i) {
        double left = time_slots[i] - slot_duration / 2;
        double right = time_slots[i] + slot_duration / 2;
        for (unsigned int j = 0; j < ph_time[i].size(); ++j) {

            double t = ph_time[i][j];
//            qDebug() << "Импульс:" << ph_time[i][j];

            if (t < left)
                continue;;

            if (t > right)
                break;

            if (t >= last_registered_time + detector_dead){
                time_reg.push_back(t);
                last_registered_time = t;
                num_reg_pulses++;
//                qDebug() << "Регистрация импульса:" << t;
                break;
            }
        }
    }
    return num_reg_pulses;
}
//            if (ph_time[i][j] >= time_slots[i] - (detector.time_slot / 2) &&
//                    ph_time[i][j] <= time_slots[i] + (detector.time_slot / 2) &&
//                    ph_time[i][j] > time_slots[i] + detector.dead_time){
//                time_reg.push_back(ph_time[i][j]);
//                num_reg_pulses++;

//                qDebug() << "Регистрация импульса:" << time_reg[i];
//                qDebug() << "dead time:" << detector.dead_time;
//                qDebug() << "Регистрация импульса + dead time:" << time_reg[i] + detector.dead_time;
//                break;
//            }
