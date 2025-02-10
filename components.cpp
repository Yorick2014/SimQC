#include "components.h"
#include "subsystemparam.h"
#include <algorithm>

Components::Components()
{

}

Components::~Components()
{

}

void Components::spectrum(const Laser &laser, std::vector<double> &frequency, std::vector<double> &spectrum) {

    double omega = (SPEED_LIGHT * 2 * M_PI / laser.centralWavelength); // Угловая частота
    double F = SPEED_LIGHT / laser.centralWavelength;
    double T = (1/(laser.frequencyResolution * F));
    double dt = T / laser.numberPoints;

    // Создаём сигнал
    std::vector<double> signal(laser.numberPoints);
    for (int i = 0; i < laser.numberPoints; ++i) {
        double t = i * T - (laser.numberPoints / 2) * T;
        signal[i] = 10.0 * cos(2 * M_PI * omega * t) * exp(-pow(t, 2) / (2 * pow(laser.pulseDuration, 2)));
    }

    // FFTW: входной и выходной массивы
    fftw_complex *in, *out;
    fftw_plan planFFT;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * laser.numberPoints);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * laser.numberPoints);

    for (int i = 0; i < laser.numberPoints; ++i) {
        in[i][0] = signal[i]; // Действительная часть
        in[i][1] = 0.0;       // Мнимая часть (нулевая)
    }

    // Выполняем FFT
    planFFT = fftw_plan_dft_1d(laser.numberPoints, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(planFFT);

    // Заполняем выходные массивы частоты и спектра
    frequency.clear();
    spectrum.clear();
    for (int i = 0; i < laser.numberPoints / 2; ++i) { // Только положительные частоты
        double freq = i / (laser.numberPoints * dt);
        double magnitude = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);

        frequency.push_back(freq);
        spectrum.push_back(magnitude);
    }

    // Освобождаем память
    fftw_destroy_plan(planFFT);
    fftw_free(in);
    fftw_free(out);
}
