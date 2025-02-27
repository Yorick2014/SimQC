#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <vector>
#include <cmath>
#include <fftw3.h>
#include <QVector>

struct Laser{
    double centralWavelength;
    double phase;
    double pulseDuration;
    double energy;
    double averageCountPhotons; //среднее число фотонов
    double numberPoints; //число точек для спектра
    double frequencyResolution; //множитель для частоты дискретизации
};

struct SpectrumData {
    QVector<double> frequency;
    QVector<double> intensity;
};

struct TimeDomainData {
    QVector<double> time;
    QVector<double> intensity; // огибающая временного сигнала
};

struct QuantumChannel{
    double channelLength;
    double chromaticDispersion;
    double channelAttenuation;
};

class Components
{
private:
    double gaussian_spectrum(double nu, double nu0, double sigma_nu);
public:
    Components();
    ~Components();

    SpectrumData get_spectrum(const Laser &laser);
    TimeDomainData get_time_domain(const SpectrumData &spectrum, const Laser &laser);
};

#endif // COMPONENTS_H
