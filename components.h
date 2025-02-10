#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <vector>
#include <cmath>
#include <fftw3.h>

struct Laser{
    double centralWavelength;
    double phase;
    double pulseDuration;
    double energy;
    double averageCountPhotons; //среднее число фотонов
    double numberPoints; //число точек для спектра
    double frequencyResolution; //множитель для частоты дискретизации
};

struct QuantumChannel{
    double channelLength;
    double chromaticDispersion;
    double channelAttenuation;
};

class Components
{
public:
    Components();
    ~Components();

    const double SPEED_LIGHT = 3e8;

    void spectrum(const Laser &laser, std::vector<double> &frequency, std::vector<double> &spectrum);
};

#endif // COMPONENTS_H
