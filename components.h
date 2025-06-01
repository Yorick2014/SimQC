#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <vector>
#include <cmath>
#include <QVector>

struct Laser{
    double centralWavelength;
    double phase;
    double pulseDuration;
    double energy;
    double averageCountPhotons; //среднее число фотонов
    double numberPoints; //число точек для спектра
    double repRate; // частота повторения импульсов
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
    bool isAtt;
    bool isCromDisp; //хроматическая дисперсия
};

struct Photodetector{
    double quantum_efficiency;
    double dead_time;
    double time_slot;
};

class Components
{
private:
    double gaussian_spectrum(double nu, double nu0, double sigma_nu);
public:
    Components();
    ~Components();

    SpectrumData get_spectrum(const Laser &laser);
    TimeDomainData get_time_domain(const SpectrumData &spectrum, const Laser &laser, const QuantumChannel &quantumChannel);
    int get_photons (const Laser &laser, const QuantumChannel &quantumChannel);
    TimeDomainData generateCompositePulse(const TimeDomainData &singlePulse, const Laser &laser, int numPulses,
                                          double dt, const QuantumChannel &quantumChannel);
    void gen_ph_timelabel(unsigned int numPulses, const std::vector<unsigned int>& numPhotons,
                                      std::vector<std::vector<double>>& ph_time, const Photodetector &detector);
};

#endif // COMPONENTS_H
