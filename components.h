#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <vector>
#include <cmath>
#include <QVector>

struct Laser{
    double central_wavelength;
    double phase;
    double pulse_duration;
    double energy;
    double avg_count_photons;
    double number_points; //число точек для спектра
    double rep_rate; // частота повторения импульсов
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
    double channel_length;
    double chromatic_dispersion;
    double channel_attenuation;
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
    TimeDomainData spectrum_to_time_domain(const SpectrumData &spectrum, const Laser &laser, const QuantumChannel &quantumChannel);
    int get_photons (const Laser &laser, const QuantumChannel &quantumChannel);
    TimeDomainData gen_composite_pulse(const TimeDomainData &singlePulse, const Laser &laser, int numPulses,
                                          double dt, const QuantumChannel &quantumChannel);

    void gen_ph_timelabel(unsigned int numPulses, const std::vector<unsigned int>& numPhotons,
                                      std::vector<std::vector<double>>& ph_time, const Photodetector &detector,
                                      std::vector<double> &time_slots, double &time);

    void get_time_slot(unsigned int num_pulses, const Laser &laser, std::vector<double>& time_slots, const Photodetector &detector);

    int reg_pulses(std::vector<std::vector<double>>& ph_time,
                                const Photodetector &detector, std::vector<double> &time_slots, std::vector<double> &time_reg);
};

#endif // COMPONENTS_H
