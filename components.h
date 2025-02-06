#ifndef COMPONENTS_H
#define COMPONENTS_H

struct Laser{
    double centralWavelength;
    double phase;
    double pulseDuration;
    double timeSlot;
    double energy;
    int countPhotons;
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
};
#endif // COMPONENTS_H
