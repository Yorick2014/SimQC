#ifndef COMPONENTS_H
#define COMPONENTS_H

#include <cmath>

struct Laser{
    double centralWavelength;
    double phase;
    double pulseDuration;
    double energy;
    double averageCountPhotons; //среднее число фотонов
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
