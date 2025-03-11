/*
  ==============================================================================

    Oscillators.h

  ==============================================================================
*/

#pragma once
#ifndef Oscillators_h
#define Oscillators_h

#include <cmath>
#include <JuceHeader.h>

class Phasor {
public:
    float process() {
        phase += phaseDelta;

        if (phase > 1.0f)
            phase -= 1.0f;

        return output(phase);
    }

    virtual float output(float p) {
        return p;
    }

    void setSampleRate(float SR)
    {
        sampleRate = SR;
    }

    void setFrequency(float freq) 
    {
        frequency = freq;
        phaseDelta = frequency / sampleRate;
    }

    void setPhase(float ph)
    {
        phase = ph;
    }

private:
    float frequency;
    float sampleRate;
    float phase = 0.0f;
    float phaseDelta;
};

//==================================================

//   CHILD Class
class TriOsc : public Phasor 
{
    float output(float p) override 
    {
        return fabsf(p - 0.5f) - 0.5f;
    }
};

//   CHILD Class
//==================================================
class SinOsc : public Phasor 
{
    float output(float p) override 
    {
        return std::sin(p * 2.0 * 3.14159);
    }
};

//   CHILD Class
//==================================================
class SquareOsc : public Phasor 
{
public:
    float output(float p) override 
    {
        float outVal = 0.5;
        if (p > pulseWidth)
            outVal = -0.5;
        return outVal;
    }
    void setPulseWidth(float pw) 
    {
        pulseWidth = pw;
    }
    float pulseWidth = 0.5f;
private:
};

//   CHILD Class
//==================================================
class SawOsc : public Phasor
{
public:
    /// get saw output
    /// 
    /// @param phase
    /// @return saw output between -1 - 1
    float output(float p) override
    {
        return 2.0f * p - 1.0f;
    }
};
#endif // Oscillators_h