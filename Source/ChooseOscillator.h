/*
  ==============================================================================

    ChooseOscillator.h
    Created: 27 Apr 2024 12:53:33pm
    Author:  Junting

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Oscillators.h"

/// <summary>
/// use to choose different wave shape
/// </summary>
class ChooseOscillator
{
public:
    void setSampleRate(float _sampleRate)
    {
        sampleRate = _sampleRate;
    }

    void setFrequency(float _frequency)
    {
        frequency = _frequency;
    }

    /// get different wave shape and return its output
    /// 
    /// @param  oscillator type (0 - nothing, 1 - sine wave, 2 - triangle wave, 3 - square wave, 4 - saw wave)
    float setType(int _oscType)
    {
        oscType = _oscType;
        float output;
        if (oscType == 1)
        {
            sinOsc.setSampleRate(sampleRate);
            sinOsc.setFrequency(frequency);
            output = sinOsc.process();
            return output;
        }
        else if (oscType == 2)
        {
            triOsc.setSampleRate(sampleRate);
            triOsc.setFrequency(frequency);
            output = triOsc.process();
            return output;
        }
        else if (oscType == 3)
        {
            squareOsc.setSampleRate(sampleRate);
            squareOsc.setFrequency(frequency);
            output = squareOsc.process();
            return output;
        }
        else if (oscType == 4)
        {
            sawOsc.setSampleRate(sampleRate);
            sawOsc.setFrequency(frequency);
            output = sawOsc.process();
            return output;
        }
    }

private:
    float sampleRate = 0.0f;
    float frequency;
    int oscType;

    // get the wave instance
    SinOsc sinOsc;
    TriOsc triOsc;
    SquareOsc squareOsc;
    SawOsc sawOsc;
};