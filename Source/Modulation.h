/*
  ==============================================================================

    Modulation.h
    Created: 1 May 2024 11:27:24am
    Author:  Junting

  ==============================================================================
*/

#pragma once
#include"Oscillators.h"

/// <summary>
/// for more obvious modulation, different from lfo
/// </summary>
class Modulation
{
public:
    void setModulatorSampleRate(float _moduSampleRate)
    {
        moduSampleRate = _moduSampleRate;
        smoothedModuValue.reset(_moduSampleRate, 0.01);
    }

    void setCarrierSampleRate(float _carrSampleRate)
    {
        carrSampleRate = _carrSampleRate;
    }

    void setModulatorFrequency(float _moduFreq)
    {
        moduFreq = _moduFreq;
    }

    void setCarrierFrequency(float _carrFreq)
    {
        carrFreq = _carrFreq;
    }

    void setModuDepth(float _moduDepth)
    {
        moduDepth = _moduDepth;
    }

    void setFreqBase(float _freqBase)
    {
        freqBase = _freqBase;
    }

    /// initialize modulator and carrier wave
    /// 
    /// @param sampleRate and frequency of modulator and carrier, modulation depth 
    void moduStartNote(float _moduSampleRate, float _carrSampleRate, float _moduFreq, float _carrFreq, float _moduDepth)
    {
        setCarrierSampleRate(_carrSampleRate);
        setCarrierFrequency(_carrFreq);
        setModulatorSampleRate(_moduSampleRate);
        setModulatorFrequency(_moduFreq);
        setModuDepth(_moduDepth);
        setFreqBase(_carrFreq);
        moduOsc.setSampleRate(moduSampleRate);
        moduOsc.setFrequency(moduFreq);
        smoothedModuValue.setCurrentAndTargetValue(0.0f);
    }

    /// use modulator on carrier wave and get modulation output
    /// 
    /// @return frequency modulation out
    float fmProcess()
    {
        float moduSample = moduDepth * moduOsc.process();
        smoothedModuValue.setTargetValue(moduSample);
        float smoothedModuSample = smoothedModuValue.getNextValue();
        float carrOutFreq = freqBase + (carrSampleRate * smoothedModuSample / (2 * 3.14159 * moduFreq));
        return carrOutFreq;
    }

    /// use modulator on carrier wave and get modulation output
    /// 
    /// @return phase modulation out
    float pmProcess()
    {
        float moduSample = moduDepth * 0.1f * moduOsc.process();
        smoothedModuValue.setTargetValue(moduSample);
        float smoothedModuSample = smoothedModuValue.getNextValue();
        float carrOutFreq = freqBase + (carrSampleRate * smoothedModuSample / (2 * 3.14159));
        return carrOutFreq;
    }

private:
    float moduFreq;
    float carrFreq;
    float moduSampleRate;
    float carrSampleRate;
    float moduDepth;
    float freqBase;

    SinOsc moduOsc;
    juce::SmoothedValue<float> smoothedModuValue;
};