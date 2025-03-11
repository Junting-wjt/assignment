/*
  ==============================================================================

    Filter.h
    Created: 26 Apr 2024 6:55:18pm
    Author:  Junting

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h> 

/// <summary>
/// Filter class used to set filter type and parameters
/// </summary>
class Filter
{
public:

    /// get the sample processed by filter
    /// 
    /// @param input sample and filter type
    /// @return filter processed sample
    float process(float _inSample, int _filterType)
    {
        float freq = cutoffBase + lfoCutoff;
        setFrequency(freq);
        setType(_filterType);
        resetLfoModu();
        return filter.processSingleSampleRaw(_inSample);
    }

    void setSampleRate(float _sampleRate)
    {
        sampleRate = _sampleRate;
    }

    /// set filter frequency
    /// 
    /// @param frequency
    void setFrequency(float _frequency)
    {
        cutoff = _frequency;
    }

    /// set filter resonance
    /// 
    /// @param resonance
    void setResonance(float _resonance)
    {
        Q = _resonance;
    }

    /// set filter type
    /// 
    /// @param filter type (0 - low pass, 1 - high pass, 2 - band pass)
    void setType(int _filterType)
    {
        switch (_filterType)
        {
        case 0:
            filter.setCoefficients(juce::IIRCoefficients::makeLowPass(sampleRate, cutoff, Q));
            break;

        case 1:
            filter.setCoefficients(juce::IIRCoefficients::makeHighPass(sampleRate, cutoff, Q));
            break;

        case 2:
            filter.setCoefficients(juce::IIRCoefficients::makeBandPass(sampleRate, cutoff, Q));
            break;
        }
    }

    /// <summary>
    /// initialize filter
    /// </summary>
    /// <param name="_sampleRate"></param>
    /// <param name="_frequency"></param>
    /// <param name="_resonance"></param>
    /// <param name="_filterType"></param>
    void startNote(float _sampleRate, float _frequency, float _resonance, int _filterType)
    {
        filter.reset();
        setSampleRate(_sampleRate);
        setFrequency(_frequency);
        setResonance(_resonance);
        setCutoffBase(_frequency);
        setType(_filterType);
    }
    
    void setLfoCutoff(float _lfosample)
    {
        lfoCutoff = _lfosample;
    }

    /// set cutoff frequency without lfo
    /// 
    /// @param cutoff without lfo
    void setCutoffBase(float _cutoff)
    {
        cutoffBase = _cutoff;
    }

private:
    float sampleRate = 0.0f;

    juce::IIRFilter filter;
    float  cutoff; // frequency
    float Q; // resonance
    float lfoCutoff = 0.0f;
    float cutoffBase = 4000.0f;

    void resetLfoModu()
    {
        lfoCutoff = 0.0f;
    }
};
