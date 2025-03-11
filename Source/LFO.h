/*
  ==============================================================================

    LFO.h
    Created: 28 Apr 2024 10:55:51am
    Author:  Junting

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include"ChooseOscillator.h" // for set LFO type and parameters

/// <summary>
/// LFO class used to set LFO and get the output of it
/// </summary>
class LFO
{
public:

    /// set sampleRate
    /// 
    /// @param sampleRate
    void setSampleRate(float _sampleRate)
    {
        sampleRate = _sampleRate;
        lfo.setSampleRate(sampleRate);
        smoothedLFOValue.reset(_sampleRate, 0.01);
    }

    /// set lfo frequency
    /// 
    /// @param frequency
    void setFrequency(float _frequency)
    {
        frequency = _frequency;
        lfo.setFrequency(frequency);
    }

    /// set lfo type and get output sample without amplitude
    /// 
    /// @param filter type
    /// @return smoothed lfo out
    float setLfoType(int _lfoType)
    {
        float lfoOut = lfo.setType(_lfoType);
        smoothedLFOValue.setTargetValue(lfoOut);
        float smoothedLfoOut = smoothedLFOValue.getNextValue();
        return smoothedLfoOut;
    }

    /// set lfo amplitude
    /// 
    /// @param amplitude
    /// @return smoothed lfo out
    float setAmplitude(float _lfoAmplitude)
    {
        lfoAmplitude = _lfoAmplitude;
        return lfoAmplitude;
    }

    /// initialize lfo
    /// 
    /// @param sampleRate and frequency
    void lfoStartNote(float _sampleRate, float _frequency)
    {
        setSampleRate(_sampleRate);
        setFrequency(_frequency);
        smoothedLFOValue.setCurrentAndTargetValue(0.0f); // reset lfo
    }

    /// get the output of lfo with amplitude
    /// 
    /// @param lfo type and lfo amplitude
    /// @return final lfo out
    float lfoProcess( int _lfoType, float _lfoAmplitude)
    {
       
        float lfoSample = setLfoType(_lfoType);
        float lfoAm = setAmplitude(_lfoAmplitude);
        float lfoOutSample = lfoSample * lfoAm;
        return lfoOutSample;
    }

    /// when choose osc1, return true
    /// 
    /// @param lfo target
    /// @return bool value
    bool lfoApplyToOscillator1(int _lfoTarget)
    {
        if (_lfoTarget == 1)
            return true;
        else
            return false;
    }

    /// when choose osc2, return true
    /// 
    /// @param lfo target
    /// @return bool value
    bool lfoApplyToOscillator2(int _lfoTarget)
    {
        if (_lfoTarget == 2)
            return true;
        else
            return false;
    }

    /// when choose cutoff, return true
   /// 
   /// @param lfo target
   /// @return bool value
    bool lfoApplyToFilter(int _lfoTarget)
    {
        if (_lfoTarget == 3)
            return true;
        else
            return false;
    }

private:
    ChooseOscillator lfo;
    juce::SmoothedValue<float> smoothedLFOValue;

    float sampleRate = 0.0f;
    float frequency = 0.0f;
    float lfoAmplitude = 0.0f;
};