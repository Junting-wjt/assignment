/*
  ==============================================================================

    This file contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Synth.h"

//==============================================================================
/**
*/
class PolyphonicSynthAudioProcessor  : public juce::AudioProcessor
{
public:
    //==============================================================================
    PolyphonicSynthAudioProcessor();
    ~PolyphonicSynthAudioProcessor() override;

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const juce::String getProgramName (int index) override;
    void changeProgramName (int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation (juce::MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

private:
    juce::Synthesiser synth;
    int voicecount = 4;

    //reverb
    juce::Reverb reverb;
    std::atomic<float>* reverbOnParam;

    //UI
    juce::AudioProcessorValueTreeState apvts;
    juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout()
    {
        juce::AudioProcessorValueTreeState::ParameterLayout layout;

        //Reverb
        layout.add(std::make_unique<juce::AudioParameterBool>("reverbOn", "Reverb: on", false));

        //oscillator1
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("osc1", 1), "Oscillator1",
            juce::StringArray({ "none","sine","triangle","square","saw" }), 0));

        //harmonic1
        layout.add(std::make_unique<juce::AudioParameterBool>("harmonicOn", "Harmonic: on", false));
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("harmonicCalculate", 1), "Harmonic Type",
            juce::StringArray({ "pow","exp","log"}), 0));

        //Gain1
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("gain1", 1), "Gain 1", 0, 3, 1));

        //env1
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("attack", 1), "Attack 1", 0.1, 5, 1));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("decay", 1), "Decay 1", 0.0, 5, 0.5));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("sustain", 1), "Sustain 1", 0.0, 1, 0.2));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("release", 1), "Release 1", 0.0, 5, 1));
        
        //oscillator2
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("osc2", 1), "Oscillator2",
            juce::StringArray({ "none","sine","triangle","square","saw" }), 0));

        //harmonic2
        layout.add(std::make_unique<juce::AudioParameterBool>("harmonicOn2", "Harmonic: on", false));
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("harmonicCalculate2", 1), "Harmonic Type",
            juce::StringArray({ "pow","exp","log" }), 0));

        //Gain2
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("gain2", 1), "Gain 2", 0, 3, 1));

        //env2
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("attack2", 1), "Attack 2", 0.1, 5, 1));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("decay2", 1), "Decay 2", 0.0, 5, 0.5));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("sustain2", 1), "Sustain 2", 0.0, 1, 0.2));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("release2", 1), "Release 2", 0.0, 5, 1));

        //filter
        layout.add(std::make_unique<juce::AudioParameterBool>("filterOn", "Filter: on", false));
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("filterType", 1), "Filter",
            juce::StringArray({"LowPass","HighPass","BandPass" }), 0));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("cutOff", 1), "Cutoff", 20, 8000, 120));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("Q", 1), "Q", 0.0, 1, 0.1));

        //LFO1
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("lfo1Apply", 1), "LFO1 apply to",
            juce::StringArray({ "none", "oscillator1","oscillator2","cutoff frequency"}), 0));
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("lfo1", 1), "LFO1 Shape",
            juce::StringArray({ "none", "sine","triangle","square","saw" }), 0));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("lfo1Freq", 1), "LFO1 Frequency", 0, 10, 0));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("lfo1Amplitude", 1), "LFO1 Amplitude (%)", 0, 100, 0));
        //LFO2
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("lfo2Apply", 1), "LFO2 apply to",
            juce::StringArray({ "none", "oscillator1","oscillator2","cutoff frequency" }), 0));
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("lfo2", 1), "LFO2 Shape",
            juce::StringArray({ "none", "sine","triangle","square","saw" }), 0));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("lfo2Freq", 1), "LFO2 Frequency", 0, 10, 0));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("lfo2Amplitude", 1), "LFO2 Amplitude (%)", 0, 100, 0));

        //Modulation1
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("moduTarget", 1), "Modulation Target 1",
            juce::StringArray({ "none", "oscillator1","oscillator2"}), 0));
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("moduType", 1), "Modulation Type 1",
            juce::StringArray({ "FM", "PM"}), 0));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("moduFreq", 1), "Modulator Frequency 1 (%)", 1, 100, 50));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("moduDepth", 1), "Modulation Depth 1", 0, 10, 0));
        //Modulation2
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("moduTarget2", 1), "Modulation Target 2",
            juce::StringArray({ "none", "oscillator1","oscillator2" }), 0));
        layout.add(std::make_unique<juce::AudioParameterChoice>(juce::ParameterID("moduType2", 1), "Modulation Type 2",
            juce::StringArray({ "FM", "PM" }), 0));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("moduFreq2", 1), "Modulator Frequency 2 (%)", 1, 100, 50));
        layout.add(std::make_unique<juce::AudioParameterFloat>(juce::ParameterID("moduDepth2", 1), "Modulation Depth 2", 0, 10, 0));

        return layout;
    }
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (PolyphonicSynthAudioProcessor)
};
