/*
  ==============================================================================

    Synth.h
    Created: 27 Mar 2024 12:50:39pm
    Author:  Junting

  ==============================================================================
*/

#pragma once
#include <JuceHeader.h>
#include "Oscillators.h"
#include "Filter.h"
#include "ChooseOscillator.h"
#include "LFO.h"
#include "Modulation.h"

class synthSound : public juce::SynthesiserSound
{
public:
    bool appliesToNote(int midiNoteNumber) override
    {
        return true;
    }
    bool appliesToChannel(int midiChannel) override
    {
        return true;
    }
    /** The class is reference-counted, so this is a handy pointer class for it. */
};

class synthVoice : public juce::SynthesiserVoice
{
public:

    void setParametersFromApvts(juce::AudioProcessorValueTreeState& apvts) 
    {
        //env
        attackParam[0] = apvts.getRawParameterValue("attack");
        decayParam[0] = apvts.getRawParameterValue("decay");
        sustainParam[0] = apvts.getRawParameterValue("sustain");
        releaseParam[0] = apvts.getRawParameterValue("release");

        attackParam[1] = apvts.getRawParameterValue("attack2");
        decayParam[1] = apvts.getRawParameterValue("decay2");
        sustainParam[1] = apvts.getRawParameterValue("sustain2");
        releaseParam[1] = apvts.getRawParameterValue("release2");

        //osc
        oscTypeParam[0] = apvts.getRawParameterValue("osc1");
        oscTypeParam[1] = apvts.getRawParameterValue("osc2");

        //harmonic
        harmonicOn = apvts.getRawParameterValue("harmonicOn");
        harmType = apvts.getRawParameterValue("harmonicCalculate");
        harmonicOn2 = apvts.getRawParameterValue("harmonicOn2");
        harmType2 = apvts.getRawParameterValue("harmonicCalculate2");

        //gain
        gainParam1 = apvts.getRawParameterValue("gain1");
        gainParam2 = apvts.getRawParameterValue("gain2");

        //filter
        filterOnParam = apvts.getRawParameterValue("filterOn");
        filterTypeParam = apvts.getRawParameterValue("filterType");
        cutOffParam = apvts.getRawParameterValue("cutOff");
        QParam = apvts.getRawParameterValue("Q");

        //lfo
         lfoApplyParam[0] = apvts.getRawParameterValue("lfo1Apply");
         lfoParam[0] = apvts.getRawParameterValue("lfo1");
         lfoFreqParam[0] = apvts.getRawParameterValue("lfo1Freq");
         lfoAmplitudeParam[0] = apvts.getRawParameterValue("lfo1Amplitude");

         lfoApplyParam[1] = apvts.getRawParameterValue("lfo2Apply");
         lfoParam[1] = apvts.getRawParameterValue("lfo2");
         lfoFreqParam[1] = apvts.getRawParameterValue("lfo2Freq");
         lfoAmplitudeParam[1] = apvts.getRawParameterValue("lfo2Amplitude");

         //Modulation
         moduTypeParam[0] = apvts.getRawParameterValue("moduType");
         moduTargetParam[0] = apvts.getRawParameterValue("moduTarget");
         moduFreqParam[0] = apvts.getRawParameterValue("moduFreq");
         moduDepthParam[0] = apvts.getRawParameterValue("moduDepth");

         moduTypeParam[1] = apvts.getRawParameterValue("moduType2");
         moduTargetParam[1] = apvts.getRawParameterValue("moduTarget2");
         moduFreqParam[1] = apvts.getRawParameterValue("moduFreq2");
         moduDepthParam[1] = apvts.getRawParameterValue("moduDepth2");
    }

    void startNote(int midiNoteNumber,
        float velocity,
        juce::SynthesiserSound* sound,
        int currentPitchWheelPosition) override 
    {
        playing = true;
       
        // set basic osc 
        chooseOscillator1.setSampleRate(getSampleRate());
        chooseOscillator2.setSampleRate(getSampleRate());
        // get basic osc freq
        freq = juce::MidiMessage::getMidiNoteInHertz(midiNoteNumber);
        chooseOscillator1.setFrequency(freq);
        chooseOscillator2.setFrequency(freq);

        // set harmonic series
        for (int i = 0; i < 4; i++)
        {
            harmonic[i].setSampleRate(getSampleRate());
            harmonic[i].setFrequency(freq * (i + 2));
        }
        for (int i = 0; i < 4; i++)
        {
            harmonic2[i].setSampleRate(getSampleRate());
            harmonic2[i].setFrequency(freq * (i + 2));
        }

        // set filter
        userfilter.startNote(getSampleRate(), *cutOffParam, *QParam, *filterTypeParam);

        // set lfo
        userlfo1.lfoStartNote(getSampleRate(), *lfoFreqParam[0]);
        userlfo2.lfoStartNote(getSampleRate(), *lfoFreqParam[1]);
         
        //set modulation
        modulation.moduStartNote(getSampleRate(), getSampleRate(), *moduFreqParam[0], freq, *moduDepthParam[0]);
        modulation2.moduStartNote(getSampleRate(), getSampleRate(), *moduFreqParam[0], freq, *moduDepthParam[0]);

        // set env
        env.setSampleRate(getSampleRate());
        juce::ADSR::Parameters envPara;
        envPara.attack = *attackParam[0];
        envPara.decay = *decayParam[0];
        envPara.sustain = *sustainParam[0];
        envPara.release = *releaseParam[0];
        env.setParameters(envPara);
        env.noteOn();

        env2.setSampleRate(getSampleRate());
        juce::ADSR::Parameters envPara2;
        envPara2.attack = *attackParam[1];
        envPara2.decay = *decayParam[1];
        envPara2.sustain = *sustainParam[1];
        envPara2.release = *releaseParam[1];
        env2.setParameters(envPara2);
        env2.noteOn();
    }

    void stopNote(float velocity, bool allowTailOff) override
    {
        env.noteOff();
        env2.noteOff();
     }

    void renderNextBlock(juce::AudioBuffer<float>& outputBuffer,
        int startSample,
        int numSamples) override
    {
        if (playing)
        {
            // DSP LOOP 
            // for each sample        
            for (int i = startSample; i < startSample+numSamples; i++)
            {
                // reset harmonic out
                harmonicSample = 0;
                harmonicSample2 = 0;

                // use lfo1 to osc1
                if (userlfo1.lfoApplyToOscillator1(*lfoApplyParam[0]))
                {
                    float lfoModu = userlfo1.lfoProcess(*lfoParam[0], *lfoAmplitudeParam[0]); // get lfo out
                    chooseOscillator1.setFrequency(freq + lfoModu); // set new frequency
                }
                outputSample1 = chooseOscillator1.setType(*oscTypeParam[0]); // get output of oscillator

                // harmonic for osc1
                if (*harmonicOn == true)
                {
                    if (*harmType == 0) // use different methods to attenuate harmonic amplitudes 
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            float amplitude = 1 / pow(i + 2, 2); // pow type
                            float sampleValue = harmonic[i].setType(*oscTypeParam[0]);
                            sampleValue = amplitude * sampleValue;
                            harmonicSample += sampleValue;
                        }
                    }
                    else if (*harmType == 1)
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            float amplitude = exp(-0.5 * (i + 2)); // exp type
                            float sampleValue = harmonic[i].setType(*oscTypeParam[0]);
                            sampleValue = amplitude * sampleValue;
                            harmonicSample += sampleValue;
                        }
                    }
                    else if (*harmType == 2)
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            float amplitude = 1 / log(i + 2); // log type
                            float sampleValue = harmonic[i].setType(*oscTypeParam[0]);
                            sampleValue = amplitude * sampleValue;
                            harmonicSample += sampleValue;
                        }
                    }
                    outputSample1 = chooseOscillator1.setType(*oscTypeParam[0]) + harmonicSample*0.25;
                }

                // use lfo1 to osc2
                if (userlfo1.lfoApplyToOscillator2(*lfoApplyParam[0]))
                {
                    float lfoModu = userlfo1.lfoProcess( *lfoParam[0], *lfoAmplitudeParam[0]);
                    chooseOscillator2.setFrequency(freq + lfoModu);
                }
                outputSample2 = chooseOscillator2.setType(*oscTypeParam[1]);

                // harmonic for osc2
                if (*harmonicOn2 == true)
                {
                    if (*harmType2 == 0)
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            float amplitude = 1 / pow(i + 2, 2);
                            float sampleValue = harmonic2[i].setType(*oscTypeParam[1]);
                            sampleValue = amplitude * sampleValue;
                            harmonicSample2 += sampleValue;
                        }
                    }
                    else if (*harmType2 == 1)
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            float amplitude = exp(-0.5 * (i + 2));
                            float sampleValue = harmonic2[i].setType(*oscTypeParam[1]);
                            sampleValue = amplitude * sampleValue;
                            harmonicSample2 += sampleValue;
                        }
                    }
                    else if (*harmType2 == 2)
                    {
                        for (int i = 0; i < 4; i++)
                        {
                            float amplitude = 1 / log(i + 2);
                            float sampleValue = harmonic2[i].setType(*oscTypeParam[1]);
                            sampleValue = amplitude * sampleValue;
                            harmonicSample2 += sampleValue;
                        }
                    }
                    outputSample2 = chooseOscillator2.setType(*oscTypeParam[1]) + harmonicSample2*0.25;
                }

                // use lfo2 to osc1
                if (userlfo2.lfoApplyToOscillator1(*lfoApplyParam[1]))
                {
                    float lfoModu = userlfo2.lfoProcess(*lfoParam[1], *lfoAmplitudeParam[1]);
                    chooseOscillator1.setFrequency(freq + lfoModu);
                    outputSample1 = chooseOscillator1.setType(*oscTypeParam[0]);
                }

                // use lfo2 to osc2
                if (userlfo2.lfoApplyToOscillator2(*lfoApplyParam[1]))
                {
                    float lfoModu = userlfo2.lfoProcess(*lfoParam[1], *lfoAmplitudeParam[1]);
                    chooseOscillator2.setFrequency(freq + lfoModu);
                    outputSample2 = chooseOscillator2.setType(*oscTypeParam[1]);
                }
          
                // FM osc1
                if (*moduTargetParam[0] == 1)
                {
                    if (*moduTypeParam[0] == 0)
                        chooseOscillator1.setFrequency(modulation.fmProcess());
                    outputSample1 = chooseOscillator1.setType(*oscTypeParam[0]);
                }

                // FM osc2
                if (*moduTargetParam[0] == 2)
                {
                    if (*moduTypeParam[0] == 0)
                        chooseOscillator2.setFrequency(modulation.fmProcess());
                    outputSample2 = chooseOscillator2.setType(*oscTypeParam[1]);
                }

                // PM osc1
                if (*moduTargetParam[0] == 1)
                {
                    if (*moduTypeParam[0] == 1)
                        chooseOscillator1.setFrequency((modulation.pmProcess())/20.0f);
                    outputSample1 = chooseOscillator1.setType(*oscTypeParam[0]);
                }

                // PM osc2
                if (*moduTargetParam[0] == 2)
                {
                    if (*moduTypeParam[0] == 1)
                        chooseOscillator2.setFrequency((modulation.pmProcess())/20.0f);
                    outputSample2 = chooseOscillator2.setType(*oscTypeParam[1]);
                }

                // FM2 osc1
                if (*moduTargetParam[1] == 1)
                {
                    if (*moduTypeParam[1] == 0)
                        chooseOscillator1.setFrequency(modulation2.fmProcess());
                    outputSample1 = chooseOscillator1.setType(*oscTypeParam[0]);
                }

                // FM2 osc2
                if (*moduTargetParam[1] == 2)
                {
                    if (*moduTypeParam[1] == 0)
                        chooseOscillator2.setFrequency(modulation2.fmProcess());
                    outputSample2 = chooseOscillator2.setType(*oscTypeParam[1]);
                }

                // PM2 osc1
                if (*moduTargetParam[1] == 1)
                {
                    if (*moduTypeParam[1] == 1)
                        chooseOscillator1.setFrequency((modulation2.pmProcess()) / 20.0f); // make the modulator frequency vary between 0-5
                    outputSample1 = chooseOscillator1.setType(*oscTypeParam[0]);
                }

                // PM2 osc2
                if (*moduTargetParam[1] == 2)
                {
                    if (*moduTypeParam[1] == 1)
                        chooseOscillator2.setFrequency((modulation2.pmProcess()) / 20.0f);
                    outputSample2 = chooseOscillator2.setType(*oscTypeParam[1]);
                }

                // use env
                float envvalue = env.getNextSample();
                float envvalue2 = env2.getNextSample();

                // add two basic osc
                outputSample = outputSample1 * envvalue * (*gainParam1) + outputSample2 * envvalue2 * (*gainParam2);

                // use filter
                if (*filterOnParam == true)
                {
                    if (userlfo1.lfoApplyToFilter(*lfoApplyParam[0]))
                    {
                       float lfoModu = userlfo1.lfoProcess(*lfoParam[0], *lfoAmplitudeParam[0]);
                        userfilter.setLfoCutoff(lfoModu*10.0f);
                    }

                    if (userlfo2.lfoApplyToFilter(*lfoApplyParam[1]))
                    {
                        float lfoModu = userlfo2.lfoProcess(*lfoParam[1], *lfoAmplitudeParam[1]);
                        userfilter.setLfoCutoff(lfoModu*10.0f); 
                    }
                    outputSample = userfilter.process(outputSample, *filterTypeParam);
                }

                //for each channel
                for (int chan = 0; chan < outputBuffer.getNumChannels(); chan++)
                {
                    outputBuffer.addSample(chan, i, outputSample * 0.5);
                }

                if ( ! env.isActive() && ! env2.isActive() )
                {
                    playing = false;
                    clearCurrentNote();
                }
            }
        }
    }

    bool canPlaySound(juce::SynthesiserSound*)override
    {
        return true;
     }

    void pitchWheelMoved(int newPitchWheelValue) override
    {

     }

    void controllerMoved(int controllerNumber, int newControllerValue) override
    {

     }

private:
    bool playing = true;

    float freq; // note frequency

    // basic osc output
    float outputSample;
    float outputSample1;
    float outputSample2;
    float harmonicSample;
    float harmonicSample2;
    
    // basic oscillator
    ChooseOscillator chooseOscillator1;
    ChooseOscillator chooseOscillator2;

    // harmonic series
    ChooseOscillator harmonic[4];
    ChooseOscillator harmonic2[4];

    // filter
    Filter userfilter;

    // LFO
    LFO userlfo1;
    LFO userlfo2;

    // modulation (FM, PM)
    Modulation modulation;
    Modulation modulation2;

    //env
    juce::ADSR env;
    juce::ADSR env2;

    //env parameters
    std::atomic<float>* attackParam[2];
    std::atomic<float>* decayParam[2];
    std::atomic<float>* sustainParam[2];
    std::atomic<float>* releaseParam[2];

    //osc param
    std::atomic<float>* oscTypeParam[2];

    //harmonic param
    std::atomic<float>* harmonicOn;
    std::atomic<float>* harmType;
    std::atomic<float>* harmonicOn2;
    std::atomic<float>* harmType2;

    //gain
    std::atomic<float>* gainParam1;
    std::atomic<float>* gainParam2;

    //filter param
    std::atomic<float>* filterOnParam;
    std::atomic<float>* filterTypeParam;
    std::atomic<float>* cutOffParam;
    std::atomic<float>* QParam;

    //lfo param
    std::atomic<float>* lfoApplyParam[2];
    std::atomic<float>* lfoParam[2];
    std::atomic<float>* lfoFreqParam[2];
    std::atomic<float>* lfoAmplitudeParam[2];

    //modulation param
    std::atomic<float>* moduTargetParam[2];
    std::atomic<float>* moduTypeParam[2];
    std::atomic<float>* moduFreqParam[2];
    std::atomic<float>* moduDepthParam[2];
};
