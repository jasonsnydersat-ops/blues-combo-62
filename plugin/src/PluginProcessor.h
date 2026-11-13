#pragma once

#include <juce_audio_formats/juce_audio_formats.h>
#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "AmpCore.h"

class ComboProcessor : public juce::AudioProcessor, private juce::Timer {
public:
  ComboProcessor();
  ~ComboProcessor() override;

  void prepareToPlay(double sampleRate, int samplesPerBlock) override;
  void releaseResources() override {}
  bool isBusesLayoutSupported(const BusesLayout& layouts) const override;
  void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
  using juce::AudioProcessor::processBlock;

  juce::AudioProcessorEditor* createEditor() override;
  bool hasEditor() const override { return true; }

  const juce::String getName() const override { return JucePlugin_Name; }
  bool acceptsMidi() const override { return false; }
  bool producesMidi() const override { return false; }
  bool isMidiEffect() const override { return false; }
  double getTailLengthSeconds() const override { return 1.0; }

  int getNumPrograms() override { return 1; }
  int getCurrentProgram() override { return 0; }
  void setCurrentProgram(int) override {}
  const juce::String getProgramName(int) override { return {}; }
  void changeProgramName(int, const juce::String&) override {}

  void getStateInformation(juce::MemoryBlock& destData) override;
  void setStateInformation(const void* data, int sizeInBytes) override;

  juce::AudioProcessorValueTreeState apvts;

  /* Cabinet impulse response chosen by the user. Call from the message thread. */
  bool loadUserIR(const juce::File& file);
  void clearUserIR();
  juce::String irDescription() const;

  std::atomic<float> htVolts { 445.0f };

private:
  static juce::AudioProcessorValueTreeState::ParameterLayout createLayout();
  void timerCallback() override;
  void processChunk(juce::AudioBuffer<float>& buffer, int start, int len);
  void applyCabIR();          // built-in or user IR -> convolver (needs irLock)
  bool readIRFile(const juce::File& file);

  std::atomic<float>* pKnob[bb::AmpCore::NUM_KNOBS] = {};
  std::atomic<float> *pInput = nullptr, *pSens = nullptr, *pTrem = nullptr, *pStandby = nullptr, *pTrim = nullptr;
  std::atomic<float> *pCab = nullptr, *pMic = nullptr, *pPos = nullptr, *pDist = nullptr, *pRoom = nullptr, *pLevel = nullptr, *pOS = nullptr;

  juce::CriticalSection coreLock;   // audio thread try-locks; held briefly to swap the engine
  std::unique_ptr<bb::AmpCore> core;
  int coreOS = 0;

  juce::CriticalSection irLock;     // prepareToPlay vs. timer vs. editor (never the audio thread)
  juce::dsp::Convolution cab, room;
  double rate = 0;
  float lastPos = -1, lastDist = -1; int lastMic = -1;
  juce::String irPath, irName, irInfo;
  std::vector<std::vector<float>> userRaw;   // decoded file, at userRawRate
  double userRawRate = 0;

  juce::AudioBuffer<float> work, wet;
  juce::SmoothedValue<float> roomGain, outGain, cabMix;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ComboProcessor)
};
