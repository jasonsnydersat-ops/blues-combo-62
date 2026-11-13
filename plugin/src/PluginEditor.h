#pragma once

#include <juce_audio_utils/juce_audio_utils.h>

#include "PluginProcessor.h"

/* The host-style parameter list, plus a bar for choosing a cabinet impulse response. */
class ComboEditor : public juce::AudioProcessorEditor, private juce::Timer {
public:
  explicit ComboEditor(ComboProcessor&);
  ~ComboEditor() override;
  void paint(juce::Graphics&) override;
  void resized() override;

private:
  void timerCallback() override;

  ComboProcessor& proc;
  juce::GenericAudioProcessorEditor params;
  juce::TextButton loadButton { "Load IR..." }, builtInButton { "Use built-in cab" };
  juce::Label irLabel, htLabel;
  std::unique_ptr<juce::FileChooser> chooser;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(ComboEditor)
};
