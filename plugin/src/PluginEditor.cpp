#include "PluginEditor.h"

namespace { const int barHeight = 64; }

ComboEditor::ComboEditor(ComboProcessor& p) : juce::AudioProcessorEditor(&p), proc(p), params(p) {
  addAndMakeVisible(params);
  addAndMakeVisible(loadButton); addAndMakeVisible(builtInButton);
  addAndMakeVisible(irLabel); addAndMakeVisible(htLabel);
  irLabel.setMinimumHorizontalScale(0.6f);
  htLabel.setJustificationType(juce::Justification::centredRight);

  loadButton.onClick = [this] {
    chooser = std::make_unique<juce::FileChooser>("Choose a cabinet impulse response", juce::File(), "*.wav;*.aif;*.aiff;*.flac");
    const auto flags = juce::FileBrowserComponent::openMode | juce::FileBrowserComponent::canSelectFiles;
    chooser->launchAsync(flags, [this](const juce::FileChooser& fc) {
      const juce::File f = fc.getResult();
      if (f == juce::File()) return;
      if (!proc.loadUserIR(f))
        juce::AlertWindow::showMessageBoxAsync(juce::MessageBoxIconType::WarningIcon, "Impulse response",
                                               "That file could not be read as an impulse response.");
      timerCallback();
    });
  };
  builtInButton.onClick = [this] { proc.clearUserIR(); timerCallback(); };

  setSize(juce::jmax(420, params.getWidth()), juce::jmin(720, params.getHeight()) + barHeight);
  timerCallback();
  startTimerHz(8);
}

ComboEditor::~ComboEditor() { stopTimer(); }

void ComboEditor::paint(juce::Graphics& g) {
  g.fillAll(getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));
}

void ComboEditor::resized() {
  auto area = getLocalBounds();
  auto bar = area.removeFromTop(barHeight).reduced(8, 6);
  auto top = bar.removeFromTop(bar.getHeight() / 2);
  loadButton.setBounds(top.removeFromLeft(100).reduced(0, 2));
  top.removeFromLeft(6);
  builtInButton.setBounds(top.removeFromLeft(130).reduced(0, 2));
  htLabel.setBounds(top);
  irLabel.setBounds(bar);
  params.setBounds(area);
}

void ComboEditor::timerCallback() {
  irLabel.setText("Cabinet IR: " + proc.irDescription(), juce::dontSendNotification);
  htLabel.setText("HT " + juce::String(juce::roundToInt(proc.htVolts.load())) + " V", juce::dontSendNotification);
}
