#include "PluginProcessor.h"

#include "CabIR.h"
#include "IRPrep.h"
#include "PluginEditor.h"

using Knob = bb::AmpCore::Knob;

namespace {
const char* const knobIds[bb::AmpCore::NUM_KNOBS] = { "vol1", "vol2", "treble", "mid", "bass", "presence", "speed", "depth" };
const char* const knobNames[bb::AmpCore::NUM_KNOBS] = { "Loudness I", "Loudness II", "Treble", "Middle", "Bass", "Presence", "Tremolo Speed", "Tremolo Intensity" };
const float knobDefaults[bb::AmpCore::NUM_KNOBS] = { 5.5f, 3.0f, 6.0f, 6.0f, 5.0f, 5.0f, 5.0f, 0.0f };

juce::AudioBuffer<float> toBuffer(const std::vector<float>* ch, int numCh) {
  juce::AudioBuffer<float> b(numCh, (int)ch[0].size());
  for (int c = 0; c < numCh; c++) b.copyFrom(c, 0, ch[c].data(), (int)ch[c].size());
  return b;
}
}  // namespace

juce::AudioProcessorValueTreeState::ParameterLayout ComboProcessor::createLayout() {
  using namespace juce;
  AudioProcessorValueTreeState::ParameterLayout l;
  auto f = [&](const char* id, const String& name, float lo, float hi, float step, float def, const String& unit = {}) {
    l.add(std::make_unique<AudioParameterFloat>(ParameterID { id, 1 }, name, NormalisableRange<float>(lo, hi, step), def,
                                                AudioParameterFloatAttributes().withLabel(unit)));
  };
  for (int i = 0; i < 6; i++) f(knobIds[i], knobNames[i], 0, 10, 0.1f, knobDefaults[i]);
  l.add(std::make_unique<AudioParameterChoice>(ParameterID { "input", 1 }, "Input", StringArray { "I - High Treble", "II - Normal", "Jumpered" }, 0));
  l.add(std::make_unique<AudioParameterChoice>(ParameterID { "sens", 1 }, "Input Jack", StringArray { "High", "Low" }, 0));
  l.add(std::make_unique<AudioParameterBool>(ParameterID { "trem", 1 }, "Tremolo", false));
  for (int i = 6; i < 8; i++) f(knobIds[i], knobNames[i], 0, 10, 0.1f, knobDefaults[i]);
  l.add(std::make_unique<AudioParameterBool>(ParameterID { "standby", 1 }, "Standby", false));
  f("trim", "Input Trim", -18, 18, 0.5f, 0, "dB");
  l.add(std::make_unique<AudioParameterBool>(ParameterID { "cab", 1 }, "Cabinet", true));
  l.add(std::make_unique<AudioParameterChoice>(ParameterID { "mic", 1 }, "Mic", StringArray { "Dynamic", "Ribbon" }, 0));
  f("pos", "Mic Position", 0, 1, 0.01f, 0.35f);
  f("dist", "Mic Distance", 0, 1, 0.01f, 0.15f);
  f("room", "Room", 0, 1, 0.01f, 0.12f);
  f("level", "Output Level", -40, 6, 0.5f, -12, "dB");
  l.add(std::make_unique<AudioParameterChoice>(ParameterID { "os", 1 }, "Oversampling", StringArray { "2x", "4x", "8x" }, 1));
  return l;
}

ComboProcessor::ComboProcessor()
  : AudioProcessor(BusesProperties()
                     .withInput("Input", juce::AudioChannelSet::mono(), true)
                     .withOutput("Output", juce::AudioChannelSet::stereo(), true)),
    apvts(*this, nullptr, "state", createLayout()) {
  for (int i = 0; i < bb::AmpCore::NUM_KNOBS; i++) pKnob[i] = apvts.getRawParameterValue(knobIds[i]);
  pInput = apvts.getRawParameterValue("input"); pSens = apvts.getRawParameterValue("sens");
  pTrem = apvts.getRawParameterValue("trem"); pStandby = apvts.getRawParameterValue("standby");
  pTrim = apvts.getRawParameterValue("trim"); pCab = apvts.getRawParameterValue("cab");
  pMic = apvts.getRawParameterValue("mic"); pPos = apvts.getRawParameterValue("pos");
  pDist = apvts.getRawParameterValue("dist"); pRoom = apvts.getRawParameterValue("room");
  pLevel = apvts.getRawParameterValue("level"); pOS = apvts.getRawParameterValue("os");
  startTimerHz(15);
}

ComboProcessor::~ComboProcessor() { stopTimer(); }

bool ComboProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const {
  const auto in = layouts.getMainInputChannelSet(), out = layouts.getMainOutputChannelSet();
  const bool inOk = in == juce::AudioChannelSet::mono() || in == juce::AudioChannelSet::stereo();
  const bool outOk = out == juce::AudioChannelSet::mono() || out == juce::AudioChannelSet::stereo();
  return inOk && outOk;
}

void ComboProcessor::prepareToPlay(double sampleRate, int samplesPerBlock) {
  const juce::ScopedLock il(irLock);
  const int os = 2 << juce::jlimit(0, 2, (int)pOS->load());
  {
    auto fresh = std::make_unique<bb::AmpCore>(sampleRate, os);
    const juce::ScopedLock cl(coreLock);
    core = std::move(fresh); coreOS = os;
  }
  rate = sampleRate;
  setLatencySamples(core->latencySamples());

  const int block = juce::jmax(64, samplesPerBlock);
  work.setSize(2, block); wet.setSize(2, block);
  const juce::dsp::ProcessSpec spec { sampleRate, (juce::uint32)block, 2 };

  // queue the impulse responses first so they are in place for the first block
  std::vector<float> rir[2];
  bb::makeRoomIR(sampleRate, rir);
  room.loadImpulseResponse(toBuffer(rir, 2), sampleRate, juce::dsp::Convolution::Stereo::yes,
                           juce::dsp::Convolution::Trim::no, juce::dsp::Convolution::Normalise::no);
  lastPos = -1;   // force a rebuild of the built-in cabinet at this rate
  applyCabIR();
  cab.prepare(spec); room.prepare(spec);
  cab.reset(); room.reset();

  roomGain.reset(sampleRate, 0.05); outGain.reset(sampleRate, 0.05); cabMix.reset(sampleRate, 0.03);
  const float r = pRoom->load();
  roomGain.setCurrentAndTargetValue(r * r * 1.6f);
  outGain.setCurrentAndTargetValue(juce::Decibels::decibelsToGain(pLevel->load()));
  cabMix.setCurrentAndTargetValue(pCab->load() > 0.5f ? 1.0f : 0.0f);
}

/* Loads the user's IR if there is one, otherwise the built-in cabinet for the current mic settings. */
void ComboProcessor::applyCabIR() {
  if (rate <= 0) return;
  using C = juce::dsp::Convolution;
  if (!userRaw.empty()) {
    std::vector<std::vector<float>> ch;
    for (const auto& c : userRaw) ch.push_back(bb::resampleIR(c, userRawRate, rate));
    if (bb::trimAndLevelIR(ch, rate)) {
      const int n = (int)ch.size();
      irInfo = juce::String(juce::roundToInt(1000.0 * (double)ch[0].size() / rate)) + " ms, " + (n == 2 ? "stereo" : "mono");
      cab.loadImpulseResponse(toBuffer(ch.data(), n), rate, n == 2 ? C::Stereo::yes : C::Stereo::no, C::Trim::no, C::Normalise::no);
      return;
    }
    userRaw.clear(); irInfo = "could not be used";
  }
  lastPos = pPos->load(); lastDist = pDist->load(); lastMic = (int)pMic->load();
  std::vector<float> ir[2];
  bb::makeCabIR(rate, lastPos, lastDist, lastMic == 1, ir);
  cab.loadImpulseResponse(toBuffer(ir, 2), rate, C::Stereo::yes, C::Trim::no, C::Normalise::no);
}

bool ComboProcessor::readIRFile(const juce::File& file) {
  juce::AudioFormatManager fm;
  fm.registerBasicFormats();
  std::unique_ptr<juce::AudioFormatReader> reader(fm.createReaderFor(file));
  if (reader == nullptr || reader->sampleRate <= 0 || reader->lengthInSamples <= 0) return false;
  const int nch = juce::jmin(2, (int)reader->numChannels);
  const int len = (int)std::min<juce::int64>(reader->lengthInSamples, (juce::int64)(reader->sampleRate * 4.0));
  juce::AudioBuffer<float> buf(nch, len);
  if (!reader->read(&buf, 0, len, 0, true, nch > 1)) return false;
  if (buf.getMagnitude(0, len) <= 0.0f) return false;
  userRaw.clear();
  for (int c = 0; c < nch; c++) userRaw.emplace_back(buf.getReadPointer(c), buf.getReadPointer(c) + len);
  userRawRate = reader->sampleRate;
  return true;
}

bool ComboProcessor::loadUserIR(const juce::File& file) {
  const juce::ScopedLock il(irLock);
  if (!readIRFile(file)) return false;
  irPath = file.getFullPathName(); irName = file.getFileName(); irInfo = {};
  applyCabIR();
  return !userRaw.empty();
}

void ComboProcessor::clearUserIR() {
  const juce::ScopedLock il(irLock);
  userRaw.clear(); irPath = {}; irName = {}; irInfo = {};
  applyCabIR();
}

juce::String ComboProcessor::irDescription() const {
  const juce::ScopedLock il(irLock);
  if (irName.isEmpty()) return "Built-in 2x12 (modelled)";
  return irName + (irInfo.isNotEmpty() ? "  (" + irInfo + ")" : juce::String());
}

/* Message thread: follow the oversampling choice and the mic controls. */
void ComboProcessor::timerCallback() {
  const juce::ScopedLock il(irLock);
  if (rate <= 0) return;
  const int os = 2 << juce::jlimit(0, 2, (int)pOS->load());
  if (os != coreOS) {
    auto fresh = std::make_unique<bb::AmpCore>(rate, os);
    {
      const juce::ScopedLock cl(coreLock);
      core.swap(fresh); coreOS = os;
    }
  }
  if (userRaw.empty() && (pPos->load() != lastPos || pDist->load() != lastDist || (int)pMic->load() != lastMic)) applyCabIR();
}

void ComboProcessor::processBlock(juce::AudioBuffer<float>& buffer, juce::MidiBuffer&) {
  juce::ScopedNoDenormals noDenormals;
  const int n = buffer.getNumSamples();
  const juce::ScopedTryLock cl(coreLock);
  if (!cl.isLocked() || core == nullptr || work.getNumSamples() == 0) { buffer.clear(); return; }

  using K = bb::AmpCore;
  for (int i = 0; i < K::NUM_KNOBS; i++) core->setKnob((Knob)i, pKnob[i]->load() * 0.1);
  core->mode = juce::jlimit(0, 2, (int)pInput->load());
  core->sens = pSens->load() > 0.5f ? 0.5 : 1.0;
  core->inScale = 2.0 * std::pow(10.0, pTrim->load() / 20.0);
  core->setTremolo(pTrem->load() > 0.5f);
  const float r = pRoom->load();
  roomGain.setTargetValue(r * r * 1.6f);
  outGain.setTargetValue(juce::Decibels::decibelsToGain(pLevel->load()));
  cabMix.setTargetValue(pCab->load() > 0.5f ? 1.0f : 0.0f);

  const int maxChunk = work.getNumSamples();
  for (int start = 0; start < n; start += maxChunk) processChunk(buffer, start, juce::jmin(maxChunk, n - start));
  htVolts.store((float)core->Vb, std::memory_order_relaxed);
}

void ComboProcessor::processChunk(juce::AudioBuffer<float>& buffer, int start, int len) {
  const int inCh = getTotalNumInputChannels(), outCh = buffer.getNumChannels();
  float* l = work.getWritePointer(0); float* rr = work.getWritePointer(1);

  // mono guitar in
  if (pStandby->load() > 0.5f || inCh == 0) {
    juce::FloatVectorOperations::clear(l, len);
  } else {
    juce::FloatVectorOperations::copy(l, buffer.getReadPointer(0, start), len);
    if (inCh > 1 && buffer.getNumChannels() > 1) {
      juce::FloatVectorOperations::add(l, buffer.getReadPointer(1, start), len);
      juce::FloatVectorOperations::multiply(l, 0.5f, len);
    }
  }

  core->process(l, l, len);
  juce::FloatVectorOperations::copy(rr, l, len);

  // cabinet (or the bare power-stage signal, attenuated, when it is switched off)
  float* wl = wet.getWritePointer(0); float* wr = wet.getWritePointer(1);
  if (cabMix.isSmoothing() || cabMix.getTargetValue() > 0.0f) {
    juce::FloatVectorOperations::copy(wl, l, len); juce::FloatVectorOperations::copy(wr, rr, len);
    float* chans[2] = { wl, wr };
    juce::dsp::AudioBlock<float> blk(chans, 2, (size_t)len);
    cab.process(juce::dsp::ProcessContextReplacing<float>(blk));
    for (int i = 0; i < len; i++) {
      const float m = cabMix.getNextValue(), d = 0.5f * (1.0f - m);
      l[i] = wl[i] * m + l[i] * d; rr[i] = wr[i] * m + rr[i] * d;
    }
  } else {
    juce::FloatVectorOperations::multiply(l, 0.5f, len); juce::FloatVectorOperations::multiply(rr, 0.5f, len);
  }

  // room
  if (roomGain.isSmoothing() || roomGain.getTargetValue() > 0.0f) {
    juce::FloatVectorOperations::copy(wl, l, len); juce::FloatVectorOperations::copy(wr, rr, len);
    float* chans[2] = { wl, wr };
    juce::dsp::AudioBlock<float> blk(chans, 2, (size_t)len);
    room.process(juce::dsp::ProcessContextReplacing<float>(blk));
    for (int i = 0; i < len; i++) { const float g = roomGain.getNextValue(); l[i] += wl[i] * g; rr[i] += wr[i] * g; }
  }

  // output level
  for (int i = 0; i < len; i++) { const float g = outGain.getNextValue(); l[i] *= g; rr[i] *= g; }

  if (outCh >= 2) {
    buffer.copyFrom(0, start, l, len); buffer.copyFrom(1, start, rr, len);
    for (int c = 2; c < outCh; c++) buffer.clear(c, start, len);
  } else if (outCh == 1) {
    juce::FloatVectorOperations::add(l, rr, len);
    buffer.copyFrom(0, start, l, len, 0.5f);
  }
}

juce::AudioProcessorEditor* ComboProcessor::createEditor() { return new ComboEditor(*this); }

void ComboProcessor::getStateInformation(juce::MemoryBlock& destData) {
  auto state = apvts.copyState();
  {
    const juce::ScopedLock il(irLock);
    state.setProperty("irPath", irPath, nullptr);
  }
  if (auto xml = state.createXml()) copyXmlToBinary(*xml, destData);
}

void ComboProcessor::setStateInformation(const void* data, int sizeInBytes) {
  auto xml = getXmlFromBinary(data, sizeInBytes);
  if (xml == nullptr || !xml->hasTagName(apvts.state.getType())) return;
  auto state = juce::ValueTree::fromXml(*xml);
  const juce::String path = state.getProperty("irPath").toString();
  apvts.replaceState(state);
  const juce::ScopedLock il(irLock);
  userRaw.clear(); irPath = {}; irName = {}; irInfo = {};
  if (path.isNotEmpty()) {
    const juce::File f(path);
    irPath = path; irName = f.getFileName();
    if (!readIRFile(f)) irInfo = "file not found, using built-in";
  }
  lastPos = -1;
  applyCabIR();
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() { return new ComboProcessor(); }
