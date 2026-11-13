cmake_minimum_required(VERSION 3.22)

set(CMAKE_OSX_DEPLOYMENT_TARGET "11.0" CACHE STRING "Minimum macOS version")
project(BluesCombo62 VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(MSVC)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")   # no VC++ runtime install needed
endif()

# ---- engine self-test: no JUCE needed (cmake -DBUILD_PLUGIN=OFF builds only this)
enable_testing()
add_executable(engine_test test/engine_test.cpp)
add_test(NAME engine_test COMMAND engine_test)

option(BUILD_PLUGIN "Build the audio plugin (downloads JUCE)" ON)
if(NOT BUILD_PLUGIN)
  return()
endif()

include(FetchContent)
FetchContent_Declare(JUCE
  GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
  GIT_TAG 8.0.15
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(JUCE)

if(APPLE)
  set(PLUGIN_FORMATS AU VST3 Standalone)
elseif(WIN32)
  set(PLUGIN_FORMATS VST3 Standalone)
else()
  set(PLUGIN_FORMATS VST3 LV2 Standalone)
endif()

juce_add_plugin(BluesCombo62
  PRODUCT_NAME "Blues Combo 62"
  COMPANY_NAME "Combo62"
  BUNDLE_ID "com.combo62.bluescombo62"
  PLUGIN_MANUFACTURER_CODE Cb62
  PLUGIN_CODE Bc62
  FORMATS ${PLUGIN_FORMATS}
  IS_SYNTH FALSE
  NEEDS_MIDI_INPUT FALSE
  NEEDS_MIDI_OUTPUT FALSE
  IS_MIDI_EFFECT FALSE
  EDITOR_WANTS_KEYBOARD_FOCUS FALSE
  COPY_PLUGIN_AFTER_BUILD FALSE
  AU_MAIN_TYPE kAudioUnitType_Effect
  VST3_CATEGORIES Fx Distortion
  LV2URI "urn:combo62:bluescombo62")

target_sources(BluesCombo62 PRIVATE
  src/PluginProcessor.cpp
  src/PluginEditor.cpp)

target_compile_definitions(BluesCombo62 PUBLIC
  JUCE_WEB_BROWSER=0
  JUCE_USE_CURL=0
  JUCE_VST3_CAN_REPLACE_VST2=0
  JUCE_DISPLAY_SPLASH_SCREEN=0)

target_link_libraries(BluesCombo62
  PRIVATE
    juce::juce_audio_utils
    juce::juce_dsp
  PUBLIC
    juce::juce_recommended_config_flags
    juce::juce_recommended_lto_flags)
