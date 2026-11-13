# Blues Combo 62

A circuit-level model of a 1962-style 2x12 valve combo: two-channel ECC83 preamp, tone
stack, long-tail phase inverter, push-pull KT66 power stage with output transformer and
reactive speaker load, GZ34 supply sag, valve tremolo, and a 2x12 cabinet with an impulse
response loader.

It exists in two forms that share the same engine design:

- **Web version** - `index.html` + `amp-processor.js`. Serve the folder over HTTP
  (for example `python3 -m http.server`) and open it in a browser.
- **Plugin** - `plugin/`, a C++ port built with JUCE: Audio Unit and VST3 on macOS,
  VST3 on Windows, VST3 and LV2 on Linux.

## Getting the plugin

Every push to `main` builds all three platforms on GitHub Actions. Open the latest run
under the **Actions** tab and download the zip for your system from **Artifacts**, then
follow `INSTALL.txt` inside it (also at `plugin/INSTALL.txt`).

## Building it yourself

    cmake -S plugin -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release --parallel

JUCE 8 is downloaded automatically. On Linux, install the packages listed in
`.github/workflows/build.yml` first. Results land in `build/BluesCombo62_artefacts/Release/`.

To build and run only the engine self-test (no JUCE needed):

    cmake -S plugin -B build -DBUILD_PLUGIN=OFF && cmake --build build && ctest --test-dir build

## Layout

| Path | What |
|---|---|
| `plugin/src/AmpCore.h` | The amp engine (standard C++ only), a port of `amp-processor.js` |
| `plugin/src/CabIR.h` | Built-in cabinet and room impulse responses |
| `plugin/src/IRPrep.h` | Resampling, trimming and level-matching of loaded IRs |
| `plugin/src/Plugin*.cpp` | JUCE plugin wrapper and editor |
| `plugin/test/engine_test.cpp` | Engine self-test |

The plugin is built with [JUCE](https://juce.com), which is used here under its own licence terms.
