# Stimulated-echo T2 LUT correction

This toolbox creates a LUT table based on sequence parameters and slice profile using EPG simulations for a range of plausible T2 values in the human brain at 3T. 

This toolbox was developped to correct for T2 maps reconstructed with the GRAPPATINI Siemens sequence: the raw signal is not available and the T2 map is contaminated by stimulated echo times that cannot be corrected with a standard StimFit approach ([Lebel and Wilman, MRM 2010](doi.org/10.1002/mrm.22487)). 

Note on assumptions: 
- This approach assumes a negligible effect from B1         inhomogeneity. Further integration may be necessary to optimize the model for regions subject to strong B1 inhomogeneities. With additional measurements required. 
- The range of T2 values simulated for this LUT table is 30 to 130ms. This range may not match all field strength and tissue and may require adjustments.
- The "GRAPPATINI" mode assumes that the monoexponential fit applied during reconstruction skipped the first echo of the multi echo train length and was reconstructed from complex data, hence the simulation do not account for a noise floor. Applying the model to other sequence schemes may require further adjustments.
- The slice profile needs to be provided or simulated as per [StimFit](https://github.com/usc-mrel/StimFit/tree/main).

## Dependencies
- [StimFit](https://github.com/usc-mrel/StimFit/tree/main)

## Inputs
A LUT is created for a given protocol and must be updated if the protocol changes. Three input files are required to build the LUT: due to origins and sharing restrictions of proprietary data from scanner manufacturers, these files are not merged. 

| # | File | Origin | Shareable |
|---|------|--------|-----------|
| 1 | Acquisition parameters JSON | Protocol extracted from the scanner / sequence data | Yes |
| 2 | Pulse properties JSON | Manufacturer / might be provided with dicom metadata | No/Yes |
| 3 | RF waveform files (`.pta`) | Manufacturer | No |

### 1. Acquisition parameters

Describes the sequence protocol. A BIDS sidecar produced by
`dcm2bids`/`dcm2niix` already uses these names, so a sidecar can be used directly, provided it carries the required fields:

| Field | Type | Units | Meaning |
|---|---|---|---|
| `EchoTime` | array of numbers | seconds | The echo train. A single number is not enough — the LUT needs every echo. |
| `FlipAngle` | number | degrees | By DICOM convention: nominal **refocusing** flip angle (180 for a standard spin echo). |
| `SliceThickness` | number | millimetres | Sets the bounds of the simulated slice profile. |
| `ProtocolName` | string | — | Identifies the protocol. |

**Recommended** — not used by the simulation, but informs on *when the LUT must be rebuilt*. 

| Field | Type | Units | Meaning |
|---|---|---|---|
| `SoftwareVersions` | string | — | Scanner software line, e.g. `syngo MR E11`. |
| `SequenceName` | string | — | e.g. `WIP_se2d16`. |
| `PulseSequenceDetails` | string | — | e.g. `%CustomerSeq%\accT2map_wip899c`. |
| `Manufacturer`, `ManufacturersModelName` | string | — | Scanner make and model. |
| `MagneticFieldStrength` | number | tesla | |
| `EchoTrainLength` | integer | — | Number of echoes. Must equal `len(EchoTime)`. |
| `RepetitionTime` | number | seconds | Defines the sequence TR. |

Any other keys are carried through and ignored.

**Example** (trimmed; a real dcm2bids sidecar has many more keys):

```json
{
    "ProtocolName": "T2MAP_adapted_from_GRAPPATINI",
    "Manufacturer": "Siemens",
    "ManufacturersModelName": "Biograph_mMR",
    "MagneticFieldStrength": 3,
    "SoftwareVersions": "syngo MR E11",
    "SequenceName": "WIP_se2d16",
    "PulseSequenceDetails": "%CustomerSeq%\\accT2map_wip899c",
    "EchoTime": [0.0111, 0.0222, 0.0333, 0.0444, 0.0555, 0.0666, 0.0777, 0.0888,
                 0.0999, 0.111, 0.1221, 0.1332, 0.1443, 0.1554, 0.1665, 0.1776],
    "EchoTrainLength": 16,
    "FlipAngle": 180,
    "SliceThickness": 2,
    "RepetitionTime": 4.43
}
```

**Producing the acquisition JSON**

May be written out by hand from the table above, or using a dcm2bids sidecar. In the event that the echo train is missing, the sidecare may be merged with a protocol printout (fonctionality only tested with Siemens protocol PDF printout):

```bash
python utils/populate_json_sidecar_from_PDF.py <sidecar.json> <protocol.pdf> -o <out.json>
```

The sidecar is authoritative: fields present in both are kept from the sidecar and any disagreement is reported, not silently resolved. `EchoTime` is the exception and is taken from the printout (a sidecar written for a single echo holds only that echo's TE). Use `--force` to change which keys may be overwritten, `-i` to edit in place, and no flag at all to print to stdout without writing anything.

The tool will throw an error if the printout and the sidecar name different protocols.

**Scope:** PDF parsing is tested on **Siemens** protocol printouts only, on two host generations (syngo MR VB/VE and Numaris/X VA30). Other manufacturers' exports are not supported, and the tool fails rather than guessing at an unfamiliar layout; `--no_strict` downgrades this to a warning. 

Note also that a printout is a hand-made export — this route exists for protocols whose DICOMs are no longer available, not as a routine pipeline step.

### 2. Pulse properties JSON

RF pulse durations and the excitation flip angle may not be provided in dicom metadata or protocol printout but may be obtained from the manufacturer and supplied through in JSON format.

| Field | Type | Units | Meaning |
|---|---|---|---|
| `Excitation.FlipAngle` | number | degrees | Nominal excitation angle (90 for a spin echo). |
| `Excitation.FlipAngleUnit` | string | — | `"deg"`. |
| `Excitation.PulseDuration` | number | seconds | Duration of the excitation pulse. |
| `Excitation.PulseDurationUnits` | string | — | `"s"`. |
| `Refocusing.*` | | | The same four fields for the refocusing pulse. |

`Refocusing.FlipAngle` duplicates `FlipAngle` in the acquisition JSON on purpose the two are cross-checked, and a mismatch means the wrong pulse file has been paired with the protocol.

**Example.** The durations below are placeholders; substitute the values for your
own protocol.

```json
{
    "Excitation": {
        "FlipAngle": 90,
        "FlipAngleUnit": "deg",
        "PulseDuration": 0.0025,
        "PulseDurationUnits": "s"
    },
    "Refocusing": {
        "FlipAngle": 180,
        "FlipAngleUnit": "deg",
        "PulseDuration": 0.0038,
        "PulseDurationUnits": "s"
    }
}
```

### 3. RF waveform files

One waveform per pulse, giving the real slice profile. Siemens exports these as `.pta`; another vendor format may be substituted provided the LUT builder can read it. These are read with StimFit's `readRF_Siemens`.

## Generating pulse profiles

The function `build_opt_from_protocol.m` generates a slice profile for given input files as detailed above. 

**Usage**
```matlab
P = 'proprietary_data_folder';
opt = build_opt_from_protocol( ...
        fullfile(P,'acquisition_parameters.json'), ...
        fullfile(P,'pulse_properties.json'), ...
        fullfile(P,'excitation_pulse.pta'), ...
        fullfile(P,'refocusing_pulse.pta'), true);
```

Setting verbose to `true` will plot slice profiles obtained for a prescribed slice thickness and slice select gradient generated using StimFit bloch simulations. 

**Example slice profile for 2mm slice thickness**

<p align="center">
<img src="docs/figures/slice_profiles.png" alt="Slice profiles" width="50%"/>
</p>