# Data provenance

## `bidmc_47_Signals.csv`

Record **47** of the **BIDMC PPG and Respiration Dataset**, obtained from PhysioNet.

An 8-minute recording sampled at 125 Hz from an adult intensive-care patient, carrying
simultaneously:

| Column | Signal |
|---|---|
| `Time [s]` | elapsed time |
| `RESP` | impedance respiration |
| `ABP` | invasive arterial blood pressure (mmHg) — **ground truth only** |
| `PLETH` | finger photoplethysmogram |
| `II`, `III`, `V` | ECG leads |

Only `II`, `PLETH` and `ABP` are used by this project. The `ABP` channel is used
exclusively as reference truth and never as a model input.

The dataset was extracted from the MIMIC II Waveform Database Matched Subset.

## Citation

If you use this data, cite both the dataset and PhysioNet:

> Pimentel, M.A.F., Johnson, A.E.W., Charlton, P.H., Birrenkott, D., Watkinson, P.J.,
> Tarassenko, L., and Clifton, D.A. "Toward a robust estimation of respiratory rate from
> pulse oximeters." *IEEE Transactions on Biomedical Engineering*, 64(8), 1914–1923, 2017.

> Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P.C., Mark, R., Mietus,
> J.E., Moody, G.B., Peng, C.K., and Stanley, H.E. "PhysioBank, PhysioToolkit, and
> PhysioNet: Components of a New Research Resource for Complex Physiologic Signals."
> *Circulation* 101(23), e215–e220, 2000.

## Licence

The BIDMC PPG and Respiration Dataset is published on PhysioNet under the
**Open Data Commons Attribution License v1.0 (ODC-BY 1.0)**, which permits redistribution
with attribution. This single record is redistributed here under that licence so the
project is runnable without a separate download.

> **Before making this repository public, confirm the licence on the dataset's PhysioNet
> page.** PhysioNet licence terms are per-dataset and can change between versions. If the
> terms do not permit redistribution, delete this CSV and replace it with a download step
> in the top-level README.

Source: <https://physionet.org/content/bidmc/>
