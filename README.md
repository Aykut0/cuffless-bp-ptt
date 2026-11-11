# Cuffless Blood Pressure from ECG + PPG

> **An experimental project.** I built this to learn how cuffless blood-pressure estimation works and to see, honestly, how far a single pulse-transit-time measurement can go on real data.

Beat-to-beat systolic blood pressure estimated from **pulse transit time** alone, graded
against a **real invasive arterial line** on an intensive-care recording.

No cuff. No arterial pressure as an input. Calibrate on 90 seconds, predict the next
seven minutes from two ordinary wearable-grade signals.

| | |
|---|---|
| **Held-out MAE** | **3.0 mmHg** over 377 beats |
| **Correlation** | r = 0.45 |
| **Calibration** | first 90 s only (73 beats) |
| **Data** | BIDMC record 47, PhysioNet — real ICU patient, 125 Hz, 8 min |

![Cuffless estimate against the arterial line](figures/sbp_tracking.png)

---

## What this is

Pulse transit time is the delay between the heart firing and the pulse arriving at the
finger. Stiffer arteries carry the pulse faster, and arteries stiffen as pressure rises —
so the delay shortens when blood pressure goes up. That relationship is the basis of every
cuffless blood-pressure wearable.

This project measures that delay on a real patient and tests how far it actually gets you.

**The chain:**

1. Detect each heartbeat on the ECG (R-peak, 5–20 Hz band-pass).
2. Detect the pulse arriving at the finger — the point of **maximum upslope** on the PPG.
3. The delay between them is the pulse transit time.
4. Read true systolic pressure from the arterial line for the same beat.
5. Reject motion and line-flush artifacts with a MAD outlier test.
6. Fit the Bramwell–Hill relation `SBP = a/PTT² + b` on the **first 90 seconds only**.
7. Predict pressure for the remaining seven minutes from ECG + PPG, and score it.

![ECG, PPG and arterial pressure with per-beat markers](figures/ptt_alignment.png)

The top trace is the ECG with detected R-peaks; the middle is the finger PPG with the
maximum-upslope point marked on each beat; the bottom is the arterial line the estimate is
graded against. The horizontal gap between a red triangle and the following green dot *is*
the pulse transit time.

---

## Honest reading of the result

**3.0 mmHg mean absolute error is a good number. r = 0.45 is a modest one.** Both are real,
and the gap between them is the interesting part of this project.

The error is small mostly because this patient's pressure only moves across a ~20 mmHg
range, so any estimate near the mean scores well on absolute error. The correlation is what
actually measures whether the method tracks the *changes*, and 0.45 says it partly does.
That is the honest state of single-site PTT, and reporting it is more useful than tuning
until a number looks impressive.

### This is not an AAMI validation

The ANSI/AAMI/ISO 81060-2 standard for blood-pressure devices requires mean error ≤ 5 mmHg
and standard deviation ≤ 8 mmHg **measured across at least 85 subjects**. This is one
record from one patient. The 3.0 mmHg figure is a single-record result and should not be
read as meeting that standard.

### Sampling resolution was tested and ruled out

At 125 Hz, one sample is 8 ms. A sample-indexed PTT can therefore only take **9 distinct
values** across the entire record, which makes the raw estimate a visible staircase — an
obvious thing to blame for the weak correlation.

So both the R-peak and the PPG upslope point are refined to sub-sample precision by fitting
a parabola through the three points around each extremum. That takes PTT from 9 distinct
values to 450:

| | sample-indexed | sub-sample refined |
|---|---|---|
| distinct PTT values | 9 | 450 |
| held-out MAE | 3.00 mmHg | 2.99 mmHg |
| held-out r | 0.452 | 0.453 |

**It changes nothing.** The quantisation was real but was never the limiting factor — the
ceiling is the PTT–pressure relationship itself. The refinement stays in the code because it
is the correct thing to do and it removes the explanation; the table stays in this README
because ruling a hypothesis out is a result.

Run `python cuffless_bp_ptt.py --no-interp` to reproduce the quantised column.

### Other known limits

- **Single-site PTT ignores pre-ejection period.** The delay from R-peak to mechanical
  ejection varies with contractility and is folded into the measurement here.
- **Calibration does not persist.** The `a` and `b` coefficients drift over hours in
  reality; this record is eight minutes.
- **Maximum upslope, not foot.** Two conventions exist for pulse arrival. Maximum upslope
  is used here because it resists the baseline wander that dominates a finger PPG. They
  differ by a near-constant offset the calibration absorbs.

![Pulse transit time against arterial pressure](figures/ptt_vs_sbp.png)

---

## Running it

```bash
pip install -r requirements.txt
python cuffless_bp_ptt.py
```

Writes the three figures into `figures/` and the numbers into `metrics.txt`.

A MATLAB implementation is included as `cuffless_bp_ptt.m` and reproduces the same
numbers. Run it from this directory:

```matlab
cuffless_bp_ptt
```

---

## Data

`data/bidmc_47_Signals.csv` — record 47 of the **BIDMC PPG and Respiration Dataset**
(PhysioNet), an 8-minute recording at 125 Hz carrying ECG lead II, finger PPG and an
invasive arterial pressure waveform simultaneously. See [`data/README.md`](data/README.md)
for full citation and licence.

---

## Why this project exists

I work on wearable ultrasound blood-pressure monitoring ([FIWUS preprint](https://doi.org/10.21203/rs.3.rs-9000131/v1))
at the Bio-Integrated Microdevices Laboratory, Koç University, where pressure is recovered
from arterial wall motion measured by ultrasound. This project approaches the same clinical
problem from the opposite direction — optical and electrical signals a smartwatch already
has — which makes the trade-off concrete: PTT needs no transducer and no acoustic coupling,
and pays for it in exactly the correlation ceiling shown above.

---

## Licence

Code is MIT (see [LICENSE](LICENSE)). The BIDMC data carries its own licence and citation
requirements — see [`data/README.md`](data/README.md).

---

**Author:** Yaşar Aykut Altuntop · MSc student, Biomedical Sciences and Engineering, Koç University  
[LinkedIn](https://linkedin.com/in/aykutaltuntop) · [GitHub](https://github.com/Aykut0)
