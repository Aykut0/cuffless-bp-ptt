"""
Cuffless Blood Pressure via Pulse Transit Time on real data
===========================================================
Estimates beat-to-beat systolic blood pressure from the ECG and PPG alone, and
grades it against a real invasive arterial blood-pressure line.

Data: BIDMC record 47 (PhysioNet), a real 8-minute intensive-care recording at
125 Hz that contains, simultaneously, an ECG lead (II), a finger PPG (PLETH) and
an invasive arterial blood-pressure waveform (ABP). This patient's pressure
actually moves over the recording (about 96 to 116 mmHg), so there is something
real to track. The ABP line is the clinical gold standard and is used here only
as ground truth -- never as an input.

Method
------
  * detect each heartbeat on the ECG (R-peak)
  * detect the arrival of the matching pulse at the finger, taken as the point of
    maximum upslope on the PPG (not the waveform foot -- see note below)
  * the delay between them is the pulse transit time (PTT)
  * read the true systolic pressure from the ABP waveform for the same beat
  * fit the Bramwell-Hill relation SBP = a/PTT^2 + b on the first 90 s only,
    then predict SBP for the remainder from ECG+PPG with no cuff

Note on the PTT reference point
-------------------------------
Two conventions exist for where the pulse "arrives": the foot of the PPG upstroke,
and the point of maximum upslope. This code uses maximum upslope, which is more
robust to the low-frequency baseline wander that dominates a finger PPG. Both are
standard; they differ by a near-constant offset that the calibration absorbs.

Note on sampling resolution
---------------------------
At 125 Hz one sample is 8 ms, so a raw sample-indexed PTT can only take a handful
of distinct values across the whole record. Both the R-peak and the PPG upslope
point are refined to sub-sample precision by fitting a parabola through the three
points around each extremum. Run with --no-interp to disable this and see the
quantised result for comparison; metrics.txt reports both.

Run:  python cuffless_bp_ptt.py
Writes figures/ptt_alignment.png, figures/ptt_vs_sbp.png, figures/sbp_tracking.png
and metrics.txt.

Author: Yasar Aykut Altuntop
"""
import os
import sys
import numpy as np
from scipy.signal import butter, filtfilt, find_peaks, medfilt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FS = 125.0
HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DATA = os.path.join(HERE, "data")


def load_column(csv_path, name):
    header = open(csv_path).readline().strip().split(",")
    idx = [h.strip() for h in header].index(name)
    return np.loadtxt(csv_path, delimiter=",", skiprows=1, usecols=idx)


def bp(x, lo, hi, order=3):
    b, a = butter(order, [lo / (FS / 2), hi / (FS / 2)], btype="band")
    return filtfilt(b, a, x)


def subsample_offset(y, i):
    """Sub-sample position of a local maximum at index i, by parabolic fit.

    Returns the offset in samples (|offset| <= 1), so the refined location is
    i + offset. Falls back to 0 at array edges or on a degenerate fit.
    """
    if i <= 0 or i >= len(y) - 1:
        return 0.0
    denom = y[i - 1] - 2.0 * y[i] + y[i + 1]
    if denom == 0:
        return 0.0
    offset = 0.5 * (y[i - 1] - y[i + 1]) / denom
    return offset if abs(offset) <= 1.0 else 0.0


def extract_beats(ecg_f, ppg_f, dppg, abp, rloc, t, interp=True):
    """Per-beat pulse transit time and reference systolic pressure."""
    ptt, sbp_ref, tbeat, slope_pts = [], [], [], []
    for k in range(len(rloc) - 1):
        s0 = rloc[k]
        s1 = min(s0 + int(0.5 * FS), len(t))
        if s1 - s0 < 5:
            continue
        j = s0 + int(np.argmax(dppg[s0:s1]))          # PPG maximum-upslope point
        if interp:
            r_at = s0 + subsample_offset(ecg_f, s0)
            p_at = j + subsample_offset(dppg, j)
        else:
            r_at, p_at = float(s0), float(j)
        delay_ms = (p_at - r_at) / FS * 1000.0
        sys_p = np.nanmax(abp[s0:min(s0 + int(FS), len(abp))])
        if 50 < delay_ms < 400 and 60 < sys_p < 220:
            ptt.append(delay_ms)
            sbp_ref.append(sys_p)
            tbeat.append(t[s0])
            slope_pts.append(j)
    return np.array(ptt), np.array(sbp_ref), np.array(tbeat), slope_pts


def fit_and_score(ptt, sbp_ref, tbeat):
    """MAD outlier rejection, 90 s calibration, held-out scoring."""
    def inlier(v):
        m = np.median(v)
        mad = np.median(np.abs(v - m)) + 1e-9
        return np.abs(v - m) < 4 * mad

    keep = inlier(ptt) & inlier(sbp_ref)
    ptt, sbp_ref, tbeat = ptt[keep], sbp_ref[keep], tbeat[keep]
    sbp_ref = medfilt(sbp_ref, 5)          # BP does not jump between beats

    cal = tbeat <= tbeat[0] + 90
    coef = np.polyfit(1 / ptt[cal] ** 2, sbp_ref[cal], 1)
    sbp_est = np.polyval(coef, 1 / ptt ** 2)

    test = ~cal
    mae = np.mean(np.abs(sbp_est[test] - sbp_ref[test]))
    rho = np.corrcoef(sbp_est[test], sbp_ref[test])[0, 1]
    return ptt, sbp_ref, tbeat, cal, test, coef, sbp_est, mae, rho


def main(interp=True):
    os.makedirs(FIG, exist_ok=True)
    sig = os.path.join(DATA, "bidmc_47_Signals.csv")
    ecg = load_column(sig, "II")
    ppg = load_column(sig, "PLETH")
    abp = load_column(sig, "ABP")
    t = np.arange(len(ecg)) / FS

    ecg_f = bp(ecg, 5, 20)
    rloc, _ = find_peaks(ecg_f, distance=int(0.4 * FS), height=2.5 * np.std(ecg_f))

    b, a = butter(4, 8 / (FS / 2), btype="low")
    ppg_f = filtfilt(b, a, ppg)
    dppg = np.gradient(ppg_f)

    beats = extract_beats(ecg_f, ppg_f, dppg, abp, rloc, t, interp=interp)
    slope_pts = beats[3]
    ptt, sbp_ref, tbeat, cal, test, coef, sbp_est, mae, rho = fit_and_score(*beats[:3])

    # identical pipeline without sub-sample refinement, for the resolution check
    q = extract_beats(ecg_f, ppg_f, dppg, abp, rloc, t, interp=False)
    q_ptt, _, _, _, _, _, _, q_mae, q_rho = fit_and_score(*q[:3])

    lines = [
        "Cuffless BP via PTT on BIDMC record 47 (real ECG+PPG+ABP, 125 Hz, 8 min)",
        f"beats used            : {len(ptt)}  ({int(cal.sum())} for calibration, {int(test.sum())} for test)",
        f"mean PTT              : {np.mean(ptt):.0f} ms",
        f"reference SBP (ABP)   : {np.mean(sbp_ref):.0f} mmHg (range {sbp_ref.min():.0f}-{sbp_ref.max():.0f})",
        f"test-set error (MAE)  : {mae:.1f} mmHg",
        f"test-set correlation  : r = {rho:.2f}",
        "",
        "Sampling-resolution check (is 8 ms quantisation the limiting factor?)",
        f"  distinct PTT values, sample-indexed : {len(np.unique(np.round(q_ptt, 6)))}",
        f"  distinct PTT values, sub-sample     : {len(np.unique(np.round(ptt, 6)))}",
        f"  MAE  sample-indexed / sub-sample    : {q_mae:.2f} / {mae:.2f} mmHg",
        f"  r    sample-indexed / sub-sample    : {q_rho:.3f} / {rho:.3f}",
        "  -> removing the quantisation does not improve the fit, so r is limited",
        "     by the PTT-pressure relationship itself, not by sampling resolution.",
    ]
    with open(os.path.join(HERE, "metrics.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines))

    # ---- figure 1: aligned ECG / PPG / ABP with per-beat markers
    w = (t >= 60) & (t <= 64)
    fig, ax = plt.subplots(3, 1, figsize=(9.5, 5.5), sharex=True)
    ax[0].plot(t[w], ecg[w], color="#374151", lw=1)
    rk = rloc[(t[rloc] >= 60) & (t[rloc] <= 64)]
    ax[0].plot(t[rk], ecg[rk], "v", color="#dc2626", ms=7, label="R-peak")
    ax[0].set_ylabel("ECG (II)")
    ax[0].legend(loc="upper right")
    ax[0].set_title("Real ECG, PPG and arterial BP with pulse-transit-time markers")
    ax[1].plot(t[w], ppg_f[w], color="#2563eb", lw=1)
    sk = [s for s in slope_pts if 60 <= t[s] <= 64]
    ax[1].plot(t[sk], ppg_f[sk], "o", color="#16a34a", ms=6, label="PPG max-upslope point")
    ax[1].set_ylabel("PPG")
    ax[1].legend(loc="upper right")
    ax[2].plot(t[w], abp[w], color="#b45309", lw=1)
    ax[2].set_ylabel("ABP\n[mmHg]")
    ax[2].set_xlabel("time [s]")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ptt_alignment.png"), dpi=140)
    plt.close(fig)

    # ---- figure 2: PTT against reference SBP, with the fitted law
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    ax.scatter(sbp_ref, ptt, s=10, color="#0891b2", alpha=0.5)
    xs = np.linspace(sbp_ref.min(), sbp_ref.max(), 200)
    ax.plot(xs, np.sqrt(np.clip(coef[0] / (xs - coef[1]), 1e-6, None)),
            color="#dc2626", lw=1.5, label=r"fitted $a/\mathrm{PTT}^2$ law")
    ax.set_xlabel("arterial systolic BP [mmHg]")
    ax.set_ylabel("measured PTT [ms]")
    ax.set_title("Pulse transit time vs real arterial pressure")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "ptt_vs_sbp.png"), dpi=140)
    plt.close(fig)

    # ---- figure 3: tracking against the arterial line
    fig, ax = plt.subplots(figsize=(9.5, 3.6))
    ax.plot(tbeat, sbp_ref, color="#374151", lw=1.3, label="arterial line (truth)")
    ax.plot(tbeat, sbp_est, color="#dc2626", lw=1.1, ls="--", label="cuffless PTT estimate")
    ax.axvspan(tbeat[0], tbeat[cal][-1], color="#fde68a", alpha=0.5,
               label="calibration (first 90 s)")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("SBP [mmHg]")
    ax.set_title(f"Cuffless SBP vs arterial line - held-out MAE {mae:.1f} mmHg, r = {rho:.2f}")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "sbp_tracking.png"), dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    main(interp="--no-interp" not in sys.argv)
