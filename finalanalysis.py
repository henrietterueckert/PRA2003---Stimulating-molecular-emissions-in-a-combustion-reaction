# Output per pair: counts per event, overall asymmetry, asymmetry vs momentum P (with a straight-line fit),
# and the momentum spectrum dN/dP.
# Method (sub-sampling): each data file is one sub-sample. The result is the weighted mean over the
# sub-samples, and the uncertainty is the standard deviation of the sub-sample results.
# Put the data files in the same directory as this script!

import numpy as np
import matplotlib.pyplot as plt

# Setting
FILES = [f"output-Set{k}.txt" for k in range(1, 11)]       # output-Set1.txt ... output-Set10.txt

# (normal ID, variant ID)
PAIRS = [(211, -211), (321, -321), (2212, -2212), (3122, -3122), (3312, -3312), (3334, -3334)]

# ID -> (name, formula). The formula uses matplotlib notation: $_2$ = subscript, $^-$ = superscript
NAMES = {
    211:  ("Carbon monoxide", "CO"),        -211:  ("Carbon-13 monoxide", "$^{13}$CO"),
    321:  ("Nitric oxide", "NO"),           -321:  ("Ionised NO", ""),
    2212: ("Water", "H$_2$O"),              -2212: ("Heavy water", "D$_2$O"),
    3122: ("Methane", "CH$_4$"),            -3122: ("Methyl ion", "CH$_3^-$"),
    3312: ("Ethylene", "C$_2$H$_4$"),       -3312: ("Ionised ethylene", "C$_2$H$_3^-$"),
    3334: ("Ozone", "O$_3$"),               -3334: ("Superoxide anion", "O$_2^-$"),
}

# Count events with 0 molecules in the number of events? Only affects counts per event and dN/dP.
INCLUDE_EMPTY_EVENTS = False

# Momentum bin edges. The last bin is everything above 8 
BINS = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, np.inf])
N_BINS = len(BINS) - 1
MIN_COUNTS = 20     # a bin is marked "low statistics" if a sub-sample has fewer N + N' than this

# For plotting, draw the open last bin from 8 to 10
EDGES = BINS.copy()
EDGES[-1] = 10.0
CENTERS = (EDGES[:-1] + EDGES[1:]) / 2
WIDTHS = np.diff(EDGES)
P_LABEL = "P ($10^{-23}$ kg m/s)"

# Functions
def read_file(filename):
    """Read one data file. Returns (number of events, momentum of every molecule, ID of every molecule)."""
    n_events = 0
    momenta, ids = [], []
    with open(filename) as infile:
        for header in infile:  # header line of an event: "event ID, number of molecules"
            if not header.split(): # blank line at the end of the file
                break
            n_molecules = int(header.split()[1])
            if n_molecules > 0 or INCLUDE_EMPTY_EVENTS:
                n_events += 1
            for _ in range(n_molecules):  # one line per molecule: "px py pz ID"
                px, py, pz, code = next(infile).split()
                momenta.append(np.sqrt(float(px)**2 + float(py)**2 + float(pz)**2))
                ids.append(int(code))
    return n_events, np.array(momenta), np.array(ids)


def asymmetry(n_normal, n_variant):
    """A = (N - N') / (N + N'), and nan (not a number) where N + N' = 0."""
    n_normal = np.asarray(n_normal, float)
    n_variant = np.asarray(n_variant, float)
    total = n_normal + n_variant
    with np.errstate(invalid="ignore"):  
        return np.where(total > 0, (n_normal - n_variant) / total, np.nan)


def combine(values, weights):
    """Combine the sub-sample results: weighted mean, and std of the sub-samples (n - 1 in the denominator).
    Sub-samples with no result (nan) or weight 0 are skipped. Needs at least 2 sub-samples, otherwise nan."""
    ok = np.isfinite(values) & (weights > 0)
    if ok.sum() < 2:
        return np.nan, np.nan
    return np.average(values[ok], weights=weights[ok]), np.std(values[ok], ddof=1)


def n_sigma(value, error):
    """How many standard deviations the value is from zero."""
    return abs(value) / error if error > 0 else np.nan


# Reading the data
# counts[pair, 0 = normal / 1 = variant, file, bin] = number of molecules
events = np.zeros(len(FILES))
counts = np.zeros((len(PAIRS), 2, len(FILES), N_BINS))

for f, filename in enumerate(FILES):
    print("Reading", filename)
    events[f], momenta, ids = read_file(filename)
    for i, pair in enumerate(PAIRS):
        for j, code in enumerate(pair):
            # np.histogram counts how many momenta fall in each bin
            counts[i, j, f], _ = np.histogram(momenta[ids == code], bins=BINS)

print("\nEvents per sub-sample:", events.astype(int))
print(f"(max/min = {events.max() / events.min():.3f}; close to 1 means similar-sized sub-samples)")
print("Empty events", "INCLUDED" if INCLUDE_EMPTY_EVENTS else "EXCLUDED", "in the number of events")


# ANALYSIS PER PAIR
results = [] # one entry per pair

for i, (normal, variant) in enumerate(PAIRS):
    n_norm, n_var = counts[i, 0], counts[i, 1] # [file, bin]
    total_norm, total_var = n_norm.sum(axis=1), n_var.sum(axis=1)   # all momenta together, per file

    print("\n" + "=" * 70)
    print(f"Pair {normal}/{variant}: {NAMES[normal][0]} / {NAMES[variant][0]}")
    print("=" * 70)

    # Counts per event (weight = events per file)
    for code, total in ((normal, total_norm), (variant, total_var)):
        m, s = combine(total / events, events)
        print(f"ID {code:>6} ({NAMES[code][0]}): {m:.5f} +- {s:.5f} counts per event")

    # Overall asymmetry (weight = N + N' per file)
    a_all, s_all = combine(asymmetry(total_norm, total_var), total_norm + total_var)
    print(f"Overall asymmetry: {a_all:+.5f} +- {s_all:.5f}  ({n_sigma(a_all, s_all):.1f} sigma from zero)")

    # Asymmetry in each momentum bin (weight = N + N' in that bin, per file)
    A = np.full(N_BINS, np.nan)
    S = np.full(N_BINS, np.nan)
    print("\n  P range      asymmetry +- std     |A|/std  |A-A_overall|/std")
    for b in range(N_BINS):
        A[b], S[b] = combine(asymmetry(n_norm[:, b], n_var[:, b]), n_norm[:, b] + n_var[:, b])
        low = "  <- low statistics" if (n_norm[:, b] + n_var[:, b]).min() < MIN_COUNTS else ""
        print(f"  {BINS[b]:>3g} - {BINS[b + 1]:<4g}  {A[b]:+.5f} +- {S[b]:.5f}   "
              f"{n_sigma(A[b], S[b]):5.1f}   {n_sigma(A[b] - a_all, S[b]):5.1f}{low}")

    # Straight-line fit A(P) = intercept + slope * P, weights 1/std.
    use = np.isfinite(S) & (S > 0)
    use[-1] = False
    fit = None
    if use.sum() >= 3:      # a line has 2 parameters, so at least 3 points are needed
        (slope, intercept), cov = np.polyfit(CENTERS[use], A[use], 1, w=1 / S[use], cov="unscaled")
        slope_err = np.sqrt(cov[0, 0])
        fit = (slope, intercept, slope_err)
        print(f"\n  Slope dA/dP: {slope:+.5f} +- {slope_err:.5f} per 1e-23 kg m/s "
              f"({n_sigma(slope, slope_err):.1f} sigma from 0)")
    else:
        print("\n  Slope dA/dP: not enough bins for a fit")

    # Bins where the asymmetry is more than 3 sigma from zero
    strong = [f"{BINS[b]:g}-{BINS[b + 1]:g}" for b in range(N_BINS) if n_sigma(A[b], S[b]) > 3]
    print("  Bins with |A|/std > 3:", ", ".join(strong) if strong else "none")

    results.append(dict(A=A, S=S, a_all=a_all, s_all=s_all, fit=fit))



# PLOTS: one panel per pair, 3 panels per row
def title(code):
    name, formula = NAMES[code]
    return f"{name} ({formula})" if formula else name


n_cols = min(3, len(PAIRS))
n_rows = int(np.ceil(len(PAIRS) / n_cols))

# Asymmetry vs P 
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.5 * n_cols, 6 * n_rows), squeeze=False)
axes = axes.ravel()
for ax in axes[len(PAIRS):]:
    ax.set_visible(False) # hide empty panels

for ax, (normal, variant), r in zip(axes, PAIRS, results):
    # Leave out bins without a real error bar (too few molecules)
    shown = np.isfinite(r["S"]) & (r["S"] > 0)
    ax.errorbar(CENTERS[shown], r["A"][shown], yerr=r["S"][shown], xerr=WIDTHS[shown] / 2,
                fmt="o", color="tab:blue", capsize=3,
                label="Weighted mean of sub-samples (vertical bar = std, horizontal bar = bin width)")
    ax.axhline(0, color="gray", ls="--", lw=0.8)
    ax.axhspan(r["a_all"] - r["s_all"], r["a_all"] + r["s_all"], color="tab:orange", alpha=0.25,
               label="Overall asymmetry $\\pm$ std")
    if r["fit"] is not None:
        slope, intercept, slope_err = r["fit"]
        p = np.array([0, BINS[-2]]) # draw the line only over the fitted range 0-8
        ax.plot(p, intercept + slope * p, color="tab:green",
                label=f"Linear fit: slope = {slope:+.4f} $\\pm$ {slope_err:.4f} "
                      f"({n_sigma(slope, slope_err):.1f}$\\sigma$ from 0)")

    # Each panel its own y-range (rare pairs have much larger error bars)
    reach = np.abs(r["A"][shown]) + r["S"][shown]
    ymax = min(1.2 * reach.max(), 1.2) if shown.any() else 1.2
    ax.set_ylim(-ymax, ymax)
    ax.set_xlim(-0.25, 10.25)
    ax.set_xticks(EDGES[:-1])
    ax.set_xticklabels([f"{e:g}" for e in BINS[:-1]])
    ax.text(0.98, 0.02, "last bin: P > 8", transform=ax.transAxes, ha="right", fontsize=8, color="gray")
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel("(N - N') / (N + N')")
    ax.set_title(f"{title(normal)} vs {title(variant)}\n(IDs {normal} / {variant})", fontsize=10)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14))

plt.tight_layout()
plt.savefig("asymmetry_vs_P_1-9.png", dpi=150, bbox_inches="tight")
plt.show()

# Momentum spectrum dN/dP (counts per event per unit P)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(5.5 * n_cols, 4.5 * n_rows), squeeze=False)
axes = axes.ravel()
for ax in axes[len(PAIRS):]:
    ax.set_visible(False)

for i, (ax, pair) in enumerate(zip(axes, PAIRS)):
    for j, (code, colour) in enumerate(zip(pair, ("tab:blue", "tab:red"))):
        density = counts[i, j] / events[:, None] / WIDTHS # [file, bin]
        m, s = np.array([combine(density[:, b], events) for b in range(N_BINS - 1)]).T
        keep = m > 0 # a log scale cannot show 0
        ax.errorbar(CENTERS[:-1][keep], m[keep], yerr=s[keep], xerr=WIDTHS[:-1][keep] / 2,
                    fmt="o", color=colour, capsize=3, label=f"{title(code)}, ID {code}")
    ax.set_yscale("log")
    ax.set_xlabel(P_LABEL)
    ax.set_ylabel("dN/dP (counts per event per unit P)")
    ax.set_title(f"Momentum spectrum: {title(pair[0])} vs {title(pair[1])}\n(IDs {pair[0]} / {pair[1]})",
                 fontsize=10)
    ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("spectrum_dNdP_1-9.png", dpi=150)
plt.show()