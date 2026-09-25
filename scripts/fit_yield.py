import json
import glob

import numpy as np
from scipy.optimize import curve_fit
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

results = {}
for f in glob.glob("runs/lammps/campaign/result_e*.json"):
    d = json.load(open(f))
    results[d["energy_ev"]] = d

energies = np.array(sorted(results.keys()))
yields = np.array([results[e]["yield_mean"] for e in energies])
n_shots = np.array([results[e]["n_shots"] for e in energies])
totals = np.array([results[e]["total"] for e in energies])
# Binomial-ish standard error on a rate estimated from n_shots trials
yerr = np.sqrt(np.maximum(totals, 1)) / n_shots

def yield_fn(E, A, Eth):
    return np.maximum(A * (np.sqrt(np.maximum(E - Eth, 0))), 0.0)

try:
    popt, pcov = curve_fit(yield_fn, energies, yields, p0=[0.03, 20.0],
                           sigma=np.maximum(yerr, 0.03), bounds=([0.0, 0.0], [1.0, 60.0]))
    A, Eth = popt
    perr = np.sqrt(np.diag(pcov))
except Exception as ex:
    A, Eth, perr = None, None, None
    print("fit failed:", ex)

print(f"Data: E={list(energies)}  Y={list(yields)}  n_shots={list(n_shots)}  totals={list(totals)}")
if A is not None:
    print(f"Fit: A={A:.4f} +/- {perr[0]:.4f}   Eth={Eth:.1f} +/- {perr[1]:.1f} eV")

# Literature anchors from theory.md sec 7
anchor_E = [100, 500, 1000]
anchor_Y = [0.07, 0.65, 0.93]

fig, ax = plt.subplots(figsize=(7,5))
ax.errorbar(energies, yields, yerr=yerr, fmt='o', color='#d62728', capsize=3,
           label=f'this MD (10 impacts/energy, reduced statistics)')
for e, y, n in zip(energies, yields, n_shots):
    ax.annotate(f'n={n}', (e, y), textcoords='offset points', xytext=(6,6), fontsize=8)
if A is not None:
    Efine = np.linspace(0, 550, 300)
    ax.plot(Efine, yield_fn(Efine, A, Eth), 'r--', lw=1.5,
           label=f'fit: A={A:.3f}, $E_{{th}}$={Eth:.0f} eV')
ax.plot(anchor_E, anchor_Y, 'k^', ms=10, label='literature anchors (theory.md sec.7)')
ax.set_xlabel('Ar+ energy (eV)')
ax.set_ylabel('sputter yield Y (atoms/ion)')
ax.set_title('Ar+ -> Si(100) sputter yield: this MD vs literature')
ax.legend(fontsize=9, frameon=False)
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('figures/sputter_yield.png', dpi=150)

np.savez('runs/lammps/campaign/yield_summary.npz', energies=energies, yields=yields,
        n_shots=n_shots, totals=totals, A=A, Eth=Eth)
print("saved figures/sputter_yield.png")
