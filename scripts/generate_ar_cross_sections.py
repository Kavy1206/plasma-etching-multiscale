"""Generate data/cross_sections/ar/*.csv from the analytic Phelps & Petrovic
(1999) / Phelps (1994) argon formulas, transcribed from the eduPIC reference
code (Donko et al. 2021, PSST 30 095017, github.com/donkozoltan/eduPIC,
functions set_electron_cross_sections_ar / set_ion_cross_sections_ar).

See data/cross_sections/ar/PROVENANCE.md for the full citation and for notes
on the lumped-excitation and derived-backscatter simplifications, both of
which are eduPIC's own modeling choices, carried over here unchanged.

Run: python -m scripts.generate_ar_cross_sections
"""
import numpy as np

E_EXC_TH = 11.5   # eV, lumped effective excitation threshold
E_ION_TH = 15.8   # eV, ionization threshold


def qmel(en):
    en = np.asarray(en, dtype=float)
    return 1e-20 * (np.abs(6.0 / (1.0 + (en / 0.1) + (en / 0.6) ** 2.0) ** 3.3
                    - 1.1 * en ** 1.4 / (1.0 + (en / 15.0) ** 1.2)
                    / np.sqrt(1.0 + (en / 5.5) ** 2.5 + (en / 60.0) ** 4.1))
                    + 0.05 / (1.0 + en / 10.0) ** 2.0
                    + 0.01 * en ** 3.0 / (1.0 + (en / 12.0) ** 6.0))


def qexc(en):
    en = np.asarray(en, dtype=float)
    out = np.zeros_like(en)
    m = en > E_EXC_TH
    x = en[m]
    out[m] = 1e-20 * (0.034 * (x - 11.5) ** 1.1 * (1.0 + (x / 15.0) ** 2.8) / (1.0 + (x / 23.0) ** 5.5)
                      + 0.023 * (x - 11.5) / (1.0 + x / 80.0) ** 1.9)
    return out


def qion(en):
    en = np.asarray(en, dtype=float)
    out = np.zeros_like(en)
    m = en > E_ION_TH
    x = en[m]
    out[m] = 1e-20 * (970.0 * (x - 15.8) / (70.0 + x) ** 2.0 + 0.06 * (x - 15.8) ** 2.0 * np.exp(-x / 9.0))
    return out


def qiso(e_lab):
    e_lab = np.asarray(e_lab, dtype=float)
    return 2e-19 * e_lab ** -0.5 / (1.0 + e_lab) + 3e-19 * e_lab / (1.0 + e_lab / 3.0) ** 2.0


def qmom(e_lab):
    e_lab = np.asarray(e_lab, dtype=float)
    return 1.15e-18 * e_lab ** -0.1 * (1.0 + 0.015 / e_lab) ** 0.6


def qback(e_lab):
    # Phelps (1994)'s own decomposition of the momentum-transfer cross section
    # into isotropic + backward parts. Goes slightly negative at a few energies
    # from fit/extrapolation edge effects -- floored at zero below.
    return (qmom(e_lab) - qiso(e_lab)) / 2.0


def main():
    import os
    outdir = "data/cross_sections/ar"
    os.makedirs(outdir, exist_ok=True)

    e_elec = np.concatenate([[1e-4], np.linspace(0.01, 1000.0, 20000)])
    e_ion = np.concatenate([[1e-4], np.linspace(0.01, 2000.0, 20000)])

    tables = {
        "Elastic_Ar.csv": (e_elec, qmel(e_elec)),
        "Excitation_Ar.csv": (e_elec, qexc(e_elec)),
        "Ionization_Ar.csv": (e_elec, qion(e_elec)),
        "Isotropic_Ar.csv": (e_ion, qiso(e_ion)),
        "Backscattering_Ar.csv": (e_ion, np.maximum(qback(e_ion), 0.0)),
    }
    for name, (e, sigma) in tables.items():
        np.savetxt(f"{outdir}/{name}", np.column_stack([e, sigma]), delimiter=";", fmt="%.6e")
        print(f"wrote {outdir}/{name}  ({len(e)} points, max sigma={sigma.max():.3e} m^2)")


if __name__ == "__main__":
    main()
