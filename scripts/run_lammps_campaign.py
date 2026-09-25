import json
import sys
import time

sys.path.insert(0, ".")
from src.md_lammps.run_impact_series import run_series

energy = float(sys.argv[1])
seed = int(sys.argv[2])
t0 = time.time()
res = run_series(energy_ev=energy, angle_deg=0.0, n_shots=10, seed=seed,
                 workdir="runs/lammps/campaign", n_impact_steps=1200, n_rethermalize_steps=300)
out = dict(energy_ev=energy, sputtered_per_shot=res, total=sum(res), n_shots=len(res),
          yield_mean=sum(res) / len(res), elapsed_s=round(time.time() - t0, 1))
with open(f"runs/lammps/campaign/result_e{energy:g}.json", "w") as f:
    json.dump(out, f, indent=2)
print(json.dumps(out))
