"""Advance the Turner Case 1 run by ONE wall-clock-bounded chunk and exit.
Call this repeatedly (each call is a separate foreground process) until it
prints DONE. Safe against this environment's ~300s per-command limit."""
import json
import sys
import time

from src.pic_mcc.config import turner_case
from src.pic_mcc.simulation import Simulation

budget = float(sys.argv[1]) if len(sys.argv) > 1 else 250.0

cfg = turner_case(1)
sim = Simulation(cfg, "runs/case1.npz")
t0 = time.time()
res = sim.run(max_wall_seconds=budget, progress_every=4000, checkpoint_every=4000)
dt = time.time() - t0
print(json.dumps({**res, "chunk_wall_s": round(dt, 1),
                  "pct": round(100 * res["step"] / res["n_steps"], 2)}))
if res["step"] >= res["n_steps"]:
    print("DONE")
