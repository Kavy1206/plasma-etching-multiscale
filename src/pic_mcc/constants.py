"""Physical constants.

CODATA 2006 values, as required by Turner et al. (2013) sec. II ("Further
fundamental constants are to be represented by the 2006 CODATA values,
expressed to at least four places of decimals").

NOTE: the benchmark PRESCRIBES m_e = 9.109e-31 kg and m_i = 6.67e-27 kg
exactly, which are rounded versions of the CODATA/isotope values. For the
benchmark run we use the prescribed numbers, not these. See config.py.
"""

E_CHARGE = 1.602176487e-19   # C
EPS0 = 8.854187817e-12       # F/m
M_E = 9.10938215e-31         # kg
K_B = 1.3806504e-23          # J/K
AMU = 1.660538782e-27        # kg

# Convenience
EV = E_CHARGE                # 1 eV in joules
TORR = 133.322368            # Pa per Torr
MTORR = TORR * 1e-3
