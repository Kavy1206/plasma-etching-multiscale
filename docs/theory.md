# Theory foundation — *Plasma to Profile*

**Phase 0 deliverable.** Everything here is either derived from scratch or sourced
with a citation. Where a number is a literature value that this project has *not*
yet independently reproduced, it is marked **[lit]**. Where it is a prediction that
Phase 1 or Phase 2 is supposed to test, it is marked **[predict]**.

Section 9 lists four places where the original project plan does not match the
benchmark paper. Read that section before writing any solver code.

---

## 0. Notation and constants

| Symbol | Meaning | Value / units |
|---|---|---|
| $e$ | elementary charge | $1.602177\times10^{-19}$ C |
| $\varepsilon_0$ | vacuum permittivity | $8.854188\times10^{-12}$ F/m |
| $m_e$ | electron mass | $9.109\times10^{-31}$ kg |
| $M$ | ion mass | He⁺: $6.67\times10^{-27}$ kg; Ar⁺: $6.63\times10^{-26}$ kg |
| $T_e$ | electron temperature | quoted in **volts** unless written $k_BT_e$ |
| $n_0$ | bulk (mid-plane) plasma density | m⁻³ |
| $n_s$ | density at the sheath edge | m⁻³ |
| $\lambda_D$ | Debye length | $\sqrt{\varepsilon_0 T_e / (e n)}$ |
| $\omega_{pe}$ | electron plasma frequency | $\sqrt{e^2 n_e/(\varepsilon_0 m_e)}$ |
| $u_B$ | Bohm speed | $\sqrt{e T_e/M}$ |
| $s$ | sheath thickness | m |

Two handy engineering forms, used throughout and worth memorising:

$$\lambda_D\,[\mathrm{m}] = 7434\sqrt{\frac{T_e\,[\mathrm{eV}]}{n\,[\mathrm{m^{-3}}]}},
\qquad
\omega_{pe}\,[\mathrm{rad/s}] = 5.64\times10^{4}\sqrt{n_e\,[\mathrm{cm^{-3}}]}$$

### The benchmark operating point (Turner et al. 2013, Case 1)

All numbers below are from Table I and Table III of the paper. These are the
targets Phase 1 must hit.

| Quantity | Case 1 value |
|---|---|
| Gas | **helium** (see §9.1) |
| Gap $L$ | 6.7 cm |
| Neutral density $N$ | $9.64\times10^{20}$ m⁻³ at 300 K $\Rightarrow$ 4.00 Pa $\approx$ **30 mTorr** |
| Drive | $V(t) = 450\sin(2\pi\cdot13.56\times10^6 t)$ V, zero at $t=0$ |
| Secondary emission | **none** — all fluxes fully absorbed, $\gamma=0$ |
| Initial $n_0$ / $T_e$ / $T_i$ | $2.56\times10^{14}$ m⁻³ / 30 000 K / 300 K |
| Cells | $\Delta x = L/128 = 523.4\ \mu$m, 512 particles/cell initially |
| Time step | $\Delta t = (400f)^{-1} = 184.4$ ps |
| Run length | 512 000 steps = **1280 RF cycles**, averaged over the last **32 cycles** |
| **Result: $n_i$ at mid-plane** | $\mathbf{1.40\times10^{14}}$ **m⁻³** |
| **Result: $k_BT_e$ at mid-plane** | $\mathbf{9.36}$ **eV** |
| Result: $J_i$ at either electrode | 0.219 A/m² |
| Result: $\omega_{pe}\Delta t$ / $\lambda_D/\Delta x$ | 0.121 / 3.72 |

Sanity check on our own formulas, before we trust them anywhere else:
at $n_e=1.363\times10^{14}$ m⁻³ and $T_e=9.36$ eV,
$\lambda_D = 7434\sqrt{9.36/1.363\times10^{14}} = 1.95$ mm, so
$\lambda_D/\Delta x = 1.95/0.523 = 3.72$ ✓, and
$\omega_{pe} = 5.64\times10^4\sqrt{1.363\times10^{8}} = 6.59\times10^{8}$ rad/s, so
$\omega_{pe}\Delta t = 0.121$ ✓. Both reproduce the paper's Table III exactly. Our
unit conventions are right.

---

## 1. Why a sheath exists at all

Put a plasma in contact with a wall. Electron thermal speed
$\bar v_e=\sqrt{8k_BT_e/\pi m_e}$ exceeds the ion thermal speed by
$\sqrt{M T_e/(m_e T_i)}$ — for He at $T_e=9.36$ eV, $T_i=0.026$ eV that is a factor
of $\sim\!2500$. Electrons therefore hit the wall far faster than ions. The wall
charges negative, the plasma floats positive, and a thin positive-space-charge layer
(the **sheath**) forms that reflects most electrons and accelerates ions. In steady
state the sheath potential self-adjusts until the two wall fluxes are equal.

The sheath is where essentially all of the interesting etch physics lives: it sets
the ion energy hitting the wafer and, because the field there is normal to the
surface, the ion directionality that makes anisotropic etching possible at all.

---

## 2. The Bohm criterion

### Setup

1D, planar, collisionless sheath. Let $x=0$ be the sheath edge, $\phi(0)=0$,
$\phi(x)<0$ inside. Ions are cold and enter monoenergetically with speed $u_s$.
Electrons are in Boltzmann equilibrium with the potential:

$$n_e(x) = n_s\exp\!\left(\frac{e\phi}{k_BT_e}\right)$$

Ion flux conservation ($n_i u = n_s u_s$) plus energy conservation
($\tfrac12 M u^2 = \tfrac12 M u_s^2 - e\phi$) give

$$u(x) = u_s\sqrt{1 - \frac{2e\phi}{Mu_s^2}}
\qquad\Longrightarrow\qquad
n_i(x) = n_s\left(1 - \frac{2e\phi}{Mu_s^2}\right)^{-1/2}$$

Poisson's equation:

$$\frac{d^2\phi}{dx^2} = \frac{e n_s}{\varepsilon_0}
\left[\exp\!\left(\frac{e\phi}{k_BT_e}\right)
- \left(1-\frac{2e\phi}{Mu_s^2}\right)^{-1/2}\right]$$

### First integral

Multiply by $d\phi/dx$ and integrate from $0$ to $x$, using $\phi(0)=0$ and
$\phi'(0)\simeq0$ (the sheath edge is where the field is still weak):

$$\frac12\left(\frac{d\phi}{dx}\right)^{2}
= \frac{e n_s}{\varepsilon_0}\left\{
\frac{k_BT_e}{e}\left[e^{e\phi/k_BT_e}-1\right]
+ \frac{Mu_s^2}{e}\left[\sqrt{1-\frac{2e\phi}{Mu_s^2}}-1\right]\right\}$$

### The inequality

The left side is a square, so the brace must be $\ge 0$ for every $\phi<0$ in the
sheath. Expand to second order in small $|\phi|$, with
$X \equiv e\phi/k_BT_e$ and $Y \equiv 2e\phi/(Mu_s^2)$:

- $e^{X}-1 \approx X + X^2/2$, so the first term $\approx \phi + \dfrac{e\phi^2}{2k_BT_e/e}$
- $\sqrt{1-Y}-1 \approx -Y/2 - Y^2/8$, so the second term $\approx -\phi - \dfrac{e\phi^2}{2Mu_s^2}$

The linear terms cancel exactly — that is the whole point, and it is why the
criterion is a *second-order* statement. What is left is

$$\frac{e\phi^{2}}{2}\left[\frac{e}{k_BT_e} - \frac{1}{Mu_s^{2}}\right]\ \ge\ 0
\qquad\Longrightarrow\qquad
M u_s^{2} \ge k_BT_e$$

$$\boxed{\;u_s \ \ge\ u_B \equiv \sqrt{\frac{k_BT_e}{M}} = \sqrt{\frac{eT_e}{M}}\;}$$

### What this means physically

Inside the sheath, moving toward the wall, ions speed up (so $n_i$ falls) and
electrons are repelled (so $n_e$ falls). For net positive space charge to be
sustained, $n_e$ must fall *faster* than $n_i$. Electrons fall off exponentially
with a scale set by $T_e$; ions fall off algebraically with a scale set by their
kinetic energy. The Bohm condition is exactly the statement that the ion kinetic
energy at entry is large enough to make the ion falloff the slower of the two.
If ions entered too slowly, the solution would oscillate instead of decaying
monotonically, and no monotonic sheath could exist.

### Presheath and the 0.61 factor

Ions must be accelerated to $u_B$ *before* the sheath, which requires a weak
quasi-neutral **presheath** with a potential drop of $T_e/2$:
$\tfrac12 M u_B^2 = \tfrac12 k_BT_e = e\cdot(T_e/2)$. Through the presheath the
Boltzmann relation still holds, so

$$n_s = n_0\,e^{-1/2} \approx 0.61\,n_0$$

### **[predict]** Test this against the benchmark

For Case 1, $T_e=9.36$ eV, He⁺:
$u_B=\sqrt{1.602\times10^{-19}\times9.36/6.67\times10^{-27}} = 1.50\times10^{4}$ m/s.
With $n_s = 0.61\times1.40\times10^{14} = 8.55\times10^{13}$ m⁻³:

$$J_i = e n_s u_B = 0.205\ \mathrm{A/m^2}$$

The benchmark reports **0.219 A/m²**. That is 6.4% low, which is about what you
expect from a formula that assumes a sharp sheath edge and a flat presheath in a
discharge where the sheath occupies a good fraction of the gap. **This is a real,
cheap, first-day validation check** — if your PIC gives a wall ion flux that is
off by a factor of two from $0.61\,n_0 e u_B$, the bug is in the solver, not in
the physics.

Note that the Bohm criterion is *emergent* in a PIC simulation. You do not impose
it. Plotting the ion drift velocity vs. $x$ and confirming it crosses $u_B$ right
where quasi-neutrality breaks is one of the strongest "this code is correct"
plots you can put in a README.

---

## 3. Sheath thickness: matrix → Child–Langmuir

### 3a. Matrix sheath (the crude version)

Assume the sheath is a slab of uniform ion density $n_s$ with no electrons at all.
$d^2\phi/dx^2 = -e n_s/\varepsilon_0$, integrate twice with $\phi(0)=0$,
$\phi'(0)=0$, $\phi(s)=-V_0$:

$$V_0 = \frac{e n_s s^2}{2\varepsilon_0}
\qquad\Longrightarrow\qquad
s = \lambda_{Ds}\sqrt{\frac{2V_0}{T_e}}$$

Useful for order-of-magnitude, wrong in detail because it ignores that ions
accelerate (and therefore thin out) as they cross.

### 3b. Child–Langmuir sheath

Now keep ion acceleration. Current density $J = e n_i u$ is constant across the
sheath. Ions enter with negligible energy compared to $eV_0$, so
$u = \sqrt{-2e\phi/M}$ and

$$\frac{d^2\phi}{dx^2} = -\frac{e n_i}{\varepsilon_0}
= -\frac{J}{\varepsilon_0}\sqrt{\frac{M}{2e}}\,(-\phi)^{-1/2}$$

Multiply by $\phi'$ and integrate, using $\phi'(0)=0$:

$$\frac12(\phi')^2
= -\frac{J}{\varepsilon_0}\sqrt{\frac{M}{2e}}\int_0^{\phi}(-\phi')^{-1/2}d\phi'
= \frac{2J}{\varepsilon_0}\sqrt{\frac{M}{2e}}\,(-\phi)^{1/2}$$

Write $K \equiv (J/\varepsilon_0)\sqrt{M/2e}$, so $d(-\phi)/dx = 2\sqrt{K}(-\phi)^{1/4}$.
Separate and integrate from $0$ to $s$:

$$\frac{4}{3}(-\phi)^{3/4}\Big|_0^{s} = 2\sqrt{K}\,s
\qquad\Longrightarrow\qquad
\frac43 V_0^{3/4} = 2\sqrt{K}\,s$$

$$\boxed{\;J = \frac{4}{9}\varepsilon_0\left(\frac{2e}{M}\right)^{1/2}
\frac{V_0^{3/2}}{s^{2}}\;}$$

the **Child–Langmuir law**. The $V^{3/2}$ scaling is the signature: space charge
limits the current, and the limit stiffens with voltage more slowly than Ohm's law
would.

### 3c. Matched form — sheath width in terms of Debye length

The ion current crossing the sheath is set upstream by the Bohm flux,
$J = e n_s u_B$. Substituting:

$$s^2 = \frac{4}{9}\cdot\frac{\varepsilon_0 V_0^{3/2}\sqrt2}{e n_s \sqrt{T_e}}
= \frac{4\sqrt2}{9}\,\lambda_{Ds}^2\left(\frac{V_0}{T_e}\right)^{3/2}$$

$$\boxed{\;s = \frac{\sqrt2}{3}\,\lambda_{Ds}\left(\frac{2V_0}{T_e}\right)^{3/4}\;}
\qquad\text{(Lieberman \& Lichtenberg eq. 6.3.13)}$$

So a sheath is "many Debye lengths thick" with the excess set by the voltage
overdrive to the $3/4$ power. At $V_0/T_e \sim 20$ you get $s\sim 7\lambda_{Ds}$;
at $V_0/T_e\sim 200$, $s\sim40\lambda_{Ds}$.

### 3d. Numbers for Case 1 **[predict]**

$\lambda_{Ds} = \lambda_D/\sqrt{0.61} = 1.95/0.781 = 2.49$ mm. Taking a
representative time-averaged sheath voltage $\bar V \approx 180$ V (roughly
$0.4V_0$ for a symmetric discharge):

$$s \approx \frac{\sqrt2}{3}(2.49\ \mathrm{mm})\left(\frac{360}{9.36}\right)^{3/4}
\approx 1.8\ \mathrm{cm}$$

That is **over a quarter of the 6.7 cm gap on each side**. Case 1 is a
sheath-dominated, low-density, high-voltage discharge — which is exactly why
Turner chose it as a hard test. Expect the time-averaged density profile to be
broad and rounded, not a flat-top with thin edges. If your PIC gives thin sheaths
here, something is wrong.

### 3e. Collisional correction

The above assumed ions cross ballistically. That holds when the ion mean free path
$\lambda_i \gtrsim s$. At 30 mTorr in He, and worse in Ar where resonant charge
exchange gives a very large cross section, this can fail. In the mobility-limited
(collisional) limit the scaling softens to $J \propto V^{2}/s^{3}$
(the Warren/collisional Child law). Practically: ions arrive with less than the
full sheath voltage, and the IEDF develops a large low-energy population — see §6.

---

## 4. RF sheaths and the DC self-bias

### Why a self-bias appears

Drive one electrode through a **blocking capacitor**. No DC current can flow
through a capacitor, so averaged over an RF cycle the net charge collected by that
electrode must be exactly zero:

$$\langle I_e \rangle + \langle I_i \rangle = 0$$

Ions are heavy; their flux is essentially constant over the cycle at the Bohm
value. Electrons, being light, reach the electrode only during the brief window
each cycle when the sheath collapses. If the electrode sat at zero DC potential,
that window would deliver far more electron charge than the ions supply. So charge
accumulates on the blocking capacitor, pushing the electrode DC potential negative,
until the collapse window is short enough that the two fluxes balance. That
steady-state offset is the **DC self-bias** $V_{dc}$.

In a *symmetric* discharge both sheaths behave identically and $V_{dc}=0$. A
self-bias requires an asymmetry.

### The area-ratio scaling

Take a capacitive discharge with a powered electrode of area $A_p$ and a grounded
one of area $A_g$, sheath thicknesses $s_p$, $s_g$, sustained sheath voltages
$V_p$, $V_g$. Two facts:

**(i) The sheaths are in series and carry the same RF current.** Treating each as a
capacitor $C = \varepsilon_0 A/s$:

$$I = \omega C_p V_p = \omega C_g V_g
\qquad\Longrightarrow\qquad
\frac{V_p}{V_g} = \frac{C_g}{C_p} = \frac{A_g/s_g}{A_p/s_p}$$

**(ii) Both sheaths obey Child's law at the same bulk density**, so the same ion
current density flows into each, and $J \propto V^{3/2}/s^2$ gives
$s \propto V^{3/4}$.

Combine:

$$\frac{V_p}{V_g} = \frac{A_g}{A_p}\cdot\frac{s_p}{s_g}
= \frac{A_g}{A_p}\left(\frac{V_p}{V_g}\right)^{3/4}
\qquad\Longrightarrow\qquad
\left(\frac{V_p}{V_g}\right)^{1/4} = \frac{A_g}{A_p}$$

$$\boxed{\;\frac{V_p}{V_g} = \left(\frac{A_g}{A_p}\right)^{q},\quad q=4\ \text{(ideal)}\;}$$

Measured exponents are **[lit]** typically $q\approx1$–$2.5$, not 4, because real
sheaths are collisional, the density is not identical at both electrodes, and
"grounded area" is ill-defined in a real chamber. The *direction* is robust though:
**the smaller electrode takes the larger sheath voltage and goes DC-negative.**
That is why in a real etcher the wafer sits on the small powered electrode — you
want the high-energy directional ions to land on the wafer, not on the chamber
walls.

For a sinusoidal drive with the powered electrode much smaller than ground,
$V_{dc} \to -V_0$ in the limit; for comparable areas it is a fraction of that.

### **[flag]** This project's solver cannot show a geometric self-bias

A 1D Cartesian model has no area ratio — both electrodes have the same (notional)
area by construction. So a 1D PIC of a symmetric gap gives $V_{dc}=0$ up to
statistical noise, and you should expect that. Do not "find" a self-bias in
Phase 1 and report it.

If you want a self-bias in a 1D code, the honest routes are:
- **Circuit asymmetry:** add a series blocking capacitor and an explicitly
  prescribed effective area ratio in the external circuit equation.
- **Electrical asymmetry effect (EAE):** drive with $V_0[\sin(\omega t) +
  \sin(2\omega t + \theta)]$. A waveform whose positive and negative excursions are
  unequal produces a self-bias even in perfectly symmetric geometry, tunable by
  $\theta$. This is a real and currently-fashionable technique and would be a
  strong bonus result.
- **Material asymmetry:** different $\gamma_{SEE}$ on the two electrodes.

Pick one deliberately, or state in the README that the reactor model is symmetric
and the self-bias discussion is analytic.

---

## 5. Why ICP decouples ion flux from ion energy and CCP does not

### CCP: one knob

In a CCP, the same RF voltage does two jobs at once:

1. **Sustains the plasma.** Power is delivered to electrons mainly by (a) ohmic
   heating in the bulk and (b) stochastic/sheath heating — the oscillating sheath
   edge acts as a moving wall that Fermi-accelerates electrons. Both scale with the
   RF voltage/current amplitude.
2. **Accelerates ions into the wafer.** The sheath voltage is set by the same
   amplitude.

Raise $V_0$ to get more ion flux and the ion energy goes up too. Lower $V_0$ to
soften the ions and the density collapses. Roughly, $n_e$ and $\langle E_i\rangle$
both increase monotonically with $V_0$ — **one knob, two outputs.** For etch, this
is a real process limitation: you cannot independently dial "lots of ions" and
"gentle ions", which is precisely what selective or damage-free etching needs.

### ICP: two knobs

An ICP couples power **inductively**. A coil outside a dielectric window drives an
azimuthal RF current; the resulting time-varying $B_z$ induces an azimuthal
$E_\theta$ inside the plasma, which drives electron current and deposits power. The
coupling is transformer-like: the plasma is the one-turn secondary.

The crucial point is that this power transfer **does not require a large sheath
voltage**. The induced $E_\theta$ is parallel to the wafer surface and sits in the
bulk; the wafer sheath is only the small floating sheath ($\sim5T_e$, i.e. tens of
volts). So:

- **Source (coil) power** $\Rightarrow$ ionisation rate $\Rightarrow$ $n_e$ and
  hence **ion flux**. Typical ICP densities $10^{17}$–$10^{18}$ m⁻³, one to two
  orders above a CCP at the same pressure **[lit]**.
- **A separate RF bias supply on the wafer chuck** $\Rightarrow$ sheath voltage
  $\Rightarrow$ **ion energy**, with only weak effect on $n_e$ because the bias
  power is small compared with the source power.

Two knobs. This is why every high-aspect-ratio dielectric etcher and every ALE tool
is an ICP or a dual/triple-frequency CCP rather than a single-frequency CCP.

### The honest caveats, because an interviewer will push here

- Decoupling is **not perfect**. Bias power still adds to total absorbed power, so
  at high bias the density does rise. The usual rule of thumb is that decoupling
  holds while $P_{bias} \ll P_{source}$.
- At low coil currents an ICP runs in **E-mode**: the coil couples capacitively
  through its own voltage, density is low, and you have effectively a CCP. The
  E→H transition is a sharp, hysteretic jump.
- **Multi-frequency CCPs** recover much of the same freedom: a high frequency
  (e.g. 60–100 MHz) dominantly sets density because sheath heating efficiency rises
  with $\omega$, and a low frequency (2 MHz) dominantly sets ion energy because the
  sheath voltage at low $\omega$ divides preferentially there. Not fully
  independent, but two useful knobs.

---

## 6. The ion energy distribution at the wafer (IEDF)

This is the output of Phase 1 that Phase 3 consumes, so it matters that we get the
*mechanism* right and not just the picture.

### The controlling parameter is transit time, not pressure

An ion crossing the sheath takes

$$\tau_i \approx \frac{3s}{\bar v} = 3s\sqrt{\frac{M}{2e\bar V}}$$

(the factor 3 comes from integrating across the Child-law field profile;
$2s/v_f$ from a uniform field is the cruder version). Compare with the RF period
$\tau_{RF}=1/f$:

- **$\tau_i \ll \tau_{RF}$** — the ion crosses within one phase of the RF. It
  samples the *instantaneous* sheath voltage, which runs over the full range from
  near-zero to $V_{max}$. Result: a **broad, strongly bimodal** IEDF with peaks at
  the extremes (the distribution piles up at the turning points of a sinusoid).
- **$\tau_i \gg \tau_{RF}$** — the ion averages over many cycles and arrives with
  the *time-averaged* sheath energy. Result: a **single narrow peak** near
  $e\bar V$, with a small residual splitting.

The splitting width in the second limit scales as

$$\Delta E \ \sim\ C\,\frac{eV_1}{\omega s}\sqrt{\frac{2e\bar V}{M}}$$

with $V_1$ the modulation amplitude of the sheath voltage and $C$ of order unity
($C\approx2$ from a uniform-field estimate; Lieberman & Lichtenberg §11.2 and the
original Benoit-Cattin & Bernard (1968) treatment give the Child-law prefactor —
**check the source before quoting a prefactor in an interview**, the scalings below
are the defensible part):

$$\Delta E \ \propto\ \frac{1}{f}\,,\qquad
\Delta E \ \propto\ M^{-1/2}\,,\qquad
\Delta E \ \propto\ \frac{1}{s}$$

Lighter ions and lower frequency $\Rightarrow$ more splitting. **Pressure does not
appear in this expression at all.**

### **[predict]** What Case 1 should actually give

$\bar V \approx 180$ V, $s\approx1.8$ cm, He⁺:
$v_f=\sqrt{2e\bar V/M}=9.3\times10^{4}$ m/s, so $\tau_i \approx 3(0.018)/9.3\times10^4
= 5.8\times10^{-7}$ s $= 580$ ns, versus $\tau_{RF}=73.7$ ns.

$$\tau_i/\tau_{RF} \approx 8$$

So **Case 1 is deep in the time-averaged regime and should give a narrow,
single-peaked high-energy IEDF**, not a textbook bimodal one. For a representative
Ar case (10 mTorr, $n=10^{16}$ m⁻³, $T_e=3$ eV, $\bar V=100$ V) the same estimate
gives $s\approx1.8$ mm, $\tau_i\approx250$ ns, ratio $\approx3.4$ — modest
splitting, visible but not dramatic.

### **[flag]** The plan's "bimodal collapses with pressure" claim needs restating

What pressure actually changes is **collisionality in the sheath**, chiefly
**resonant charge exchange** $\mathrm{Ar^+ + Ar \to Ar + Ar^+}$, which has a huge
cross section ($\sim5\times10^{-19}$ m² at a few eV **[lit]**) and effectively
stops an ion dead partway across the sheath. Each CX event converts a
would-be-full-energy ion into a slow one that then falls through only the remaining
potential. So as pressure rises:

- the high-energy peak **loses population** and eventually disappears;
- a **large low-energy population** grows, with structure (discrete secondary peaks
  in the collisionless-ish regime, smearing into a monotonic low-energy ramp when
  $s/\lambda_{CX}\gg1$);
- the mean ion energy drops well below $e\bar V$.

Both statements — "bimodal splitting is set by $\omega\tau_i$" and "the low-energy
peak at high pressure is charge-exchange" — are true, but they are *different
mechanisms*. Writing "bimodal at low pressure collapses to single peak at high
pressure" as though pressure controlled the splitting is the kind of sentence that
falls apart under one follow-up question. Report them separately:
$\omega\tau_i$ for the splitting, $s/\lambda_{CX}$ for the low-energy structure.

**If you want a clearly bimodal case to show**, force $\omega\tau_i<1$
deliberately: lower the frequency (a 2 MHz run), raise the density (thinner
sheath), or lower the voltage. Do it as an explicit sweep in $\omega\tau_i$ and
label the axis $\omega\tau_i$, not pressure. That is a much stronger figure anyway
because it shows you know what the controlling parameter is.

### IADF

Ions enter the sheath with a thermal spread $\sim T_i$ transverse, then get
accelerated to $eV_{sh}$ normal. The angular half-width is roughly

$$\theta_{1/2} \sim \sqrt{\frac{T_i}{eV_{sh}}}$$

For $T_i\approx0.05$ eV and $V_{sh}=100$ V that is $\sim0.02$ rad $\approx1.3°$ —
near-perfectly collimated, which is the entire basis of anisotropic etching.
Sheath collisions broaden this substantially, and the broadening is what Phase 3
will feel as sidewall attack. Extract the *joint* $f(E,\theta)$, not two separate
1D histograms — the correlation matters (low-energy ions from CX are also the wide-
angle ones).

---

## 7. From ions to removed material

### Sputter yield

Physical sputtering near threshold is well described **[lit]** by the
Steinbrüchel (1985) form used in Phase 2:

$$Y(E) = A\left(\sqrt{E}-\sqrt{E_{th}}\right)$$

#### Threshold energy — the semi-empirical formulas overshoot for Si

An earlier draft of this document quoted only the Bohdansky estimate, for
$M_1/M_2>0.3$: $E_{th}\approx 8U_s(M_1/M_2)^{2/5}$, giving
$8(4.7)(1.4286)^{0.4}\approx 43$ eV for Ar on Si at $U_s=4.7$ eV. Yamamura's
form ($E_{th}/U_s = 6.7/\gamma$ for $M_1\ge M_2$, with
$\gamma = 4M_1M_2/(M_1{+}M_2)^2 = 0.969$) gives 33 eV, and the
$7.0A^{-0.54}+0.15A^{1.12}$ fit gives 40 eV. **All of these are roughly twice
the measured value.** Experiment:

| Source | Method | $E_{th}$ (Ar⁺→Si) |
|---|---|---|
| Wolsky & Zdanuk, via Kuschel et al. (2025) | extrapolation of 34–800 eV yields on c-Si | 15–20 eV |
| Wu et al., *J. Appl. Phys.* **106**, 054902 (2009) | in-situ weight loss + QCM, near-threshold | 18–20 eV |
| Kuschel et al., arXiv:2509.01171 | ICP etch rate, tailored-waveform bias, a-Si, linear fit 30–45 eV | 23 eV |
| Barone & Graves | MD, physical sputtering of (fluorinated) Si | 20 eV |
| COMSOL default, Ar on polysilicon | commercial-tool default parameter | 50 eV |

*(For contrast: SiO₂ is a genuinely higher-threshold material — 36–37 eV from
three independent monoenergetic-IEDF studies. Don't mix the two up.)*

So the honest range for **Ar⁺ → Si physical sputtering is ~15–25 eV**, with the
outliers at 50 eV coming from polysilicon and from tool defaults rather than from
near-threshold measurements.

**Why the formulas are wrong here**, which is the part worth being able to say out
loud:

1. $U_s=4.7$ eV is the **bulk cohesive energy**. The quantity the formula wants is
   a surface binding energy for a bombardment-amorphised Si(100) surface with
   dangling bonds, which is lower. Back-solving the Bohdansky prefactor from
   $E_{th}=20$ eV gives $U_s\approx2.2$ eV. Note that the competing rule of thumb
   $E_{th}=4U_s$ gives 18.8 eV with the *same* $U_s=4.7$ eV — the two rules
   disagree by a factor of 2 about which $U_s$ to use, which is the tell that
   they're being extrapolated outside the (mostly metallic) data they were fitted
   to. Si is covalent; these fits are not calibrated for it.
2. The **$E_{th}$ in $Y=A(\sqrt E-\sqrt{E_{th}})$ is a fit parameter, not a
   measured onset.** It's the x-intercept of a $\sqrt E$ extrapolation, which
   systematically lands below the true physical onset. Steinbrüchel's own framing
   is that $n=0.5$ is a universal *representation* of the data, with $A$ and
   $E_{th}$ as characteristic constants — not that $E_{th}$ is a physical
   threshold energy.
3. Real etched Si is already **amorphised** by prior bombardment, and a-Si sputters
   at lower energy than c-Si (23 eV vs higher in the ICP study above).

**[predict]** Revised Phase 2 acceptance band: **$E_{th}$ anywhere in 15–45 eV is
consistent with the literature** and should be reported with the spread, not as a
single number. Outside 10–60 eV, suspect the MD setup (thermostat bleeding energy,
slab too thin, ZBL splice mis-joined) rather than the physics.

One more caveat specific to this project: the MD threshold is set by the
**interatomic potential's** surface binding energy, not by nature's. Stillinger–Weber
and Tersoff have different Si cohesive energies (≈4.34 and ≈4.63 eV/atom), so the
two will not give the same $E_{th}$. Report which potential produced the number.

#### A better acceptance criterion than $E_{th}$

Because $E_{th}$ is a fit parameter with a factor-of-two literature spread, it is a
weak validation target. **Absolute yield at well-measured energies is much
stronger** — these are the numbers to check Phase 2 against first:

| Energy | Experimental $Y$ (Ar⁺→Si, normal incidence) **[lit]** |
|---|---|
| 100 eV | ≈ 0.07 |
| 500 eV | ≈ 0.6–0.7 |
| 1 keV | ≈ 0.93 |

For reference, an MD study using a DFT-derived Ar–Si repulsive potential obtained
0.59 at 500 eV and 0.88 at 1 keV, the latter about 23% below the experimental 0.93 —
which is a realistic expectation for how close SW/Tersoff + ZBL will land. If our
Phase 2 yields sit within a few tens of percent of the table above, that is a
genuine result. If they are off by 10×, the setup is wrong.

### Angular dependence

$Y(\theta)$ rises from normal incidence, peaks around $50$–$70°$, then falls at
grazing incidence as ions reflect instead of depositing energy **[lit]**. The rise
is because oblique impacts deposit their collision cascade closer to the surface;
the fall is reflection. **This non-monotonic shape is the direct cause of
microtrenching in Phase 3**: ions specularly reflected off a sloped sidewall land
at the base corner, arriving at an angle near the yield maximum, and cut a groove
there. When Phase 3 shows microtrenching, prove it by turning specular reflection
off and confirming the trench disappears.

---

## 8. Atomic layer etching and the synergy window

### The idea

ALE splits etching into two self-limiting half-cycles:

- **A — modification.** Dose the surface with a reactant (e.g. Cl₂ or Cl radicals)
  that chemisorbs and forms a thin, weakly-bound modified layer (SiCl$_x$). The
  reaction stops when the surface is saturated: self-limiting.
- **B — removal.** Bombard with low-energy ions (Ar⁺) at an energy chosen to
  desorb the modified layer but *not* to sputter the underlying bulk. When the
  modified layer is gone, removal stops: self-limiting.

Etch per cycle (EPC) is then set by the thickness of the modified layer, not by
time — which is what buys you atomic-scale control and excellent uniformity.

### The synergy window, defined

Physically, the window is the ion-energy interval

$$\boxed{\;E_{th}^{\text{modified}} \;<\; E_{ion} \;<\; E_{th}^{\text{bulk}}\;}$$

Below the lower bound nothing is removed even after modification. Above the upper
bound you sputter unmodified material, the process stops being self-limiting, and
you have simply reinvented continuous sputter etching with extra steps.

**[lit]** For Cl/Si the window is roughly 10–40 eV. **The whole point of the
optional Phase 2 Cl-passivation run is to produce our own two numbers for those
bounds instead of quoting someone else's.** If that run happens, the README claim
becomes "we measured a window of X–Y eV," which is genuinely differentiating.

### Synergy, quantified

Kanarik et al. (*J. Vac. Sci. Technol. A* **33**, 020802 (2015)) define the ALE
synergy as

$$S = \frac{\mathrm{EPC}_{AB} - \mathrm{EPC}_{A} - \mathrm{EPC}_{B}}{\mathrm{EPC}_{AB}}\times100\%$$

where $\mathrm{EPC}_A$ is what the modification step alone removes per cycle,
$\mathrm{EPC}_B$ what the ion step alone removes per cycle, and
$\mathrm{EPC}_{AB}$ the full two-step cycle. $S\to100\%$ means neither step does
anything on its own and all removal comes from their combination — ideal ALE.
Low $S$ means one of the steps is doing continuous etching and you have lost
self-limitation.

This is a clean, quotable definition and worth having cold: it is the standard
figure of merit in the field and Applied Materials interviewers in the etch org
will recognise it.

---

## 9. **Corrections to the project plan** — read this before coding

Four mismatches between the kickoff plan and the benchmark paper. The first is
fundamental.

### 9.1 The Turner benchmark is **helium**, not argon

From §II of the paper: the gap is filled with **helium**, the ion mass is
$6.67\times10^{-27}$ kg (He⁺), electron cross sections are the **Biagi v7.1 He**
set (elastic + two excitations + ionization), and ion–neutral scattering uses
**Phelps' He⁺–He** analytic model split into isotropic and backscattering
components.

You cannot validate an argon simulation against a helium benchmark. The ion mass
differs by a factor of 10, the ionization threshold differs (24.59 eV vs 15.76 eV),
and Ar has strong resonant charge exchange that He does not.

**Recommendation — run both gases, deliberately:**

1. **Validation run: helium**, exactly Case 1 parameters. This is the number the
   README's validation claim rests on. The cross sections and the reference profile
   are already in `data/` (see §10).
2. **Production runs: argon**, for the pressure/voltage sweep and the IEDF/IADF
   that Phase 3 consumes, because argon is what an etch reactor actually uses.

This costs one extra cross-section set and a `gas` config switch, and it makes the
validation claim honest. The alternative — running Ar and eyeballing it against a
He curve — is not validation and would not survive a technical interview.

### 9.2 The benchmark quantity is **ion** density, not electron density

The paper says explicitly that they focus on ion density because it has an
unambiguous definition in any simulation procedure and is highly sensitive to both
numerical effects and implementation details. All comparison figures (7–10) are
$n_i(x)$. Overlay ion density.

### 9.3 Secondary electron emission must be **off** for the validation run

The benchmark assumes all charged-particle fluxes at the electrodes are completely
absorbed with **no secondary emission**. The plan's $\gamma\approx0.07$ would
change the discharge and break the comparison. Implement $\gamma$ as a config
parameter, set it to 0 for the Turner run, and turn it on for the argon production
runs (where 0.05–0.1 for Ar⁺ on a metal electrode is reasonable **[lit]**).

### 9.4 Run length: 1280 RF cycles, not 300–500

Case 1 is 512 000 steps at 400 steps/cycle = 1280 cycles, averaged over the final
32. The higher-pressure cases need 5120 and 15 360 cycles. Budget for this: Case 1
is the cheap one.

### 9.5 State the agreement metric before you measure it

"Agrees to within X%" is meaningless unless X is defined. Turner's own framing:
implementation uncertainty between the five codes is $<0.5\%$; residual numerical
error in the base-case parameters is a few percent at the mid-plane. So define, in
advance, in `PROGRESS.md`:

- **Primary:** relative error in mid-plane $n_i$ vs. $1.40\times10^{14}$ m⁻³.
- **Secondary:** relative error in mid-plane $k_BT_e$ vs. 9.36 eV.
- **Profile:** normalised RMS error
  $\sqrt{\langle(n_i^{ours}-n_i^{ref})^2\rangle}/\max(n_i^{ref})$ over the full
  gap, reported separately for the bulk and for the sheath regions (the sheath
  edge has steep gradients and will dominate any global metric — the paper makes
  the same point about its own refinement study).

A realistic target for a first independent implementation is **a few percent at
the mid-plane**. Matching to 1% would be excellent. If you land at 20%, that is a
bug, not a tolerance.

---

## 10. Cross-section data: what is here and where it came from

### Already in the repo

`data/cross_sections/turner_benchmark_he/` — the exact set the benchmark
prescribes:

| File | Process | Threshold |
|---|---|---|
| `Elastic_He.csv` | e⁻ + He elastic momentum transfer | — |
| `Excitation1_He.csv` | e⁻ + He → He* | 19.82 eV |
| `Excitation2_He.csv` | e⁻ + He → He** | 20.61 eV |
| `Ionization_He.csv` | e⁻ + He → He⁺ + 2e⁻ | 24.59 eV |
| `Isotropic_He.csv` | He⁺ + He isotropic scattering | — (CM energy) |
| `Backscattering_He.csv` | He⁺ + He backward scattering | — (CM energy) |

Format: semicolon-delimited, `energy_eV;sigma_m2`. Electron energies are in the lab
frame; **ion energies are in the centre-of-mass frame** — a classic bug source, get
this right in the loader. The benchmark prescribes **linear interpolation** between
tabulated points, and clamping to the last value above the maximum tabulated
energy.

`data/benchmark/turner2013_case1_reference.csv` — space-delimited, column 0 is
$x$ in metres on the $L/128$ grid, **column 1 is electron density**, **column 4 is
ion density** (m⁻³). Its mid-plane ion density is $1.4046\times10^{14}$ m⁻³, which
reproduces the paper's Table III value of $0.140\times10^{15}$ m⁻³ to three digits —
that is the check that confirms it is a faithful Case 1 reference.

### **[flag]** Provenance caveat — this is a secondary source

I obtained these from the public reimplementation
`github.com/lase-unb/ccp-benchmark`, which carries them with the LXCat/Phelps
attribution reproduced in `LICENSE.txt` alongside them. That is a *secondary*
source. Before the repo goes on a resume:

1. Download the official supplementary material attached to the Phys. Plasmas paper
   (doi:10.1063/1.4775084) and diff it against these files.
2. Register on lxcat.net (free) and pull Biagi v7.1 He yourself.
3. Keep the LXCat attribution block. LXCat's terms require citing the database and
   the originating contributor — Biagi for the electron set, Phelps for the ion
   set. That attribution is already preserved in `LICENSE.txt`; do not strip it.

### Still needed: the argon set

I could not reach LXCat from the environment these files were generated in, so the
argon cross sections are **not** here yet. For the production runs you need, from
LXCat (Biagi v7.1 or the Phelps Ar set — pick one and use it consistently, do not
mix compilations):

- e⁻ + Ar: elastic momentum transfer, effective excitation (lumped is fine at this
  level), ionization ($E_{th}=15.76$ eV)
- Ar⁺ + Ar: elastic/isotropic **and resonant charge exchange** — the CX channel is
  the one that produces the low-energy IEDF structure in §6, so it is not optional

A note on consistency that is worth knowing: momentum-transfer cross sections are
only consistent with transport data when used with the scattering model they were
derived for. Biagi's sets assume isotropic scattering in the CM frame, and the
benchmark makes isotropic scattering a *requirement* for that reason. If you later
swap in an anisotropic scattering model, you must swap the cross sections too.

---

## 11. References

1. M. M. Turner, A. Derzsi, Z. Donkó, D. Eremin, S. J. Kelly, T. Lafleur,
   T. Mussenbrock, "Simulation benchmarks for low-pressure plasmas: Capacitive
   discharges," *Phys. Plasmas* **20**, 013507 (2013). doi:10.1063/1.4775084;
   preprint arXiv:1211.5246.
2. M. A. Lieberman & A. J. Lichtenberg, *Principles of Plasma Discharges and
   Materials Processing*, 2nd ed., Wiley (2005). Ch. 6 (sheaths), Ch. 11 (RF
   sheaths and IEDFs), Ch. 12 (ICP).
3. S. F. Biagi, cross-section compilation v7.1 (2004), via www.lxcat.net.
4. A. V. Phelps, *J. Appl. Phys.* **76**, 747 (1994); analytic compilation at
   jila.colorado.edu/~avp/.
5. C. K. Birdsall & A. B. Langdon, *Plasma Physics via Computer Simulation*,
   Adam Hilger (1991).
6. V. Vahedi & M. Surendra, "A Monte Carlo collision model for the particle-in-cell
   method," *Comput. Phys. Commun.* **87**, 179 (1995). — the null-collision method
   for Phase 1 step (5).
7. C. Steinbrüchel, "Universal energy dependence of physical and ion-enhanced
   chemical etch yields at low ion energy," *Appl. Phys. Lett.* **55**, 1960 (1989);
   J. Muri & Ch. Steinbrüchel, "Universal energy dependence of sputtering yields at
   low ion energy," *MRS Proc.* (1991).
7a. Wu et al., "Sputtering yields of Ru, Mo and Si under low energy Ar+
   bombardment," *J. Appl. Phys.* **106**, 054902 (2009). — near-threshold Si data.
7b. Kuschel et al., "Multi-diagnostic characterization of inductively coupled
   discharges with tailored waveform substrate bias," arXiv:2509.01171 (2025).
   — a-Si threshold 23 eV; SiO2 37 eV; review of prior Si threshold values.
8. J. Bohdansky, "A universal relation for the sputtering yield of monatomic solids
   at normal ion incidence," *Nucl. Instrum. Methods B* **2**, 587 (1984).
9. K. J. Kanarik et al., "Overview of atomic layer etching in the semiconductor
   industry," *J. Vac. Sci. Technol. A* **33**, 020802 (2015).
10. V. A. Godyak, R. B. Piejak, B. M. Alexandrovich, *Plasma Sources Sci. Technol.*
    **1**, 36 (1992). — the experiment the benchmark conditions were modelled on.
