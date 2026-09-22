# Study notes — plain-English explainer

This file exists so the person who built this repo can actually defend it in
an interview. Opus: after each phase, write a real, specific answer here —
tied to what this repo's code actually does and what its actual results
showed, not a generic textbook answer. If a question can't yet be answered
because that phase isn't built, leave it explicitly marked TODO rather than
answering from general knowledge.

1. **Why does a sheath form at all, and what does the Bohm criterion require
   of ions entering it?**
   *(answer)*

2. **Why does an asymmetric CCP develop a negative DC self-bias?**
   *(answer)*

3. **Why does ICP decouple plasma density from ion energy, and why can't
   CCP do that?**
   *(answer)*

4. **Why is a low-pressure IEDF bimodal, and what collapses it at high
   pressure? What did YOUR simulation show at each pressure you swept?**
   *(answer — reference your actual figures)*

5. **What makes plasma etching anisotropic when wet etching isn't?**
   *(answer)*

6. **What is the ion–neutral synergy in ion-enhanced etching?**
   *(answer)*

7. **What makes ALE self-limiting, and what defines the energy window on
   both sides? What did your yield-threshold fit give for E_th?**
   *(answer — reference your actual fitted numbers)*

8. **What causes ARDE (aspect-ratio-dependent etching)? Name two distinct
   mechanisms. Which one dominates in YOUR model?**
   *(answer)*

9. **What causes microtrenching and sidewall bowing? Did your feature-scale
   model reproduce them, and if so, from what mechanism in your code?**
   *(answer)*

10. **Why is ALD self-limiting — what are GPC, precursor saturation, and the
    ALD temperature window?** *(from coursework, not this project — but you
    should be fluent since the JD leads with deposition)*
    *(answer)*

11. **How does PECVD let you deposit at lower temperature than thermal CVD,
    and what do you trade away?**
    *(answer)*

12. **Why is step coverage better in ALD than in sputtering?**
    *(answer)*

---

## What I'd say if asked "walk me through this project" (60 seconds)

*(Opus: draft this once Phase 4 is done. It should be something the actual
author could say out loud without sounding like they're reciting the README.)*

## Honest weak points I should be ready to defend

*(Opus: list the real simplifications — e.g. "I used a fitted empirical yield
law instead of full reactive MD for the etch step because reactive MD at this
scale wasn't tractable in the time I had" — and a defensible answer for each.)*
