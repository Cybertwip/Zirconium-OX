# Zirconium-90 breakup, neutron beta decay and proton trap (educational simulation)

This project calculates conservation-law bookkeeping and the rest-mass energy deficit for the hypothetical breakup

`90Zr -> 16O + 32p + 42n`.

**This is not a demonstrated nuclear reaction.** A negative Q value gives only the minimum mass-energy balance. It does not predict attainable yield, incident-beam energy, production rate, power requirements, radiation shielding needs, or a feasible device. There is no known conversion of detached protons directly to radiation through the Standard Model, and proton decay has not been observed. The toy proton-capture percentage is a user-set assumption, *not* a computed capture probability.

## Run

```bash
python -m pip install -r requirements.txt
python zirconium_proton_sim.py
```

Keyboard: `B` animate hypothetical breakup, `M` cycle starting Zr-90 mass, `C` change assumed capture percentage, `T` change accelerated elapsed time (neutron beta decay), `Space` pause, `R` reset, `Esc` quit.

For exact numbers, including in an environment without Pygame:

```bash
python zirconium_proton_sim.py --headless --mass-mg 1 --capture-efficiency 0.75 --elapsed-seconds 3600
python zirconium_proton_sim.py --self-test
# Optional photon energy budget from *assumed* proton kinetic energy, NOT proton annihilation:
python zirconium_proton_sim.py --headless --proton-kinetic-mev 2 --stopping-photon-fraction 0.1
```

## Core formulas

Atomic masses in unified mass units `u` are tabulated for neutral Zr-90, O-16, H-1 and the neutron. Hydrogen **atomic** mass is used instead of bare proton mass to preserve the electron mass count in the first reaction.

```text
Q_breakup = [m(Zr-90) - m(O-16) - 32 m(H-1) - 42 m(n)] * 931.49410242 MeV/u
          = -656.279 MeV per nucleus
E_min     = -Q_breakup * (mass_mg * 1e-3 / 89.9046977) * N_A * MeV_to_J
S_p(Zr-90 -> Y-89 + p) = [m(Y-89)+m(H-1)-m(Zr-90)]c² = 8.353 MeV

free neutron beta decay: n -> p + e- + anti_nu_e
Q_beta    = [m(n)-m(H-1)] c² = 0.782347 MeV
f_decay   = 1 - exp(-t/878.4 s)
protons   = N_nuclei * [32 + 42*f_decay]
remaining_free_neutrons = N_nuclei * 42*(1-f_decay)
trapped_protons_idealized = protons * capture_efficiency

Penning trap: F = q(E + v x B)
Phi(r,z) = V0/(2d²)*(z²-r²/2)
omega_c = q B/m_p
omega_z² = q V0/(m_p d²)
omega_± = (omega_c ± sqrt(omega_c²-2omega_z²))/2
stability in linear idealized trap: omega_c² >= 2 omega_z²

QCD: F^a_{mu nu} = d_mu A^a_nu - d_nu A^a_mu + g_s f^abc A^b_mu A^c_nu
The non-Abelian f^abc describes gluon self-interactions. It does not give a nuclear transmutation channel.
V_confinement(r) ~ sigma*r: at r=1 fm, sigma=0.19 GeV² -> V~0.963 GeV.
This estimates the growth of the quark-antiquark potential, not the nuclear proton-removal energy.

Proton decay p -> e+ pi0 is *unobserved*. A lower limit on its PARTIAL lifetime is
  tau/B > 2.4e34 years at 90% confidence (Super-Kamiokande).
Optional kinetic-energy-to-photon bookkeeping (only kinetic energy, protons survive):
  E_prompt_kin = 32 N_nuclei K_p
  E_photons    = E_prompt_kin * user_assumed_photon_fraction
  E_nonphoton  = E_prompt_kin - E_photons
These are not predicted yields. The extra kinetic energy has to come from somewhere.

The code prints (N_initial_protons * observation_years / 2.4e34 years)
only as a *reference upper scale* from that limit. It is not the actual expected proton-decay count or a predicted rate.
```

## Example at 1 mg of isotopically pure Zr-90

- About 6.698e18 nuclei.
- Minimum breakup mass deficit: 656.279 MeV per nucleus = **7.043e8 J** = **195.64 kWh** per mg.
- 2.143e20 immediately free protons; 2.813e20 immediately free neutrons.
- If all free neutrons eventually decay, 74 free protons (including the 32 prompt ones) per original nucleus and 42 emitted electrons + 42 antineutrinos per nucleus.
- Decay rest-energy release (full neutron conversion): about **35.26 MJ** for 1 mg, not entirely electromagnetic radiation.
- At 75% *assumed* capture, total ultimately trapped protons: about 3.718e20 (if all protons were appropriately cooled and loaded).

Electrons themselves are stable in the Standard Model; neutron beta decay is a **weak** interaction. Quark confinement and string breaking create hadrons, not free quarks or a pathway for transmuting a bulk zirconium sample.

The animation displays idealized Penning motion (rescaled in time) and projectile icons; it does not solve actual post-reaction ion optics, decay trajectories or cooling. Trapping enormous quantities of like-charged particles would involve severe space-charge limits and is not demonstrated.

## Source data

- NIST Zr isotope masses: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Zr
- NIST O isotope masses: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii&ele=O
- NIST H isotope masses: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii&ele=H
- NIST Y isotope masses: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ascii=ascii&isotype=some
- Particle Data Group neutron properties: https://pdg.lbl.gov/2025/AtomicNuclearProperties/neutron.html
- Super-K proton-decay limit (p -> e+ pi0): https://escholarship.org/uc/item/9gz4c03p
- CERN BASE Penning trap overview: https://base.web.cern.ch/base-experiment

## GUI font-cache fix (2026-10-25)

The first GUI version raised `KeyError: 18` because its prebuilt fonts omitted
size 18. The updated `get_font(size, bold)` loads fonts on demand and caches
them by `(size, bold)`, preventing the same error for other sizes.

To replace an older script, overwrite `zirconium_proton_sim.py` with the new
version, or extract the updated ZIP and run:

```
python3 zirconium_proton_sim.py --self-test
python3 zirconium_proton_sim.py
```
