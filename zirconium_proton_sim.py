#!/usr/bin/env python3
"""Zr-90 -> O-16 + nucleons: conservation/energy and illustrative proton trap.

Usage:
  python zirconium_proton_sim.py --headless --mass-mg 1 --capture-efficiency 0.75 --elapsed-seconds 3600
  python zirconium_proton_sim.py                  # Pygame GUI
  python zirconium_proton_sim.py --self-test

IMPORTANT: This is NOT a validated nuclear reaction or reactor design.
The mass-deficit calculation is exact to the tabulated masses used, while the
conversion, neutron release and capture efficiency are counterfactual inputs.
Proton decay is not simulated: the quoted bound is a limit for ONE decay mode.
"""

import argparse
from dataclasses import asdict, dataclass
import json
import math
import random
from typing import Optional

# Atomic masses (u); neutral atoms keep the electron mass count balanced
# between 90Zr and 16O + 32(1H) + 42n.
M_ZR90_U = 89.9046977    # NIST
M_O16_U = 15.99491461957 # NIST
M_H1_U = 1.00782503223   # NIST (p + bound electron)
M_N_U = 1.00866491595    # PDG
M_Y89_U = 88.9058403     # NIST; for 90Zr -> 89Y + 1H
MEV_PER_U = 931.49410242
MEV_TO_J = 1.602176634e-13
N_AVOGADRO = 6.02214076e23
N_LIFETIME_S = 878.4
P_DECAY_PARTIAL_LIFETIME_BOUND_Y = 2.4e34  # p -> e+ pi0, 90% CL, Super-K
HBAR_C_GEV_FM = 0.1973269804
SIGMA_GEV2 = 0.19  # QCD string tension illustrative
CHARGE_C = 1.602176634e-19
PROTON_KG = 1.67262192369e-27

# Net energy required for the artificial breakup. Q < 0 means endothermic.
Q_BREAKUP_MEV = (M_ZR90_U - M_O16_U - 32*M_H1_U - 42*M_N_U) * MEV_PER_U
MIN_BREAKUP_MEV = -Q_BREAKUP_MEV
Q_BETA_MEV = (M_N_U - M_H1_U) * MEV_PER_U
FIRST_PROTON_SEPARATION_MEV = (M_Y89_U + M_H1_U - M_ZR90_U) * MEV_PER_U
STRING_ENERGY_1FM_GEV = SIGMA_GEV2 * (1.0/HBAR_C_GEV_FM)


@dataclass
class Inputs:
    mass_mg: float = 1.0
    capture_efficiency: float = 0.75  # assumption: cooled protons already admitted to trap
    elapsed_seconds: float = 3600.0
    observation_years: float = 1.0
    field_tesla: float = 5.0
    trap_voltage: float = 10.0
    trap_length_m: float = 0.01
    prompt_proton_kinetic_MeV: float = 0.0  # extra assumed kinetic energy per prompt proton
    stopping_photon_fraction: float = 0.0 # assumed fraction of that kinetic energy radiated as photons


# f^{abc} for eight Gell-Mann SU(3) generators, using standard signs.
# The commutator is [T^a,T^b] = i f^{abc} T^c.
SU3_INDEPENDENT_F = {
    (1,2,3):1.0,
    (1,4,7):0.5, (1,5,6):-0.5, (2,4,6):0.5, (2,5,7):0.5,
    (3,4,5):0.5, (3,6,7):-0.5,
    (4,5,8):math.sqrt(3)/2, (6,7,8):math.sqrt(3)/2
}


def f_abc(a: int, b: int, c: int) -> float:
    """Antisymmetric SU(3) constants; a,b,c label generators from 1 to 8."""
    if min(a,b,c)<1 or max(a,b,c)>8 or len({a,b,c})<3:
        return 0.0
    for triple,value in SU3_INDEPENDENT_F.items():
        if set((a,b,c)) == set(triple):
            inversions=sum((a,b,c)[i]>(a,b,c)[j] for i in range(3) for j in range(i+1,3))
            return value * (-1 if inversions % 2 else 1)
    return 0.0


def gluon_nonlinear_term(a: int, A_mu: dict, A_nu: dict, g_s: float=1.0) -> float:
    """Illustrative g_s sum_(b,c) f^(abc) A_mu^b A_nu^c (field units suppressed)."""
    return g_s*sum(f_abc(a,b,c)*A_mu.get(b,0.0)*A_nu.get(c,0.0)
                   for b in range(1,9) for c in range(1,9))


def calculate(params: Inputs) -> dict:
    if not math.isfinite(params.mass_mg) or params.mass_mg < 0:
        raise ValueError("mass_mg must be a non-negative finite number")
    if not 0 <= params.capture_efficiency <= 1:
        raise ValueError("capture_efficiency must be from 0 to 1")
    if params.elapsed_seconds < 0 or params.observation_years < 0:
        raise ValueError("time parameters must be nonnegative")
    if params.field_tesla <= 0 or params.trap_length_m <= 0 or params.trap_voltage <= 0:
        raise ValueError("trap field, length, and voltage must be positive")
    if params.prompt_proton_kinetic_MeV < 0 or not 0 <= params.stopping_photon_fraction <= 1:
        raise ValueError("kinetic energy must be nonnegative; photon fraction within [0,1]")

    # mg -> g; 90Zr molar mass approximately its relative atomic mass in g/mol.
    nuclei = (params.mass_mg * 1e-3 / M_ZR90_U) * N_AVOGADRO
    beta_frac = -math.expm1(-params.elapsed_seconds / N_LIFETIME_S)
    beta_per_nucleus = 42 * beta_frac
    protons_per_nucleus = 32 + beta_per_nucleus

    # Penning trap frequencies from Phi=V0/(2d^2) (z^2-r^2/2).
    omega_c = CHARGE_C * params.field_tesla / PROTON_KG
    omega_z2 = CHARGE_C * params.trap_voltage / (PROTON_KG * params.trap_length_m**2)
    discriminant = omega_c**2 - 2 * omega_z2
    trap_stable = discriminant >= 0
    omega_plus = (omega_c + math.sqrt(discriminant)) / 2 if trap_stable else None
    omega_minus = (omega_c - math.sqrt(discriminant)) / 2 if trap_stable else None

    return {
        "input": asdict(params),
        "nuclide": "zirconium-90 (isotopically pure)",
        "hypothetical_reaction": "90Zr -> 16O + 32p + 42n",
        "conservation": {"initial_protons": 40, "initial_neutrons": 50,
                         "final_initial_protons": 8+32, "final_initial_neutrons": 8+42},
        "initial_nuclei": nuclei,
        "oxygen16_nuclei": nuclei,
        "q_breakup_MeV_per_nucleus": Q_BREAKUP_MEV,
        "minimum_breakup_MeV_per_nucleus": MIN_BREAKUP_MEV,
        "minimum_breakup_energy_J": MIN_BREAKUP_MEV * nuclei * MEV_TO_J,
        "minimum_breakup_energy_kWh": MIN_BREAKUP_MEV * nuclei * MEV_TO_J / 3.6e6,
        # Kinetic energy is separate from Q and must be supplied by an unspecified mechanism.
        "assumed_prompt_proton_kinetic_energy_J": 32*nuclei*params.prompt_proton_kinetic_MeV*MEV_TO_J,
        "assumed_photon_energy_from_proton_stopping_J": 32*nuclei*params.prompt_proton_kinetic_MeV*MEV_TO_J*params.stopping_photon_fraction,
        "assumed_nonphotonic_stopping_energy_J": 32*nuclei*params.prompt_proton_kinetic_MeV*MEV_TO_J*(1-params.stopping_photon_fraction),
        "first_proton_separation_MeV": FIRST_PROTON_SEPARATION_MEV,
        "prompt_free_protons": 32*nuclei,
        "initial_free_neutrons": 42*nuclei,
        "neutron_decay_fraction": beta_frac,
        "free_neutrons_remaining": (42-beta_per_nucleus)*nuclei,
        "beta_decay_protons_produced": beta_per_nucleus*nuclei,
        "beta_electrons_produced": beta_per_nucleus*nuclei,
        "beta_antineutrinos_produced": beta_per_nucleus*nuclei,
        "available_free_protons": protons_per_nucleus*nuclei,
        "captured_protons_expected_idealized": protons_per_nucleus*nuclei*params.capture_efficiency,
        "untrapped_protons_expected_idealized": protons_per_nucleus*nuclei*(1-params.capture_efficiency),
        "beta_Q_MeV_per_neutron": Q_BETA_MEV,
        "beta_decay_total_energy_J": beta_per_nucleus*nuclei*Q_BETA_MEV*MEV_TO_J,
        "net_mass_energy_input_J_after_beta": (MIN_BREAKUP_MEV - beta_per_nucleus*Q_BETA_MEV)*nuclei*MEV_TO_J,
        "proton_decay_mode": "p -> e+ + pi0 (not seen; NOT neutron beta decay)",
        "upper_expected_proton_mode_events_in_observation_time": (
            40*nuclei*params.observation_years/P_DECAY_PARTIAL_LIFETIME_BOUND_Y),
        "qcd_linear_potential_GeV_at_1fm": STRING_ENERGY_1FM_GEV,
        "su3_example": {"f123": f_abc(1,2,3), "f213":f_abc(2,1,3),
                        "gluon_Fxy3_nonlinear_only_gs1_Ax1_Ay2_1":
                         gluon_nonlinear_term(3,{1:1},{2:1})},
        "trap": {"stable_linear_penning_model": trap_stable,
                 "cyclotron_Hz": omega_c/(2*math.pi),
                 "axial_Hz": math.sqrt(omega_z2)/(2*math.pi),
                 "modified_cyclotron_Hz": (omega_plus/(2*math.pi) if trap_stable else None),
                 "magnetron_Hz": (omega_minus/(2*math.pi) if trap_stable else None)},
        "notes": [
            "Breakup is counterfactual; mass balance is not an achievable reaction rate.",
            "Minimum breakup energy excludes inefficiencies, kinetic barriers and device costs.",
            "Neutron beta decay makes protons, electrons and antineutrinos; protons do not become light.",
            "Optional photons from slowing prompt protons consume assumed kinetic energy only, not proton rest mass.",
            "Capture efficiency is an assumption; a Penning trap cannot trap high-energy protons without prior cooling/acceptance.",
            "Super-K bound is a 90%-CL partial lifetime lower limit for one proton decay channel, not a measured decay rate."
        ]
    }


def self_test() -> None:
    z = calculate(Inputs(mass_mg=1.0, capture_efficiency=1.0, elapsed_seconds=0))
    assert abs(z["minimum_breakup_MeV_per_nucleus"]-656.27897) < 0.02
    assert z["conservation"]["initial_protons"] == z["conservation"]["final_initial_protons"]
    assert z["conservation"]["initial_neutrons"] == z["conservation"]["final_initial_neutrons"]
    assert abs(z["available_free_protons"] - 32*z["initial_nuclei"]) < 1e-8*32*z["initial_nuclei"]
    late = calculate(Inputs(mass_mg=1.0, capture_efficiency=0.5, elapsed_seconds=2e5))
    assert abs(late["available_free_protons"] / late["initial_nuclei"] - 74) < 1e-10
    assert abs(late["captured_protons_expected_idealized"] / late["initial_nuclei"] - 37) < 1e-10
    assert late["trap"]["stable_linear_penning_model"]
    assert math.isclose(late["minimum_breakup_energy_J"], z["minimum_breakup_energy_J"])
    assert late["beta_decay_total_energy_J"] > 0
    assert f_abc(1,2,3)==1 and f_abc(2,1,3)==-1
    assert f_abc(4,5,8)==math.sqrt(3)/2
    assert gluon_nonlinear_term(3,{1:1.0},{2:1.0})==1
    scenario=calculate(Inputs(mass_mg=1, prompt_proton_kinetic_MeV=1, stopping_photon_fraction=0.25))
    assert math.isclose(scenario["assumed_photon_energy_from_proton_stopping_J"],
                        0.25*scenario["assumed_prompt_proton_kinetic_energy_J"])
    print("Self-test passed: mass-energy, nucleon and charge counts, decay, and trap stability.")


# ---------------------------- optional Pygame GUI ---------------------------
W, H = 1440, 850
BG=(11,17,31)
WHITE=(233,242,255)
MUTED=(140,158,184)
BORDER=(50,69,94)
CYAN=(91,221,236)
BLUE=(91,149,232)
ORANGE=(255,164,100)
GREEN=(104,226,178)
PURPLE=(172,157,253)
RED=(250,104,113)


def gui() -> None:
    try:
        import pygame
    except ImportError as ex:
        raise SystemExit("Pygame is not installed. Run: python -m pip install pygame\n"
                         "Or use --headless for calculations without graphics.") from ex
    pygame.init()
    screen = pygame.display.set_mode((W, H), pygame.RESIZABLE)
    pygame.display.set_caption("Zirconium-90 nuclear breakup | Conservation & Penning trap toy")
    canvas = pygame.Surface((W,H))
    clock=pygame.time.Clock()
    # Fonts are loaded lazily. A label can request any size, including 18,
    # without depending on a hard-coded list of font dictionary keys.
    font_cache = {}

    def get_font(size, bold=False):
        key = (int(size), bool(bold))
        if key not in font_cache:
            try:
                font_cache[key] = pygame.font.SysFont("dejavusans", key[0], bold=key[1])
            except Exception:
                font_cache[key] = pygame.font.Font(None, key[0] + 5)
        return font_cache[key]

    def t(s, x, y, size=16, color=WHITE, bold=False):
        surf=get_font(size, bold).render(str(s), True, color)
        canvas.blit(surf,(x,y))
        return surf.get_width()

    def panel(rect,fill=(20,30,48),border=BORDER,radius=17):
        pygame.draw.rect(canvas,fill,rect,border_radius=radius)
        pygame.draw.rect(canvas,border,rect,1,border_radius=radius)

    random.seed(7123)
    stars=[(random.randrange(W),random.randrange(H),random.choice([1,1,2])) for _ in range(80)]
    zr_dots=[]
    for i in range(90):
        angle=i*2.3999632
        rad=15+100*math.sqrt((i+.5)/90)
        zr_dots.append((math.cos(angle)*rad,math.sin(angle)*rad,RED if i<40 else BLUE))
    o_dots=[]
    for i in range(16):
        angle=i*2.3999632
        rad=7+43*math.sqrt((i+.5)/16)
        o_dots.append((math.cos(angle)*rad,math.sin(angle)*rad,RED if i<8 else BLUE))

    mass_levels=[0.000001,0.00001,0.0001,0.001,0.01,0.1,1.0,10.0,100.0]
    idx=6
    capture=0.75
    elapsed=0.0
    accelerated=120.0
    paused=False
    particles=[]
    demo_started=False
    rng=random.Random(49)
    seconds_visible=0.0
    
    def spawn():
        nonlocal particles, demo_started, elapsed, seconds_visible
        particles=[]
        elapsed=0.0
        seconds_visible=0.0
        demo_started=True
        for i in range(32):
            particles.append({"kind":"p","age":0.0,"start":0.1+i*.06,
                              "vy":rng.uniform(-90,90),"capture":rng.random()<capture,
                              "phase":rng.uniform(0,math.tau)})
        for i in range(42):
            particles.append({"kind":"n","age":0.0,"start":0.1+i*.06,
                              "vy":rng.uniform(-125,125),"capture":rng.random()<capture,
                              "phase":rng.uniform(0,math.tau),
                              "decay_time":rng.expovariate(1/N_LIFETIME_S)})

    def update(dt):
        nonlocal elapsed, seconds_visible
        if paused:
            return
        seconds_visible+=dt
        if demo_started:
            elapsed+=dt*accelerated
        for p in particles:
            p["age"]=max(0.0,seconds_visible-p["start"])

    def draw():
        canvas.fill(BG)
        for sx,sy,r in stars:
            pygame.draw.circle(canvas,(29,44,66),(sx,sy),r)
        t("ZIRCONIUM  /  OXYGEN  /  PROTON ACCOUNTING",32,18,27,WHITE,True)
        t("Nuclear mass deficit  |  neutron beta decay  |  Penning trap",34,56,15,MUTED)
        panel(pygame.Rect(28,86,1384,40),(54,37,36),(115,70,62),9)
        t("CONCEPTUAL MODEL ONLY: conserving nucleons is not a proven transmutation or proton-decay mechanism.",45,97,16,(255,207,170),True)
        panel(pygame.Rect(28,143,955,466))
        panel(pygame.Rect(998,143,414,466))
        panel(pygame.Rect(28,623,1384,199))
        t("PARTICLE DIAGRAM: one nucleus (not to scale)",48,161,17,WHITE,True)
        # Atom silhouette
        zx,zy=190,355
        pygame.draw.circle(canvas,(40,58,87),(zx,zy),119,2)
        for dx,dy,c in zr_dots:
            pygame.draw.circle(canvas,c,(int(zx+dx),int(zy+dy)),6)
        t("Zr-90",145,485,22,WHITE,True)
        t("40 p + 50 n",125,514,14,MUTED)
        # reaction arrow
        pygame.draw.line(canvas,CYAN,(325,356),(386,356),3)
        pygame.draw.polygon(canvas,CYAN,[(392,356),(380,348),(380,364)])
        t("Q < 0",330,331,15,ORANGE)
        ox,oy=464,355
        pygame.draw.circle(canvas,(35,59,68),(ox,oy),57,2)
        for dx,dy,c in o_dots:
            pygame.draw.circle(canvas,c,(int(ox+dx),int(oy+dy)),8)
        t("O-16",436,485,22,WHITE,True)
        t("8 p + 8 n",431,514,14,MUTED)
        # particle channel labels
        t("32 free p + 42 free n",560,231,16,WHITE,True)
        t("After neutron beta decay:",560,258,15,MUTED)
        t("+p  +electron  +anti-neutrino",560,280,14,ORANGE)
        pygame.draw.line(canvas,(60,103,129),(548,350),(642,350),2)
        pygame.draw.polygon(canvas,(60,103,129),[(648,350),(638,342),(638,358)])
        # Penning trap drawing
        cx,cy=803,377
        for radius,c in [(117,(47,74,101)),(90,(61,98,126)),(61,(40,74,94))]:
            pygame.draw.circle(canvas,c,(cx,cy),radius,2)
        pygame.draw.arc(canvas,CYAN,(cx-117,cy-117,234,234),0.42,2.74,5)
        pygame.draw.arc(canvas,CYAN,(cx-117,cy-117,234,234),3.59,5.89,5)
        t("PENNING TRAP",730,515,18,WHITE,True)
        t("electric + magnetic fields",698,543,13,MUTED)

        params=Inputs(mass_mg=mass_levels[idx],capture_efficiency=capture,elapsed_seconds=elapsed)
        res=calculate(params)
        trap=res["trap"]
        omega_plus=trap["modified_cyclotron_Hz"]*2*math.pi
        omega_minus=trap["magnetron_Hz"]*2*math.pi
        # This is a time-rescaled analytical XY Penning-orbit projection.
        sim_orbit_t=seconds_visible*1e-8
        if capture>0:
            count=max(1,round(22*capture))
            for i in range(count):
                phase=(i/count)*math.tau
                angle=omega_minus*sim_orbit_t+phase
                cycl=omega_plus*sim_orbit_t+phase*2.73
                x=cx + 84*math.cos(angle) + 14*math.cos(cycl)
                y=cy - 84*math.sin(angle) - 14*math.sin(cycl)
                pygame.draw.circle(canvas,GREEN,(int(x),int(y)),5)
        # Emission streaks, purely illustrative and not linked to exact particle trajectory.
        if demo_started:
            for p in particles:
                a=p["age"]
                if a<=0 or a>5:
                    continue
                phase=min(1.0,a/2.6)
                x=550+218*phase
                y=352+p["vy"]*(phase**1.25)
                if p["kind"]=="p":
                    if phase>.85 and p["capture"]:continue
                    color=GREEN if p["capture"] else RED
                else:
                    if phase>.90: continue
                    color=BLUE
                pygame.draw.circle(canvas,color,(int(x),int(y)),3)
        t("red: p",65,571,13,RED)
        t("blue: n",160,571,13,BLUE)
        t("green: modeled captured p",258,571,13,GREEN)

        # Metrics right column
        t("ENERGY & PARTICLES",1017,161,17,WHITE,True)
        def metric(y,label,value,unit=""):
            t(label,1017,y,13,MUTED)
            t(value,1017,y+19,22,WHITE,True)
            if unit:t(unit,1290,y+26,12,MUTED)
        metric(199,"Q value, one nucleus",f"{Q_BREAKUP_MEV:,.2f} MeV")
        metric(260,"Minimum energy, selected mass",fmt_si(res["minimum_breakup_energy_J"],"J"))
        metric(321,"Prompt free protons",fmt_num(res["prompt_free_protons"]))
        metric(382,"Free protons after neutron decay",fmt_num(res["available_free_protons"]))
        metric(443,"Captured, if precooled and accepted",fmt_num(res["captured_protons_expected_idealized"]))
        t("QCD SU(3): f123 = +1, f213 = -1",1017,552,12,CYAN)
        t(f"String potential at 1 fm: {STRING_ENERGY_1FM_GEV:.3f} GeV",1017,571,12,MUTED)
        t("No proton rest-mass -> light conversion",1017,589,12,ORANGE)

        # bottom compact controls and dynamics
        t("CONTROL PANEL",48,638,16,WHITE,True)
        t(f"[M] Mass: {mass_levels[idx]:g} mg",48,671,17,CYAN,True)
        t(f"[C] Capture: {capture:.0%} (assumed)",340,671,17,GREEN,True)
        t(f"[T] Time scale: {accelerated:g} simulated s/s",660,671,17,PURPLE,True)
        t("[B] Breakup animation  [SPACE] Pause  [R] Reset  [ESC] Quit",48,706,15,MUTED)
        t(f"Neutron decay: {res['neutron_decay_fraction']:.1%}",48,745,15,WHITE,True)
        t(f"Beta energy: {fmt_si(res['beta_decay_total_energy_J'],'J')}",350,745,15,ORANGE,True)
        t(f"Free neutrons left: {fmt_num(res['free_neutrons_remaining'])}",720,745,15,BLUE,True)
        t("Ideal trap stability: " + ("LINEARIZED STABLE" if trap["stable_linear_penning_model"] else "UNSTABLE"),48,782,13,GREEN)
        t(f"B=5 T, V=10 V, d=1 cm  |  cyclotron={trap['cyclotron_Hz']/1e6:.2f} MHz  |  axial={trap['axial_Hz']/1e3:.1f} kHz",390,782,13,MUTED)
        
        # Small neutron decay curve at center of bottom panel to indicate actual model.
        # (drawn behind/no overlap with controls if moved below)
        win=pygame.display.get_surface()
        tw,th=win.get_size()
        scaled=pygame.transform.smoothscale(canvas,(tw,th)) if (tw,th)!=(W,H) else canvas
        win.blit(scaled,(0,0))
        pygame.display.flip()

    def fmt_num(v):
        return f"{v:,.2e}" if abs(v)>=1e7 else f"{v:,.3g}"
    def fmt_si(v,unit):
        if abs(v)>=1e9:return f"{v/1e9:,.3f} G{unit}"
        if abs(v)>=1e6:return f"{v/1e6:,.3f} M{unit}"
        if abs(v)>=1e3:return f"{v/1e3:,.3f} k{unit}"
        if abs(v)>=1:return f"{v:,.3f} {unit}"
        return f"{v:.3e} {unit}"

    running=True
    while running:
        dt=min(clock.tick(60)/1000.0,0.08)
        for event in pygame.event.get():
            if event.type==pygame.QUIT:running=False
            if event.type==pygame.KEYDOWN:
                if event.key==pygame.K_ESCAPE: running=False
                elif event.key==pygame.K_m:idx=(idx+1)%len(mass_levels)
                elif event.key==pygame.K_c:
                    capture=round((capture+0.05)%1.05,2)
                    # Re-spawn shows new accepted fraction.
                    if demo_started:spawn()
                elif event.key==pygame.K_t:
                    accelerated={120.0:1200.0,1200.0:12000.0,12000.0:120.0}.get(accelerated,120.0)
                elif event.key==pygame.K_SPACE:paused=not paused
                elif event.key==pygame.K_b:spawn()
                elif event.key==pygame.K_r:
                    demo_started=False
                    particles=[]
                    elapsed=seconds_visible=0.0
                    paused=False
        update(dt)
        draw()
    pygame.quit()


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--headless",action="store_true",help="print computed model JSON; no Pygame needed")
    p.add_argument("--self-test",action="store_true",help="verify formulas without Pygame")
    p.add_argument("--mass-mg",type=float,default=1.0)
    p.add_argument("--capture-efficiency",type=float,default=0.75)
    p.add_argument("--elapsed-seconds",type=float,default=3600)
    p.add_argument("--observation-years",type=float,default=1.0)
    p.add_argument("--proton-kinetic-mev",type=float,default=0.0,
                   help="ASSUMED extra kinetic energy for each prompt proton (not predicted)")
    p.add_argument("--stopping-photon-fraction",type=float,default=0.0,
                   help="ASSUMED fraction of that kinetic energy emitted as photons, from 0 to 1")
    a=p.parse_args()
    if a.self_test:
        self_test()
        return
    if a.headless:
        res=calculate(Inputs(mass_mg=a.mass_mg,capture_efficiency=a.capture_efficiency,
                             elapsed_seconds=a.elapsed_seconds,observation_years=a.observation_years,
                             prompt_proton_kinetic_MeV=a.proton_kinetic_mev,
                             stopping_photon_fraction=a.stopping_photon_fraction))
        print(json.dumps(res,indent=2,ensure_ascii=False))
    else:gui()

if __name__=="__main__":
    main()
