"""Snap-fit, wall and moulding numbers for the SH-100 / SH-200 housing.

Plain Python, no Fusion needed:  python snap_calc.py
Geometry values match the defaults in snap_fit_housing.py. Material data is
SABIC CYCOLAC MG47F (ABS), datasheet values from ISO test methods.
"""
import math

# --- material: SABIC CYCOLAC MG47F ABS
E_t = 2370.0      # MPa tensile modulus, ISO 527 1 mm/min
E_f = 2200.0      # MPa flexural modulus, ISO 178 (used as secant modulus at <1 % strain)
sig_y = 47.0      # MPa tensile stress at yield, ISO 527 50 mm/min
eps_y = 0.026     # tensile strain at yield, ISO 527
rho = 1.04        # g/cm3
mu = 0.6          # ABS on ABS, upper end of the usual 0.5-0.6 snap-fit design range

# --- allowable strain (the rule in the Celanese and Covestro snap-fit design guides):
eps_once = 0.70 * eps_y          # one-time assembly: 70 % of yield strain
eps_rep = 0.60 * eps_once        # engaged/disengaged repeatedly: 60 % of that

# --- arm geometry (snap_fit_housing.py)
L = 12.1          # mm, arm root (cover underside) to the catch face
t = 1.2           # mm, arm thickness (0.6 x 2.0 wall)
b = 12.0          # mm, arm width at root (narrows 0.5 deg/side to 11.76 at the tip)
barb = 0.7        # mm, barb overhang off the arm face
clr = 0.10        # mm, arm-to-wall clearance at the catch
lead = 30.0       # deg, lead-in face from the insertion axis
stop_gap = 0.85   # mm, over-deflection stop behind the arm tip

# --- section and moulding inputs
tw = 2.0          # mm nominal wall
alpha = 0.12      # mm2/s ABS thermal diffusivity (typical)
T_melt, T_mould, T_eject = 240.0, 60.0, 90.0   # deg C
vol_body_cm3, vol_cover_cm3 = 15.799, 7.056     # measured from the built model


def strain(y, L=L, t=t):
    """Constant-section cantilever with a rigid root: eps = 1.5 t y / L^2."""
    return 1.5 * t * y / L ** 2


def deflection_force(eps, L=L, t=t, b=b, E=E_f):
    """Tip force for strain eps: P = (b t^2 / 6) E eps / L."""
    return b * t ** 2 / 6.0 * E * eps / L


def mating_force(P, a_deg, mu=mu):
    """Axial push to deflect the arm up a lead-in of angle a with friction mu."""
    ta = math.tan(math.radians(a_deg))
    return P * (mu + ta) / (1 - mu * ta)


def tee_inscribed(rib_t, wall=tw):
    """Diameter of the largest circle at a rib/wall T-junction (sharp corners).
    Circle tangent to the far wall face, passing the rib root corners:
    R^2 = (wall - R)^2 + (rib_t/2)^2  ->  D = (wall^2 + (rib_t/2)^2) / wall."""
    return (wall ** 2 + (rib_t / 2) ** 2) / wall


def corner_inscribed(wall=tw, r_in=1.0):
    """Floor/wall corner, inner fillet r_in: circle tangent to both outer faces and
    the fillet, from (wall + r_in - r) * sqrt2 = r_in + r."""
    return 2 * ((wall + r_in) * math.sqrt(2) - r_in) / (1 + math.sqrt(2))


def cooling_time(wall=tw):
    """Plate cooling time to ejection temperature (1-D conduction, first term)."""
    return wall ** 2 / (math.pi ** 2 * alpha) * math.log(8 / math.pi ** 2 * (T_melt - T_mould) / (T_eject - T_mould))


def calc():
    out = {'allowable': {'one_time': eps_once, 'repeated': eps_rep}}
    y_nom = barb - clr            # deflection to clear the wall, nominal
    y_max = barb                  # zero clearance (worst case of the tolerance stack)
    for name, y in (('assembly_nominal', y_nom), ('assembly_worst', y_max), ('at_stop', stop_gap)):
        e = strain(y)
        out[name] = {'y_mm': y, 'strain': e, 'SF_vs_repeated': eps_rep / e,
                     'SF_vs_once': eps_once / e, 'P_N': deflection_force(e)}
    W1 = mating_force(out['assembly_nominal']['P_N'], lead)
    out['push_on'] = {'per_arm_N': W1, 'two_arms_N': 2 * W1}
    # retention: a 90 deg catch is self-locking, so the limit is shearing the barb off
    barb_root_len = 1.51          # mm, barb root length along the arm (z 4.39 -> 5.90)
    shear_width = 0.98 * b        # mm, arm width at the barb after the 0.5 deg side draft
    tau = 0.6 * sig_y             # MPa, shear strength ~0.6 x tensile yield
    out['retention'] = {'barb_shear_N_per_arm': shear_width * barb_root_len * tau,
                        'arm_tension_N_per_arm': b * t * sig_y}
    # root thickness of each feature that meets the 2.0 wall (after draft)
    out['sections'] = {
        name: {'root_thickness': rt, 'ratio_to_wall': rt / tw,
               'inscribed': tee_inscribed(rt), 'thick_spot': tee_inscribed(rt) / tw - 1}
        for name, rt in (('cradle_rib', 1.0), ('stop_rib', 1.08), ('boss_wall', 1.25))}
    out['sections']['floor_wall_corner'] = {'inscribed': corner_inscribed(),
                                            'ratio_to_wall': corner_inscribed() / tw}
    out['moulding'] = {'cooling_s': cooling_time(), 'flow_L_over_t_body': 60 / tw,
                       'flow_L_over_t_cover': 62 / tw}
    out['mass_g'] = {'body': vol_body_cm3 * rho, 'cover': vol_cover_cm3 * rho}
    return out


def main():
    out = calc()

    def pct(x):
        return '%.2f %%' % (100 * x)

    print('allowable strain: one-time %s, repeated %s'
          % (pct(out['allowable']['one_time']), pct(out['allowable']['repeated'])))
    for k in ('assembly_nominal', 'assembly_worst', 'at_stop'):
        d = out[k]
        print('%-17s y=%.2f  eps=%s  SF(rep)=%.2f  SF(once)=%.2f  P=%.2f N' % (
            k, d['y_mm'], pct(d['strain']), d['SF_vs_repeated'], d['SF_vs_once'], d['P_N']))
    print('push-on: %.1f N per arm, %.1f N to close' % (out['push_on']['per_arm_N'], out['push_on']['two_arms_N']))
    print('retention: barb shear %.0f N/arm, arm tension %.0f N/arm' % (
        out['retention']['barb_shear_N_per_arm'], out['retention']['arm_tension_N_per_arm']))
    for k, v in out['sections'].items():
        print('%-18s inscribed %.2f mm (%.0f %% of wall)' % (k, v['inscribed'], 100 * v['inscribed'] / tw))
    print('cooling %.1f s, mass body %.1f g, cover %.1f g'
          % (out['moulding']['cooling_s'], out['mass_g']['body'], out['mass_g']['cover']))


if __name__ == '__main__':
    main()
