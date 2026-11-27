"""Independent section and moving-load calculation for Design 0, Load Case 1."""

import math


# Appendix II dimensions, measured from the underside of the bottom flange (mm).
t = 1.27
deck_width = 100
box_width = 80
deck_underside = 75
tab_width = 5
height = deck_underside + t
web_height = deck_underside - t
glue_width = 2 * tab_width

# Each tuple gives width, height and centroid elevation of a rectangular part.
rectangles = (
    (deck_width, t, deck_underside + t / 2),
    (box_width, t, t / 2),
    (2 * t, web_height, t + web_height / 2),  # Both webs
    (glue_width, t, deck_underside - t / 2),  # Both tabs
)
area = sum(b * h for b, h, y in rectangles)
yBar = sum(b * h * y for b, h, y in rectangles) / area
I = sum(b * h**3 / 12 + b * h * (y - yBar)**2 for b, h, y in rectangles)

# First moments of the bottom flange and lower webs at the neutral axis,
# and of the deck alone immediately above the glue joint (mm³).
Q_cent = box_width * t * (yBar - t / 2) + t * (yBar - t)**2
Q_glue = deck_width * t * (deck_underside + t / 2 - yBar)

P = -400  # Downward total train load (N).
span = 1200
axle_offsets = (0, 176, 340, 516, 680, 856)
train_positions = range(-axle_offsets[-1], span + 1)  # Leftmost axle, from the left support.
# Same diaphragm locations as config_0, expressed relative to the left support.
diaphragms = (-15, 15, 400, 800, 1185, 1215)


def Define_bridge(Xini, P):
    """Set the axle loads on the span and the two support reactions (N)."""
    global Forces, Ay, By
    Forces = {Xini + dx: P / len(axle_offsets) for dx in axle_offsets
              if 0 <= Xini + dx <= span}
    By = -sum(x * force for x, force in Forces.items()) / span
    Ay = -sum(Forces.values()) - By


Define_bridge(172, P)


def find_shear(mm, *, right=False):
    """Shear at an in-span cut (N). At supports, use the limit inside the span."""
    if mm == 0 or (right and mm < span):
        return Ay + sum(force for x, force in Forces.items() if x <= mm)
    return Ay + sum(force for x, force in Forces.items() if x < mm)


def find_moment(mm):
    """Internal bending moment from the forces to the left of a cut (N·mm)."""
    return Ay * mm + sum(force * (mm - x) for x, force in Forces.items() if x < mm)


def Highest_moment(Xini, P):
    Define_bridge(Xini, P)
    # Moment is linear between axles and zero at each support.
    return max((find_moment(x) for x in Forces), default=0.0)


def find_max_moments():
    return max(Highest_moment(position, P) for position in train_positions)


def find_max_shear(start=0, end=span):
    """Peak absolute shear within a span interval for the 1 mm train sweep."""
    peak = 0.0
    for position in train_positions:
        Define_bridge(position, P)
        # Downward axle loads make shear nonincreasing between the supports.
        peak = max(peak, abs(find_shear(start, right=True)), abs(find_shear(end)))
    return peak


def Calc_flex_yield_tens(MaxM, I, Ytens):
    return 30 / (MaxM * Ytens / I)


def Calc_flex_yield_comp(MaxM, I, Ycomp):
    return 6 / (MaxM * Ycomp / I)


def Calc_shear_stress_fail(MaxV, I, Bintrest, Qintrest, matType):
    strength = {"glue": 2, "matboard": 4}[matType]
    return strength / (MaxV * Qintrest / (I * Bintrest))


def Calc_plate_buckling(CaseNum, layers, b, h=0, a=0, *, panel_start=0, panel_end=span):
    plate_factor = math.pi**2 * 4000 / (12 * (1 - 0.2**2))
    thickness = layers * t
    if CaseNum == 4:
        capacity = 5 * plate_factor * ((thickness / h)**2 + (thickness / a)**2)
        shear = find_max_shear(panel_start, panel_end)
        applied = shear * Q_cent / (I * 2 * t)
    else:
        coefficient = {1: 4, 2: 0.425, 3: 6}[CaseNum]
        capacity = coefficient * plate_factor * (thickness / b)**2
        applied = find_max_moments() * (height - yBar) / I
    return capacity / applied


if __name__ == "__main__":
    moment = find_max_moments()
    shear = find_max_shear()
    # Use each panel's diaphragm spacing and the maximum shear in that same panel.
    shear_buckling = min(
        Calc_plate_buckling(4, 1, 0, web_height, right - left,
                            panel_start=max(0, left), panel_end=min(span, right))
        for left, right in zip(diaphragms, diaphragms[1:])
    )
    factors = {
        "Matboard tension": Calc_flex_yield_tens(moment, I, yBar),
        "Matboard compression": Calc_flex_yield_comp(moment, I, height - yBar),
        "Matboard shear": Calc_shear_stress_fail(shear, I, 2 * t, Q_cent, "matboard"),
        "Glue shear": Calc_shear_stress_fail(shear, I, glue_width, Q_glue, "glue"),
        "Middle flange buckling": Calc_plate_buckling(1, 1, box_width - 2 * t),
        "Side flange buckling": Calc_plate_buckling(2, 1, (deck_width - box_width) / 2),
        "Web compression buckling": Calc_plate_buckling(3, 1, deck_underside - yBar),
        "Web shear buckling": shear_buckling,
    }
    print(f"Design 0, moving {abs(P)} N train (1 mm steps)")
    print(f"Centroid from bottom = {yBar:.3f} mm, I = {I:,.2f} mm⁴")
    print(f"Maximum shear = {shear:.2f} N, maximum moment = {moment:,.2f} N·mm")
    for mode, fos in factors.items():
        print(f"{mode:<25} {fos:.3f}")
    governing = min(factors, key=factors.get)
    print(f"Overall FOS = {factors[governing]:.3f} ({governing})")
    print(f"Predicted failure load = {abs(P) * factors[governing]:.1f} N")
