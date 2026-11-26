"""Run the CIV102 bridge calculations. Use --help for cases and output options."""

import argparse
import csv
from importlib import import_module
from pathlib import Path

import numpy as np

from compute import u
from Simulation import Simulation
from strength import analyze, check_section, MODE_LABELS, SHEAR_MODES


LOAD_CASES = {
    '1': ('Load Case 1', np.full(6, 400 / 6)),
    # Use the original centered example's ordering (locomotive at the right).
    '2': ('Load Case 2, base pass', np.array([67.5, 67.5, 67.5, 67.5, 91., 91.])),
}


def peak(envelope, prefix):
    """Return the signed extreme and station with greatest absolute demand."""
    low, high = envelope[prefix + '_min'], envelope[prefix + '_max']
    record = low if -low.values.min() > high.values.max() else high
    return record, np.abs(record.values).argmax()


def centered_check(sim, data, plot):
    position = (sim.span - sim.axle_positions[-1]) / 2
    sim.loads(position)
    print(f'  Centered train: reactions L = {sim.reac_L.m_as(u.N):.3f} N, '
          f'R = {sim.reac_R.m_as(u.N):.3f} N')
    m = sim.bmd.m_as(u.N * u.mm)
    # Include the actual wheel locations even when the reporting grid is coarse.
    event_x, _, _, event_m = sim.diagram(position)
    off_grid = ~np.isin(event_x, sim.x.magnitude)
    extra = [(moment, check_section(sim.beam, x))
             for x, moment in zip(event_x[off_grid], event_m[off_grid])]
    for mode in ('tension', 'compression'):
        utilization = np.maximum(m / data[mode + '_capacity_pos'], m / data[mode + '_capacity_neg'])
        worst = utilization.max()
        for moment, (_, positive, negative) in extra:
            limit = positive[mode] if moment >= 0 else negative[mode]
            worst = max(worst, abs(moment) / limit.capacity)
        print(f'    {MODE_LABELS[mode]} FOS = {1 / worst:.4f}')
    if plot:
        sim.plot(show=False)


def summarize(sim, envelope, data):
    x = data['x_mm']
    local = x.size // 2
    print(f'  Section at x = {x[local]:g} mm:')
    for key, values in data.items():
        if key in ('area_mm2', 'ybar_mm', 'Iz_mm4', 'Q_NA_mm3', 'b_NA_mm') or key.startswith('Q_glue_'):
            print(f'    {key} = {values[local]:.3f}')
    for prefix, scale, unit in [('V', 1, 'N'), ('M', 1000, 'N·m')]:
        record, i = peak(envelope, prefix)
        pose = record.first_axle_mm[i]
        side = {-1: 'just below', 0: 'at', 1: 'just above'}[record.approach[i]]
        order = 'reversed' if record.reversed[i] else 'as supplied'
        print(f'  Peak |{prefix}| = {abs(record.values[i]) / scale:.3f} {unit} at x = {x[i]:g} mm'
              f' (leftmost axle {side} {pose:g} mm from left support, order {order})')
    print(f'  {"Failure mode":27} {"Min FOS":>9} {"x (mm)":>9} {"Stress (MPa)":>14} {"Limit (MPa)":>13}')
    for mode, label in MODE_LABELS.items():
        i = data[mode + '_fos'].argmin()
        print(f'  {label:27} {data[mode + "_fos"][i]:9.4f} {x[i]:9.2f}'
              f' {data[mode + "_stress_MPa"][i]:14.4f} {data[mode + "_limit_MPa"][i]:13.4f}')
    mode = min(MODE_LABELS, key=lambda name: data[name + '_fos'].min())
    fos = data[mode + '_fos'].min()
    print(f'  Governing: {MODE_LABELS[mode]}. Predicted total train load at these axle ratios: '
          f'{fos * sim.loadcase.magnitude.sum():.2f} N.')


def write_csv(path, data):
    headers = []
    for key in data:
        if key.startswith('V_') or key.startswith('M_'):
            unit = 'N' if key.startswith('V_') else 'Nmm'
            key = f'{key}_{unit}'
        elif '_capacity_' in key:
            unit = 'N' if key.split('_capacity_')[0] in SHEAR_MODES else 'Nmm'
            key = f'{key}_{unit}'
        headers.append(key)
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        writer.writerows(zip(*data.values()))


def plot_analysis(sim, envelope, data, title):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True, layout='constrained')
    x = data['x_mm']
    for j, (prefix, scale, unit) in enumerate([('V', 1, 'N'), ('M', 1000, 'N·m')]):
        ax = axes[0, j]
        ax.fill_between(x, data[prefix + '_min'] / scale, data[prefix + '_max'] / scale,
                        color='0.8', label='Moving train')
        for mode, label in MODE_LABELS.items():
            if (mode in SHEAR_MODES) != (prefix == 'V'):
                continue
            for sign in ('pos', 'neg'):
                cap = data[f'{mode}_capacity_{sign}'] / scale
                ax.plot(x, np.where(np.isfinite(cap), cap, np.nan),
                        color=f'C{list(MODE_LABELS).index(mode)}', lw=1.4,
                        label=label if sign == 'pos' else None)
        ax.set_title('Shear envelope and capacities' if j == 0 else 'Moment envelope and capacities')
        ax.set_ylabel(unit)
        record, i = peak(envelope, prefix)
        xx, vl, vr, moment = sim.diagram(record.first_axle_mm[i], record.reversed[i], record.approach[i])
        label = f'Train at peak |{prefix}|'
        axes[1, 0].plot(np.repeat(xx, 2), np.column_stack((vl, vr)).ravel(), label=label)
        axes[1, 1].plot(xx, moment / 1000, label=label)
    for ax, label in zip(axes[1], ['Shear (N)', 'Moment (N·m)']):
        ax.set_ylabel(label)
        ax.set_xlabel('Bridge x (mm)')
    for ax in axes.flat:
        ax.axhline(0, color='0.3', lw=0.7)
        ax.grid(alpha=0.25)
        ax.legend(fontsize='small')
    fig.suptitle(title)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['1', '2', 'both', 'later'], default='both',
                        help='1 = 400 N, 2 = 452 N base pass, later = subsequent LC2 ratios')
    parser.add_argument('--weight', type=float, metavar='N', help='total train load for --case later')
    parser.add_argument('--design', default='config_0:beam_0', metavar='MODULE:BEAM',
                        help='Python configuration and Beam object (default: config_0:beam_0)')
    parser.add_argument('--stations', type=int, default=1251, help='number of bridge stations before adding supports')
    parser.add_argument('--csv-dir', type=Path, help='write per-station geometry, stresses, capacities and FOS')
    parser.add_argument('--no-plot', action='store_true', help='print results without opening figures')
    args = parser.parse_args()
    if args.case == 'later':
        if args.weight is None or not np.isfinite(args.weight) or args.weight <= 0:
            parser.error('--case later needs a positive, finite --weight in newtons')
        # Far freight : near freight : locomotive = 1.10 : 1 : (1.39 × 1.10).
        ratio = np.repeat([1.10, 1., 1.39 * 1.10], 2)
        cases = {'later': ('Load Case 2, subsequent pass', args.weight * ratio / ratio.sum())}
    else:
        if args.weight is not None:
            parser.error('--weight applies only to --case later')
        cases = LOAD_CASES if args.case == 'both' else {args.case: LOAD_CASES[args.case]}
    module, separator, name = args.design.partition(':')
    if not separator or not name:
        parser.error('--design requires MODULE:BEAM')
    beam = getattr(import_module(module), name)
    if args.csv_dir:
        args.csv_dir.mkdir(parents=True, exist_ok=True)
    for key, (name, weights) in cases.items():
        sim = Simulation(weights * u.N, beam, args.stations)
        title = f'{args.design} — {name} ({weights.sum():g} N)'
        print(title)
        envelope = sim.envelopes()
        data = analyze(beam, sim.x, envelope)
        if key == '2':
            centered_check(sim, data, not args.no_plot)
        summarize(sim, envelope, data)
        if args.csv_dir:
            path = args.csv_dir / f'case_{key}.csv'
            write_csv(path, data)
            print(f'  Per-station results: {path}')
        if not args.no_plot:
            plot_analysis(sim, envelope, data, title)
    if not args.no_plot:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == '__main__':
    main()
