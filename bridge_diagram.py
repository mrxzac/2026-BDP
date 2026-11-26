"""Show the bridge and compare geometric section properties at two x stations."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon
from matplotlib.widgets import TextBox
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from compute import d, u
from config_0 import beam_0


INK, MUTED = '#203346', '#526676'
BOARD, PLATE, CUT = '#a8bbc5', '#e5a055', '#126596'


def _board_faces(board):
    """Extrude the polygon about its midplane by the actual matboard thickness."""
    offset = np.zeros(3)
    offset[d[board.dirc]] = board.t.m_as(u.mm) / 2
    first, last = board.vertices - offset, board.vertices + offset
    walls = [[first[i-1], first[i], last[i], last[i-1]] for i in range(len(first))]
    # Display x along the girder, z across it, and y vertically.
    return [np.asarray(f)[:, [0, 2, 1]] for f in [first, last, *walls]]


def _draw_bridge(ax, beam):
    boards = [(b, _board_faces(b)) for b in beam.matboards]
    points = np.concatenate([f for _, faces in boards for f in faces])
    lo, hi = points.min(axis=0), points.max(axis=0)
    mid = (lo + hi) / 2
    for board, faces in boards:
        transverse = board.dirc == 'x'
        # Open the deck and near web in the rendering so the diaphragms are visible.
        ghost = ((board.dirc == 'y' and board.vertices[:, 1].mean() > mid[2])
                 or (board.dirc == 'z' and board.vertices[:, 2].mean() < mid[1]))
        if ghost:
            outline = np.vstack((faces[0], faces[0][0]))
            ax.plot(*outline.T, color=BOARD, lw=0.65, alpha=0.65, zorder=2)
        else:
            surface = Poly3DCollection(
                faces, facecolors=PLATE if transverse else BOARD, alpha=0.94 if transverse else 0.35,
                edgecolors='#92602d' if transverse else MUTED, linewidths=0.45)
            ax.add_collection3d(surface)
    ax.set(xlim=(lo[0], hi[0]), ylim=(lo[1], hi[1]), zlim=(lo[2], hi[2]))
    ax.set_box_aspect(hi - lo, zoom=1.13)
    ax.set_proj_type('ortho')
    ax.view_init(elev=25, azim=-133)
    ax.set_axis_off()
    for x in [lo[0], hi[0]]:
        ax.text(x, lo[1], lo[2] - 22, f'{x:,.0f} mm', ha='left', va='top', color=MUTED, fontsize=9)


def _draw_section(ax, beam, x, label, color):
    ax.clear()
    outlines = []
    has_diaphragm = False
    for board in beam.matboards:
        points = board.get_sec_mm('x', x)
        if len(points) < 3:
            continue
        transverse = board.dirc == 'x'
        has_diaphragm |= transverse
        outlines.append(points)
        ax.add_patch(Polygon(points[:, [2, 1]], facecolor=PLATE if transverse else BOARD,
                             edgecolor='#97632f' if transverse else MUTED, linewidth=0.7))
    if not outlines:
        raise ValueError(f'No cross-section at x = {x:g} mm.')
    points = np.concatenate(outlines)
    low, high = points.min(axis=0), points.max(axis=0)
    cy = beam.get_centroid('y', x).m_as(u.mm)
    cz = beam.get_centroid('z', x).m_as(u.mm)
    iz = beam.get_2nd_moment('z', x).m_as(u.mm**4)
    iy = beam.get_2nd_moment('y', x).m_as(u.mm**4)
    ax.axhline(cy, color=MUTED, lw=0.7, ls='--')
    ax.axvline(cz, color=MUTED, lw=0.7, ls='--')
    ax.plot(cz, cy, '.', color=color, ms=7)
    ax.set(xlim=(low[2]-8, high[2]+8), ylim=(low[1]-6, high[1]+6), aspect='equal')
    ax.set_axis_off()
    ax.set_title(f'{label}  ·  x = {x:g} mm', loc='left', color=color, fontsize=12, pad=8)
    ax.text(1.04, 0.86, 'Diaphragm' if has_diaphragm else 'Open box',
            transform=ax.transAxes, color=INK, fontsize=11, weight='bold')
    ax.text(1.04, 0.64, f'$I_z$   {iz:,.0f} mm⁴\n$I_y$   {iy:,.0f} mm⁴\n$ȳ$     {cy:.2f} mm',
            transform=ax.transAxes, va='top', color=INK, fontsize=10, linespacing=1.65)
    margin = (high - low) * 0.05
    z0, z1 = low[2] - margin[2], high[2] + margin[2]
    y0, y1 = low[1] - margin[1], high[1] + margin[1]
    return [x]*5, [z0, z1, z1, z0, z0], [y0, y0, y1, y1, y0]


def plot_bridge(beam, x, reference=None, interactive=True):
    """Return a rotatable 3D figure. Coordinates are millimetres in the beam frame."""
    start, end = (p.m_as(u.mm) for p in beam.get_endpoints())
    if reference is None:
        reference = (start + end) / 2
    if not (start <= x <= end and start <= reference <= end):
        raise ValueError(f'Section positions must be between {start:g} and {end:g} mm.')

    fig = plt.figure(figsize=(11.4, 5.2), facecolor='white')
    fig.canvas.manager.set_window_title('Bridge section properties')
    fig.text(0.045, 0.92, 'Bridge geometry and section inertia', size=19, color=INK, weight='bold')
    ax3 = fig.add_axes([0.01, 0.065, 0.57, 0.90], projection='3d', computed_zorder=False)
    _draw_bridge(ax3, beam)
    selected_ax = fig.add_axes([0.615, 0.54, 0.16, 0.25])
    reference_ax = fig.add_axes([0.615, 0.17, 0.16, 0.25])
    ref_ring = _draw_section(reference_ax, beam, reference, 'Reference', MUTED)
    ax3.plot(*ref_ring, color=MUTED, lw=1.0, ls=':', zorder=10)
    ring, = ax3.plot(*_draw_section(selected_ax, beam, x, 'Selected section', CUT),
                     color=CUT, lw=1.5, ls='--', zorder=11)
    fig.legend(handles=[Patch(facecolor=PLATE, label='Diaphragm'),
                        Line2D([], [], color=CUT, ls='--', label='Selected cut'),
                        Line2D([], [], color=MUTED, ls=':', label='Reference cut')],
               loc='center', bbox_to_anchor=(.30, .195), ncol=3, frameon=False, fontsize=9)

    if interactive:
        entry = TextBox(fig.add_axes([0.21, 0.055, 0.10, 0.055]), 'Section x (mm)   ', initial=f'{x:g}')
        feedback = fig.text(.33, .067, 'Enter a position and press Return.', color=MUTED, size=9)

        def move_section(value):
            try:
                pos = float(value)
            except ValueError:
                pos = np.nan
            if not start <= pos <= end:
                feedback.set_text(f'Enter a value from {start:g} to {end:g} mm.')
            else:
                ring.set_data_3d(*_draw_section(selected_ax, beam, pos, 'Selected section', CUT))
                feedback.set_text('Enter a position and press Return.')
            fig.canvas.draw_idle()

        entry.on_submit(move_section)
        fig.section_entry = entry  # Matplotlib widgets require a live reference.
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--x', type=float, default=425, help='selected section in mm (default: 425)')
    parser.add_argument('--reference', type=float, help='comparison section in mm (default: midspan)')
    parser.add_argument('--save', metavar='FILE', help='save the figure instead of opening a window')
    args = parser.parse_args()
    try:
        fig = plot_bridge(beam_0, args.x, args.reference, interactive=args.save is None)
    except ValueError as error:
        parser.error(str(error))
    if args.save:
        fig.savefig(args.save, dpi=160, bbox_inches='tight', pad_inches=0.18)
        plt.close(fig)
    else:
        plt.show()


if __name__ == '__main__':
    main()
