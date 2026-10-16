"""Figures and an animation from saved prototype runs (run_prototype.py outputs)."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.animation import FuncAnimation, PillowWriter
from model import history, read_slice

# The predeclared dt=0.05 pocket run stopped on its frequency guard at t=25.5;
# the dt=0.025 'timestep' run of the same initial state reached t=100.
LONG = 'timestep'

SURFACE, INK, INK2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
SERIES = dict(pocket='#2a78d6', leaky='#eb6834', control='#1baf7a')
BLUES = LinearSegmentedColormap.from_list('blues', ['#fcfcfb', '#cde2fb', '#86b6ef', '#3987e5', '#1c5cab', '#0d366b'])
ORANGES = LinearSegmentedColormap.from_list('oranges', ['#fcfcfb', '#f9d7c6', '#f2a07c', '#eb6834', '#a8401a', '#5c220c'])
plt.rcParams.update({'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
                     'text.color': INK, 'axes.labelcolor': INK2, 'axes.edgecolor': AXIS, 'xtick.color': MUTED,
                     'ytick.color': MUTED, 'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': .6,
                     'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 2,
                     'font.family': 'sans-serif', 'font.size': 10, 'legend.frameon': False})


def end_label(ax, x, y, text, color):
    ax.annotate(text, (x[-1], y[-1]), xytext=(4, 0), textcoords='offset points', va='center', color=INK2, fontsize=9)
    ax.plot(x[-1], y[-1], 'o', ms=4, color=color)


def axion_density(s, box, mass_tv):
    """Axion energy density on a slice, with in-plane gradients only."""
    dx = box/len(s);phi, a, pa = s[..., 0], s[..., 1], s[..., 3]
    grad = ((np.roll(a, -1, 0)-a)**2+(np.roll(a, -1, 1)-a)**2)/(2*dx*dx)
    return .5*pa**2+grad+.5*mass_tv**2*phi**2*a**2


def slices(directory):
    return [read_slice(p) for p in sorted((Path(directory)/'native').glob('slice_*.bin'))]


def nearest(frames, t):
    return min(frames, key=lambda f: abs(f[0]['u']-t))


def radius_figure(source, control, output):
    pocket = json.loads((source/'pocket/summary.json').read_text())
    fig, ax = plt.subplots(figsize=(8, 4.2))
    tw = np.genfromtxt(source/'pocket/thin_wall.csv', delimiter=',', names=True)
    ax.plot(tw['u'], tw['radius'], '--', lw=1.2, color=SERIES['pocket'], alpha=.7, label='pocket, thin-wall adiabatic ODE')
    radial = np.genfromtxt(control/'control/radial_reference.csv', delimiter=',', names=True)
    # Only up to the reference's own collapse; its later crossings are post-implosion radiation.
    keep = np.arange(len(radial)) < np.argmax(radial['radius'] < 2.)
    ax.plot(radial['u'][keep], radial['radius'][keep], '--', lw=1.2, color=SERIES['control'], alpha=.8,
            label='control, radial PDE reference')
    for name, directory, suite in (('pocket', source, LONG), ('leaky', source, 'leaky'), ('control', control, 'control')):
        h = history(directory/suite)
        ax.plot(h['u'], h['fv_radius'], color=SERIES[name], label=f'{name} (3D lattice)')
        end_label(ax, h['u'], h['fv_radius'], name, SERIES[name])
    stop = history(source/'pocket')
    ax.plot(stop['u'][-1], stop['fv_radius'][-1], 'x', ms=7, color=INK2)
    ax.annotate('dt=0.05 run stopped\n(frequency guard)', (stop['u'][-1], stop['fv_radius'][-1]), xytext=(6, -22),
                textcoords='offset points', color=INK2, fontsize=8)
    for value, text, va in ((pocket['brackets']['adiabatic'], 'adiabatic equilibrium', 'top'),
                            (pocket['brackets']['energy_conserving'], 'energy-conserving equilibrium', 'bottom'),
                            (8., 'survival threshold', 'bottom')):
        ax.axhline(value, color=MUTED, lw=.8, ls=':')
        ax.text(100, value, text, color=MUTED, fontsize=8, ha='right', va=va)
    ax.set(xlabel='time t  (1/m_phi)', ylabel='false-vacuum radius R  (1/m_phi)', xlim=(0, 112), ylim=(0, 21))
    ax.set_title('False-vacuum radius: trapped gas versus controls', loc='left', fontsize=11)
    ax.legend(loc='lower right', fontsize=8, ncol=2, bbox_to_anchor=(1, .08))
    fig.tight_layout();fig.savefig(output/'radius.png', dpi=160);plt.close(fig)


def checks_figure(source, control, output):
    gates = {s: json.loads((source/s/'summary.json').read_text()) for s in ('pocket', 'timestep', 'spatial')}
    fig, axes = plt.subplots(1, 4, figsize=(14, 3.2))
    for ax, suite, base_suite, tol, title in (
            (axes[0], 'timestep', 'pocket', .05, 'R(dt=0.025) - R(dt=0.05), N=128'),
            (axes[1], 'spatial', LONG, .5, 'post hoc: R(N=96) - R(N=128, dt=0.025)')):
        h, base = history(source/suite), history(source/base_suite)
        shared, i, j = np.intersect1d(np.round(h['u'], 9), np.round(base['u'], 9), return_indices=True)
        ax.axhspan(-tol, tol, color=GRID, alpha=.6, lw=0)
        ax.plot(shared, h['fv_radius'][i]-base['fv_radius'][j], color=SERIES['pocket'])
        ax.set_title(title, loc='left', fontsize=10);ax.set_xlim(0, 100)
        ax.set_xlabel('t');ax.text(.02, .04, f'tolerance +/-{tol}', transform=ax.transAxes, color=MUTED, fontsize=8)
    for name, directory, suite in (('pocket', source, LONG), ('leaky', source, 'leaky'), ('control', control, 'control')):
        h = history(directory/suite)
        axes[2].plot(h['u'], (h['total']-h['total'][0])/h['total'][0], color=SERIES[name], lw=1.4, label=name)
    axes[2].axhspan(-1e-4, 1e-4, color=GRID, alpha=.6, lw=0)
    axes[2].set_title('energy drift (E-E0)/E0', loc='left', fontsize=10);axes[2].set_xlabel('t')
    axes[2].set_yscale('symlog', linthresh=1e-5);axes[2].legend(fontsize=8)
    h = history(control/'control');radial = np.genfromtxt(control/'control/radial_reference.csv', delimiter=',', names=True)
    shared, i, j = np.intersect1d(np.round(h['u'], 9), np.round(radial['u'], 9), return_indices=True)
    keep = radial['radius'][j] >= 8
    axes[3].axhspan(-.25, .25, color=GRID, alpha=.6, lw=0)
    axes[3].plot(shared[keep], h['fv_radius'][i][keep]-radial['radius'][j][keep], color=SERIES['control'])
    axes[3].set_title('control: R(3D) - R(radial PDE)', loc='left', fontsize=10);axes[3].set_xlabel('t')
    axes[3].text(.02, .04, 'gate +/-0.25 while R>=8', transform=axes[3].transAxes, color=MUTED, fontsize=8)
    fig.tight_layout();fig.savefig(output/'checks.png', dpi=160);plt.close(fig)


def energy_figure(source, output):
    h = history(source/LONG);axion0 = h['axion_kinetic'][0]+h['axion_gradient'][0]+h['coupling'][0]
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.6))
    for key, label, color in (('axion_fv', 'in false vacuum (phi<1/2)', SERIES['pocket']),
                              ('axion_wall', 'in wall (1/2<=phi<=0.99)', '#eb6834'),
                              ('axion_tv', 'deep true vacuum (phi>0.99)', '#1baf7a')):
        axes[0].plot(h['u'], h[key]/axion0, color=color, lw=1.6, label=label)
    axes[0].set_yscale('log');axes[0].set_ylim(1e-8, 2)
    axes[0].set_title('pocket (dt=0.025): axion energy by region / initial', loc='left', fontsize=10)
    axes[0].set_xlabel('t');axes[0].legend(fontsize=8, loc='center right')
    axes[1].plot(h['u'], h['dilaton_kinetic'], color=SERIES['pocket'], lw=1.6)
    axes[1].set_title('pocket (dt=0.025): dilaton kinetic energy', loc='left', fontsize=10)
    axes[1].set_xlabel('t')
    fig.tight_layout();fig.savefig(output/'energy.png', dpi=160);plt.close(fig)


def slice_panel(ax, s, box, kind, mass_tv, vmax=None):
    extent = (0, box, 0, box)
    if kind == 'phi':
        return ax.imshow(s[..., 0].T, origin='lower', extent=extent, cmap=BLUES, vmin=0, vmax=1)
    rho = axion_density(s, box, mass_tv)
    return ax.imshow(np.maximum(rho, vmax*1e-4).T, origin='lower', extent=extent, cmap=ORANGES,
                     norm=LogNorm(vmax*1e-4, vmax))


def slices_figure(source, control, output, mass_tv=2., leaky_mass=.3):
    frames = slices(source/LONG);box = frames[0][0]['box']
    times = (0, 18, 36, 60, 80, 100)
    vmax = max(axion_density(nearest(frames, t)[1], box, mass_tv).max() for t in times)
    fig, axes = plt.subplots(2, len(times), figsize=(2.6*len(times), 5.6))
    for col, t in enumerate(times):
        meta, s = nearest(frames, t)
        im0 = slice_panel(axes[0, col], s, box, 'phi', mass_tv)
        im1 = slice_panel(axes[1, col], s, box, 'axion', mass_tv, vmax)
        axes[0, col].set_title(f't = {meta["u"]:g}', fontsize=10)
        for ax in axes[:, col]:
            ax.set_xticks([]);ax.set_yticks([]);ax.grid(False)
    axes[0, 0].set_ylabel('dilaton phi');axes[1, 0].set_ylabel('axion energy density')
    fig.colorbar(im0, ax=axes[0], fraction=.015, pad=.01, label='phi (0 false, 1 true vacuum)')
    fig.colorbar(im1, ax=axes[1], fraction=.015, pad=.01, label='axion energy density (log)')
    fig.suptitle('Pocket run (dt=0.025): central z-slice', x=.02, ha='left', fontsize=11)
    fig.savefig(output/'slices_pocket.png', dpi=140, bbox_inches='tight');plt.close(fig)
    fig, axes = plt.subplots(3, 3, figsize=(8, 8))
    ctrl = slices(control/'control');leaky = slices(source/'leaky')
    vleak = axion_density(leaky[0][1], box, leaky_mass).max()
    for col, t in enumerate((0, 12, 24)):
        meta, s = nearest(ctrl, t);slice_panel(axes[0, col], s, box, 'phi', 0.)
        axes[0, col].set_title(f'control t = {meta["u"]:g}', fontsize=9)
    for col, t in enumerate((0, 50, 100)):
        meta, s = nearest(leaky, t)
        slice_panel(axes[1, col], s, box, 'phi', leaky_mass)
        slice_panel(axes[2, col], s, box, 'axion', leaky_mass, vleak)
        axes[1, col].set_title(f'leaky t = {meta["u"]:g}', fontsize=9)
    for ax in axes.ravel():
        ax.set_xticks([]);ax.set_yticks([]);ax.grid(False)
    axes[0, 0].set_ylabel('control phi');axes[1, 0].set_ylabel('leaky phi');axes[2, 0].set_ylabel('leaky axion energy')
    fig.suptitle('Controls: no gas (top), m_TV=0.3 leaky walls (bottom)', x=.02, ha='left', fontsize=11)
    fig.tight_layout();fig.savefig(output/'slices_controls.png', dpi=140);plt.close(fig)


def animation(source, output, mass_tv=2.):
    frames = slices(source/LONG);box = frames[0][0]['box']
    h = history(source/LONG)
    vmax = max(axion_density(s, box, mass_tv).max() for _, s in frames[::10])
    fig = plt.figure(figsize=(10, 3.6))
    ax0 = fig.add_axes([.01, .05, .28, .8]);ax1 = fig.add_axes([.31, .05, .28, .8]);ax2 = fig.add_axes([.67, .17, .31, .68])
    im0 = slice_panel(ax0, frames[0][1], box, 'phi', mass_tv);im1 = slice_panel(ax1, frames[0][1], box, 'axion', mass_tv, vmax)
    for ax, title in ((ax0, 'dilaton phi'), (ax1, 'axion energy density (log)')):
        ax.set_xticks([]);ax.set_yticks([]);ax.grid(False);ax.set_title(title, fontsize=9, loc='left')
    ax2.plot(h['u'], h['fv_radius'], color=SERIES['pocket'], lw=1.4)
    marker, = ax2.plot([0], [h['fv_radius'][0]], 'o', ms=6, color=SERIES['pocket'])
    ax2.set(xlabel='t', ylabel='R', ylim=(0, 21));ax2.set_title('false-vacuum radius', fontsize=9, loc='left')
    label = fig.text(.01, .93, '', fontsize=10)

    def update(k):
        meta, s = frames[k];im0.set_data(s[..., 0].T)
        im1.set_data(np.maximum(axion_density(s, box, mass_tv), vmax*1e-4).T)
        i = np.argmin(np.abs(h['u']-meta['u']));marker.set_data([h['u'][i]], [h['fv_radius'][i]])
        label.set_text(f'relic pocket prototype   t = {meta["u"]:.0f}')
        return im0, im1, marker, label
    FuncAnimation(fig, update, frames=len(frames), blit=False).save(output/'pocket.gif', writer=PillowWriter(fps=12), dpi=80)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True);p.add_argument('--control', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True);args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    radius_figure(args.source, args.control, args.output)
    checks_figure(args.source, args.control, args.output)
    energy_figure(args.source, args.output)
    slices_figure(args.source, args.control, args.output)
    animation(args.source, args.output)


if __name__ == '__main__':
    main()
