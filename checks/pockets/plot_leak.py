"""Figures for the lifetime and leak-rate experiment (run_leak.py outputs)."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from model import history
from plot import INK2, MUTED

# Ordinal single-hue ramps: blue for the m_TV scan, orange for the resolution ladder.
SCAN = {'m1.5_N128': ('m_TV=1.5', '#86b6ef'), 'm2_N128': ('m_TV=2', '#3987e5'),
        'm3_N128': ('m_TV=3', '#1c5cab'), 'm4_N128': ('m_TV=4', '#0d366b')}
LADDER = {'m2_N96': ('N=96, dx=0.67', '#f2a07c'), 'm2_N128': ('N=128, dx=0.5', '#eb6834'),
          'm2_N192': ('N=192, dx=0.33', '#a8401a')}
OMEGA_MAX = .6


def control_line(ax, source, key):
    h = history(source/'control')
    if key == 'fv_radius':
        ax.plot(h['u'], h['fv_radius'], '--', lw=1.2, color=MUTED, label='no gas (control)')


def series(ax, source, runs, key, scale=None):
    for name, (label, color) in runs.items():
        h = history(source/name)
        y = h[key] if scale is None else h[key]/scale(h)
        if key == 'axion_far':
            # The far-field sphere is defined around the false vacuum; undefined once it is gone.
            y = np.where(h['fv_radius'] >= 2., y, np.nan)
        ax.plot(h['u'], y, color=color, lw=1.6, label=label)


def radius_figure(source, output):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4), sharey=True)
    for ax, runs, title in ((axes[0], SCAN, 'false-vacuum radius, m_TV scan (N=128)'),
                            (axes[1], LADDER, 'false-vacuum radius, resolution ladder (m_TV=2)')):
        control_line(ax, source, 'fv_radius');series(ax, source, runs, 'fv_radius')
        ax.axhline(8, color=MUTED, lw=.8, ls=':');ax.text(395, 8, 'death threshold R=8', color=MUTED, fontsize=8, ha='right', va='bottom')
        ax.set(xlabel='t', xlim=(0, 410), ylim=(0, 21));ax.set_title(title, loc='left', fontsize=10)
        ax.legend(fontsize=8, loc='upper right')
    axes[0].set_ylabel('R')
    fig.tight_layout();fig.savefig(output/'leak_radius.png', dpi=160);plt.close(fig)


def fraction_figure(source, output):
    axion0 = lambda h: h['axion_kinetic'][0]+h['axion_gradient'][0]+h['coupling'][0]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4), sharey=True)
    for ax, runs, title in ((axes[0], SCAN, 'escaped axion energy f_far, m_TV scan'),
                            (axes[1], LADDER, 'escaped axion energy f_far, resolution ladder')):
        ax.axvspan(25, 125, color='#e1e0d9', alpha=.5, lw=0)
        ax.text(75, 1.2, 'leak-rate window', color=MUTED, fontsize=8, ha='center')
        series(ax, source, runs, 'axion_far', axion0)
        ax.set(xlabel='t', xlim=(0, 410), yscale='log', ylim=(1e-6, 2));ax.set_title(title, loc='left', fontsize=10)
        ax.legend(fontsize=8, loc='lower right')
    axes[0].set_ylabel('far-field axion energy / initial axion energy')
    fig.tight_layout();fig.savefig(output/'leak_fraction.png', dpi=160);plt.close(fig)


def summary_figure(source, output):
    reports = {name: json.loads((source/name/'summary.json').read_text()) for name in set(SCAN) | set(LADDER)}
    tc = json.loads((source/'report.json').read_text())['control_collapse_time']
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8))
    for name, (label, color) in SCAN.items():
        r = reports[name];x = r['spec']['mass_tv']/OMEGA_MAX
        if r['leak_rate'] is not None and r['leak_rate'] > 0:
            axes[0].plot(x, r['leak_rate'], 'o', ms=8, color=color)
            axes[0].annotate(label, (x, r['leak_rate']), xytext=(6, 4), textcoords='offset points', fontsize=8, color=INK2)
        life = r['lifetime'] if r['lifetime'] is not None else r['completed_until']
        axes[1].plot(x, life/tc, 'o' if not r['censored'] else '^', ms=8, color=color,
                     mfc=color if not r['censored'] else 'none', mew=2)
    axes[0].set(yscale='log', xlabel='m_TV / omega_max', ylabel='leak rate  d f_far / dt')
    axes[0].set_title('leak rate versus trapping margin', loc='left', fontsize=10)
    axes[1].axhline(1, color=MUTED, lw=.8, ls=':')
    axes[1].set(xlabel='m_TV / omega_max', ylabel='lifetime / no-gas collapse time', ylim=(0, None))
    axes[1].set_title('lifetime (open triangle: survived the run)', loc='left', fontsize=10)
    for name, (label, color) in LADDER.items():
        r = reports[name];dx = 64/r['spec']['n']
        if r['leak_rate'] is not None:
            axes[2].plot(dx, r['leak_rate'], 'o', ms=8, color=color)
            axes[2].annotate(label, (dx, r['leak_rate']), xytext=(0, -16), textcoords='offset points', fontsize=8,
                             color=INK2, ha='center')
    axes[2].set(xlabel='grid spacing dx', ylabel='leak rate  d f_far / dt', xlim=(0, .8), ylim=(0, .004))
    axes[2].set_title('leak rate versus resolution (m_TV=2)', loc='left', fontsize=10)
    fig.tight_layout();fig.savefig(output/'leak_summary.png', dpi=160);plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--source', type=Path, required=True);p.add_argument('--output', type=Path, required=True)
    args = p.parse_args();args.output.mkdir(parents=True, exist_ok=False)
    radius_figure(args.source, args.output);fraction_figure(args.source, args.output);summary_figure(args.source, args.output)


if __name__ == '__main__':
    main()
