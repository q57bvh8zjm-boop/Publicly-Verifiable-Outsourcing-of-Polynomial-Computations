#!/usr/bin/env python3
"""Export measured CSVs and Figures 2--5 using the current manuscript numbering."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
from run_vm_experiments import jobs_for
from vm_protocol import SCHEMES, save_json


def write_csv(path, rows):
    if rows:
        with path.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def summarize(run):
    summary = {}
    for path in sorted((run / 'jobs').glob('*/status.json')):
        if json.loads(path.read_text()).get('status') != 'complete':
            continue
        with (path.parent / 'raw.csv').open() as f:
            rows = list(csv.DictReader(f))
        if not rows or any(r['correct'] != '1' or r['mode'] != 'tcp_full_protocol' for r in rows):
            raise ValueError(f'invalid measurement set: {path.parent}')
        head = rows[0]
        key = head['scheme'], int(head['m']), int(head['d']), int(head['t'])
        values = {'scheme': key[0], 'm': key[1], 'd': key[2], 't': key[3], 'samples': len(rows)}
        for name in ('client_ms', 'native_ms', 'probgen_ms', 'verify_ms', 'decode_ms', 'server_max_ms',
                     'server_sum_ms', 'tcp_roundtrip_ms', 'online_wall_ms', 'keygen_ms', 'remote_setup_ms'):
            samples = [float(r[name]) for r in rows]
            values[name + '_mean'] = statistics.mean(samples)
            values[name + '_median'] = statistics.median(samples)
            values[name + '_stdev'] = statistics.stdev(samples) if len(samples) > 1 else 0
        values['client_ratio'] = values['client_ms_mean'] / values['native_ms_mean']
        values['log10_ratio'] = math.log10(values['client_ratio'])
        summary[key] = values
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('run', type=Path)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--figures', nargs='+', choices=['2', '3', '4', '5'], default=['2', '3', '4', '5'])
    ap.add_argument('--allow-partial', action='store_true')
    ap.add_argument('--summary-only', action='store_true')
    args = ap.parse_args()
    data = summarize(args.run)
    if not data:
        ap.error('no completed measured jobs')
    args.output.mkdir(parents=True, exist_ok=True)
    write_csv(args.output / 'measurements.csv', list(data.values()))
    if args.summary_only:
        return
    coverage = {f: [list(j) for j in jobs_for('figure' + f, SCHEMES) if j not in data] for f in args.figures}
    save_json(args.output / 'coverage.json', {f: {'complete': not missing, 'missing': missing}
                                            for f, missing in coverage.items()})
    if any(coverage.values()) and not args.allow_partial:
        ap.error('selected figures have incomplete measured grids; see coverage.json or use --allow-partial')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    labels = {'Gamma_FE_Sh': r'$\Gamma_{\mathrm{FE}}^{\mathrm{Sh}}$',
              'Gamma_OV_Sh': r'$\Gamma_{\mathrm{OV}}^{\mathrm{Sh}}$',
              'Gamma_OV_HSS': r'$\Gamma_{\mathrm{OV}}^{\mathrm{HSS}}$'}
    styles = {'Gamma_FE_Sh': ('#333333', 's'), 'Gamma_OV_Sh': ('#c73838', 'o'),
              'Gamma_OV_HSS': ('#2463a7', '^')}

    def series(ax, s, d, t, ms):
        ys = [data[(s, m, d, t)]['log10_ratio'] if (s, m, d, t) in data else math.nan for m in ms]
        color, marker = styles[s]
        ax.plot(list(ms), ys, label=labels[s] + (f', t={t}' if args_current == '3' else ''),
                color=color, marker=marker, linestyle='--' if t == 3 else '-')

    for args_current in args.figures:
        if args_current == '5':
            fig, axes = plt.subplots(1, 3, figsize=(13, 3.7), sharey=True)
            for ax, d in zip(axes, [4, 8, 10]):
                rows = []
                for m in range(1, 16):
                    row = {'m': m}
                    for s, short in zip(SCHEMES[:3], ['fe_sh', 'ov_sh', 'ov_hss']):
                        value = data.get((s, m, d, 1), {})
                        for key in ('client_ms_mean', 'native_ms_mean', 'client_ratio', 'log10_ratio'):
                            row[short + '_' + key.replace('client_ratio', 'ratio')] = value.get(key, math.nan)
                    rows.append(row)
                write_csv(args.output / f'figure5_d{d}.csv', rows)
                for s in SCHEMES[:3]:
                    series(ax, s, d, 1, range(1, 16))
                ax.set_title(f'd={d}, t=1')
        else:
            fig, ax = plt.subplots(figsize=(6, 4))
            axes = [ax]
            ms = range(200, 2001, 200)
            if args_current == '4':
                rows = []
                for d in (2, 3):
                    for t in (1, 2):
                        ys = []
                        for m in ms:
                            fe, pi = data.get(('Gamma_FE_Sh', m, d, t)), data.get(('Pi3', m, d, t))
                            ratio = fe['client_ms_mean'] / pi['client_ms_mean'] if fe and pi else math.nan
                            ys.append(ratio)
                            rows.append({'m': m, 'd': d, 't': t, 'client_ratio_fe_over_pi3': ratio})
                        ax.plot(list(ms), ys, marker='o', label=f'd={d}, t={t}')
                write_csv(args.output / 'figure4_fe_vs_pi3.csv', rows)
            else:
                names = SCHEMES[:3] if args_current == '2' else SCHEMES[:2]
                thresholds = [1] if args_current == '2' else [2, 3]
                for s in names:
                    for t in thresholds:
                        series(ax, s, 2, t, ms)
                selected = [v for k, v in data.items() if k in jobs_for('figure' + args_current, SCHEMES)]
                write_csv(args.output / f'figure{args_current}.csv', selected)
        for ax in axes:
            ax.set_xlabel('m')
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
        axes[0].set_ylabel(r'$T_c(\Gamma_{\mathrm{FE}}^{\mathrm{Sh}})/T_c(\Pi_3)$' if args_current == '4'
                           else r'$\log_{10}(T_c/T_n)$')
        fig.suptitle(f'Figure {args_current}' + (' — partial measured grid' if coverage[args_current] else ''))
        fig.tight_layout()
        fig.savefig(args.output / f'figure{args_current}.pdf')
        fig.savefig(args.output / f'figure{args_current}.png', dpi=180)
        plt.close(fig)


if __name__ == '__main__':
    main()
