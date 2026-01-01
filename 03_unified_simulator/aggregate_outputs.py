#!/usr/bin/env python3
import argparse
import glob
import json
import os
import sys

def compute_metrics(records):
    st = [x.get('coherence_factor', 0.0) for x in records]
    p = [x.get('p700_concentration', 0.0) for x in records]
    li = [x.get('light_intensity', 0.0) for x in records]
    return {
        'steps': len(records),
        'coherence_avg': sum(st)/len(st) if st else 0.0,
        'p700_avg': sum(p)/len(p) if p else 0.0,
        'light_avg': sum(li)/len(li) if li else 0.0,
    }

def load_json(path):
    with open(path) as f:
        return json.load(f)

def schedule_to_arbol(light_schedule, target="plant"):
    lines = []
    lines.append("STIMULUS light_pulse ( intensity: 0.0, duration: 0.0 );")
    lines.append("RUN {")
    for i in range(len(light_schedule) - 1):
        t0, i0 = light_schedule[i]
        t1, i1 = light_schedule[i+1]
        if i0 > 0.0 and i1 == 0.0:
            duration = float(t1) - float(t0)
            lines.append(f"    STEP({int(t0)});")
            args = f"intensity: {float(i0)}, duration: {duration}"
            lines.append(f"    APPLY light_pulse ({args}) ON \"{target}\";")
    lines.append("}")
    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser(description='Aggregate output_*.json metrics into a consolidated table or emit Arbol script from BSIM schedule.')
    parser.add_argument('--pattern', type=str, default='03_unified_simulator/output_*.json')
    parser.add_argument('--format', type=str, choices=['json', 'csv'], default='json')
    parser.add_argument('--out', type=str, default='')
    parser.add_argument('--score', action='store_true')
    parser.add_argument('--sort_by', type=str, default='')
    parser.add_argument('--desc', action='store_true')
    parser.add_argument('--max_light', type=float, default=None)
    parser.add_argument('--min_p700', type=float, default=None)
    parser.add_argument('--emit_arbol', action='store_true', help='Emit an Arbol program from a BSIM file schedule')
    parser.add_argument('--bsim', type=str, default='', help='Path to BSIM JSON file containing environment_config.light_schedule')
    parser.add_argument('--target', type=str, default='plant', help='Target identifier for stimulus application')
    args = parser.parse_args()

    # Emit Arbol program from BSIM schedule if requested
    if args.emit_arbol:
        if not args.bsim:
            print(json.dumps({'error': 'Missing --bsim path'}, indent=2))
            sys.exit(1)
        bsim = load_json(args.bsim)
        env = bsim.get('environment_config') or bsim.get('environment') or {}
        sched = env.get('light_schedule', [])
        arbol_text = schedule_to_arbol(sched, target=args.target)
        if args.out:
            with open(args.out, 'w') as f:
                f.write(arbol_text)
        print(arbol_text)
        sys.exit(0)

    files = sorted(glob.glob(args.pattern))
    if not files:
        print(json.dumps({'error': f'No files matched pattern {args.pattern}'}, indent=2))
        sys.exit(1)

    rows = []
    for path in files:
        data = load_json(path)
        # Accept both list logs and dict with 'results'
        records = data['results'] if isinstance(data, dict) and 'results' in data else data
        m = compute_metrics(records)
        row = {'file': os.path.basename(path), **m}
        if args.score:
            row['score'] = row['coherence_avg'] + row['p700_avg'] - 0.01 * row['light_avg']
        rows.append(row)

    if args.max_light is not None:
        rows = [r for r in rows if r.get('light_avg', 0.0) <= args.max_light]
    if args.min_p700 is not None:
        rows = [r for r in rows if r.get('p700_avg', 0.0) >= args.min_p700]

    if args.sort_by:
        rows.sort(key=lambda r: r.get(args.sort_by, 0.0), reverse=args.desc)

    if args.format == 'json':
        report = {'files': rows}
        output = json.dumps(report, indent=2)
    else:
        header = ['file', 'steps', 'coherence_avg', 'p700_avg', 'light_avg']
        if args.score:
            header.append('score')
        lines = [','.join(header)]
        for r in rows:
            line = ','.join(str(r.get(k, '')) for k in header)
            lines.append(line)
        output = '\n'.join(lines)

    print(output)
    if args.out:
        with open(args.out, 'w') as f:
            f.write(output)

if __name__ == '__main__':
    main()
