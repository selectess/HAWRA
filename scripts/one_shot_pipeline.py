import os
import sys
import argparse
import importlib

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
SRC = os.path.join(ROOT, '03_unified_simulator', 'src')
ARBOl_PKG = os.path.join(ROOT)

def add_src():
    if SRC not in sys.path:
        sys.path.append(SRC)

def compile_arbol(arbol_path):
    if ARBOl_PKG not in sys.path:
        sys.path.append(ARBOl_PKG)
    mod = importlib.import_module('arbol.compiler.compiler')
    # mimic CLI behavior
    sys.argv = ['compiler', arbol_path]
    try:
        mod.main()
    except SystemExit:
        pass
    return arbol_path.replace('.arbol', '.bsim.json')

def run_bsim(bsim_path, output_path):
    add_src()
    Simulator = importlib.import_module('hawra_simulator.simulator').Simulator
    sim = Simulator(bsim_script=bsim_path)
    res = sim.run()
    sim.save_results(output_path, res)
    return output_path

def aggregate(pattern, out_csv):
    # Fallback: spawn subprocess by importing runpy
    import runpy
    script_path = os.path.join(ROOT, '03_unified_simulator', 'aggregate_outputs.py')
    sys.argv = ['aggregate_outputs.py', '--pattern', pattern, '--format', 'csv', '--out', out_csv, '--score', '--sort_by', 'score', '--desc']
    runpy.run_path(script_path, run_name='__main__')

def generate_figures():
    import runpy
    script_path = os.path.join(ROOT, '03_unified_simulator', 'src', 'plot_comparisons.py')
    try:
        globals_dict = runpy.run_path(script_path, run_name='__module__')
        out_dir = os.path.join(ROOT, '03_unified_simulator', 'results')
        func = globals_dict.get('plot_spectral_sweep')
        if func:
            func(out_dir)
    except Exception:
        pass

def generate_plasmid():
    import runpy
    # GenBank
    try:
        runpy.run_path(os.path.join(ROOT, 'scripts', 'create_genbank.py'), run_name='__main__')
    except Exception:
        pass
    # Carte PNG
    try:
        vp = os.path.join(ROOT, '01_genomics', 'genome_analysis_scripts', 'visualize_plasmid.py')
        runpy.run_path(vp, run_name='__main__')
    except Exception:
        pass

def build_bundles():
    import runpy
    for p in ['scripts/build_preprint_bundle.py', 'scripts/build_full_archive.py']:
        try:
            runpy.run_path(os.path.join(ROOT, p), run_name='__main__')
        except Exception:
            pass

def main():
    parser = argparse.ArgumentParser(description='One-shot pipeline: compile .arbol, run bsim, aggregate, figures, bundles')
    parser.add_argument('--arbol', type=str, help='Path to .arbol file')
    parser.add_argument('--bsim', type=str, help='Path to .bsim.json file')
    parser.add_argument('--out', type=str, default=os.path.join(ROOT, '03_unified_simulator', 'results', 'one_shot_output.json'))
    parser.add_argument('--full', action='store_true', help='Run full pipeline with default demo')
    args = parser.parse_args()

    if args.full:
        args.arbol = os.path.join(ROOT, 'arbol', 'phytoqmmml_demo.arbol')

    bsim_path = args.bsim
    if args.arbol:
        bsim_path = compile_arbol(args.arbol)
    if not bsim_path or not os.path.exists(bsim_path):
        print(f'BSIM not found: {bsim_path}')
        sys.exit(1)
    out_json = run_bsim(bsim_path, args.out)
    aggregate(os.path.join(ROOT, '03_unified_simulator', 'output_*.json'), os.path.join(ROOT, '03_unified_simulator', 'results', 'consolidated_metrics.csv'))
    generate_figures()
    generate_plasmid()
    build_bundles()
    print('DONE:', out_json)

if __name__ == '__main__':
    main()
