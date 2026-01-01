
import argparse
from hawra_simulator.simulator import Simulator

def main():
    parser = argparse.ArgumentParser(description='HAWRA Simulator')
    parser.add_argument('--bsim', type=str, required=True, help='Path to the bsim script file')
    parser.add_argument('--output', type=str, required=True, help='Path to the output file')
    args = parser.parse_args()

    simulator = Simulator(bsim_script=args.bsim)
    results = simulator.run()
    simulator.save_results(args.output, results)

if __name__ == '__main__':
    main()
