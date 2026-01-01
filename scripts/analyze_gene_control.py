import json
import matplotlib.pyplot as plt

def analyze_gene_expression(results_path, output_path):
    with open(results_path, 'r') as f:
        data = json.load(f)

    times = [i for i, _ in enumerate(data['results'])]
    gene_a_expression = [entry['genes']['gene_a']['expression_level'] for entry in data['results']]

    plt.figure(figsize=(10, 6))
    plt.plot(times, gene_a_expression, label='Gene A Expression')
    plt.axvline(x=100, color='r', linestyle='--', label='Inhibitor Applied')
    plt.xlabel('Time')
    plt.ylabel('Expression Level')
    plt.title('Gene A Expression Over Time')
    plt.legend()
    plt.grid(True)
    plt.savefig(output_path)
    print(f"Graphique de l'expression génique sauvegardé dans {output_path}")

if __name__ == "__main__":
    analyze_gene_expression('gene_control_results.json', 'gene_a_expression.png')
