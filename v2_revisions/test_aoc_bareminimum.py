from aoc import main
from click.testing import CliRunner

def test_aoc_runs():
    """Checking if AOC runs under bare minimum flags"""
    runner = CliRunner()
    result1 = runner.invoke(main, ['-c', 'example_data/simple/clusters.txt', '-g', 'example_data/simple/graph.txt', '--inclusion-criterion', 'k', '--candidate-file', 'example_data/simple/candidates.txt', '--min-k-core', '1', '-o', '../scratch/new_cluster.txt'])
    assert result1.exit_code == 0
    result2 = runner.invoke(main, ['-c', 'example_data/simple/clusters.txt', '-g', 'example_data/simple/graph.txt', '--inclusion-criterion', 'mcd', '--candidate-file', 'example_data/simple/candidates.txt', '-o', '../scratch/new_cluster.txt'])
    assert result2.exit_code == 0