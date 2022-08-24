from aoc import main
from click.testing import CliRunner


def test_aoc_runs():
    """Checking if AOC runs under bare minimum flags"""
    runner = CliRunner()
    clusters_path = "example_data/simple/clusters.txt"
    network_path = "example_data/simple/graph.txt"
    candidates_path = "example_data/simple/candidates.txt"
    aoc_k_with_candidates = runner.invoke(
        main,
        [
            "-c",
            clusters_path,
            "-g",
            network_path,
            "--inclusion-criterion",
            "k",
            "--candidate-file",
            candidates_path,
            "--min-k-core",
            "1",
            "-o",
            "../scratch/new_cluster.txt",
        ],
    )
    assert aoc_k_with_candidates.exit_code == 0
    aoc_m_with_candidates = runner.invoke(
        main,
        [
            "-c",
            clusters_path,
            "-g",
            network_path,
            "--inclusion-criterion",
            "mcd",
            "--candidate-file",
            candidates_path,
            "-o",
            "../scratch/new_cluster.txt",
        ],
    )
    assert aoc_m_with_candidates.exit_code == 0
    aoc_k_wo_candidates = runner.invoke(
        main,
        [
            "-c",
            clusters_path,
            "-g",
            network_path,
            "--inclusion-criterion",
            "k",
            "--min-k-core",
            "1",
            "-o",
            "../scratch/new_cluster.txt",
        ],
    )
    assert aoc_k_wo_candidates.exit_code == 0
    aoc_k_bare_minimum = runner.invoke(
        main,
        [
            "-c",
            clusters_path,
            "-g",
            network_path,
            "--inclusion-criterion",
            "k",
            "--min-k-core",
            "1",
        ],
    )
    assert aoc_k_bare_minimum.exit_code == 0
