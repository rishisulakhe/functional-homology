import numpy as np

from find_pairs import (
    build_random_control,
    l2_normalize,
    parse_diamond_hits,
    top_k_cosine,
)


def test_l2_normalize_rows_and_zero_row():
    mat = np.array([[3.0, 4.0], [0.0, 0.0]], dtype=np.float32)
    out = l2_normalize(mat)
    assert np.allclose(np.linalg.norm(out[0]), 1.0)
    assert np.allclose(out[0], [0.6, 0.8])
    assert np.allclose(out[1], [0.0, 0.0])


def test_top_k_cosine_returns_sorted_matches():
    emb_a = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    emb_b = np.array([[1.0, 0.0], [0.9, 0.1], [0.0, 1.0]], dtype=np.float32)
    idx, scores = top_k_cosine(emb_a, emb_b, k=2)
    assert idx[0, 0] == 0
    assert idx[1, 0] == 2
    assert scores[0, 0] >= scores[0, 1]
    assert scores[1, 0] >= scores[1, 1]


def test_top_k_k_larger_than_b():
    emb_a = np.ones((2, 3), dtype=np.float32)
    emb_b = np.ones((1, 3), dtype=np.float32)
    idx, scores = top_k_cosine(emb_a, emb_b, k=5)
    assert idx.shape == (2, 1)
    assert scores.shape == (2, 1)


def test_parse_diamond_hits_keeps_best_and_normalizes(tmp_path):
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "sp|P1|A_ECOLI\tsp|P2|B_BACSU\t40.0\t1e-10\t100\n"
        "sp|P1|A_ECOLI\tsp|P2|B_BACSU\t50.0\t1e-20\t150\n"
        "sp|P1|A_ECOLI\tsp|P3|C_BACSU\t30.0\t1e-2\t20\n"
        "sp|P1|A_ECOLI\tsp|P4|D_BACSU\t30.0\t0.5\t10\n"
    )
    best, stats = parse_diamond_hits(hits, 1e-3, {"P1"}, {"P2", "P3", "P4"})
    assert set(best) == {("P1", "P2")}
    assert best[("P1", "P2")]["pident"] == 50.0
    assert best[("P1", "P2")]["evalue"] == 1e-20
    assert stats["non_significant"] == 2


def test_parse_diamond_hits_skips_unknown_ids(tmp_path):
    hits = tmp_path / "hits.tsv"
    hits.write_text("sp|P1|A_ECOLI\tsp|ZZZ|B_BACSU\t40.0\t1e-10\t100\n")
    best, stats = parse_diamond_hits(hits, 1e-3, {"P1"}, {"P2"})
    assert best == {}
    assert stats["unknown_id"] == 1


def test_build_random_control_deterministic_and_excludes():
    ids_a = [f"A{i}" for i in range(20)]
    ids_b = [f"B{i}" for i in range(20)]
    excluded = {("A0", "B0"), ("A1", "B1")}
    first = build_random_control(ids_a, ids_b, 50, excluded, seed=7)
    second = build_random_control(ids_a, ids_b, 50, excluded, seed=7)
    assert first == second
    assert len(first) == 50
    assert len(set(first)) == 50
    assert not (set(first) & excluded)


def test_build_random_control_zero_size():
    assert build_random_control(["A"], ["B"], 0, set(), seed=1) == []
