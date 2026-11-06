import math

from fh_common import ancestor_closure, parse_obo
from go_eval import (
    build_ic,
    jaccard,
    read_annotations,
    resnik_similarity,
)


OBO = """format-version: 1.2

[Term]
id: GO:0000001
name: root

[Term]
id: GO:0000002
name: child
is_a: GO:0000001 ! root

[Term]
id: GO:0000003
name: grandchild
is_a: GO:0000002 ! root

[Term]
id: GO:0000004
name: part child
relationship: part_of GO:0000002 ! child
"""


def write_obo(tmp_path):
    path = tmp_path / "mini.obo"
    path.write_text(OBO)
    return path


def test_parse_obo_and_ancestor_closure(tmp_path):
    blocks, parents = parse_obo(write_obo(tmp_path))
    assert set(blocks) == {"GO:0000001", "GO:0000002", "GO:0000003", "GO:0000004"}
    assert parents["GO:0000003"] == {"GO:0000002"}
    assert parents["GO:0000004"] == {"GO:0000002"}
    assert ancestor_closure({"GO:0000003"}, parents) == {
        "GO:0000001", "GO:0000002", "GO:0000003"}
    assert ancestor_closure({"GO:0000004"}, parents) == {
        "GO:0000001", "GO:0000002", "GO:0000004"}


def test_jaccard():
    assert jaccard({"GO:1", "GO:2"}, {"GO:2", "GO:3"}) == 1 / 3
    assert jaccard(set(), set()) == 0.0
    assert jaccard({"GO:1"}, {"GO:1"}) == 1.0


def test_build_ic_and_resnik():
    propagated = {
        "p1": {"GO:0000001", "GO:0000002"},
        "p2": {"GO:0000001"},
        "p3": {"GO:0000001", "GO:0000002", "GO:0000003"},
    }
    ic = build_ic(propagated)
    assert ic["GO:0000001"] == 0.0
    assert math.isclose(ic["GO:0000002"], math.log(3 / 2))
    assert math.isclose(ic["GO:0000003"], math.log(3))
    sim = resnik_similarity({"GO:0000002", "GO:0000003"}, {"GO:0000002"}, ic)
    assert math.isclose(sim, math.log(3 / 2))
    assert resnik_similarity({"GO:0000003"}, set(), ic) == 0.0


def test_read_annotations_filters_evidence_and_not(tmp_path):
    ann = tmp_path / "ann.tsv"
    ann.write_text(
        "accession\tgo_id\tevidence_code\tqualifier\n"
        "P1\tGO:0000001\tIDA\t\n"
        "P1\tGO:0000001\tIEA\t\n"
        "P2\tGO:0000002\tIDA\tNOT|involved_in\n"
        "P2\tGO:0000003\tIMP\t\n"
    )
    out, stats = read_annotations(ann)
    assert out == {"P1": {"GO:0000001"}, "P2": {"GO:0000003"}}
    assert stats["bad_evidence"] == 1
    assert stats["not_qualifier"] == 1
