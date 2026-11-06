import torch

from embed import mean_pool_masked
from fh_common import accession_from_header, fasta_names, read_fasta


def test_accession_from_uniprot_headers():
    assert accession_from_header(">sp|P12345|NAME_ECOLI Description") == "P12345"
    assert accession_from_header(">tr|A0A123|X_BACSU") == "A0A123"
    assert accession_from_header(">plain_id some text") == "plain_id"


def test_mean_pool_excludes_special_tokens():
    hidden = torch.tensor([[[100.0, 100.0], [1.0, 1.0], [3.0, 3.0], [99.0, 99.0]]])
    token_ids = torch.tensor([[0, 5, 6, 2]])
    mask = torch.tensor([[1, 1, 1, 1]])
    pooled = mean_pool_masked(hidden, token_ids, mask, special_ids={0, 2})
    assert torch.allclose(pooled, torch.tensor([[2.0, 2.0]]))


def test_mean_pool_excludes_padding():
    hidden = torch.tensor([[[50.0, 50.0], [1.0, 1.0], [500.0, 500.0], [500.0, 500.0]]])
    token_ids = torch.tensor([[0, 5, 1, 1]])
    mask = torch.tensor([[1, 1, 0, 0]])
    pooled = mean_pool_masked(hidden, token_ids, mask, special_ids={0, 1, 2})
    assert torch.allclose(pooled, torch.tensor([[1.0, 1.0]]))


def test_read_fasta_and_names(tmp_path):
    fa = tmp_path / "x.fasta"
    fa.write_text(
        ">sp|P1|A_ECOLI Alpha protein OS=Escherichia coli\n"
        "MKV\n"
        ">sp|P2|B_BACSU Beta protein\n"
        "MAA\n"
    )
    records = read_fasta(fa)
    assert records == [("P1", "MKV"), ("P2", "MAA")]
    names = fasta_names(fa)
    assert "A_ECOLI" in names["P1"]
    assert "Alpha protein" in names["P1"]
