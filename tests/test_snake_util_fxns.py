from os.path import join as pj
import pandas as pd
import pytest
from workflow.snake_utils.snake_functions import (
    get_lines_to_skip, picard_calculate_strandedness, concat_picard_insert_size,
    concat_star_log, comb_filepaths, make_fp_dict, specified_strandedness
)

## TO RUN: pytest -v tests/test_snake_util_fxns.py

## tests for picard_calculate_strandedness() function
@pytest.fixture
def make_picard_file(tmp_path):
    def _make(name, r1, r2):
        content = (
            "## header line\n"
            "## METRICS CLASS picard.analysis.RnaSeqMetrics\n"
            "PF_BASES\tPCT_R1_TRANSCRIPT_STRAND_READS\tPCT_R2_TRANSCRIPT_STRAND_READS\n"
            f"1000000\t{r1}\t{r2}\n"
        )
        file_path = pj(tmp_path, f"{name}.picard.metrics.txt")
        with open(file_path, "w") as f:
            f.write(content)
        return str(file_path)
    ## this allows me to return the function in the fixture to be used in my tests
    return _make

@pytest.mark.parametrize(
    "name,r1,r2,expected_strandedness,expected_key",
    [
        ("sample1", 0.701, 0.30, "R1", 1),
        ("sample2", 0.30, 0.701, "R2", 0),
        ("sample3", 0.69, 0.31, "unstranded", 0.5),
    ]
)
def test_picard_calculate_strandedness_variants(tmp_path, make_picard_file, name, r1, r2, expected_strandedness, expected_key):
    make_picard_file(name, r1, r2)
    out_file = pj(tmp_path, "output.tsv")
    pattern = str(pj(tmp_path, "*.picard.metrics.txt"))
    picard_calculate_strandedness(file_pattern=pattern, out_file=str(out_file))
    df = pd.read_csv(out_file, sep="\t")
    row = df[df["sampleid"] == name].iloc[0]
    assert row["strandedness"] == expected_strandedness
    assert row["rsem_strand_key"] == expected_key

def test_picard_calculate_strandedness_multiple_files(tmp_path, make_picard_file):
    for name, r1, r2 in [("sample1", 0.701, 0.30), ("sample2", 0.30, 0.701), ("sample3", 0.69, 0.31)]:
        make_picard_file(name, r1, r2)
    out_file = pj(tmp_path, "output.tsv")
    pattern = str(pj(tmp_path, "*.picard.metrics.txt"))
    picard_calculate_strandedness(file_pattern=pattern, out_file=str(out_file))
    df = pd.read_csv(out_file, sep="\t")
    assert set(df["sampleid"]) == {"sample1", "sample2", "sample3"}
    assert len(df) == 3

def test_get_lines_to_skip_not_found(tmp_path):
    file_path = pj(tmp_path, "empty.txt")
    with open(file_path, "w") as f:
        f.write("no metrics here\n")
    with pytest.raises(ValueError):
        get_lines_to_skip(file_path, r"^## METRICS")

@pytest.mark.parametrize("content,exc", [
    ("## header line\n## METRICS CLASS picard.analysis.RnaSeqMetrics\nPF_BASES\n1000000\n", KeyError),
    ("", ValueError),
    ("## header line\n## METRICS CLASS picard.analysis.RnaSeqMetrics\nPF_BASES\tPCT_R1_TRANSCRIPT_STRAND_READS\tPCT_R2_TRANSCRIPT_STRAND_READS\n1000000\tfoo\tbar\n", Exception),
])
def test_picard_calculate_strandedness_edge_cases(tmp_path, content, exc):
    file_path = pj(tmp_path, "bad.picard.metrics.txt")
    with open(file_path, "w") as f:
        f.write(content)
    out_file = pj(tmp_path, "output.tsv")
    pattern = str(pj(tmp_path, "*.picard.metrics.txt"))
    with pytest.raises(exc):
        picard_calculate_strandedness(file_pattern=pattern, out_file=str(out_file))

## test for concat_picard_insert_size() function
@pytest.fixture
def picard_insert_size_file(tmp_path):
    content = (
        "## header line\n"
        "## METRICS CLASS picard.analysis.InsertSizeMetrics\n"
        "MEDIAN_INSERT_SIZE\tWIDTH_OF_10_PERCENT\tWIDTH_OF_20_PERCENT\n"
        "150\t20\t40\n"
    )
    file_path = pj(tmp_path, "sample1.picard.insertSize.txt")
    with open(file_path, "w") as f:
        f.write(content)
    return str(file_path)

def test_concat_picard_insert_size_basic(tmp_path, picard_insert_size_file):
    out_file = pj(tmp_path, "insert_size_out.tsv")
    pattern = str(pj(tmp_path, "*.picard.insertSize.txt"))
    concat_picard_insert_size(file_pattern=pattern, out_file=str(out_file))
    df = pd.read_csv(out_file, sep="\t")
    assert "sampleid" in df.columns
    assert df.loc[0, "sampleid"] == "sample1"
    assert df.loc[0, "MEDIAN_INSERT_SIZE"] == 150

## test for concat_star_log() function
@pytest.fixture
def star_log_file(tmp_path):
    content = (
        "Number of input reads |	1000000\n"
        "Uniquely mapped reads number |	900000\n"
        "Uniquely mapped reads % |	90%\n"
    )
    file_path = pj(tmp_path, "sample1.Log.final.out")
    with open(file_path, "w") as f:
        f.write(content)
    return str(file_path)

def test_concat_star_log_basic(tmp_path, star_log_file):
    out_file = pj(tmp_path, "star_log_out.tsv")
    pattern = str(pj(tmp_path, "*.Log.final.out"))
    concat_star_log(file_pattern=pattern, out_file=str(out_file))
    df = pd.read_csv(out_file, sep="\t")
    assert "sampleid" in df.columns
    assert df.loc[0, "sampleid"] == "sample1"
    assert df.loc[0, "number_of_input_reads"] == 1000000
    assert df.loc[0, "uniquely_mapped_reads_number"] == 900000
    assert df.loc[0, "uniquely_mapped_reads_pct"] == 90

## testing combining filepaths function
def test_comb_filepaths_basic():
    assert comb_filepaths("/dir", "file.txt") == "/dir/file.txt"
    assert comb_filepaths("dir1", "dir2/file.txt") == "dir1/dir2/file.txt"

## testing sample fp dictionary creation 
def test_make_fp_dict_basic(tmp_path):
    df = pd.DataFrame({
        'sampleid': ['sampleA', 'sampleB'],
        'forward': ['A_R1.fq.gz', 'B_R1.fq.gz'],
        'reverse': ['A_R2.fq.gz', 'B_R2.fq.gz']
    })
    dataset_dir = str(tmp_path)
    result = make_fp_dict(df, dataset_dir)
    assert set(result.keys()) == {'sampleA', 'sampleB'}
    assert result['sampleA'] == [f"{dataset_dir}/A_R1.fq.gz", f"{dataset_dir}/A_R2.fq.gz"]
    assert result['sampleB'] == [f"{dataset_dir}/B_R1.fq.gz", f"{dataset_dir}/B_R2.fq.gz"]

@pytest.mark.parametrize(
    "sample,expected_key",
    [("s1", 1), ("s2", 0), ("s3", 0.5)]
)

def test_specified_strandedness_rsem_key(sample, expected_key):
    df = pd.DataFrame({'sampleid': [sample], 'strandedness': [expected_key]})
    assert specified_strandedness(df, sample).iloc[0] == expected_key

def test_specified_strandedness_missing():
    df = pd.DataFrame({'sampleid': ['s1', 's2']})
    assert specified_strandedness(df, 's1') == ""
