import csv
import json
from pathlib import Path

from piranha.analysis import clean_gaps, get_co_occurrence, get_haplotypes, variation_functions


def test_remove_gaps_replaces_deletions_with_n(tmp_path: Path):
    aln = tmp_path / "aln.fasta"
    aln.write_text(">ref\n--ATGC--\n>cns1\n--AT-C--\n")

    cns_id, cns = clean_gaps.remove_gaps(str(aln))
    assert cns_id == "cns1"
    assert cns == "ATNC"


def test_clean_cns_mask_handles_snp_insertion_deletion():
    in_cns = "ATGCCCAA"
    to_mask = {
        2: "3:CT",      # SNP mask
        4: "5:ins2",    # skip two bases
        7: "8:del2",    # mask two bases
    }

    masked = clean_gaps.clean_cns_mask(in_cns, to_mask)
    assert len(masked) <= len(in_cns)
    assert "N" in masked


def test_clean_medaka_cns_writes_output_and_returns_mask_dict(tmp_path: Path):
    aln = tmp_path / "medaka_aln.fasta"
    out = tmp_path / "cleaned.fasta"
    aln.write_text(">ref\nATGC\n>cns\nATGC\n")

    to_mask = clean_gaps.clean_medaka_cns("sample1", str(aln), str(out))
    assert isinstance(to_mask, dict)
    rendered = out.read_text()
    assert rendered.startswith(">sample1")
    assert "ATGC" in rendered


def test_get_combinations_counts_variant_patterns(tmp_path: Path):
    reads = tmp_path / "reads.fasta"
    reads.write_text(">r1\nATGC\n>r2\nATGT\n>r3\nATGC\n")

    combos = get_co_occurrence.get_combinations(
        variants="2:AG;4:CT",
        read_fasta_file=str(reads),
        reference="RefA",
        barcode="barcode01",
        threshold=1,
    )

    assert isinstance(combos, dict)
    assert any(key.startswith("2") for key in combos)


def test_get_variation_pcent_reports_site_percentages(tmp_path: Path):
    ref = tmp_path / "ref.fasta"
    reads = tmp_path / "reads.fasta"
    ref.write_text(">ref\nATGC\n")
    reads.write_text(">r1\nATGC\n>r2\nATGT\n")

    info = get_haplotypes.get_variation_pcent(str(ref), str(reads))
    site4 = next(row for row in info if row["Position"] == 4)
    assert site4["Percentage"] == 50.0


def test_parse_vcf_and_write_haplotype_outputs(tmp_path: Path):
    reads_fasta = tmp_path / "reads.fasta"
    reads_fastq = tmp_path / "reads.fastq"
    ref = tmp_path / "ref.fasta"
    vcf = tmp_path / "vars.vcf"
    hap_csv = tmp_path / "haps.csv"
    outdir = tmp_path / "hap_out"

    ref.write_text(">ref\nATGC\n")
    reads_fasta.write_text(">r1\nATGC\n>r2\nATGT\n>r3\nATGT\n")
    reads_fastq.write_text("@r1\nATGC\n+\n!!!!\n@r2\nATGT\n+\n!!!!\n@r3\nATGT\n+\n!!!!\n")
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "RefA\t4\t.\tC\tT\t50\tPASS\t.\tGT:GQ\t1:50\n"
    )

    haps = get_haplotypes.get_haplotypes(
        fasta=str(reads_fasta),
        vcf=str(vcf),
        reads=str(reads_fastq),
        ref=str(ref),
        out_haplotypes=str(hap_csv),
        outdir=str(outdir),
        taxon="RefA",
        min_reads=1,
        min_pcent=10,
    )

    assert haps
    assert hap_csv.exists()
    assert any(path.name.endswith(".fastq") for path in outdir.iterdir())
    assert any(path.name.endswith(".reference.fasta") for path in outdir.iterdir())


def test_gather_haplotype_data_writes_combined_csv_and_config(tmp_path: Path):
    hap1 = tmp_path / "h1.csv"
    out = tmp_path / "combined.csv"
    cfg = tmp_path / "config.yaml"

    hap1.write_text("taxon,sites,haplotype,num_reads,make_cns\nRefA,4,4T,5,True\n")

    get_haplotypes.gather_haplotype_data([str(hap1)], str(out), str(cfg), {"runname": "x"})
    assert out.exists()
    assert cfg.exists()


def test_parse_variant_file_and_non_ref_percent(tmp_path: Path):
    var_csv = tmp_path / "vars.csv"
    var_csv.write_text("reference,variants\nRefA,4:CT;6:ins2\n")

    var_dict = variation_functions.parse_variant_file(str(var_csv))
    assert 4 in var_dict["RefA"]
    assert 6 in var_dict["RefA"]

    pileup = {"Position": 4, "A reads": 0, "C reads": 8, "T reads": 2, "G reads": 0, "- reads": 0}
    ref_dict = {3: "C"}
    pcent, total = variation_functions.non_ref_prcnt_calc(3, pileup, ref_dict)
    assert total == 10
    assert pcent == 20.0


def test_calculate_coocc_json_returns_records():
    var_dict = {
        1: ["A", "G"],
        2: ["T", "C"],
    }
    read_vars = {
        "r1": {1: "A", 2: "T"},
        "r2": {1: "G", 2: "C"},
    }

    coocc = variation_functions.calculate_coocc_json(var_dict, read_vars)
    assert isinstance(coocc, list)
    assert coocc
    assert {"SNP1", "SNP2", "Alt", "Ref", "Total"}.issubset(set(coocc[0].keys()))
