import csv
import json
from pathlib import Path

from piranha.analysis import haplo_functions, haplotyping_functions, phylo_functions
from piranha.utils.config import (
    KEY_BARCODE,
    KEY_LOCAL_DATABASE_THRESHOLD,
    KEY_PHYLO_METADATA_COLUMNS,
    KEY_REFERENCE_GROUP,
    KEY_REFERENCE_GROUP_FIELD,
    KEY_REFERENCES_FOR_CNS,
    KEY_SUPPLEMENTARY_METADATA_COLUMNS,
    KEY_SUPPLEMENTARY_METADATA_ID_COLUMN,
    KEY_TREE_ANNOTATIONS,
    KEY_VARIANT_COUNT,
)


def test_parse_partition_file_and_support_col(tmp_path: Path):
    partition = tmp_path / "parts.txt"
    partition.write_text("#0\nr1\tX\n#1\nr2\tX\n")

    parsed = haplo_functions.parse_partition_file(str(partition))
    assert parsed[0] == {"r1"}
    assert parsed[1] == {"r2"}

    total, support = haplo_functions.parse_support_col("0:3|1:2")
    assert total == 5
    assert support == {"0": 3, "1": 2}



def test_haplo_sum_base_support_and_sd():
    merged = haplo_functions.sum_base_support({"0": 2}, {"0": 3, "1": 4})
    assert merged == {"0": 5, "1": 4}

    haps = [
        {1: {"assigned allele": "1", "base support counts": {"1": 10}}, 2: {"assigned allele": "1", "base support counts": {"1": 10}}},
    ]
    even = haplo_functions.calc_sd(haps)
    assert 1 in even
    assert even[1] == 0.0



def test_haplo_parse_vcf_and_flopp_reads(tmp_path: Path):
    vcf = tmp_path / "vars.vcf"
    flopp = tmp_path / "flopp.txt"

    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "Ref\t5\t.\tC\tT\t50\tPASS\tRO=8;AO=2\n"
    )
    flopp.write_text(
        "header\n"
        "1:5\t1\t0:3|1:2\n"
        "*****\n"
    )

    parsed_vcf = haplo_functions.parse_VCF(str(vcf))
    assert parsed_vcf[5]["0"] == "C"
    assert parsed_vcf[5]["1"] == "T"

    parsed_flopp = haplo_functions.parse_flopp_reads(str(flopp))
    assert len(parsed_flopp) == 1
    assert parsed_flopp[0][5]["assigned allele"] == "1"



def test_haplo_collapse_close_empty_flopp_returns_single_partition(tmp_path: Path):
    vcf = tmp_path / "vars.vcf"
    flopp = tmp_path / "empty_flopp.txt"
    vcf.write_text("##fileformat=VCFv4.2\n")
    flopp.write_text("header\n")

    collapsed = haplo_functions.collapse_close(str(flopp), 1, str(vcf))
    assert collapsed == [[0]]



def test_haplotyping_parse_vcf_calls_and_support(tmp_path: Path):
    v1 = tmp_path / "h1.vcf"
    v1.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "Ref\t5\t.\tC\tT\t50\tPASS\tRO=8;AO=2\n"
        "Ref\t8\t.\tAA\tTT\t50\tPASS\tRO=8;AO=2\n"
    )

    calls = haplotyping_functions.parseVCFCalls([str(v1)])
    assert len(calls) == 1
    assert (5, "T") in calls[0]

    total, support = haplotyping_functions.parseSupportCol("0:2|1:3")
    assert total == 5
    assert support == {"0": 2, "1": 3}



def test_haplotyping_sum_and_sd():
    merged = haplotyping_functions.sumBaseSupport({"0": 1}, {"0": 2})
    assert merged == {"0": 3}

    haps = [{1: {"assigned allele": "1", "base support counts": {"1": 10}}}]
    even = haplotyping_functions.calcSD(haps)
    assert even[1] == 0.0



def test_phylo_get_seqs_and_clusters_minimal(tmp_path: Path):
    sample_seqs = tmp_path / "sample.fasta"
    sample_info = tmp_path / "sample_info.json"
    supplementary_sequences = tmp_path / "supp.fasta"
    reference_sequences = tmp_path / "refs.fasta"
    outgroup_sequences = tmp_path / "outgroups.fasta"
    barcodes = tmp_path / "barcodes.csv"
    phylo_outdir = tmp_path / "phylo"
    phylo_outdir.mkdir()

    record_id = "S1|Sabin2-related|h1|E1|2026-03-26"
    sample_seqs.write_text(f">{record_id}\nATGC\n")
    sample_info.write_text(json.dumps({record_id: {KEY_BARCODE: "barcode01", KEY_VARIANT_COUNT: "1"}}))
    supplementary_sequences.write_text("")
    reference_sequences.write_text(">RefA reference_group=Sabin2-related\nATGC\n")
    outgroup_sequences.write_text(">Out reference_group=Sabin2-related\nATGC\n")
    barcodes.write_text("sample,barcode\nS1,barcode01\n")

    config = {
        KEY_REFERENCES_FOR_CNS: ["Sabin2-related"],
        KEY_SUPPLEMENTARY_METADATA_COLUMNS: [],
        KEY_SUPPLEMENTARY_METADATA_ID_COLUMN: "name",
        KEY_PHYLO_METADATA_COLUMNS: ["sample"],
        KEY_TREE_ANNOTATIONS: "",
    }

    clusters, tree_annotations = phylo_functions.get_seqs_and_clusters(
        sample_seqs=str(sample_seqs),
        sample_info=str(sample_info),
        supplementary_sequences=str(supplementary_sequences),
        reference_sequences=str(reference_sequences),
        outgroup_sequences=str(outgroup_sequences),
        barcodes_csv=str(barcodes),
        supplementary_metadata="",
        phylo_outdir=str(phylo_outdir),
        config=config,
    )

    assert "Sabin2-related" in clusters
    assert (phylo_outdir / "Sabin2-related.fasta").exists()
    assert (phylo_outdir / "annotations.csv").exists()
    assert "reference_group" in tree_annotations
    assert "call" in tree_annotations


def test_phylo_get_seqs_and_clusters_passes_metadata_fields(tmp_path: Path):
    sample_seqs = tmp_path / "sample.fasta"
    sample_info = tmp_path / "sample_info.json"
    supplementary_sequences = tmp_path / "supp.fasta"
    supplementary_metadata = tmp_path / "supp.csv"
    reference_sequences = tmp_path / "refs.fasta"
    outgroup_sequences = tmp_path / "outgroups.fasta"
    barcodes = tmp_path / "barcodes.csv"
    phylo_outdir = tmp_path / "phylo"
    phylo_outdir.mkdir()

    sample_id = "S1|Sabin2-related|h1|E1|2026-03-26"
    sample_seqs.write_text(f">{sample_id}\nATGC\n")
    sample_info.write_text(json.dumps({sample_id: {KEY_BARCODE: "barcode01", KEY_VARIANT_COUNT: "1"}}))
    supplementary_sequences.write_text(">BG1 reference_group=Sabin2-related\nATGC\n")
    supplementary_metadata.write_text("name,country\nBG1,UK\n")
    reference_sequences.write_text(">RefA reference_group=Sabin2-related\nATGC\n")
    outgroup_sequences.write_text(">Out reference_group=Sabin2-related\nATGC\n")
    barcodes.write_text("sample,barcode,institute\nS1,barcode01,NIC\n")

    config = {
        KEY_REFERENCES_FOR_CNS: ["Sabin2-related"],
        KEY_SUPPLEMENTARY_METADATA_COLUMNS: ["country"],
        KEY_SUPPLEMENTARY_METADATA_ID_COLUMN: "name",
        KEY_PHYLO_METADATA_COLUMNS: ["institute"],
        KEY_TREE_ANNOTATIONS: "",
    }

    phylo_functions.get_seqs_and_clusters(
        sample_seqs=str(sample_seqs),
        sample_info=str(sample_info),
        supplementary_sequences=str(supplementary_sequences),
        reference_sequences=str(reference_sequences),
        outgroup_sequences=str(outgroup_sequences),
        barcodes_csv=str(barcodes),
        supplementary_metadata=str(supplementary_metadata),
        phylo_outdir=str(phylo_outdir),
        config=config,
    )

    with open(phylo_outdir / "annotations.csv") as handle:
        rows = list(csv.DictReader(handle))

    sample_row = next(row for row in rows if row["name"] == sample_id)
    bg_row = next(row for row in rows if row["name"] == "BG1")
    assert sample_row["institute"] == "NIC"
    assert bg_row["country"] == "UK"



def test_phylo_update_local_database_filters_sabin_threshold(tmp_path: Path):
    sample_sequences = tmp_path / "sample.fasta"
    sample_info = tmp_path / "sample_info.json"
    detailed_csv = tmp_path / "detailed.csv"
    new_db_seqs = tmp_path / "new.fasta"
    new_db_metadata = tmp_path / "new.csv"

    record_keep = "S_keep|Sabin2-related|h1|E1|2026-03-26"
    record_drop = "S_drop|Sabin2-related|h1|E2|2026-03-26"
    sample_sequences.write_text(
        f">{record_keep} barcode=barcode01\nATGC\n>{record_drop} barcode=barcode02\nATGC\n"
    )
    sample_info.write_text(
        json.dumps(
            {
                record_keep: {KEY_VARIANT_COUNT: "8", "ddns_group": "Sabin2-related"},
                record_drop: {KEY_VARIANT_COUNT: "1", "ddns_group": "Sabin2-related"},
            }
        )
    )
    detailed_csv.write_text("sample,barcode\nS_keep,barcode01\nS_drop,barcode02\n")

    config = {
        KEY_REFERENCE_GROUP_FIELD: "ddns_group",
        KEY_LOCAL_DATABASE_THRESHOLD: 5,
        KEY_SUPPLEMENTARY_METADATA_ID_COLUMN: "name",
    }

    phylo_functions.update_local_database(
        sample_sequences=str(sample_sequences),
        sample_info=str(sample_info),
        detailed_csv=str(detailed_csv),
        new_db_seqs=str(new_db_seqs),
        new_db_metadata=str(new_db_metadata),
        config=config,
    )

    seqs = new_db_seqs.read_text()
    assert "S_keep" in seqs
    assert "S_drop" not in seqs
    assert new_db_metadata.exists()


def test_phylo_update_local_database_custom_group_field(tmp_path: Path):
    sample_sequences = tmp_path / "sample_custom.fasta"
    sample_info = tmp_path / "sample_custom.json"
    detailed_csv = tmp_path / "detailed_custom.csv"
    new_db_seqs = tmp_path / "new_custom.fasta"
    new_db_metadata = tmp_path / "new_custom.csv"

    record_keep = "S_custom|Sabin2-related|h1|E1|2026-03-26"
    sample_sequences.write_text(f">{record_keep} barcode=barcode03\nATGC\n")
    sample_info.write_text(
        json.dumps(
            {
                record_keep: {
                    KEY_VARIANT_COUNT: "9",
                    "my_group": "Sabin2-related",
                }
            }
        )
    )
    detailed_csv.write_text("sample,barcode\nS_custom,barcode03\n")

    config = {
        KEY_REFERENCE_GROUP_FIELD: "my_group",
        KEY_LOCAL_DATABASE_THRESHOLD: 5,
        KEY_SUPPLEMENTARY_METADATA_ID_COLUMN: "name",
    }

    phylo_functions.update_local_database(
        sample_sequences=str(sample_sequences),
        sample_info=str(sample_info),
        detailed_csv=str(detailed_csv),
        new_db_seqs=str(new_db_seqs),
        new_db_metadata=str(new_db_metadata),
        config=config,
    )

    assert "S_custom" in new_db_seqs.read_text()
