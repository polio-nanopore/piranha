import csv
import itertools
import json
from pathlib import Path

import pytest

from piranha.analysis import haplo_functions, phylo_functions
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



def test_haplo_collapse_close_merges_identical_and_keeps_distant(tmp_path: Path):
    # 3 sites at 50/50 in the full read set. Haplotypes 0 and 1 both carry the
    # alt allele at every site (distance 0), haplotype 2 is ref at every site
    # (distance 3 from the others), so only 0 and 1 should merge.
    vcf = tmp_path / "vars.vcf"
    flopp = tmp_path / "flopp.txt"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "Ref\t5\t.\tC\tT\t50\tPASS\tRO=20;AO=20\n"
        "Ref\t10\t.\tG\tA\t50\tPASS\tRO=20;AO=20\n"
        "Ref\t15\t.\tA\tG\t50\tPASS\tRO=20;AO=20\n"
    )
    alt, ref = "0:0|1:20", "0:20|1:0"
    flopp.write_text(
        "header\n"
        + "".join(f"{i}:{pos}\t1\t1\t0\t{alt}\t{alt}\t{ref}\n" for i, pos in enumerate([5, 10, 15], start=1))
        + "*****\n"
    )

    assert haplo_functions.collapse_close(str(flopp), 2, str(vcf)) == [[0, 1], [2]]



def test_haplo_collapse_close_compares_past_already_merged_haplotypes(monkeypatch):
    # Regression: once the inner loop hit an already-merged haplotype it stopped
    # advancing its index, so later haplotypes were never compared.
    # Here 0~2 and 1~3 are close; everything else is far apart.
    haps = [{} for _ in range(4)]
    index = {id(h): i for i, h in enumerate(haps)}
    close_pairs = {(0, 2), (1, 3)}

    monkeypatch.setattr(haplo_functions, "parse_VCF", lambda _: {})
    monkeypatch.setattr(haplo_functions, "parse_flopp_reads", lambda _: haps)
    monkeypatch.setattr(
        haplo_functions,
        "get_haplo_dist",
        lambda a, b, _: 0 if (index[id(a)], index[id(b)]) in close_pairs else 10,
    )

    assert haplo_functions.collapse_close("flopp", 2, "vcf") == [[0, 2], [1, 3]]



def _write_haplotype_curation_inputs(tmp_path: Path):
    """
    Four flopp haplotypes over five 50/50 sites. Alt alleles carried:
        HAP0: {1,2}        HAP1: {1,2,3}      HAP2: {1,2,3,4}      HAP3: {5}
    SNP distances: 0-1=1, 0-2=2, 1-2=1, 0-3=3, 1-3=4, 2-3=5
    Each haplotype has 3 reads in the flopp partition file.
    """
    alt_sets = [{1, 2}, {1, 2, 3}, {1, 2, 3, 4}, {5}]
    sites = [1, 2, 3, 4, 5]
    alt, ref = "0:0|1:10", "0:10|1:0"

    vcf = tmp_path / "freebayes.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"Ref\t{pos}\t.\tC\tT\t50\tPASS\tRO=20;AO=20\n" for pos in sites)
    )

    flopp_lines = ["header\n"]
    for i, pos in enumerate(sites, start=1):
        hap_cols = ["1" if pos in s else "0" for s in alt_sets]
        support_cols = [alt if pos in s else ref for s in alt_sets]
        flopp_lines.append(f"{i}:{pos}\t" + "\t".join(hap_cols + support_cols) + "\n")
    flopp_lines.append("*****\n")
    flopp = tmp_path / "flopp_out.txt"
    flopp.write_text("".join(flopp_lines))

    partition = tmp_path / "part.txt"
    partition.write_text(
        "".join(f"#{h}\n" + "".join(f"h{h}_r{r}\tX\n" for r in range(3)) for h in range(4))
    )
    return flopp, vcf, partition



def test_haplotype_curation_merges_close_haplotypes_reads(tmp_path: Path):
    # Default cutoff of 2: HAP1 (1 SNP) and HAP2 (2 SNPs, on the boundary) merge
    # into HAP0. HAP3 (3+ SNPs from all) stays separate. 4 flopp haplotypes -> 2 reported.
    flopp, vcf, partition = _write_haplotype_curation_inputs(tmp_path)

    merge_info = haplo_functions.collapse_close(str(flopp), 2, str(vcf))
    assert merge_info == [[0, 1, 2], [3]]

    partitions = haplo_functions.parse_partition_file(str(partition))
    all_reads = set().union(*partitions.values())
    merged_reads = haplo_functions.merge_haplotype_reads(merge_info, partitions, all_reads)

    assert merged_reads == [
        partitions[0] | partitions[1] | partitions[2],
        partitions[3],
    ]
    # no read lost, and no read reported in more than one haplotype
    assert set().union(*merged_reads) == all_reads
    assert sum(len(r) for r in merged_reads) == len(all_reads)



def test_haplotype_curation_respects_distance_cutoff(tmp_path: Path):
    # Cutoff of 1: HAP2 is 2 SNPs from HAP0 but 1 from HAP1, and HAP1 is 1 from HAP0,
    # so the chain HAP0-HAP1-HAP2 merges into one group. HAP3 is 3+ from all.
    flopp, vcf, _ = _write_haplotype_curation_inputs(tmp_path)

    assert haplo_functions.collapse_close(str(flopp), 1, str(vcf)) == [[0, 1, 2], [3]]
    assert haplo_functions.collapse_close(str(flopp), 0, str(vcf)) == [[0], [1], [2], [3]]



@pytest.mark.parametrize("cutoff", [0, 1, 2, 3])
def test_haplotype_curation_is_independent_of_flopp_order(monkeypatch, cutoff):
    # Same four haplotypes as above (alt SNP sets {1,2}, {1,2,3}, {1,2,3,4}, {5}),
    # fed to collapse_close in every possible order. The groups of haplotypes that
    # get merged must be the same whichever order flopp lists them in.
    alt_sets = {"HAP0": {1, 2}, "HAP1": {1, 2, 3}, "HAP2": {1, 2, 3, 4}, "HAP3": {5}}
    sites = range(1, 6)

    def make_haplotype(alts):
        return {
            pos: {
                "assigned allele": "1" if pos in alts else "0",
                "total assigned": 10,
                "base support counts": {"0": 0, "1": 10} if pos in alts else {"0": 10, "1": 0},
            }
            for pos in sites
        }

    monkeypatch.setattr(
        haplo_functions,
        "parse_VCF",
        lambda _: {pos: {"frequencies": {"0": 0.5, "1": 0.5}} for pos in sites},
    )

    outcomes = set()
    for order in itertools.permutations(alt_sets):
        haps = [make_haplotype(alt_sets[name]) for name in order]
        monkeypatch.setattr(haplo_functions, "parse_flopp_reads", lambda _, haps=haps: haps)
        groups = haplo_functions.collapse_close("flopp", cutoff, "vcf")
        outcomes.add(frozenset(frozenset(order[i] for i in group) for group in groups))

    assert len(outcomes) == 1



def test_haplotype_curation_no_flopp_uses_all_reads():
    all_reads = ["r1", "r2", "r3"]
    assert haplo_functions.merge_haplotype_reads([[0]], {}, all_reads) == [set(all_reads)]



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
