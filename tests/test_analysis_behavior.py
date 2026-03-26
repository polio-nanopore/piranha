import json
import csv
from pathlib import Path

import pytest

pytest.importorskip("pysam")

from piranha.analysis import consensus_functions, stool_functions, preprocessing
from piranha.utils.config import (
    KEY_ALN_BLOCK_LEN,
    KEY_BARCODE,
    KEY_DESCRIPTION,
    KEY_MAP_QUALITY,
    KEY_MAPPED,
    KEY_NUM_READS,
    KEY_PERCENT,
    KEY_READ_NAME,
    KEY_REFERENCE,
    KEY_REFERENCE_GROUP,
    KEY_REFERENCE_HIT,
    KEY_READ_HIT_START,
    KEY_READ_HIT_END,
    KEY_DIRECTION,
    KEY_STATUS,
    KEY_UNMAPPED,
    KEY_REFERENCE_GROUP,
    KEY_REFERENCE_GROUP_FIELD,
    KEY_SAMPLE,
)


def test_merge_indels_and_find_variants_basic():
    assert consensus_functions.merge_indels([5, 6, 7], "del") == ["5:del3"]

    reference = "ATGC-TA"
    query = "ATGTCTA"
    variants = consensus_functions.find_variants(reference, query)
    assert "4:CT" in variants
    assert "5:ins1" in variants


def test_find_ambiguity_percent():
    assert consensus_functions.find_ambiguity_pcent("ATGN-") == 20.0


def test_group_consecutive_sites():
    grouped = stool_functions.group_consecutive_sites([1, 2, 3, 6, 7])
    grouped_lists = [list(group) for group in grouped]
    assert grouped_lists == [[1, 2, 3], [6, 7]]


def test_get_mask_dict_and_mask_low_coverage(tmp_path: Path):
    mask_file = tmp_path / "mask.json"
    sequences = tmp_path / "input.fasta"
    output = tmp_path / "masked.fasta"

    mask_file.write_text(json.dumps({"seq1": [2, 3]}))
    sequences.write_text(">seq1|Sabin|hap|epi|2026-03-26\nATGC\n>seq2|Sabin|hap|epi|2026-03-26\nTTTT\n")

    mask_dict = stool_functions.get_mask_dict(str(mask_file))
    assert [list(r) for r in mask_dict["seq1"]] == [[2, 3]]

    stool_functions.mask_low_coverage(str(mask_file), str(sequences), str(output))
    out = output.read_text()
    assert "ANNN" in out
    assert "TTTT" in out


def test_gather_fasta_files_writes_json_key_matching_record_id(tmp_path: Path):
    summary = tmp_path / "summary.csv"
    barcodes_csv = tmp_path / "barcodes.csv"
    cns = tmp_path / "cns.fasta"
    out_fasta = tmp_path / "all.fasta"
    out_info = tmp_path / "all.json"
    publish_dir = tmp_path / "published"

    summary.write_text("barcode,reference,reference_group\nbarcode01,RefA,Sabin2-related\n")
    barcodes_csv.write_text("barcode,sample,EPID,date\nbarcode01,S1,E1,2026-03-26\n")
    cns.write_text(">RefA.h1|barcode01|1|3:AT\nATGC\n")

    config = {KEY_REFERENCE_GROUP_FIELD: KEY_REFERENCE_GROUP}

    stool_functions.gather_fasta_files(
        summary_info=str(summary),
        barcodes_csv=str(barcodes_csv),
        input_cns_list=[str(cns)],
        all_metdata=False,
        runname="run1",
        output_file=str(out_fasta),
        output_info=str(out_info),
        publish_dir=str(publish_dir),
        config=config,
    )

    payload = json.loads(out_info.read_text())
    key = next(iter(payload.keys()))
    assert " " not in key
    assert payload[key][KEY_BARCODE] == "barcode01"
    assert payload[key][KEY_SAMPLE] == "S1"
    assert payload[key][KEY_REFERENCE] == "RefA"


def test_gather_fasta_files_writes_configured_group_field(tmp_path: Path):
    summary = tmp_path / "summary.csv"
    barcodes_csv = tmp_path / "barcodes.csv"
    cns = tmp_path / "cns.fasta"
    out_fasta = tmp_path / "all.fasta"
    out_info = tmp_path / "all.json"
    publish_dir = tmp_path / "published"

    summary.write_text("barcode,reference,reference_group\nbarcode01,RefA,Sabin2-related\n")
    barcodes_csv.write_text("barcode,sample,EPID,date\nbarcode01,S1,E1,2026-03-26\n")
    cns.write_text(">RefA.h1|barcode01|1|3:AT\nATGC\n")

    config = {KEY_REFERENCE_GROUP_FIELD: "ddns_group"}

    stool_functions.gather_fasta_files(
        summary_info=str(summary),
        barcodes_csv=str(barcodes_csv),
        input_cns_list=[str(cns)],
        all_metdata=False,
        runname="run1",
        output_file=str(out_fasta),
        output_info=str(out_info),
        publish_dir=str(publish_dir),
        config=config,
    )

    payload = json.loads(out_info.read_text())
    first_record = payload[next(iter(payload.keys()))]
    assert first_record["ddns_group"] == "Sabin2-related"


def test_make_match_field_to_reference_group_map_classifies_values():
    mapping = preprocessing.make_match_field_to_reference_group_map(
        {
            "r1": "WPV1",
            "r2": "non-polio-enterovirus",
            "r3": "polio2",
            "r4": "sabin3",
            "r5": "wt-polio1",
            "r6": "p",
        }
    )

    assert mapping["r1"] == "WPV1"
    assert mapping["r2"] == "NonPolioEV"
    assert mapping["r3"] == "Sabin2-related"
    assert mapping["r4"] == "Sabin3-related"
    assert mapping["r5"] == "WPV1"
    assert mapping["r6"] == "PositiveControl"


def test_parse_line_and_add_to_hit_dict_records_mapped_hit():
    line = "read1\t100\t1\t90\t+\tRefA\t120\t2\t91\t80\t89\t60\n"
    mapping = preprocessing.parse_line(line)

    hits = {"RefA": set()}
    unmapped, status, description = preprocessing.add_to_hit_dict(
        hits=hits,
        mapping=mapping,
        min_map_len=50,
        min_map_quality=20,
        unmapped=0,
    )

    assert mapping[KEY_READ_HIT_START] == 1
    assert mapping[KEY_READ_HIT_END] == 90
    assert mapping[KEY_DIRECTION] == "+"
    assert mapping[KEY_ALN_BLOCK_LEN] == 89
    assert status == KEY_MAPPED
    assert unmapped == 0
    assert "MAPQ" in description
    assert len(hits["RefA"]) == 1


def test_group_hits_handles_ambiguous_reads(tmp_path: Path):
    paf = tmp_path / "reads.paf"
    mapping_filter = tmp_path / "mapping_filter.csv"

    paf.write_text(
        "read1\t100\t1\t90\t+\tRefA\t120\t2\t91\t80\t89\t60\n"
        "read1\t100\t1\t90\t+\tRefB\t120\t2\t91\t80\t89\t60\n"
        "read2\t100\t1\t90\t+\tRefA\t120\t2\t91\t80\t89\t60\n"
    )

    ref_group_hits, unmapped, ambiguous, total_reads, multi_hits, ref_group_ref = preprocessing.group_hits(
        paf_file=str(paf),
        ref_name_map={"RefA": "Sabin2-related", "RefB": "WPV1"},
        min_aln_block=50,
        min_map_quality=20,
        mapping_filter_file=str(mapping_filter),
    )

    assert total_reads == 2
    assert unmapped == 0
    assert ambiguous == 1
    assert multi_hits["RefA|RefB"] == 1
    assert "Sabin2-related" in ref_group_hits
    assert ref_group_ref["Sabin2-related"] == "RefA"

    with mapping_filter.open() as handle:
        rows = list(csv.DictReader(handle))
    assert rows[0][KEY_READ_NAME] == "read1"
    assert rows[0][KEY_STATUS] == KEY_MAPPED
    assert "ambiguous" in rows[0][KEY_DESCRIPTION]


def test_check_which_refs_to_write_filters_threshold_and_unmapped(tmp_path: Path):
    sample_hits = tmp_path / "sample_hits.csv"
    with sample_hits.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[KEY_BARCODE, KEY_REFERENCE, KEY_NUM_READS, KEY_PERCENT, KEY_REFERENCE_GROUP],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({KEY_BARCODE: "barcode01", KEY_REFERENCE: KEY_UNMAPPED, KEY_NUM_READS: 30, KEY_PERCENT: 50, KEY_REFERENCE_GROUP: KEY_UNMAPPED})
        writer.writerow({KEY_BARCODE: "barcode01", KEY_REFERENCE: "RefA", KEY_NUM_READS: 30, KEY_PERCENT: 50, KEY_REFERENCE_GROUP: "Sabin2-related"})
        writer.writerow({KEY_BARCODE: "barcode01", KEY_REFERENCE: "RefB", KEY_NUM_READS: 1, KEY_PERCENT: 1, KEY_REFERENCE_GROUP: "WPV1"})

    refs = preprocessing.check_which_refs_to_write(str(sample_hits), min_reads=5, min_pcent=10)
    assert refs == ["RefA"]
