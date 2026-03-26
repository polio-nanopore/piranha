import csv
import json
from pathlib import Path

from piranha.report import make_report
from piranha.utils.config import (
    KEY_BARCODE_REPORT_TEMPLATE,
    KEY_DATE,
    KEY_EPID,
    KEY_INSTITUTE,
    KEY_LANGUAGE,
    KEY_MASKED_SITES,
    KEY_MIN_PCENT,
    KEY_MIN_READS,
    KEY_NEGATIVE,
    KEY_NOTES,
    KEY_ORIENTATION,
    KEY_OUTDIR,
    KEY_OUTPUT_PREFIX,
    KEY_POSITIVE,
    KEY_READ_COUNT,
    KEY_REFERENCE_GROUP_FIELD,
    KEY_REPORT_TEMPLATE,
    KEY_RUNNAME,
    KEY_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS,
    KEY_SUMMARY_DATA,
    KEY_SNP_SITES,
    KEY_TEMPDIR,
    KEY_VARIANTS,
    KEY_VARIANT_COUNT,
    KEY_REFERENCE,
    KEY_BARCODES,
    KEY_DETAILED_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS,
    KEY_RUN_PHYLO,
    KEY_BARCODE,
    KEY_COMPOSITION_TABLE,
    KEY_PERCENT,
    KEY_REFERENCE,
    KEY_REFERENCE_GROUP,
    KEY_SAMPLE,
    KEY_SUMMARY_TABLE,
)


def test_get_snipit_returns_file_contents(tmp_path: Path):
    snipit = tmp_path / "RefA.svg"
    snipit.write_text("<svg>\nline\n</svg>\n")

    payload = make_report.get_snipit("RefA", str(snipit))
    assert payload.startswith("<svg>\n")
    assert "line\n" in payload


def test_assign_bcode_to_well_vertical_and_horizontal():
    vertical = make_report.assign_bcode_to_well("vertical")
    horizontal = make_report.assign_bcode_to_well("horizontal")

    assert vertical["A01"] == "barcode01"
    assert vertical["H12"] == "barcode96"
    assert horizontal["A01"] == "barcode01"
    assert horizontal["A02"] == "barcode02"


def test_barcode_to_well_uses_well_column_when_present(tmp_path: Path):
    barcode_csv = tmp_path / "barcodes.csv"
    barcode_csv.write_text("sample,barcode,well\nS1,barcode01,B03\n")

    mapping = make_report.barcode_to_well(str(barcode_csv), "vertical")
    assert mapping == {"B03": "barcode01"}


def test_data_for_plate_viz_marks_presence(tmp_path: Path):
    barcode_csv = tmp_path / "barcodes.csv"
    barcode_csv.write_text("sample,barcode,well\nS1,barcode01,A01\n")

    positives = {"barcode01": {"Sabin2-related": 100}}
    plate_json, positive_types = make_report.data_for_plate_viz(
        positives_for_plate_viz=positives,
        barcode_csv=str(barcode_csv),
        orientation="vertical",
        barcodes=["barcode01"],
    )

    parsed = json.loads(plate_json)
    assert len(parsed) == 96
    a01 = next(item for item in parsed if item["x"] == 1 and item["y"] == "A")
    assert a01["All"] == "Present"
    assert "Sabin2-related" in positive_types


def test_get_background_data_reads_csv(tmp_path: Path):
    metadata = tmp_path / "background.csv"
    with metadata.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "country"])
        writer.writeheader()
        writer.writerow({"name": "seqA", "country": "UK"})

    data = make_report.get_background_data(str(metadata), config={})
    parsed = json.loads(data)
    assert parsed["seqA"]["country"] == "UK"


def test_make_detailed_csv_writes_combined_fields(tmp_path: Path):
    barcodes = tmp_path / "barcodes.csv"
    detailed = tmp_path / "detailed.csv"

    barcodes.write_text("sample,barcode,EPID,institute\nS1,barcode01,EPI1,Inst1\n")

    data_for_report = {
        KEY_SUMMARY_TABLE: [
            {
                KEY_BARCODE: "barcode01",
                KEY_SAMPLE: "S1",
                KEY_REFERENCE_GROUP: "Sabin2-related",
                KEY_REFERENCE: "RefA",
                "Number of mutations": 2,
                KEY_PERCENT: 99.2,
                "Sample classification": "Sabin-like",
            }
        ],
        KEY_COMPOSITION_TABLE: [
            {
                KEY_SAMPLE: "S1",
                KEY_BARCODE: "barcode01",
                "Sabin2-related": "120",
                "NonPolioEV": "3",
                "unmapped": "0",
            }
        ],
    }

    make_report.make_detailed_csv(
        data_for_report=data_for_report,
        barcodes_csv=str(barcodes),
        epi_csv="",
        output=str(detailed),
        detailed_header_fields=[
            "Sabin2-related|closest_reference",
            "Sabin2-related|num_reads",
            "Sabin2-related|nt_diff_from_reference",
            "Sabin2-related|pcent_match",
            "Sabin2-related|classification",
            "NonPolioEV|num_reads",
            "comments",
        ],
    )

    with detailed.open() as handle:
        row = next(csv.DictReader(handle))

    assert row["Sabin2-related|closest_reference"] == "RefA"
    assert row["Sabin2-related|num_reads"] == "120"
    assert row["Sabin2-related|classification"] == "Sabin-like"


def test_make_sample_report_end_to_end(tmp_path: Path):
    tempdir = tmp_path / "temp"
    snipit_dir = tempdir / "barcode01" / "snipit"
    snipit_dir.mkdir(parents=True)
    (snipit_dir / "RefA.h1.svg").write_text("<svg>ok</svg>\n")

    variation = tmp_path / "variation.json"
    consensus = tmp_path / "consensus.fasta"
    consensus_info = tmp_path / "consensus_info.json"
    masked = tmp_path / "masked.csv"
    template = tmp_path / "barcode_report.mako"
    output = tmp_path / "barcode01_report.html"

    consensus.write_text(
        ">S1|Sabin2-related|h1|E1|2026-03-26 barcode=barcode01\nATGC\n"
    )
    consensus_info.write_text(
        json.dumps(
            {
                "S1|Sabin2-related|h1|E1|2026-03-26": {
                    KEY_BARCODE: "barcode01",
                    KEY_REFERENCE: "RefA",
                    KEY_VARIANT_COUNT: "1",
                    KEY_VARIANTS: "3:CT",
                }
            }
        )
    )
    variation.write_text(
        json.dumps(
            {
                "RefA.h1": {
                    "variation": [{"Position": 3, "Depth": 50}, {"Position": 8, "Depth": 20}],
                    "coocc": {},
                }
            }
        )
    )
    masked.write_text("reference,site\nRefA.h1,8\n")
    template.write_text("Barcode: ${barcode} Sample: ${sample} Seqs: ${sequences}")

    config = {
        KEY_TEMPDIR: str(tempdir),
        KEY_LANGUAGE: "English",
        KEY_BARCODE_REPORT_TEMPLATE: str(template),
    }

    cns_config = {
        "barcode01": ["RefA.h1"],
        "RefA.h1": 120,
    }

    make_report.make_sample_report(
        report_to_generate=str(output),
        variation_file=str(variation),
        consensus_seqs=str(consensus),
        consensus_info=str(consensus_info),
        masked_variants=str(masked),
        barcode="barcode01",
        cns_config=cns_config,
        config=config,
    )

    rendered = output.read_text()
    assert "Barcode: barcode01" in rendered
    assert "Sample: S1" in rendered


def test_make_output_report_end_to_end_minimal(tmp_path: Path):
    outdir = tmp_path / "out"
    outdir.mkdir()

    barcodes_csv = tmp_path / "barcodes.csv"
    epi_csv = tmp_path / "epi.csv"
    preprocessing_summary = tmp_path / "preprocessing_summary.csv"
    sample_composition = tmp_path / "sample_composition.csv"
    consensus = tmp_path / "consensus.fasta"
    consensus_info = tmp_path / "consensus_info.json"
    detailed_csv_out = tmp_path / "detailed.csv"
    annotations = tmp_path / "annotations.csv"
    template = tmp_path / "report_template.mako"
    report_out = tmp_path / "report.html"

    barcodes_csv.write_text("sample,barcode,EPID,institute\nS1,barcode01,E1,Inst\n")
    epi_csv.write_text("sample,extra\nS1,val\n")
    preprocessing_summary.write_text(
        "sample,barcode,Sabin2-related,NonPolioEV,unmapped\nS1,barcode01,120,10,0\n"
    )
    sample_composition.write_text("sample,barcode\nS1,barcode01\n")
    consensus.write_text(
        ">S1|Sabin2-related|h1|E1|2026-03-26 barcode=barcode01 variant_count=1 variants=3:CT reference=RefA\nATGC\n"
    )
    consensus_info.write_text(
        json.dumps(
            {
                "S1|Sabin2-related|h1|E1|2026-03-26": {
                    KEY_BARCODE: "barcode01",
                    KEY_VARIANT_COUNT: "1",
                    KEY_VARIANTS: "3:CT",
                    KEY_REFERENCE: "RefA",
                }
            }
        )
    )
    annotations.write_text("name,country\nseqA,UK\n")
    template.write_text("Run: ${run_name} Samples: ${len(data_for_report['summary_table'])}")

    config = {
        KEY_NEGATIVE: [],
        KEY_POSITIVE: [],
        KEY_MIN_READS: 5,
        KEY_MIN_PCENT: 10,
        KEY_ORIENTATION: "vertical",
        KEY_BARCODES: ["barcode01"],
        KEY_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS: ["sample", "barcode", "Sabin2-related", "NonPolioEV", "unmapped"],
        KEY_DETAILED_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS: [
            "Sabin2-related|closest_reference",
            "Sabin2-related|num_reads",
            "Sabin2-related|nt_diff_from_reference",
            "Sabin2-related|pcent_match",
            "Sabin2-related|classification",
            "NonPolioEV|num_reads",
            "comments",
        ],
        KEY_RUN_PHYLO: False,
        KEY_LANGUAGE: "English",
        KEY_REPORT_TEMPLATE: str(template),
        KEY_RUNNAME: "testrun",
        KEY_OUTDIR: str(outdir),
        KEY_OUTPUT_PREFIX: "analysis",
        KEY_NOTES: "",
        KEY_REFERENCE_GROUP_FIELD: "reference_group",
    }

    make_report.make_output_report(
        report_to_generate=str(report_out),
        barcodes_csv=str(barcodes_csv),
        epi_csv=str(epi_csv),
        preprocessing_summary=str(preprocessing_summary),
        sample_composition=str(sample_composition),
        consensus_seqs=str(consensus),
        consensus_info=str(consensus_info),
        detailed_csv_out=str(detailed_csv_out),
        annotations_file=str(annotations),
        config=config,
    )

    assert report_out.exists()
    assert detailed_csv_out.exists()
    assert "Run: testrun" in report_out.read_text()
