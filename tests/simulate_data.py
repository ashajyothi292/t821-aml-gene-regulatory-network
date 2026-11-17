"""
simulate_data.py - small fake inputs for tests/run_tests.sh.

  python tests/simulate_data.py <dir>                        MACS2 summits, blacklist, RNA counts, GMT
  python tests/simulate_data.py <dir> --homer-table          HOMER-style table over the consensus peaks
  python tests/simulate_data.py <dir> --footprints-and-chip  footprints and ChIP peaks
"""

import random
import sys
from pathlib import Path

random.seed(1)
CONTROLS = ["PBSC-1", "PBSC-2"]
PATIENTS = ["t821-1", "t821-1rel", "t821-3", "t821-4", "t821-KIT"]
SAMPLES = CONTROLS + PATIENTS
N_SITES = 1500


def site_positions():
    return [("chr1" if i % 2 else "chr2", 10_000 + i * 2_000) for i in range(N_SITES)]


def base(d: Path):
    (d / "reference").mkdir(parents=True, exist_ok=True)
    (d / "dnase" / "peaks").mkdir(parents=True, exist_ok=True)
    (d / "reference" / "hg38-blacklist.v2.bed").write_text("chr1\t4900\t5100\n")

    for s in SAMPLES:
        rows = ["chr1\t5000\t5001\tbl\t40", "chrUn_KI270302v1\t100\t101\tun\t40"]
        for chrom, pos in site_positions():
            if random.random() < 0.85:
                jitter = random.randint(-30, 30)
                rows.append(f"{chrom}\t{pos + jitter}\t{pos + jitter + 1}\tp\t{random.uniform(5, 60):.2f}")
        (d / "dnase" / "peaks" / f"{s}_summits.bed").write_text("\n".join(rows) + "\n")

    # RNA-seq: featureCounts table with real Ensembl IDs so symbol mapping works
    genes = [("ENSG00000141510", "17"), ("ENSG00000136997", "8"), ("ENSG00000105568", "19"),
             ("ENSG00000159216", "21"), ("ENSG00000102145", "X"), ("ENSG00000198804", "MT")]
    extra = [f"ENSG{n:011d}" for n in range(100_000, 102_000)]
    try:
        import subprocess
        out = subprocess.run(
            ["Rscript", "-e",
             "suppressMessages(library(org.Hs.eg.db)); k <- keys(org.Hs.eg.db, 'ENSEMBL');"
             "set.seed(1); cat(sample(k, 2000), sep='\\n')"],
            capture_output=True, text=True, check=True).stdout.split()
        extra = out
    except Exception:
        pass
    all_genes = genes + [(g, str(random.randint(1, 22))) for g in extra]
    rna_samples = ["PBSC-1", "PBSC-2", "t821-1", "t821-1R", "t821-3", "t821-4", "t821-ITDx2-5", "t821-KIT-2"]
    (d / "rnaseq" / "counts").mkdir(parents=True, exist_ok=True)
    lines = ["# Program:featureCounts v2.0.4",
             "Geneid\tChr\tStart\tEnd\tStrand\tLength\t" + "\t".join(f"{s}.bam" for s in rna_samples)]
    for i, (g, chrom) in enumerate(all_genes):
        mu = random.uniform(20, 2000)
        shift = 4.0 if i % 10 == 0 else (0.25 if i % 10 == 1 else 1.0)
        counts = [int(random.gauss(mu, mu * 0.15)) if j < 2 else int(random.gauss(mu * shift, mu * shift * 0.15))
                  for j in range(8)]
        counts = [max(c, 0) for c in counts]
        lines.append(f"{g}\t{chrom};{chrom}\t1;50\t40;90\t+;+\t80\t" + "\t".join(map(str, counts)))
    (d / "rnaseq" / "counts" / "featureCounts_gene_counts.txt").write_text("\n".join(lines) + "\n")

    # GMT with sets built from the mapped symbols happens after DESeq2; use a
    # broad placeholder here and fill it once rank files exist (see below)
    (d / "reference" / "test_reactome.gmt").write_text("")


def homer_table(d: Path):
    peaks = [l.split("\t")[:3] for l in (d / "dnase" / "consensus" / "final_consensus_peaks.bed")
             .read_text().splitlines() if l]
    head = ["PeakID (cmd=annotatePeaks.pl)", "Chr", "Start", "End", "Strand", "Peak Score",
            "Focus Ratio/Region Size", "Annotation", "Detailed Annotation", "Distance to TSS",
            "Nearest PromoterID", "Entrez ID", "Nearest Unigene", "Nearest Refseq", "Nearest Ensembl",
            "Gene Name", "Gene Alias", "Gene Description", "Gene Type"] + [f"{s} bedGraph" for s in SAMPLES]
    feats = ["intron (NM_1, intron 2 of 5)", "Intergenic", "promoter-TSS (NM_2)", "exon (NM_3)", "TTS (NM_4)"]
    rows = ["\t".join(head)]
    for i, (c, s, e) in enumerate(peaks):
        base_sig = random.uniform(20, 200)
        change = 6.0 if i % 8 == 0 else (1 / 6.0 if i % 8 == 1 else 1.0)
        sig = [random.gauss(base_sig, base_sig * 0.1) if smp in CONTROLS
               else random.gauss(base_sig * change, base_sig * change * 0.1) for smp in SAMPLES]
        ann = [f"peak_{i + 1}", c, s, e, "+", "0", "0", random.choice(feats), "", str(random.randint(-50000, 50000)),
               "", "", "", "", "", f"GENE{i % 300}", "", "", "protein-coding"]
        rows.append("\t".join(ann + [f"{max(v, 0):.1f}" for v in sig]))
    (d / "dnase" / "consensus" / "peak_signal_annotated.txt").write_text("\n".join(rows) + "\n")


def footprints_and_chip(d: Path):
    da = d / "dnase" / "differential_accessibility"
    gained = [l.split("\t") for l in (da / "gained_peaks.bed").read_text().splitlines()]
    lost = [l.split("\t") for l in (da / "lost_peaks.bed").read_text().splitlines()]

    def mid(r):
        return int(r[1]) + (int(r[2]) - int(r[1])) // 2

    def bed(rows, frac, width):
        out = []
        for r in rows:
            if random.random() < frac:
                m = mid(r)
                out.append(f"{r[0]}\t{m - width}\t{m + width}")
        return "\n".join(sorted(out)) + "\n"

    fp = d / "dnase" / "footprinting"
    for sample, rows in (("t821-1", gained), ("PBSC-1", lost)):
        (fp / sample).mkdir(parents=True, exist_ok=True)
        (fp / sample / f"{sample}.output.bam.WellingtonFootprints.FDR.0.01.bed").write_text(bed(rows, 0.4, 10))
    chip = d / "chipseq" / "peaks"
    chip.mkdir(parents=True, exist_ok=True)
    (chip / "RUNX1-ETO_peaks.bed").write_text(bed(gained, 0.6, 300))
    (chip / "RUNX1_peaks.bed").write_text(bed(lost, 0.6, 300))
    (chip / "H3K27ac_peaks.bed").write_text(bed(gained + lost, 0.7, 500))


def gmt_from_ranks(d: Path):
    ranks = d / "rnaseq" / "deseq2" / "ranks_all.rnk"
    genes = [l.split("\t")[0] for l in ranks.read_text().splitlines()]
    sets = [f"REACTOME_TOP_SET_{k}\tna\t" + "\t".join(genes[k * 20:(k + 1) * 20 + 10]) for k in range(3)]
    sets.append("REACTOME_BOTTOM_SET\tna\t" + "\t".join(genes[-40:]))
    sets += [f"REACTOME_RANDOM_{k}\tna\t" + "\t".join(random.sample(genes, 30)) for k in range(6)]
    (d / "reference" / "test_reactome.gmt").write_text("\n".join(sets) + "\n")


if __name__ == "__main__":
    target = Path(sys.argv[1])
    if "--homer-table" in sys.argv:
        homer_table(target)
    elif "--footprints-and-chip" in sys.argv:
        footprints_and_chip(target)
    elif "--gmt" in sys.argv:
        gmt_from_ranks(target)
    else:
        base(target)
