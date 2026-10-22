"""Build a small reference-only snapshot; never reads a patient data file."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
GTF_URL = 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz'
REFSEQ_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NM_001211.6&rettype=gb&retmode=text'
PROTEIN_URL = 'https://rest.uniprot.org/uniprotkb/O60566.json'
TRANSCRIPT = 'ENST00000287598.11'


def fetch(url):
    with urlopen(Request(url, headers={'User-Agent': 'BUB1B-public-reference-audit/0.1'}), timeout=45) as response:
        raw = response.read(40_000_001)
    if len(raw) > 40_000_000:
        raise ValueError('Unexpected public reference size')
    return raw


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gtf', type=Path, help='Optional local copy of the public MANE 1.5 GTF')
    args = parser.parse_args()
    gtf_raw = args.gtf.read_bytes() if args.gtf else fetch(GTF_URL)
    gb_raw = fetch(REFSEQ_URL)
    protein_raw = fetch(PROTEIN_URL)
    gb = gb_raw.decode()
    if not re.search(r'^VERSION\s+NM_001211\.6\s*$', gb, re.M):
        raise ValueError('Unexpected RefSeq version')
    cds_match = re.search(r'^     CDS\s+(\d+)\.\.(\d+)\s*$', gb, re.M)
    if not cds_match:
        raise ValueError('Expected single contiguous RefSeq CDS')
    start, end = map(int, cds_match.groups())
    transcript_sequence = re.sub('[^acgt]', '', gb.split('ORIGIN', 1)[1].split('//', 1)[0]).upper()
    cds_with_stop = transcript_sequence[start - 1:end]
    if cds_with_stop[-3:] not in ('TAA', 'TGA', 'TAG'):
        raise ValueError('Missing terminal stop')
    cds = cds_with_stop[:-3]
    protein = json.loads(protein_raw)['sequence']['value']
    # Independent, standard genetic code; CDS excludes its terminal stop.
    from itertools import product
    code = dict(zip(map(''.join, product('TCAG', repeat=3)),
                    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
    if ''.join(code[cds[i:i + 3]] for i in range(0, len(cds), 3)) != protein:
        raise ValueError('RefSeq CDS does not match independent UniProt protein')
    exons, coding = [], []
    for line in gzip.decompress(gtf_raw).decode().splitlines():
        if line.startswith('#') or f'transcript_id "{TRANSCRIPT}";' not in line:
            continue
        fields = line.split('\t')
        if fields[2] not in ('CDS', 'exon'):
            continue
        if fields[0] != 'chr15' or fields[6] != '+':
            raise ValueError('Unexpected reference contig/strand')
        interval = (int(fields[3]), int(fields[4]))
        (exons if fields[2] == 'exon' else coding).append(interval)
    exons.sort()
    coding.sort()
    if len(exons) != 23 or sum(b - a + 1 for a, b in coding) != len(cds):
        raise ValueError('Reference transcript geometry mismatch')
    rows, offset = [], 0
    for number, (a, b) in enumerate(exons, 1):
        parts = [(x, y) for x, y in coding if a <= x <= y <= b]
        if len(parts) != 1:
            raise ValueError('Expected one coding interval per exon in this transcript')
        x, y = parts[0]
        length = y - x + 1
        rows.append(dict(number=number, genomic_start1=a, genomic_end1=b,
                         coding_start1=offset + 1, coding_end1=offset + length))
        offset += length
    result = dict(gene='BUB1B', refseq='NM_001211.6', ensembl=TRANSCRIPT,
                  assembly='GRCh38', strand='+', cds_excludes_terminal_stop=True,
                  cds=cds, protein=protein, exons=rows,
                  provenance=[dict(url=url, sha256=digest(raw)) for url, raw in
                              [(GTF_URL, gtf_raw), (REFSEQ_URL, gb_raw), (PROTEIN_URL, protein_raw)]],
                  cds_sha256=digest(cds.encode()),
                  scope='Public reference sequences only; no sample sequence, calls or phenotype.')
    out = ROOT / 'data/public/bub1b_reference.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(f'Wrote {out.name}: {len(cds)} coding bases, {len(protein)} residues, {len(exons)} exons')


if __name__ == '__main__':
    main()
