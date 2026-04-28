from Bio import SeqIO
import os


def normalize_lnc_id(x):
    """
    Normalizes ENST IDs.
    Example:
    ENST00000456328.2 -> ENST00000456328
    ENST00000456328|something -> ENST00000456328
    """
    x = x.strip()
    x = x.split("|")[0]
    x = x.split(".")[0]
    return x


def normalize_mir_id(x):
    """
    Normalizes miRNA IDs.
    Keeps naming consistent.
    """
    return x.strip()


def normalize_lnc_id(x):
    x = x.strip()
    x = x.split("|")[0]
    x = x.split()[0]
    return x


def normalize_mir_id(x):
    x = x.strip()
    x = x.split()[0]
    return x


def get_paired_lnc_mirna_index(path):
    lncrna_mirna_paired = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            parts = line.split("\t")

            if len(parts) < 2:
                continue

            lnc_name = normalize_lnc_id(parts[0])
            mirna_name = normalize_mir_id(parts[1])

            if lnc_name and mirna_name:
                lncrna_mirna_paired.append((lnc_name, mirna_name))

    return lncrna_mirna_paired

def get_mirna_lncrna_seq(paired_path, mirna_path, lncrna_path):
    os.makedirs("./data", exist_ok=True)

    lncrna_mirna_paired = get_paired_lnc_mirna_index(paired_path)

    print("Pairs extracted from validated file:", len(lncrna_mirna_paired))

    mirna = {}

    for x in SeqIO.parse(mirna_path, "fasta"):
        mir_id = normalize_mir_id(str(x.id))
        seq = str(x.seq).replace("U", "T")
        mirna[mir_id] = seq

    with open("./data/mirna.list", "w") as mirna_f:
        for x in sorted(mirna.keys()):
            mirna_f.write(x + "\n")

    lncrna = {}

    for x in SeqIO.parse(lncrna_path, "fasta"):
        lnc_id = normalize_lnc_id(str(x.id))
        seq = str(x.seq).replace("U", "T")

        if len(seq) > 200:
            lncrna[lnc_id] = seq

    with open("./data/lncrna.list", "w") as lncrna_f:
        for x in sorted(lncrna.keys()):
            lncrna_f.write(x + "\n")

    lnc_mir_pairs_id_seq = []
    matched = 0
    missing_lnc = 0
    missing_mir = 0
 
    with open("./data/lnc_mir_pairs.txt", "w") as lnc_mir_pairs_f:
        for lnc, mir in lncrna_mirna_paired:
            lnc_exists = lnc in lncrna
            mir_exists = mir in mirna

            if lnc_exists and mir_exists:
                matched += 1
                lnc_mir_pairs_id_seq.append(
                    [lnc, mir, lncrna[lnc], mirna[mir]]
                )
                lnc_mir_pairs_f.write(lnc + "," + mir + "\n")
            else:
                if not lnc_exists:
                    missing_lnc += 1
                if not mir_exists:
                    missing_mir += 1

    print("lncRNA sequences loaded:", len(lncrna))
    print("miRNA sequences loaded:", len(mirna))
    print("Matched lncRNA-miRNA pairs:", matched)
    print("Missing lncRNA IDs:", missing_lnc)
    print("Missing miRNA IDs:", missing_mir)
    print("Output written to ./data/lnc_mir_pairs.txt")

    return lnc_mir_pairs_id_seq, lncrna, mirna


paired_path = "./data/mirnas_lncrnas_validated.txt"
mirna_path = "./data/homo_mature_mirna.fa"
lncrna_path = "./data/outLncRNA.fa"

get_mirna_lncrna_seq(
    paired_path=paired_path,
    mirna_path=mirna_path,
    lncrna_path=lncrna_path
)