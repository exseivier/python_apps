import click
import subprocess
import shlex

def __load_seqs(filename):
    """

    """
    FHIN = open(filename, "r")

    line = FHIN.readline()
    seqs = {}
    while line:
        line = line.strip("\n")
        if line[0] == ">":
            header = line
            seqs[header] = ""
        else:
            seqs[header] += line

        line = FHIN.readline()

    FHIN.close()
    return seqs

def __split(seqs, size, step):
    """

    """
    out_seqs = {}
    for head, seq in seqs.items():
        lseq = len(seq)
        i = 0
        start = i
        while start < lseq:
            end = start + size - 1
            out_seqs[head + "_" + str(i)] = seq[start:end]
            start = start + step
            i += 1

    return out_seqs


def __write_seqs(seqs, outfile):
    """

    """
    FHOUT = open(outfile, "w+")
    
    big_string = ""
    for head, seq in seqs.items():
        big_string += head + "\n" + seq + "\n"

    FHOUT.write(big_string)
    FHOUT.close()


def __load_headers(hd_file):
    """
        Load the headers of the sequences that hitted a subject in genomes database.
    """
    
    FHIN = open(hd_file, "r")

    line = FHIN.readline()
    heads = []
    while line:

        line = line.strip("\n")
        heads.append(line)
        line = FHIN.readline()

    FHIN.close()
    return heads

def __select_seqs(seqs, heads):
    """
        Selects those sequences which header is not in heads list.
    """

    i = 0 # iterator index for heads
    j = 0 # iterator index for seqs

    lsq = len(seqs)
    lhd = len(heads)
    list_seqs_keys = list(seqs.keys())
    out_seqs = {}

    while i < lhd:
        while j < lsq:
            if heads[i] == list_seqs_keys[j]:
                j += 1
                break
            else:
                out_seqs[list_seqs_keys[j]] = seqs[list_seqs_keys[j]]
                j += 1

        i += 1

    return out_seqs

def __catfile(data_base):
    """
        Returns a list with the lines stored in data_base file
    """
    temp = subprocess.Popen([f'cat', data_base], stdout=subprocess.PIPE)
    output = str(temp.communicate()[0])
    output = output.strip("'")
    output = output.strip("b'")
    output = output.split("\\n")
    print(output)
    return output


def __format_blast_db(output, pathway):
    """
        Preforms makeblastdb with the fasta files stored in data_base.
    """

    for item in output:
        if item != "":
            outfile = ".".join(item.split(".")[0:-1]) + ".db"
            item = item.strip("\n")

            print(outfile)
            FH_CODE = open("shuttle", "w")
            FH_CODE.write("makeblastdb \
                                -in {pathway}{item} \
                                -dbtype nucl \
                                -input_type fasta \
                                -title {item} \
                                -parse_seqids \
                                -out {pathway}{outfile}".format(pathway=pathway, \
                                                                item=item, \
                                                                outfile=outfile))
            FH_CODE.close()
            subprocess.Popen(["bash", "shuttle"]).communicate()
            subprocess.Popen(["rm", "shuttle"]).communicate()


@click.group()
def man_cli():
    pass

@man_cli.command()
@click.option('--in-file', '-i', help='Fasta input file name', required=True)
@click.option('--size', '-s', help="Genome fragment size", required=True, type=int)
@click.option('--step-size', '-S', help='Step size the slidding window takes between each cut', required=True, type=int)
@click.option('--out-file', '-o', help='file name to store split sequences', required=True)
def split_genome(in_file, size, step_size, out_file):
    """
        Split query genome in fragments of determined size and at each determined step
    """
 
    seqs = __load_seqs(in_file)
    
    seqs = __split(seqs, size, step_size)

    __write_seqs(seqs, out_file)

    print(f"Spliting genome {in_file}")

@man_cli.command()
@click.option('--in-file', '-i', help='Fasta file name of genome fragment sequences', required=True)
@click.option('--hd-file', '-hd', help='File name of fragment sequence headers that hitted on blast genomes database', required=True)
@click.option('--out-file', '-o', help='Temporary file to store the selected sequences', default='tmp.fasta', required=False)
def select_sequences(in_file, hd_file, out_file):
    """
        Selects those sequences that hit a subject in blast searching
    """
    
    seqs = __load_seqs(in_file)
    
    heads = __load_headers(hd_file)

    seqs = __select_seqs(seqs, heads)

    __write_seqs(seqs, out_file)

    print(f"Hiding matched sequences of {in_file}")


@man_cli.command()
@click.option('--data-base', '-d', help="Database file name that contains the name of the genome sequence files", required=True, type=str)
@click.option('--pathway', '-p', help='Pathway to genome sequences database', required=True, type=str)
def blast_db(data_base, pathway):
    """
        Formats fasta-formatted genome sequences to BLAST database format version 4.
    """
    print("Begins blast-db")

    output = __catfile(data_base)

    __format_blast_db(output, pathway)

    print("Ends blast-db")


@man_cli.command()

def do_blast():
    """
        Performs BLAST searching with query genome fragments.
    """

    pass

@man_cli.command()
def do_mapping():
    """
        Performs sequence mapping of fragments to reference genome.
    """
    pass

@man_cli.command()
def do_assembling():
    """
        Performs the genome-guided  assembly of the no-alike fragments.
    """
    pass

@man_cli.command()
def build_gtf():
    """
        Builds gtf file from not-alike fragments assembly and reference genome.
    """
    pass
