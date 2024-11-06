"""
Contains most of the backend for the DRAM_setup.py script, used to setup databases for each user.
"""
from os import path, mkdir
from datetime import datetime
from shutil import move, rmtree, copyfile
from glob import glob
import gzip
from collections import defaultdict
import logging
import tarfile

from skbio import read as read_sequence
from skbio import write as write_sequence

from mag_annotator.database_handler import DatabaseHandler
from mag_annotator.utils import run_process, make_mmseqs_db, download_file, merge_files, remove_prefix, remove_suffix, setup_logger
from mag_annotator.camper_kit import DRAM_SETTINGS as CAMPER_DRAM_SETTINGS
from mag_annotator.camper_kit import PROCESS_OPTIONS as CAMPER_PROCESS_OPTIONS
from mag_annotator.camper_kit import DOWNLOAD_OPTIONS as CAMPER_DOWNLOAD_OPTIONS
from mag_annotator.camper_kit import process as process_camper_tar_gz
from mag_annotator.camper_kit import download as download_camper_tar_gz

NUMBER_OF_VIRAL_FILES = 1
DEFAULT_DBCAN_RELEASE = '11'
DEFAULT_DBCAN_DATE = '08062022'
DEFAULT_UNIREF_VERSION = '90'
DEFAULT_VOGDB_VERSION = 'latest'
DFLT_OUTPUT_DIR = '.'
LOGGER = logging.getLogger("database_processing.log")
DEFAULT_MMMSPRO_DB_NAME = 'db'


KOFAM_CITATION = ("T. Aramaki, R. Blanc-Mathieu, H. Endo, K. Ohkubo, M. Kanehisa"
                  ", S. Goto, and H. Ogata, \"Kofamkoala: Kegg ortholog assignme"
                  "nt based on profile hmm and adaptive score threshold,\" Bioin"
                  "formatics, vol. 36, no. 7, pp. 2251–2252, 2020."
                  )
VIRAL_REFSEQ_CITATION = ("J. R. Brister, D. Ako-Adjei, Y. Bao, and O. Blinkova, "
                         "\"Ncbi viral genomes resource,\" Nucleic acids researc"
                         "h, vol. 43, no. D1, pp. D571–D577, 2015. [3] M. Kanehi"
                         "sa, M. Furumichi, Y. Sato, M. Ishiguro-Watanabe, and M"
                         ". Tan-abe, \"Kegg: integrating viruses and cellular or"
                         "ganisms,\" Nucleic acids research, vol. 49, no. D1, pp"
                         ". D545–D551, 2021."
                 )
KEGG_CITATION = (" M. Kanehisa, M. Furumichi, Y. Sato, M. Ishiguro-Watanabe, and"
                 " M. Tanabe, \"Kegg: integrating viruses and cellular organisms"
                 ",\" Nucleic acids research, vol. 49, no. D1, pp. D545–D551, 20"
                 "21."
                 )
PFAM_CITATION = ("J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G. A. Sal"
                 "azar, E. L. Sonnhammer, S. C. Tosatto, L. Paladin, S. Raj, L."
                 " J. Richardson et al., \"Pfam: The protein families database "
                 "in 2021,\" Nucleic acids research, vol. 49, no. D1, pp. D412–"
                 "D419, 2021."
                 )
PEPTIDASE_CITATION = ("N. D. Rawlings, A. J. Barrett, P. D. Thomas, X. Huang, A"
                      ". Bateman, and R. D. Finn, \"The merops database of prot"
                      "eolytic enzymes, their substrates and inhibitors in 2017"
                      " and a comparison with peptidases in the panther databas"
                      "e,\" Nucleic acids research, vol. 46, no. D1, pp. D624–D"
                      "632, 2018."
                      )
VOGDB_CITATION = ("J. Thannesberger, H.-J. Hellinger, I. Klymiuk, M.-T. Kastner"
                  ", F. J. Rieder, M. Schneider, S. Fister, T. Lion, K. Kosulin"
                  ", J. Laengle et al., \"Viruses comprise an extensive pool of"
                  " mobile genetic elements in eukaryote cell cultures and huma"
                  "n clinical samples,\" The FASEB Journal, vol. 31, no. 5, pp."
                  " 1987–2000, 2017."
                  )
UNIREF_CITATION = ("Y. Wang, Q. Wang, H. Huang, W. Huang, Y. Chen, P. B. McGarv"
                  "ey, C. H. Wu, C. N. Arighi, and U. Consortium, \"A crowdsour"
                  "cing open platform for literature curation in uniprot,\" PLo"
                  "S Biology, vol. 19, no. 12, p. e3001464, 2021."
                   )
DBCAN_CITATION = ("Y. Yin, X. Mao, J. Yang, X. Chen, F. Mao, and Y. Xu, \"dbcan"
                  ": a web resource for automated carbohydrate-active enzyme an"
                  "notation,\" Nucleic acids research, vol. 40, no. W1, pp. W44"
                  "5–W451, 2012."
                  )
DRAM_CITATION = ("M. Shaffer, M. A. Borton, B. B. McGivern, A. A. Zayed, S. L. "
                 "La Rosa, L. M. Solden, P. Liu, A. B. Narrowe, J. Rodríguez-Ra"
                 "mos, B. Bolduc et al., \"Dram for distilling microbial metabo"
                 "lism to automate the curation of microbiome function,\" Nucle"
                 "ic acids research, vol. 48, no. 16, pp. 8883–8900, 2020."
                 )
# check 
# if dbcapfam is down, raise appropriate error
# TODO: upgrade to pigz?


KOFAM_CITATION = ("T. Aramaki, R. Blanc-Mathieu, H. Endo, K. Ohkubo, M. Kanehisa"
                  ", S. Goto, and H. Ogata, \"Kofamkoala: Kegg ortholog assignme"
                  "nt based on profile hmm and adaptive score threshold,\" Bioin"
                  "formatics, vol. 36, no. 7, pp. 2251–2252, 2020."
                  )
VIRAL_REFSEQ_CITATION = ("J. R. Brister, D. Ako-Adjei, Y. Bao, and O. Blinkova, "
                         "\"Ncbi viral genomes resource,\" Nucleic acids researc"
                         "h, vol. 43, no. D1, pp. D571–D577, 2015. [3] M. Kanehi"
                         "sa, M. Furumichi, Y. Sato, M. Ishiguro-Watanabe, and M"
                         ". Tan-abe, \"Kegg: integrating viruses and cellular or"
                         "ganisms,\" Nucleic acids research, vol. 49, no. D1, pp"
                         ". D545–D551, 2021."
                 )
KEGG_CITATION = (" M. Kanehisa, M. Furumichi, Y. Sato, M. Ishiguro-Watanabe, and"
                 " M. Tanabe, \"Kegg: integrating viruses and cellular organisms"
                 ",\" Nucleic acids research, vol. 49, no. D1, pp. D545–D551, 20"
                 "21."
                 )
PFAM_CITATION = ("J. Mistry, S. Chuguransky, L. Williams, M. Qureshi, G. A. Sal"
                 "azar, E. L. Sonnhammer, S. C. Tosatto, L. Paladin, S. Raj, L."
                 " J. Richardson et al., \"Pfam: The protein families database "
                 "in 2021,\" Nucleic acids research, vol. 49, no. D1, pp. D412–"
                 "D419, 2021."
                 )
PEPTIDASE_CITATION = ("N. D. Rawlings, A. J. Barrett, P. D. Thomas, X. Huang, A"
                      ". Bateman, and R. D. Finn, \"The merops database of prot"
                      "eolytic enzymes, their substrates and inhibitors in 2017"
                      " and a comparison with peptidases in the panther databas"
                      "e,\" Nucleic acids research, vol. 46, no. D1, pp. D624–D"
                      "632, 2018."
                      )
VOGDB_CITATION = ("J. Thannesberger, H.-J. Hellinger, I. Klymiuk, M.-T. Kastner"
                  ", F. J. Rieder, M. Schneider, S. Fister, T. Lion, K. Kosulin"
                  ", J. Laengle et al., \"Viruses comprise an extensive pool of"
                  " mobile genetic elements in eukaryote cell cultures and huma"
                  "n clinical samples,\" The FASEB Journal, vol. 31, no. 5, pp."
                  " 1987–2000, 2017."
                  )
UNIREF_CITATION = ("Y. Wang, Q. Wang, H. Huang, W. Huang, Y. Chen, P. B. McGarv"
                  "ey, C. H. Wu, C. N. Arighi, and U. Consortium, \"A crowdsour"
                  "cing open platform for literature curation in uniprot,\" PLo"
                  "S Biology, vol. 19, no. 12, p. e3001464, 2021."
                   )
DBCAN_CITATION = ("Y. Yin, X. Mao, J. Yang, X. Chen, F. Mao, and Y. Xu, \"dbcan"
                  ": a web resource for automated carbohydrate-active enzyme an"
                  "notation,\" Nucleic acids research, vol. 40, no. W1, pp. W44"
                  "5–W451, 2012."
                  )
DRAM_CITATION = ("M. Shaffer, M. A. Borton, B. B. McGivern, A. A. Zayed, S. L. "
                 "La Rosa, L. M. Solden, P. Liu, A. B. Narrowe, J. Rodríguez-Ra"
                 "mos, B. Bolduc et al., \"Dram for distilling microbial metabo"
                 "lism to automate the curation of microbiome function,\" Nucle"
                 "ic acids research, vol. 48, no. 16, pp. 8883–8900, 2020."
                 )

def get_iso_date():
    return datetime.today().strftime('%Y%m%d')


def download_pfam_hmm(output_dir='.', logger=LOGGER, verbose=True):
    pfam_hmm = path.join(output_dir, 'Pfam-A.hmm.dat.gz')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
    url_http = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz'
    logger.debug(f"Downloading Pfam from: {url}")
    download_file(url, pfam_hmm, logger, alt_urls=[url_http], verbose=verbose)
    return pfam_hmm


def download_dbcan(output_dir='.', logger=LOGGER, dbcan_hmm=None, version=DEFAULT_DBCAN_RELEASE, verbose=True):
    dbcan_hmm = path.join(output_dir, f"dbCAN-HMMdb-V{version}.txt" )
    if int(version) < int(DEFAULT_DBCAN_RELEASE):
        link_path = f"http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V{version}.txt"
    else:
        link_path = f"http://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V{version}.txt"

    logger.debug(f"Downloading dbCAN from: {link_path}")
    download_file(link_path, dbcan_hmm, logger, verbose=verbose)
    return dbcan_hmm


def download_dbcan_fam_activities (output_dir='.', logger=LOGGER, version=DEFAULT_DBCAN_RELEASE, upload_date=DEFAULT_DBCAN_DATE,
                                verbose=True):
    dbcan_fam_activities = path.join(output_dir, f'CAZyDB.{upload_date}.fam-activities.txt')
    url = f"https://bcb.unl.edu/dbCAN2/download/Databases/V{version}/CAZyDB.{upload_date}.fam-activities.txt"
    logger.info(f"Downloading dbCAN family activities from : {url}")
    download_file(url, dbcan_fam_activities, logger, verbose=verbose)
    return dbcan_fam_activities


def download_dbcan_subfam_ec(output_dir='.', logger=LOGGER, version=DEFAULT_DBCAN_RELEASE, upload_date=DEFAULT_DBCAN_DATE, verbose=True):
    dbcan_subfam_ec = path.join(output_dir, f"CAZyDB.{upload_date}.fam.subfam.ec.txt")
    url = (f"https://bcb.unl.edu/dbCAN2/download/Databases/"
                 f"V{version}/CAZyDB.{upload_date}.fam.subfam.ec.txt")
    logger.info(f"Downloading dbCAN sub-family encumber from : {url}")
    download_file(url, dbcan_subfam_ec, logger, verbose=verbose)
    return dbcan_subfam_ec


def download_kofam_hmm(output_dir='.', logger=LOGGER, verbose=False):
    kofam_profile_tar_gz = path.join(output_dir, 'kofam_profiles.tar.gz')
    url = 'ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz'  
    url_http = 'https://www.genome.jp/ftp/db/kofam/profiles.tar.gz'
    download_file(url, kofam_profile_tar_gz, logger, alt_urls=[url_http], verbose=verbose)
    return kofam_profile_tar_gz

def generate_modified_kegg_fasta(kegg_fasta, gene_ko_link_loc=None):
    """
    Takes kegg fasta file and gene ko link file, adds kos not already in headers to headers
    Whish I knew about this, oh well I may split this out.
    """
    genes_ko_dict = defaultdict(list)
    if gene_ko_link_loc is not None:
        if gene_ko_link_loc.endswith('.gz'):
            gene_ko_link_fh = gzip.open(gene_ko_link_loc, 'rt')
        else:
            gene_ko_link_fh = open(gene_ko_link_loc)
        for line in gene_ko_link_fh:
            gene, ko = line.strip().split()
            genes_ko_dict[gene].append(remove_prefix(ko, 'ko:'))
    for seq in read_sequence(kegg_fasta, format='fasta'):
        new_description = seq.metadata['description']
        for ko in genes_ko_dict[seq.metadata['id']]:
            if ko not in new_description:
                new_description += '; %s' % ko
        seq.metadata['description'] = new_description
        yield seq



def process_kegg(kegg_loc, output_dir, logger, gene_ko_link_loc=None, download_date=None, 
                    threads=10, verbose=True):
    if download_date is None:
        download_date = get_iso_date()
    if gene_ko_link_loc is not None:
        # add KOs to end of header where KO is not already there
        kegg_mod_loc = path.join(output_dir, 'kegg.mod.fa')
        write_sequence(generate_modified_kegg_fasta(kegg_loc, gene_ko_link_loc),
                       format='fasta', into=kegg_mod_loc)
    else:
        kegg_mod_loc = kegg_loc
    # make mmseqsdb from modified kegg fasta
    kegg_mmseqs_db = path.join(output_dir, 'kegg.%s.mmsdb' % download_date)
    make_mmseqs_db(kegg_mod_loc, kegg_mmseqs_db, logger, create_index=True, threads=threads, verbose=verbose)
    LOGGER.info('KEGG database processed')
    return {'kegg': kegg_mmseqs_db}


def process_kofam_hmm(kofam_profile_tar_gz, output_dir=DFLT_OUTPUT_DIR, logger=LOGGER, threads=1, verbose=False):
    kofam_profiles = path.join(output_dir, 'kofam_profiles')
    mkdir(kofam_profiles)
    run_process(['tar', '-xzf', kofam_profile_tar_gz, '-C', kofam_profiles], logger, verbose=verbose)
    merged_kofam_profiles = path.join(output_dir, 'kofam_profiles.hmm')
    merge_files(glob(path.join(kofam_profiles, 'profiles', '*.hmm')), merged_kofam_profiles)
    run_process(['hmmpress', '-f', merged_kofam_profiles], logger, verbose=verbose)
    LOGGER.info('KOfam database processed')
    return {'kofam_hmm': merged_kofam_profiles}


def download_kofam_ko_list(output_dir='.', logger=LOGGER, verbose=False):
    kofam_ko_list_gz = path.join(output_dir, 'kofam_ko_list.tsv.gz')
    url = 'ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz'
    url_http = 'https://www.genome.jp/ftp/db/kofam/ko_list.gz'
    download_file(url, kofam_ko_list_gz, logger, alt_urls=[url_http], verbose=verbose)
    return kofam_ko_list_gz


def download_pfam(output_dir='.', logger=LOGGER, verbose=True):
    pfam_full_zipped = path.join(output_dir, 'Pfam-A.full.gz')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz'
    http_url = 'http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz'
    download_file(url, pfam_full_zipped, logger, alt_urls=[http_url], verbose=verbose)
    return pfam_full_zipped


def download_viral(output_dir='.', logger=LOGGER, viral_files=NUMBER_OF_VIRAL_FILES, verbose=True):
    """Can only download newest version"""
    # download all of the viral protein files, need to know the number of files
    # TODO: Make it so that you don't need to know number of viral files in refseq viral
    faa_base_name = 'viral.%s.protein.faa.gz'
    viral_faa_glob = path.join(output_dir, faa_base_name % '*')
    for number in range(viral_files):
        number += 1
        url = 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.%s.protein.faa.gz' % number
        url_http = 'https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.%s.protein.faa.gz' % number
        output_name= path.join(output_dir, faa_base_name % number)
        download_file(url, output_name, logger, alt_urls=[url_http], verbose=verbose)
    # then merge files from above
    merged_viral_faas = path.join(output_dir, 'viral.merged.protein.faa.gz')
    run_process(['cat %s > %s' % (' '.join(glob(viral_faa_glob)), merged_viral_faas)], logger, shell=True)
    return merged_viral_faas


def download_uniref(output_dir='.', logger=LOGGER, version=DEFAULT_UNIREF_VERSION,
                    threads=10, verbose=True):
    uniref_fasta_zipped = path.join(output_dir, 'uniref%s.fasta.gz' % version)
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' % \
          (version, version)
    url_http = 'https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref%s/uniref%s.fasta.gz' % \
                 (version, version)
    download_file(url, uniref_fasta_zipped, logger, alt_urls=[url_http], verbose=verbose)
    return uniref_fasta_zipped


def download_peptidase(output_dir='.', logger=LOGGER, verbose=True):
    save_name = path.join(output_dir, 'merops_peptidases_nr.faa')
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib'
    url_http = 'https://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib'
    download_file(url, save_name, logger, alt_urls=[url_http], verbose=verbose)
    return save_name

def download_vogdb(output_dir='.', logger=LOGGER, version=DEFAULT_VOGDB_VERSION, verbose=True):
    vog_hmm_targz = path.join(output_dir, 'vog.hmm.tar.gz')
    vogdb_url = f'http://fileshare.csb.univie.ac.at/vog/{version}/vog.hmm.tar.gz'
    download_file(vogdb_url, vog_hmm_targz, logger, verbose=verbose)
    return vog_hmm_targz


def process_kofam_ko_list(kofam_ko_list_gz, output_dir='.', logger=LOGGER, threads=1, verbose=False):
    # TODO: fix this so that it is gunzipped to the path
    kofam_ko_list = path.join(output_dir, 'kofam_ko_list.tsv')
    run_process(['gunzip', '-c', kofam_ko_list_gz], logger, save_output=kofam_ko_list, verbose=verbose)
    LOGGER.info('KOfam ko list processed')
    return {'kofam_ko_list': kofam_ko_list}


def process_uniref(uniref_fasta_zipped, output_dir='.', logger=LOGGER,
                   version=DEFAULT_UNIREF_VERSION, threads=10,
                   verbose=True):
    uniref_mmseqs_db = path.join(output_dir, 'uniref%s.%s.mmsdb' % (version, get_iso_date()))
    make_mmseqs_db(uniref_fasta_zipped, uniref_mmseqs_db, logger, create_index=True, threads=threads, verbose=verbose)
    LOGGER.info('UniRef database processed')
    return {'uniref': uniref_mmseqs_db}


def process_mmspro(full_alignment, output_dir, logger=LOGGER, db_name=DEFAULT_MMMSPRO_DB_NAME, threads=10, verbose=True):
    mmseqs_msa = path.join(output_dir, '%s.mmsmsa' % db_name)
    run_process(['mmseqs', 'convertmsa', full_alignment, mmseqs_msa], logger, verbose=verbose)
    mmseqs_profile = path.join(output_dir, '%s.mmspro' % db_name)
    run_process(['mmseqs', 'msa2profile', mmseqs_msa, mmseqs_profile, '--match-mode', '1', '--threads', str(threads)],
                logger,
                verbose=verbose)
    tmp_dir = path.join(output_dir, 'tmp')
    run_process(['mmseqs', 'createindex', mmseqs_profile, tmp_dir, '-k', '5', '-s', '7', '--threads', str(threads)],
                logger,
                verbose=verbose)
    return mmseqs_profile


def process_pfam(pfam_full_zipped, output_dir='.', logger=LOGGER, threads=10, verbose=True):
    pfam_profile = process_mmspro(pfam_full_zipped, output_dir, logger, 'pfam', threads, verbose)
    LOGGER.info('PFAM database processed')
    return {'pfam': pfam_profile}


def process_dbcan(input, output_dir, logger=LOGGER, verbose=True, threads=1):
    output = path.join(output_dir, path.basename(input))
    move(input, output)
    run_process(['hmmpress', '-f', output], logger, verbose=verbose)
    LOGGER.info('dbCAN database processed')
    return {'dbcan': output}


def process_viral(merged_viral_faas, output_dir='.', logger=LOGGER, viral_files=NUMBER_OF_VIRAL_FILES, threads=10, verbose=True):
    refseq_viral_mmseqs_db = path.join(output_dir, 'refseq_viral.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(merged_viral_faas, refseq_viral_mmseqs_db, logger, create_index=True, threads=threads, verbose=verbose)
    LOGGER.info('RefSeq viral database processed')
    return {'viral': refseq_viral_mmseqs_db}



def process_peptidase(peptidase_faa, output_dir='.', logger=LOGGER, threads=10, verbose=True):
    peptidase_mmseqs_db = path.join(output_dir, 'peptidases.%s.mmsdb' % get_iso_date())
    make_mmseqs_db(peptidase_faa, peptidase_mmseqs_db, logger, create_index=True, threads=threads, verbose=verbose)
    LOGGER.info('MEROPS database processed')
    return {'peptidase': peptidase_mmseqs_db}


def process_vogdb(vog_hmm_targz, output_dir='.', logger=LOGGER, version=DEFAULT_VOGDB_VERSION, threads=1, verbose=True):
    hmm_dir = path.join(output_dir, 'vogdb_hmms')
    mkdir(hmm_dir)
    vogdb_targz = tarfile.open(vog_hmm_targz)
    vogdb_targz.extractall(hmm_dir)
    vog_hmms = path.join(output_dir, f'vog_{version}_hmms.txt')
    merge_files(glob(path.join(hmm_dir, 'hmm/VOG*.hmm')), vog_hmms)
    run_process(['hmmpress', '-f', vog_hmms], logger, verbose=verbose)
    LOGGER.info('VOGdb database processed')
    return {'vogdb': vog_hmms}


def download_vog_annotations(output_dir, logger=LOGGER, version=DEFAULT_VOGDB_VERSION, verbose=True):
    vog_annotations = path.join(output_dir, 'vog_annotations_%s.tsv.gz' % version)
    download_file('http://fileshare.csb.univie.ac.at/vog/%s/vog.annotations.tsv.gz' % version,
                  vog_annotations, logger, verbose=verbose)
    return vog_annotations


def download_genome_summary_form(output_dir, logger=LOGGER, branch='master', verbose=True):
    genome_summary_form = path.join(output_dir, 'genome_summary_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/genome_summary_form.tsv' % branch,
                  genome_summary_form, logger, verbose=verbose)
    return genome_summary_form


def download_module_step_form(output_dir, logger=LOGGER, branch='master', verbose=True):
    function_heatmap_form = path.join(output_dir, 'module_step_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/module_step_form.tsv' % branch,
                  function_heatmap_form, logger, verbose=verbose)
    return function_heatmap_form


def download_etc_module_database(output_dir, logger=LOGGER, branch='master', verbose=True):
    etc_module_database = path.join(output_dir, 'etc_mdoule_database.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/etc_module_database.tsv' % branch,
                  etc_module_database, logger, verbose=verbose)
    return etc_module_database


def download_function_heatmap_form(output_dir, logger=LOGGER, branch='master', verbose=True):
    function_heatmap_form = path.join(output_dir, 'function_heatmap_form.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/function_heatmap_form.tsv' % branch,
                  function_heatmap_form, logger, verbose=verbose)
    return function_heatmap_form


def download_amg_database(output_dir, logger=LOGGER, branch='master', verbose=True):
    amg_database = path.join(output_dir, 'amg_database.%s.tsv' % get_iso_date())
    download_file('https://raw.githubusercontent.com/WrightonLabCSU/DRAM/%s/data/amg_database.tsv' % branch, 
                  amg_database, logger, verbose=verbose,)
    return amg_database


def check_file_exists(*paths):
    for i in paths:
        if i is None:
            continue
        elif path.isfile(i):
            continue
        else:
            raise ValueError(f"Database location does not exist: {i}")


def prepare_databases(output_dir, loggpath=None, kegg_loc=None, gene_ko_link_loc=None, kofam_hmm_loc=None,
                      kofam_ko_list_loc=None, kegg_download_date=None, uniref_loc=None,
                      uniref_version=DEFAULT_UNIREF_VERSION, pfam_loc=None, pfam_hmm_loc=None,
                      dbcan_loc=None, dbcan_fam_activities:str=None, dbcan_subfam_ec:str=None, dbcan_version=DEFAULT_DBCAN_RELEASE,
                      dbcan_date=DEFAULT_DBCAN_DATE, viral_loc=None, peptidase_loc=None,
                      vogdb_loc=None, vogdb_version=DEFAULT_VOGDB_VERSION, vog_annotations=None,
                      genome_summary_form_loc=None, module_step_form_loc=None,
                      etc_module_database_loc=None, function_heatmap_form_loc=None,
                      amg_database_loc=None,
                      number_of_viral_files=NUMBER_OF_VIRAL_FILES,
                      skip_uniref=False, keep_database_files=False,
                      branch='master', threads=10, verbose=True, select_db=None, clear_config=False):


    dram_settings = {
        'kegg':                  {'name': 'KEGG db',
                                  'description_db_updated': 'Unknown, or Never',
                                  'citation': KEGG_CITATION},
        'gene_ko_link':          {'name': 'KEGG Gene KO link', 'citation': KEGG_CITATION},
        "kofam_hmm":             {'name': 'KOfam db', 'citation': KOFAM_CITATION},
        "kofam_ko_list":         {'name': 'KOfam KO list', 'citation': KOFAM_CITATION},
        'uniref':                {'name': 'UniRef db',
                                  'description_db_updated': 'Unknown, or Never',
                                  'citation': UNIREF_CITATION},
        'pfam':                  {'name': 'Pfam db', 'citation': PFAM_CITATION},
        "pfam_hmm":          {'name': 'Pfam hmm dat',
                                  'description_db_updated': 'Unknown, or Never',
                                  'citation': PFAM_CITATION},
        "dbcan":                 {'name': 'dbCAN db', 'citation': DBCAN_CITATION},
        "dbcan_fam_activities":  {'name': 'dbCAN family activities', 'citation': DBCAN_CITATION},
        "dbcan_subfam_ec":       {'name': 'dbCAN subfamily EC numbers', 'citation': DBCAN_CITATION},
        "vogdb":                 {'name': 'VOGDB db', 'citation': VOGDB_CITATION},
        "vog_annotations":       {'name': 'VOG annotations',
                                  'description_db_updated': 'Unknown, or Never',
                                  'citation': VOGDB_CITATION},
        "viral":          {'name': 'RefSeq Viral db',
                                  'description_db_updated': 'Unknown, or Never',
                                  'citation': VIRAL_REFSEQ_CITATION},
        "peptidase":             {'name': 'MEROPS peptidase db',
                                  'description_db_updated': 'Unknown, or Never',
                                  'citation': PEPTIDASE_CITATION},
        "genome_summary_form":   {'name': 'Genome summary form'},
        "module_step_form":      {'name': 'Module step form'},
        "function_heatmap_form": {'name': 'Function heatmap form'},
        "amg_database":          {'name': 'AMG database'},
        "etc_module_database":   {'name': 'ETC module database'},
    }
    database_settings= {
        'kegg':                  {},
        "gene_ko_link":          {},
        "kofam_hmm":             {},
        "kofam_ko_list":         {},
        'uniref':                {'version': uniref_version},
        'pfam':                  {},
        "pfam_hmm":          {},
        "dbcan":                 {'version': dbcan_version},
        "dbcan_fam_activities":  {'version': dbcan_version, 'upload_date': dbcan_date},
        "dbcan_subfam_ec":       {'version': dbcan_version, 'upload_date': dbcan_date},
        "vogdb":                 {'version': vogdb_version},
        "vog_annotations":       {'version': vogdb_version},
        "viral":                 {},
        "peptidase":             {},
        "genome_summary_form":   {"branch": branch},
        "module_step_form":      {"branch": branch},
        "function_heatmap_form": {"branch": branch},
        "amg_database":          {"branch": branch},
        "etc_module_database":   {"branch": branch},
    }
    process_settings = {
        'kegg': {},
        "gene_ko_link": {},
        "kofam_hmm": {},
        "kofam_ko_list": {},
        'uniref': {"version": uniref_version},
        'pfam': {},
        "pfam_hmm": {},
        "dbcan": {},
        "dbcan_fam_activities": {},
        "dbcan_subfam_ec": {},
        "vogdb": {},
        "vog_annotations": {},
        "viral": {'viral_files': number_of_viral_files},
        "peptidase": {},
        "genome_summary_form": {},
        "module_step_form": {},
        "function_heatmap_form": {},
        "amg_database": {},
        "etc_module_database": {},
    }

    # setup temp, logging, and db_handler
    if not path.isdir(output_dir):
        mkdir(output_dir)
    temporary = path.join(output_dir, 'database_files')
    mkdir(temporary)
    main_log = path.join(output_dir, 'database_processing.log')
    setup_logger(LOGGER, *[(main_log, loggpath) if loggpath is not None else main_log])
    db_handler = DatabaseHandler(logger=LOGGER)
    if clear_config or select_db is None:
        db_handler.clear_config()

    db_handler.config['log_path'] = main_log
    db_handler.write_config()
    LOGGER.info('Starting the process of downloading data')

    if skip_uniref:
        LOGGER.info('Skipping UniRef')
        del database_settings['uniref']
    
    
    locs = {remove_suffix(i, '_loc'):j for i, j in locals().items() if i.endswith('_loc') and j is not None}
    download_functions = {i[9:]:j for i,j in globals().items() if callable(j) and i.startswith('download_')}
    process_functions = {i[8:]:j for i,j in globals().items() if callable(j) and i.startswith('process_')}
    functions = {i:j for i,j in globals().items() if callable(j) and i.startswith('download_')}

    # Check any specified paths exist
    missing_user_inputs = [f"{i} at {j}" for i, j in locs.items() if not path.exists(j)]
    if len(missing_user_inputs) > 1:
        raise ValueError(f"The fallowing user provided paths don't seem to exist: {', '.join(missing_user_inputs)}")

    un_obtainable = [i for i in database_settings if i not in locs and i not in download_functions]

    for i in un_obtainable:
        LOGGER.info(f"The {i}_loc argument was not used to specify a downloaded {i} file, and dram can not"
                     " download it its self. So it is assumed that the user wants to set up DRAM without it")
        del database_settings[i]


    # check inputs
    if skip_uniref and uniref_loc is not None:
        raise ValueError('Cannot skip UniRef processing and provide a location of UniRef.'
                         ' Skipping UniRef will cause provided UniRef file to not be used.')

    if select_db is not None:
        miss_name = [i for i in select_db if i not in database_settings]
        user_inputs = [i for i in locs if i not in select_db]
        if len(miss_name) > 0:
            LOGGER.error("Only the databases in the db list can be pased to select_db, "
                            f"you passed {miss_name} which is/are not in the list,"
                            f" select from {list(database_settings.keys())}.")
            raise ValueError('Bad user input, see log')

        if len(user_inputs) > 0:
            LOGGER.error(f"The user provided location for {user_inputs}, but required it not be used by proving"
                             f" the select_db argument for other databases. This would suggest that the"
                             " user may have made a mistake and so this error is rased.")
            raise ValueError('Bad user input, see log')

        database_settings = {i:database_settings[i] for i in select_db}


    LOGGER.info('Database preparation started')

    # Download DBs
    for i, j in database_settings.items():
        if locs.get(i) is None:
            LOGGER.info(f"Downloading {i}")
            if i in process_functions:
                locs[i] = download_functions[i](
                    temporary, LOGGER, **j, verbose=verbose)
            else:
                locs[i] = download_functions[i](
                    output_dir, LOGGER, **j, verbose=verbose)
            j['Download time'] = datetime.now().strftime("%m/%d/%Y, %H:%M:%S")
            j['Origin'] = "Downloaded by DRAM"
        else:
            j = {k: "Unknown" for k in j}
            j['Download time'] = "Unknown"
            j['Origin'] = "Provided by user"
            j['Original path'] = locs[i]
            if i not in process_functions:
                LOGGER.info(f"Copying {locs[i]} to output_dir")
                locs[i] = copyfile(locs[i], path.join(output_dir, path.basename(locs[i])))

    LOGGER.info("All raw data files were downloaded successfully")

    # Process databases
    for i in locs:
        processed_locs = {}
        if i in process_functions:
            LOGGER.info(f"Processing {i}")
            processed_locs = process_functions[i](locs[i], output_dir, LOGGER,
                                                       threads=threads, verbose=verbose, **process_settings[i])
        else:
            processed_locs = {i:locs[i]}
        for k, v in processed_locs.items():
            final_dest = path.join(output_dir, path.basename(v))
            if v != final_dest:
                for db_file in glob('%s*' % v):
                    move(db_file, path.join(output_dir, path.basename(db_file)))
                v = path.join(output_dir, path.basename(v))
            # update_dram_forms the settings per OUTPUT fill, including the process_settings
            #  and database_settings, which are per input file.
            if db_handler.config.get('setup_info') is None:
                db_handler.config['setup_info'] = {}
            db_handler.config['setup_info'][k] = {**dram_settings[k], **process_settings[i],
                                                  **database_settings[i]}
            db_handler.set_database_paths(**{f"{k}_loc":v})
            db_handler.write_config()
            LOGGER.info(f'Moved {k} to final destination, configuration updated')

    LOGGER.info('Populating the description db, this may take some time')
    db_handler.config['description_db'] = path.realpath(path.join(output_dir, 'description_db.sqlite'))
    db_handler.write_config()
    db_handler.populate_description_db(db_handler.config['description_db'], select_db, update_config=False)
    # todo make db handler such that the destruction on success write_config
    db_handler.write_config()
    LOGGER.info('DRAM description database populated')

    if not keep_database_files:
        rmtree(temporary)
    LOGGER.info('Database preparation completed')


def update_dram_forms(output_dir, branch='master'):
    if not path.isdir(output_dir):
        mkdir(output_dir)
    form_locs = dict()
    form_locs['genome_summary_form_loc'] = download_genome_summary_form(output_dir, branch)
    form_locs['module_step_form_loc'] = download_module_step_form(output_dir, branch)
    form_locs['etc_module_database_loc'] = download_etc_module_database(output_dir, branch)
    form_locs['function_heatmap_form_loc'] = download_function_heatmap_form(output_dir, branch)
    form_locs['amg_database_loc'] = download_amg_database(output_dir, branch)
    main_log = path.join(output_dir, 'database_processing.log')
    setup_logger(LOGGER, main_log)
    db_handler = DatabaseHandler(LOGGER)
    db_handler.set_database_paths(**form_locs)
