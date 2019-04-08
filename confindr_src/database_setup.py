#!/usr/bin/env python

from rauth import OAuth1Session
from Bio import SeqIO
import argparse
import datetime
import logging
import urllib.request
import tarfile
import shutil
import glob
import csv
import re
import os


class RmlstRest(object):

    def get_session_token(self):
        session_request = OAuth1Session(self.consumer_key,
                                        self.consumer_secret,
                                        access_token=self.access_token,
                                        access_token_secret=self.access_secret)
        url = self.test_rest_url + '/oauth/get_session_token'
        r = session_request.get(url)
        if r.status_code == 200:
            self.session_token = r.json()['oauth_token']
            self.session_secret = r.json()['oauth_token_secret']
        else:
            logging.error('ERROR: Couldn\'t get a session token for rMLST database download. Check that your consumer '
                          'secret and access token files have valid credentials and try again.')
            quit(code=1)

    def get_loci_and_scheme_url(self):
        session = OAuth1Session(self.consumer_key,
                                self.consumer_secret,
                                access_token=self.session_token,
                                access_token_secret=self.session_secret)
        r = session.get(self.test_rest_url)
        if r.status_code == 200 or r.status_code == 201:
            if re.search('json', r.headers['content-type'], flags=0):
                decoded = r.json()
            else:
                decoded = r.text
            # Extract the URLs from the returned data
            self.loci = decoded['loci']
            self.profile = decoded['schemes']
        else:
            logging.error('ERROR: Could not find URLs for rMLST download, they may have moved. Please open an issue '
                          'at https://github.com/OLC-Bioinformatics/ConFindr/issues and we\'ll get things sorted out.')
            quit(code=1)

    def download_loci(self):
        session = OAuth1Session(self.consumer_key,
                                self.consumer_secret,
                                access_token=self.session_token,
                                access_token_secret=self.session_secret)
        r = session.get(self.loci)
        if r.status_code == 200 or r.status_code == 201:
            if re.search('json', r.headers['content-type'], flags=0):
                decoded = r.json()
            else:
                decoded = r.text
            # Extract all the URLs in the decoded dictionary under the key 'loci'
            for locus_url in decoded['loci']:
                output_file = os.path.join(self.output_folder, '{}.tfa'.format(os.path.split(locus_url)[1]))
                logging.info('Downloading {}...'.format(os.path.split(locus_url)[1]))
                with open(output_file, 'w') as f:
                    download = session.get(locus_url + '/alleles_fasta')
                    if download.status_code == 200 or download.status_code == 201:
                        if re.search('json', download.headers['content-type'], flags=0):
                            decoded = download.json()
                        else:
                            decoded = download.text
                        with open(output_file, 'w') as locus_fasta:
                            locus_fasta.write(decoded)
                
        else:
            logging.error('ERROR: Could not find URLs for rMLST download, they may have moved. Please open an issue '
                          'at https://github.com/OLC-Bioinformatics/ConFindr/issues and we\'ll get things sorted out.')
            quit(code=1)

    def download_profile(self):
        profile_file = os.path.join(self.output_folder, 'profiles.txt')
        session = OAuth1Session(self.consumer_key,
                                self.consumer_secret,
                                access_token=self.session_token,
                                access_token_secret=self.session_secret)
        r = session.get(self.profile + '/1/profiles_csv')
        logging.info('Downloading rMLST profiles...')
        if r.status_code == 200 or r.status_code == 201:
            if re.search('json', r.headers['content-type'], flags=0):
                decoded = r.json()
            else:
                decoded = r.text
            # Write the profile file to disk
            with open(profile_file, 'w') as profile:
                profile.write(decoded)

    def get_request_token(self):
        session = OAuth1Session(consumer_key=self.consumer_key,
                                consumer_secret=self.consumer_secret)
        # Use the test URL in the GET request
        r = session.request(method='GET',
                            url=self.request_token_url,
                            params={'oauth_callback': 'oob'})
        if r.status_code == 200:
            self.request_token = r.json()['oauth_token']
            self.request_secret = r.json()['oauth_token_secret']

    def get_access_token(self):
        authorize_url = self.test_web_url + '&page=authorizeClient&oauth_token=' + self.request_token
        print('Visit this URL in your browser: ' + authorize_url)
        verifier = input('Enter oauth_verifier from browser: ')
        session_request = OAuth1Session(consumer_key=self.consumer_key,
                                        consumer_secret=self.consumer_secret,
                                        access_token=self.request_token,
                                        access_token_secret=self.request_secret)
        # Perform a GET request with the appropriate keys and tokens
        r = session_request.get(self.access_token_url,
                                params={
                                    'oauth_verifier': verifier
                                })
        # If the status code is '200' (OK), proceed
        if r.status_code == 200:
            # Save the JSON-decoded token secret and token
            self.access_token = r.json()['oauth_token']
            self.access_secret = r.json()['oauth_token_secret']

    def __init__(self, consumer_secret_file, output_folder):
        self.test_rest_url = 'http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef'
        self.test_web_url = 'http://pubmlst.org/cgi-bin/bigsdb/bigsdb.pl?db=pubmlst_rmlst_seqdef'
        self.request_token_url = self.test_rest_url + '/oauth/get_request_token'
        self.access_token_url = self.test_rest_url + '/oauth/get_access_token'
        self.authorize_url = self.test_web_url + '&page=authorizeClient'
        self.output_folder = output_folder
         
        # Get the consumer secret set up.
        if not os.path.isfile(consumer_secret_file):
            logging.error('ERROR: Could not find consumer secret file. Please make sure the file you specified ({}) exists '
                          'and try again.'.format(consumer_secret_file))
            quit(code=1)
        with open(consumer_secret_file) as f:
            lines = f.readlines()
        try:
            self.consumer_key = lines[0].rstrip()
            self.consumer_secret = lines[1].rstrip()
        except IndexError:
            logging.error('ERROR: Could not parse your consumer secret file. File should have supplied consumer key '
                          'on first line, and consumer secret on the second line.')
            quit(code=1) 
            
        self.session_secret = str()
        self.session_token = str()
        self.loci = str()
        self.profile = str()
        self.request_token = str()
        self.request_secret = str()
        self.access_token = str()
        self.access_secret = str()


def create_gene_allele_file(profiles_file, gene_allele_file):
    genus_allele_info = dict()
    with open(profiles_file) as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            genus = row['genus']
            if genus not in genus_allele_info:
                genus_allele_info[genus] = list()
            for i in range(1, 66):
                if i < 10:
                    gene = 'BACT00000' + str(i)
                else:
                    gene = 'BACT0000' + str(i)
                if gene in row:
                    allele_number = row[gene]
                    gene_allele = '{}_{}'.format(gene, allele_number)
                    if allele_number != 'N' and gene_allele not in genus_allele_info[genus]:
                        genus_allele_info[genus].append(gene_allele)
    with open(gene_allele_file, 'w') as f:
        for genus in genus_allele_info:
            f.write(str(genus) + ':')
            for allele in genus_allele_info[genus]:
                f.write(str(allele) + ',')
            f.write('\n')


def setup_confindr_database(output_folder, consumer_secret):
    # Remove previous output folder if it existed.

    # Go through the REST API in order to get profiles downloaded.
    rmlst_rest = RmlstRest(consumer_secret_file=consumer_secret,
                           output_folder=output_folder)
    rmlst_rest.get_request_token()
    rmlst_rest.get_access_token()
    rmlst_rest.get_session_token()
    rmlst_rest.get_loci_and_scheme_url()
    rmlst_rest.download_loci()
    rmlst_rest.download_profile()

    # With the sequences downloaded, make a file of all rMLST sequences combined.
    logging.info('Combining rMLST files...')
    with open(os.path.join(output_folder, 'rMLST_combined.fasta'), 'w') as f:
        locus_files = sorted(glob.glob(os.path.join(output_folder, 'BACT*.tfa')))
        for locus_file in locus_files:
            for record in SeqIO.parse(locus_file, 'fasta'):
                record.id = record.id.replace('-', '_')
                record.seq._data = record.seq._data.replace('-', '').replace('N', '')
                record.name = ''
                record.description = ''
                SeqIO.write(record, f, 'fasta')
            # Clean up individual file.
            try:
                os.remove(locus_file)
            except OSError:
                logging.warning('WARNING: Could not delete {}. This won\'t affect ConFindr performance, but '
                                ' you may want to delete it to save on disk space.'.format(locus_file))

    logging.info('Assigning alleles to genera...')
    # Parse profiles so that we know what alleles are found with each genus.
    create_gene_allele_file(profiles_file=os.path.join(output_folder, 'profiles.txt'),
                            gene_allele_file=os.path.join(output_folder, 'gene_allele.txt'))
    
    
def download_mash_sketch(output_folder):
    logging.info('Downloading mash refseq sketch...')
    urllib.request.urlretrieve('https://github.com/OLC-Bioinformatics/ConFindr/raw/master/refseq_sketch/refseq.msh',
                               os.path.join(output_folder, 'refseq.msh'))


def download_cgmlst_derived_data(output_folder):
    logging.info('Downloading cgMLST-derived data for Escherichia, Salmonella, and Listeria...')
    urllib.request.urlretrieve('https://ndownloader.figshare.com/files/14771267',
                               os.path.join(output_folder, 'confindr_db.tar.gz'))
    confindr_tar = os.path.join(output_folder, 'confindr_db.tar.gz')
    tar = tarfile.open(confindr_tar)
    tar.extractall(path=output_folder)
    tar.close()
    os.remove(confindr_tar)


def main():
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output_folder',
                        default=os.environ.get('CONFINDR_DB', os.path.expanduser('~/.confindr_db')),
                        help='Path to download databases to - if folder does not exist, will be created. If folder does '
                             'exist, will be deleted and updated sequences downloaded. Defaults to ~/.confindr_db, or '
                             'the CONFINDR_DB environmental variable.')
    parser.add_argument('-s', '--secret_file',
                        type=str,
                        help='Path to consumer secret file for rMLST database.')
    args = parser.parse_args()
    if os.path.isdir(args.output_folder):
        logging.info('Removing old databases...')
        shutil.rmtree(args.output_folder)
    os.makedirs(args.output_folder)
    download_cgmlst_derived_data(args.output_folder)
    if args.secret_file is None:
        logging.warning('WARNING: Without an rMLST secret file, data will only be downloaded for Escherichia, '
                        'Salmonella, and Listeria. See '
                        'https://olc-bioinformatics.github.io/ConFindr/install/#downloading-confindr-databases for '
                        'instructions on how to get access to rMLST databases so ConFindr can be used for other species'
                        ' as well')
    else:
        setup_confindr_database(args.output_folder,
                                args.secret_file)
    download_mash_sketch(args.output_folder)
    current_year = datetime.datetime.utcnow().year
    current_month = datetime.datetime.utcnow().month
    current_day = datetime.datetime.utcnow().day
    with open(os.path.join(args.output_folder, 'download_date.txt'), 'w') as f:
        f.write('{}-{}-{}'.format(current_year, current_month, current_day))
    logging.info('Done downloading ConFindr databases!')
    
    
if __name__ == '__main__':
    main()
