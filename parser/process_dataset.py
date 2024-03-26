import argparse
import sys
import os
import socket
import traceback

import requests
import time
import ftplib
import logging.config
import gc
import shutil
from urllib.parse import urlparse

from parser.MzIdParser import MzIdParser
import logging.config
from parser.api_writer import APIWriter
from config.config_parser import get_conn_str
from parser.database_writer import DatabaseWriter

logging_config_file = os.path.join(os.path.dirname(__file__), '../config/logging.ini')
logging.config.fileConfig(logging_config_file)
logger = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description='Process mzIdentML files in a dataset and load them into a relational database.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-p', '--pxid', nargs='+',
                       help='proteomeXchange accession, should be of the form PXDnnnnnn or numbers only', )
    group.add_argument('-f', '--ftp',
                       help='process files from specified ftp location, e.g. ftp://ftp.jpostdb.org/JPST001914/')
    group.add_argument('-d', '--dir',
                       help='process files in specified local directory, e.g. /home/user/data/JPST001914')
    parser.add_argument('-i', '--identifier',
                        help='identifier to use for dataset (if providing '
                             'proteome exchange accession these are always used instead and this arg is ignored)')
    parser.add_argument('--dontdelete', action='store_true', help='Do not delete downloaded data after processing')
    parser.add_argument('-t', '--temp', action='store_true', help='Temp folder to download data files into')
    parser.add_argument('-n', '--nopeaklist',
                        help='No peak list files available, only works in comination with --dir arg',
                        action='store_true')
    parser.add_argument('-w', '--writer', help='Save data to database(-w db) or API(-w api)')
    args = parser.parse_args()
    try:
        logger.info("process_dataset.py is running!")
        print("process_dataset.py is running!")
        if args.temp:
            temp_dir = os.path.expanduser(args.temp)
        else:
            temp_dir = os.path.expanduser('~/mzId_convertor_temp')

        if args.writer:
            writer_method = args.writer
            if not (writer_method.lower() == 'api' or writer_method.lower() == 'db'):
                raise ValueError('Writer method not supported! please use "api" or "database"')

        if args.pxid:
            px_accessions = args.pxid
            for px_accession in px_accessions:
                # convert_pxd_accession(px_accession, temp_dir, args.dontdelete)
                convert_pxd_accession_from_pride(px_accession, temp_dir, writer_method, args.dontdelete)
        elif args.ftp:
            ftp_url = args.ftp
            if args.identifier:
                project_identifier = args.identifier
            else:
                parsed_url = urlparse(ftp_url)
                project_identifier = parsed_url.path.rsplit("/", 1)[-1]
            convert_from_ftp(ftp_url, temp_dir, project_identifier, writer_method, args.dontdelete)
        else:
            local_dir = args.dir
            if args.identifier:
                project_identifier = args.identifier
            else:
                project_identifier = local_dir.rsplit("/", 1)[-1]
            convert_dir(local_dir, project_identifier, writer_method, nopeaklist=args.nopeaklist)
        sys.exit(0)
    except Exception as ex:
        logger.error(ex)
        traceback.print_stack(ex)
        sys.exit(1)


def convert_pxd_accession(px_accession, temp_dir, dont_delete=False):
    # get ftp location from PX
    px_url = 'https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=' + px_accession + '&outputMode=JSON'
    logger.info('GET request to ProteomeExchange: ' + px_url)
    px_response = requests.get(px_url)
    r = requests.get(px_url)
    if r.status_code == 200:
        logger.info('ProteomeExchange returned status code 200')
        px_json = px_response.json()
        ftp_url = None
        for dataSetLink in px_json['fullDatasetLinks']:
            # name check is necessary because some things have wrong acc, e.g. PXD006574
            if dataSetLink['accession'] == "MS:1002852" or dataSetLink['name'] == "Dataset FTP location":
                ftp_url = dataSetLink['value']
                convert_from_ftp(ftp_url, temp_dir, px_accession, dont_delete)
                break
        if not ftp_url:
            raise Exception('Error: Dataset FTP location not found in ProteomeXchange response')
    else:
        raise Exception('Error: ProteomeXchange returned status code ' + str(px_response.status_code))


def convert_pxd_accession_from_pride(px_accession, temp_dir, writer_method, dont_delete=False):
    # get ftp location from PRIDE API
    px_url = 'https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession=' + px_accession
    logger.info('GET request to PRIDE API: ' + px_url)
    pride_response = requests.get(px_url)
    r = requests.get(px_url)
    if r.status_code == 200:
        logger.info('PRIDE API returned status code 200')
        pride_json = pride_response.json()
        ftp_url = None

        if len(pride_json) > 0:
            for protocol in pride_json[0]['publicFileLocations']:
                if protocol['name'] == "FTP Protocol":

                    # Parse the FTP path
                    parsed_url = urlparse(protocol['value'])

                    # Split the path into segments
                    path_segments = parsed_url.path.split('/')

                    # Construct the parent folder path
                    parent_folder = parsed_url.scheme + '://' + parsed_url.netloc

                    for segment in path_segments[:-1]:  # Exclude the filename with extension part
                        parent_folder += segment + '/'
                    ftp_url = parent_folder

                    logger.info('PRIDE FTP path : ' + parent_folder)
                    break;
        convert_from_ftp(ftp_url, temp_dir, px_accession,writer_method, dont_delete)
        if not ftp_url:
            raise Exception('Error: Public File location not found in PRIDE API response')
    else:
        raise Exception('Error: PRIDE API returned status code ' + str(pride_response.status_code))


def convert_from_ftp(ftp_url, temp_dir, project_identifier, writer_method, dont_delete):
    if not ftp_url.startswith('ftp://'):
        raise Exception('Error: FTP location must start with ftp://')
    if not os.path.isdir(temp_dir):
        try:
            os.mkdir(temp_dir)
        except OSError as e:
            logger.error('Failed to create temp directory ' + temp_dir)
            logger.error('Error: ' + e.strerror)
            raise e
    logger.info('FTP url: ' + ftp_url)
    parsed_url = urlparse(ftp_url)
    path = os.path.join(temp_dir, project_identifier)
    try:
        os.mkdir(path)
    except OSError:
        pass
    ftp_ip = socket.getaddrinfo(parsed_url.hostname, 21)[0][4][0]
    files = get_ftp_file_list(ftp_ip, parsed_url.path)
    for f in files:
        # check file not already in temp dir
        if not (os.path.isfile(os.path.join(path, f))
                or f.lower == "generated"  # dunno what these files are but they seem to make ftp break
                or f.lower().endswith('raw')
                or f.lower().endswith('raw.gz')
                or f.lower().endswith('all.zip')
                or f.lower().endswith('csv')
                or f.lower().endswith('txt')):
            logger.info('Downloading ' + f + ' to ' + path)
            ftp = get_ftp_login(ftp_ip)
            try:
                ftp.cwd(parsed_url.path)
                ftp.retrbinary("RETR " + f, open(os.path.join(path, f), 'wb').write)
                ftp.quit()
            except ftplib.error_perm as e:
                ftp.quit()
                # error_msg = "%s: %s" % (f, e.args[0])
                # self.logger.error(error_msg)
                raise e
    convert_dir(path, project_identifier, writer_method)
    if not dont_delete:
        # remove downloaded files
        try:
            shutil.rmtree(path)
        except OSError as e:
            logger.error('Failed to delete temp directory ' + path)
            logger.error('Error: ' + e.strerror)
            raise e


def get_ftp_login(ftp_ip):
    time.sleep(10)
    try:
        ftp = ftplib.FTP(ftp_ip)
        ftp.login()  # Uses password: anonymous@
        return ftp
    except ftplib.all_errors as e:
        logger.error('FTP fail at ' + time.strftime("%c"))
        raise e


def get_ftp_file_list(ftp_ip, ftp_dir):
    ftp = get_ftp_login(ftp_ip)
    try:
        ftp.cwd(ftp_dir)
    except ftplib.error_perm as e:
        error_msg = "%s: %s" % (ftp_dir, e.args[0])
        logger.error(error_msg)
        ftp.quit()
        raise e
    try:
        filelist = ftp.nlst()
    except ftplib.error_perm as resp:
        if str(resp) == "550 No files found":
            logger.info("FTP: No files in this directory")
        else:
            error_msg = "%s: %s" % (ftp_dir, ftplib.error_perm.args[0])
            logger.error(error_msg)
        raise resp
    ftp.close()
    return filelist


def convert_dir(local_dir, project_identifier,writer_method, nopeaklist=False):
    # logging.basicConfig(level=logging.DEBUG,
    #                     format='%(asctime)s %(levelname)s %(name)s %(message)s')
    # logger = logging.getLogger(__name__)
    peaklist_dir = local_dir if not nopeaklist else None
    #  iterate over files in local_dir
    for file in os.listdir(local_dir):
        if file.endswith(".mzid") or file.endswith(".mzid.gz"):
            logger.info("Processing " + file)
            conn_str = get_conn_str()
            if writer_method.lower() == 'api':
                writer = APIWriter(conn_str, pxid=project_identifier)
            elif writer_method.lower() == 'db':
                writer = DatabaseWriter(conn_str, pxid=project_identifier)
            id_parser = MzIdParser(os.path.join(local_dir, file), local_dir, peaklist_dir, writer, logger)
            try:
                id_parser.parse()
                # logger.info(id_parser.warnings + "\n")
            except Exception as e:
                print("Error parsing " + file)
                print(type(e).__name__, e)
                raise e
            gc.collect()
        else:
            continue


if __name__ == "__main__":
   main()
