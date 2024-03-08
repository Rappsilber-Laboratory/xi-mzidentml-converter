import argparse
import ntpath
import os
import sys
import traceback

import logging.config
import gc
from shutil import copyfile

from parser import NoPeakListsCsvParser
from parser.writer import Writer
from db_config_parser import get_conn_str
import logging.config

logging.config.fileConfig("logging.ini")
logger = logging.getLogger(__name__)


def main(args):
    # logging.basicConfig(level=logging.DEBUG,
    #                     format='%(asctime)s %(levelname)s %(name)s %(message)s')
    # logger = logging.getLogger(__name__)
    # if args.temp:
    #     temp_dir = os.path.expanduser(args.temp)
    # else:
    temp_dir = os.path.expanduser('~/mzId_convertor_temp')

    # copy fasta file to tmpdir so it is being read by the parser
    copyfile(args.fasta, os.path.join(str(temp_dir), ntpath.basename(args.fasta)))

    file = args.csv
    logger.info("Processing " + file)
    conn_str = get_conn_str()
    writer = Writer(conn_str, pxid=os.path.basename(file))
    # id_parser = MzIdParser(os.path.join(local_dir, file), local_dir, peaklist_dir, writer, logger)
    # parse the mzid file
    id_parser = NoPeakListsCsvParser(file, str(temp_dir), None, writer, logger)
    id_parser.check_required_columns()

    try:
        id_parser.parse()
        # logger.info(id_parser.warnings + "\n")
    except Exception as e:
        print("Error parsing " + file)
        print(type(e).__name__, e)
        raise e
    gc.collect()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Process csv a file and fasta and load data into a relational database.')
    parser.add_argument('-c', '--csv', help='CSV file to process', required=True)
    parser.add_argument('-f', '--fasta', help='Associated FASTA file', required=True)
    # parser.add_argument('-p', '--peaklist',
    #                     help='Not implemented yet.')')
    try:
        logger.info("process_csv_dataset.py is running!")
        main(parser.parse_args())
        sys.exit(0)
    except Exception as ex:
        logger.error(ex)
        traceback.print_stack(ex)
        sys.exit(1)
