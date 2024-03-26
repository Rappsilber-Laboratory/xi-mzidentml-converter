from parser import FullCsvParser, NoPeakListsCsvParser, LinksOnlyCsvParser
from parser.api_writer import Writer
from sqlalchemy import text
from uuid import uuid4
from .db_pytest_fixtures import *


def parse_full_csv_into_postgresql(mzid_file, peaklist, tmpdir, logger, use_database, engine):
    # create temp user for user_id
    user_id = 1
    # create writer
    writer = Writer(engine.url, user_id)
    engine.dispose()

    # parse the mzid file
    id_parser = FullCsvParser(mzid_file, str(tmpdir), peaklist, writer, logger)
    id_parser.check_required_columns()
    id_parser.parse()
    return id_parser


def parse_no_peak_lists_csv_into_postgresql(mzid_file, peaklist, tmpdir, logger, use_database, engine):
    # create temp user for user_id
    user_id = 1
    # create writer
    writer = Writer(engine.url, user_id)
    engine.dispose()

    # parse the mzid file
    id_parser = NoPeakListsCsvParser(mzid_file, str(tmpdir), peaklist, writer, logger)
    id_parser.check_required_columns()
    id_parser.parse()
    return id_parser


def parse_links_only_csv_into_postgresql(mzid_file, peaklist, tmpdir, logger, use_database, engine):
    # create temp user for user_id
    user_id = 1
    # create writer
    writer = Writer(engine.url, user_id)
    engine.dispose()

    # parse the mzid file
    id_parser = LinksOnlyCsvParser(mzid_file, str(tmpdir), peaklist, writer, logger)
    id_parser.check_required_columns()
    id_parser.parse()
    return id_parser


def parse_full_csv_into_sqllite(mzid_file, peaklist, tmpdir, logger, use_database, engine):
    # create writer
    writer = Writer(engine.url)
    engine.dispose()

    # parse the mzid file
    id_parser = FullCsvParser(mzid_file, str(tmpdir), peaklist, writer, logger)
    id_parser.check_required_columns()
    id_parser.parse()
    return id_parser


def parse_no_peak_lists_csv_into_sqllite(mzid_file, peaklist, tmpdir, logger, use_database, engine):
    # create writer
    writer = Writer(engine.url)
    engine.dispose()

    # parse the mzid file
    id_parser = NoPeakListsCsvParser(mzid_file, str(tmpdir), peaklist, writer, logger)
    id_parser.check_required_columns()
    id_parser.parse()
    return id_parser


def parse_links_only_csv_into_sqllite(mzid_file, peaklist, tmpdir, logger, use_database, engine):
    # create writer
    writer = Writer(engine.url)
    engine.dispose()

    # parse the mzid file
    id_parser = LinksOnlyCsvParser(mzid_file, str(tmpdir), peaklist, writer, logger)
    id_parser.check_required_columns()
    id_parser.parse()
    return id_parser

# def parse_csv_into_sqlite_xispec(mzid_file, peaklist, tmpdir, logger, engine):
#     # create writer
#     writer = Writer(engine.url)
#     engine.dispose()
#
#     # parse the mzid file
#     id_parser = MzIdParser.xiSPEC_MzIdParser(mzid_file, str(tmpdir), peaklist, writer, logger)
#     id_parser.parse()
#     return id_parser
