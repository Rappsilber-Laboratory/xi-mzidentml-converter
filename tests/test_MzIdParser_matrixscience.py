from parser.api_writer import Writer, Table
import os
import logging
from .db_pytest_fixtures import *
from .parse_mzid import parse_mzid_into_postgresql

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logger = logging.getLogger(__name__)


def fixture_path(file):
    current_dir = os.path.dirname(__file__)
    return os.path.join(current_dir, "fixtures", file)


def compare_spectrum_identification(results):
    assert len(results) == 249


def compare_db_sequence(results):
    assert len(results) == 47
    # ToDo:


def compare_peptide_evidence(results):
    pass
    # ToDo:


def compare_modified_peptide(results):
    assert len(results) == 286
    # ToDo:


def compare_modification(results):
    assert len(results) == 2

    assert results[0].id == 0  # id from incrementing count
    assert results[0].mod_name == 'Nethylmaleimide'  # name from <SearchModification> cvParam
    assert results[0].mass == 125.047679  # massDelta from <SearchModification>
    assert results[0].residues == 'C'  # residues from <SearchModification>
    assert results[0].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[0].fixed_mod  # fixedMod from <SearchModification>
    assert results[0].accession == 'UNIMOD:108'  # accession from <SearchModification> cvParam

    assert results[1].id == 1  # id from incrementing count
    assert results[1].mod_name == 'Oxidation'  # name from <SearchModification> cvParam
    assert results[1].mass == 15.994915  # massDelta from <SearchModification>
    assert results[1].residues == 'M'  # residues from <SearchModification>
    assert results[1].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[1].fixed_mod  # fixedMod from <SearchModification>
    assert results[1].accession == 'UNIMOD:35'  # accession from <SearchModification> cvParam


def compare_enzyme(results):
    assert len(results) == 1

    assert results[0].id == "ENZ_0"  # id from Enzyme element
    assert results[0].protocol_id == "SIP"
    assert results[0].c_term_gain == "OH"
    assert results[0].min_distance is None
    assert results[0].missed_cleavages == 2
    assert results[0].n_term_gain == "H"
    assert results[0].name == "Trypsin/P"
    assert results[0].semi_specific is False
    assert results[0].site_regexp == "(?<=[KR])"
    assert results[0].accession == "MS:1001313"


def test_psql_matrixscience_mzid_parser(tmpdir, db_info, use_database, engine):
    # file paths
    fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'mzid_parser')
    mzid = os.path.join(fixtures_dir, 'F002553.mzid')
    peak_list_folder = False

    id_parser = parse_mzid_into_postgresql(mzid, peak_list_folder, tmpdir, logger, use_database,
                                           engine)

    with engine.connect() as conn:
        # SpectrumIdentification
        stmt = Table("spectrumidentification", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_spectrum_identification(rs.fetchall())

        # DBSequence
        stmt = Table("dbsequence", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_db_sequence(rs.fetchall())

        # Modification - parsed from <SearchModification>s
        stmt = Table("searchmodification", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_modification(rs.fetchall())

        # Enzyme - parsed from SpectrumIdentificationProtocols
        stmt = Table("enzyme", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_enzyme(rs.fetchall())

        # PeptideEvidence
        stmt = Table("peptideevidence", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_peptide_evidence(rs.fetchall())

        # ModifiedPeptide
        stmt = Table("modifiedpeptide", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_modified_peptide(rs.fetchall())

        # Spectrum (peak_list_folder = False)
        stmt = Table("spectrum", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        assert len(rs.fetchall()) == 0

        # ToDo: remaining Tables

    engine.dispose()
