"""
ToDo: add test that writes multiple mzids into the db and check that results are written in
    properly.
"""
import numpy as np
from numpy.testing import assert_array_equal
from parser.api_writer import Table
import os
import logging
from sqlalchemy import text
from pyteomics import mgf
from .db_pytest_fixtures import *
from .parse_mzid import parse_mzid_into_postgresql, parse_mzid_into_sqlite_xispec
import struct

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(name)s %(message)s')
logger = logging.getLogger(__name__)


def fixture_path(file):
    current_dir = os.path.dirname(__file__)
    return os.path.join(current_dir, "fixtures", file)


def compare_db_sequence(results):
    assert len(results) == 12
    assert results[0].id == "dbseq_P0C0V0_target"  # id from mzid
    assert results[0].accession == "P0C0V0"  # accession from mzid
    assert results[0].name == (  # name from mzid
        "DEGP_ECOLI Periplasmic serine endoprotease DegP OS=Escherichia coli (strain K12) "
        "OX=83333 GN=degP PE=1 SV=1")
    assert results[0].description == (  # protein description cvParam
        "DEGP_ECOLI Periplasmic serine endoprotease DegP OS=Escherichia coli (strain K12) "
        "OX=83333 GN=degP PE=1 SV=1"
    )
    assert results[0].sequence == (  # <Seq> value from mzid
        "MKKTTLALSALALSLGLALSPLSATAAETSSATTAQQMPSLAPMLEKVMPSVVSINVEGSTTVNTPRMPRNFQQFFGDDSPFCQEG"
        "SPFQSSPFCQGGQGGNGGGQQQKFMALGSGVIIDADKGYVVTNNHVVDNATVIKVQLSDGRKFDAKMVGKDPRSDIALIQIQNPKN"
        "LTAIKMADSDALRVGDYTVAIGNPFGLGETVTSGIVSALGRSGLNAENYENFIQTDAAINRGNSGGALVNLNGELIGINTAILAPD"
        "GGNIGIGFAIPSNMVKNLTSQMVEYGQVKRGELGIMGTELNSELAKAMKVDAQRGAFVSQVLPNSSAAKAGIKAGDVITSLNGKPI"
        "SSFAALRAQVGTMPVGSKLTLGLLRDGKQVNVNLELQQSSQNQVDSSSIFNGIEGAEMSNKGKDQGVVVNNVKTGTPAAQIGLKKG"
        "DVIIGANQQAVKNIAELRKVLDSKPSVLALNIQRGDSTIYLLMQ")
    # ToDo: check more rows?


def compare_peptide_evidence(results):
    assert len(results) == 38
    # peptide_ref from <PeptideEvidence>
    assert results[0].peptide_ref == '29_KVLDSKPSVLALNIQR_30_KFDAKMVGK_1_5_p1'
    # dbsequence_ref from <PeptideEvidence>
    assert results[0].dbsequence_ref == 'dbseq_P0C0V0_target'
    assert results[0].pep_start == 148  # start from <PeptideEvidence>
    assert not results[0].is_decoy  # is_decoy from <PeptideEvidence>
    # ToDo: check more rows?


def compare_modified_peptide(results):
    assert len(results) == 38

    # id from <Peptide> id
    assert results[0].id == '29_KVLDSKPSVLALNIQR_30_KFDAKMVGK_1_5_p1'
    assert results[0].base_sequence == 'KFDAKMVGK'  # value of <PeptideSequence>
    assert results[0].mod_accessions == []
    assert results[0].mod_positions == []
    # location of <Modification> with cross-link acceptor/receiver cvParam
    assert results[0].link_site1 == 5
    # monoisotopicMassDelta of <Modification> with cross-link acceptor/receiver cvParam
    assert results[0].crosslinker_modmass == 0
    # value of cross-link acceptor/receiver cvParam
    assert results[0].crosslinker_pair_id == '1.0'

    # id from <Peptide> id
    assert results[1].id == '29_KVLDSKPSVLALNIQR_30_KFDAKMVGK_1_5_p0'
    assert results[1].base_sequence == 'KVLDSKPSVLALNIQR'  # value of <PeptideSequence>
    assert results[1].mod_accessions == []
    assert results[1].mod_positions == []
    # location of <Modification> with cross-link acceptor/receiver cvParam
    assert results[1].link_site1 == 1
    # monoisotopicMassDelta of <Modification> with cross-link acceptor/receiver cvParam
    assert results[1].crosslinker_modmass == pytest.approx(158.0037644600003, abs=1e-12)
    # value of cross-link acceptor/receiver cvParam
    assert results[1].crosslinker_pair_id == '1.0'

    # id from <Peptide> id
    assert results[2].id == '19_LLAEHNLDmetASAIKGTGVGGR_20_HLAKAPAK_13_4_p1'
    assert results[2].base_sequence == 'HLAKAPAK'  # value of <PeptideSequence>
    assert results[2].mod_accessions == []
    assert results[2].mod_positions == []
    # location of <Modification> with cross-link acceptor/receiver cvParam
    assert results[2].link_site1 == 4
    # monoisotopicMassDelta of <Modification> with cross-link acceptor/receiver cvParam
    assert results[2].crosslinker_modmass == 0
    # value of cross-link acceptor/receiver cvParam
    assert results[2].crosslinker_pair_id == '2.0'

    # id from <Peptide> id
    assert results[3].id == '19_LLAEHNLDmetASAIKGTGVGGR_20_HLAKAPAK_13_4_p0'
    assert results[3].base_sequence == 'LLAEHNLDASAIKGTGVGGR'  # value of <PeptideSequence>
    assert results[3].mod_accessions == [{'UNIMOD:34': 'Methyl'}]
    assert results[3].mod_positions == [8]
    # location of <Modification> with cross-link acceptor/receiver cvParam
    assert results[3].link_site1 == 13
    # monoisotopicMassDelta of <Modification> with cross-link acceptor/receiver cvParam
    assert results[3].crosslinker_modmass == pytest.approx(158.0037644600003, abs=1e-12)
    # value of cross-link acceptor/receiver cvParam
    assert results[3].crosslinker_pair_id == '2.0'
    # ToDo: check more rows?


def compare_modification(results):
    assert len(results) == 14

    assert results[0].id == 0  # id from incrementing count
    assert results[0].mod_name == '(158.00)'  # name from <SearchModification> cvParam / mod mass in brackets if unknown
    assert results[0].mass == 158.00377  # massDelta from <SearchModification>
    assert results[0].residues == 'STYK'  # residues from <SearchModification>
    assert results[0].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[0].fixed_mod  # fixedMod from <SearchModification>
    assert results[0].accession == 'MS:1001460'  # accession from <SearchModification> cvParam
    assert results[0].crosslinker_id == '0.0'  # value from cl donor / acceptor cv term (is a string)

    assert results[1].id == 1  # id from incrementing count
    assert results[1].mod_name == 'crosslink acceptor'  # name from <SearchModification> cvParam
    assert results[1].mass == 0  # massDelta from <SearchModification>
    assert results[1].residues == 'STYK'  # residues from <SearchModification>
    assert results[1].specificity_rules == []  # parsed from child <SpecificityRules>
    assert results[1].fixed_mod  # fixedMod from <SearchModification> (is mistake in xml file)
    assert results[1].accession == 'MS:1002510'  # accession from <SearchModification> cvParam
    assert results[1].crosslinker_id == '0.0'  # value from cl donor  / acceptor cv term (is a string)

    assert results[2].id == 2  # id from incrementing count
    assert results[2].mod_name == '(158.00)'  # name from <SearchModification> cvParam / mod mass in brackets if unknown
    assert results[2].mass == 158.00377  # massDelta from <SearchModification>
    assert results[2].residues == '.'  # residues from <SearchModification>
    assert results[2].specificity_rules == ["MS:1002057"]  # parsed from child <SpecificityRules>
    assert results[2].fixed_mod  # fixedMod from <SearchModification>
    assert results[2].accession == 'MS:1001460'  # accession from <SearchModification> cvParam
    assert results[2].crosslinker_id == '0.0'  # value from cl donor  / acceptor cv term (is a string)

    assert results[3].id == 3  # id from incrementing count
    assert results[3].mod_name == 'crosslink acceptor'  # name from <SearchModification> cvParam
    assert results[3].mass == 158.00377  # massDelta from <SearchModification> (mistake in xml?)
    assert results[3].residues == '.'  # residues from <SearchModification>
    assert results[3].specificity_rules == ["MS:1002057"]  # parsed from child <SpecificityRules>
    assert results[3].fixed_mod  # fixedMod from <SearchModification>
    assert results[3].accession == 'MS:1002510'  # accession from <SearchModification> cvParam
    assert results[3].crosslinker_id == '0.0'  # value from cl donor  / acceptor cv term (is a string)

    assert results[4].id == 4  # id from incrementing count
    assert results[4].mod_name == '(0.00)'  # name from <SearchModification> cvParam / mod mass in brackets if unknown
    assert results[4].mass == 0  # massDelta from <SearchModification>
    assert results[4].residues == '.'  # residues from <SearchModification>
    assert results[4].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[4].fixed_mod  # fixedMod from <SearchModification>
    assert results[4].accession == 'MS:1001460'  # accession from <SearchModification> cvParam
    assert results[4].crosslinker_id == '1.0'  # value from cl donor  / acceptor cv term (is a string)

    assert results[5].id == 5  # id from incrementing count
    assert results[5].mod_name == 'crosslink acceptor'  # name from <SearchModification> cvParam
    assert results[5].mass == 0  # massDelta from <SearchModification>
    assert results[5].residues == '.'  # residues from <SearchModification>
    assert results[5].specificity_rules == []  # parsed from child <SpecificityRules>
    assert results[5].fixed_mod  # fixedMod from <SearchModification>
    assert results[5].accession == 'MS:1002510'  # accession from <SearchModification> cvParam
    assert results[5].crosslinker_id == '1.0'  # value from cl donor  / acceptor cv term (is a string)

    assert results[6].id == 6  # id from incrementing count
    assert results[6].mod_name == 'Oxidation'  # name from <SearchModification> cvParam
    assert results[6].mass == 15.99491  # massDelta from <SearchModification>
    assert results[6].residues == 'M'  # residues from <SearchModification>
    assert results[6].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[6].fixed_mod  # fixedMod from <SearchModification>
    assert results[6].accession == 'UNIMOD:35'  # accession from <SearchModification> cvParam
    assert results[6].crosslinker_id is None

    assert results[7].id == 7  # id from incrementing count
    assert results[7].mod_name == '(175.03)'  # unknown modification -> name from mass
    assert results[7].mass == 175.03032  # massDelta from <SearchModification>
    assert results[7].residues == 'K'  # residues from <SearchModification>
    assert results[7].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[7].fixed_mod  # fixedMod from <SearchModification>
    assert results[7].accession == 'MS:1001460'  # accession from <SearchModification> cvParam
    assert results[7].crosslinker_id is None

    assert results[8].id == 8  # id from incrementing count
    assert results[8].mod_name == '(176.01)'  # unknown modification -> name from mass
    assert results[8].mass == 176.0143295  # massDelta from <SearchModification>
    assert results[8].residues == 'K'  # residues from <SearchModification>
    assert results[8].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[8].fixed_mod  # fixedMod from <SearchModification>
    assert results[8].accession == 'MS:1001460'  # accession from <SearchModification> cvParam
    assert results[8].crosslinker_id is None

    assert results[9].id == 9  # id from incrementing count
    assert results[9].mod_name == '(175.03)'  # unknown modification -> name from mass
    assert results[9].mass == 175.03032  # massDelta from <SearchModification>
    assert results[9].residues == '.'  # residues from <SearchModification>
    assert results[9].specificity_rules == ['MS:1002057']  # parsed from child <SpecificityRules>
    assert not results[9].fixed_mod  # fixedMod from <SearchModification>
    assert results[9].accession == 'MS:1001460'  # accession from <SearchModification> cvParam
    assert results[9].crosslinker_id is None

    assert results[10].id == 10  # id from incrementing count
    assert results[10].mod_name == '(176.01)'  # unknown modification -> name from mass
    assert results[10].mass == 176.0143295  # massDelta from <SearchModification>
    assert results[10].residues == '.'  # residues from <SearchModification>
    assert results[10].specificity_rules == ['MS:1002057']  # parsed from child <SpecificityRules>
    assert not results[10].fixed_mod  # fixedMod from <SearchModification>
    assert results[10].accession == 'MS:1001460'  # accession from <SearchModification> cvParam
    assert results[10].crosslinker_id is None

    assert results[11].id == 11  # id from incrementing count
    assert results[11].mod_name == 'Deamidated'  # name from <SearchModification> cvParam
    assert results[11].mass == 0.984016  # massDelta from <SearchModification>
    assert results[11].residues == 'NQ'  # residues from <SearchModification>
    assert results[11].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[11].fixed_mod  # fixedMod from <SearchModification>
    assert results[11].accession == 'UNIMOD:7'  # accession from <SearchModification> cvParam
    assert results[11].crosslinker_id is None

    assert results[12].id == 12  # id from incrementing count
    assert results[12].mod_name == 'Methyl'  # name from <SearchModification> cvParam
    assert results[12].mass == 14.01565  # massDelta from <SearchModification>
    assert results[12].residues == 'DE'  # residues from <SearchModification>
    assert results[12].specificity_rules == []  # parsed from child <SpecificityRules>
    assert not results[12].fixed_mod  # fixedMod from <SearchModification>
    assert results[12].accession == 'UNIMOD:34'  # accession from <SearchModification> cvParam
    assert results[12].crosslinker_id is None

    assert results[13].id == 13  # id from incrementing count
    assert results[13].mod_name == 'Carbamidomethyl'  # name from <SearchModification> cvParam
    assert results[13].mass == 57.021465  # massDelta from <SearchModification>
    assert results[13].residues == 'C'  # residues from <SearchModification>
    assert results[13].specificity_rules == []  # parsed from child <SpecificityRules>
    assert results[13].fixed_mod  # fixedMod from <SearchModification>
    assert results[13].accession == 'UNIMOD:4'  # accession from <SearchModification> cvParam
    assert results[13].crosslinker_id is None


def compare_enzyme(results):
    assert len(results) == 1
    assert results[0].id == "Trypsin_0"  # id from Enzyme element
    assert results[0].protocol_id == "SearchProtocol_1_0"
    assert results[0].c_term_gain == "OH"
    assert results[0].min_distance is None
    assert results[0].missed_cleavages == 2
    assert results[0].n_term_gain == "H"
    assert results[0].name == "Trypsin"
    assert results[0].semi_specific is False
    assert results[0].site_regexp == '(?<=[KR])(?\\!P)'
    assert results[0].accession == "MS:1001251"


def compare_spectrum_identification_protocol(results):
    assert len(results) == 1
    # parsed from <FragmentTolerance> in <SpectrumIdentificationProtocol>
    assert results[0].id == 'SearchProtocol_1_0'  # id from <SpectrumIdentificationProtocol>
    assert results[0].frag_tol == 5.0
    assert results[0].frag_tol_unit == 'ppm'
    # cvParams from <AdditionalSearchParams> 'ion series considered in search' (MS:1002473)

    assert results[0].additional_search_params == {'MS:1001211': 'parent mass type mono',
                                                   'MS:1002494': 'cross-linking search',
                                                   'MS:1001256': 'fragment mass type mono',
                                                   'MS:1001118': 'param: b ion',
                                                   'MS:1001262': 'param: y ion'}

    assert results[0].analysis_software['id'] == "xiFDR_id"


def compare_analysis_collection_mgf(results):
    assert len(results) == 2
    assert results[
               0].spectrum_identification_list_ref == 'SII_LIST_1_1_0_recal_B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mgf'
    assert results[0].spectrum_identification_protocol_ref == 'SearchProtocol_1_0'
    assert results[0].spectra_data_ref == 'SD_0_recal_B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mgf'
    assert results[
               1].spectrum_identification_list_ref == 'SII_LIST_1_1_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mgf'
    assert results[1].spectrum_identification_protocol_ref == 'SearchProtocol_1_0'
    assert results[1].spectra_data_ref == 'SD_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mgf'


def compare_analysis_collection_mzml(results):
    assert len(results) == 2
    assert results[
               0].spectrum_identification_list_ref == 'SII_LIST_1_1_0_recal_B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mzML'
    assert results[0].spectrum_identification_protocol_ref == 'SearchProtocol_1_0'
    assert results[0].spectra_data_ref == 'SD_0_recal_B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mzML'
    assert results[
               1].spectrum_identification_list_ref == 'SII_LIST_1_1_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mzML'
    assert results[1].spectrum_identification_protocol_ref == 'SearchProtocol_1_0'
    assert results[1].spectra_data_ref == 'SD_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mzML'


def compare_spectrum_mgf(conn, peak_list_folder):
    peaklists = [
        'recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mgf',
        'recal_B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mgf'
    ]
    for pl in peaklists:
        rs = conn.execute(text(f"SELECT * FROM Spectrum WHERE peak_list_file_name = '{pl}'"))
        results = rs.fetchall()
        assert len(results) == 11
        reader = mgf.read(os.path.join(peak_list_folder, pl))
        for r in results:
            # For MGF the index is encoded as e.g. index=3
            reader_idx = int(r.id.replace('index=', ''))
            spectrum = reader[reader_idx]
            # spectrumID from <SpectrumIdentificationResult>
            assert r.id == f'index={reader_idx}'  # a bit circular here
            # spectraData_ref from <SpectrumIdentificationResult>
            assert r.spectra_data_ref == f'SD_0_{pl}'
            assert r.peak_list_file_name == pl
            assert r.precursor_mz == spectrum['params']['pepmass'][0]
            assert r.precursor_charge == spectrum['params']['charge'][0]
            # check that mz and intensity values are as expected
            # 1. unpacking the blob
            assert_array_equal(struct.unpack('%sd' % (len(r.mz) // 8), r.mz), spectrum['m/z array'])
            assert_array_equal(struct.unpack('%sd' % (len(r.intensity) // 8), r.intensity), spectrum['intensity array'])
            # 2. using np.frombuffer (this works because np assumes double precision as default)
            assert_array_equal(np.frombuffer(r.mz), spectrum['m/z array'])
            assert_array_equal(np.frombuffer(r.intensity), spectrum['intensity array'])


def test_psql_db_cleared_each_test(use_database, engine):
    """Check that the database is empty."""
    with engine.connect() as conn:
        rs = conn.execute(text("SELECT * FROM Upload"))
        assert 0 == rs.rowcount
    engine.dispose()


def test_psql_mgf_mzid_parser(tmpdir, use_database, engine):
    # file paths
    fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'mzid_parser')
    mzid = os.path.join(fixtures_dir, 'mgf_ecoli_dsso.mzid')
    peak_list_folder = os.path.join(fixtures_dir, 'peaklist')

    id_parser = parse_mzid_into_postgresql(mzid, peak_list_folder, tmpdir, logger,
                                           use_database, engine)

    with engine.connect() as conn:
        # DBSequence
        stmt = Table("dbsequence", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_db_sequence(rs.fetchall())

        # SearchModification - parsed from <SearchModification>s
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
        stmt = Table("peptideevidence", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_peptide_evidence(rs.fetchall())

        # ModifiedPeptide
        stmt = Table("modifiedpeptide", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_modified_peptide(rs.fetchall())

        # Spectrum
        compare_spectrum_mgf(conn, peak_list_folder)

        # SpectrumIdentification
        stmt = Table("spectrumidentification", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        assert 22 == rs.rowcount
        results = rs.fetchall()
        assert results[0].id == 'SII_3_1'  # id from first <SpectrumIdentificationItem>
        assert results[0].spectrum_id == 'index=3'  # spectrumID from <SpectrumIdentificationResult>
        # spectraData_ref from <SpectrumIdentificationResult>
        assert results[0].spectra_data_ref == \
               'SD_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mgf'
        # peptide_ref from <SpectrumIdentificationItem>
        assert results[0].pep1_id == \
               '6_VAEmetETPHLIHKVALDPLTGPMPYQGR_11_MGHAGAIIAGGKGTADEK_11_12_p1'
        # peptide_ref from matched <SpectrumIdentificationItem> by crosslink_identification_id
        assert results[0].pep2_id == \
               '6_VAEmetETPHLIHKVALDPLTGPMPYQGR_11_MGHAGAIIAGGKGTADEK_11_12_p0'
        assert results[0].charge_state == 5  # chargeState from <SpectrumIdentificationItem>
        assert results[0].pass_threshold  # passThreshold from <SpectrumIdentificationItem>
        assert results[0].rank == 1  # rank from <SpectrumIdentificationItem>
        # scores parsed from score related cvParams in <SpectrumIdentificationItem>
        assert results[0].scores == {'xi:score': 33.814201}
        # experimentalMassToCharge from <SpectrumIdentificationItem>
        assert results[0].exp_mz == 945.677359
        # calculatedMassToCharge from <SpectrumIdentificationItem>
        assert results[0].calc_mz == pytest.approx(945.6784858667701, abs=1e-12)

        # SpectrumIdentificationProtocol
        stmt = Table("spectrumidentificationprotocol", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_spectrum_identification_protocol(rs.fetchall())

        # AnalysisCollection
        stmt = Table("analysiscollection", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_analysis_collection_mgf(rs.fetchall())

        # Upload
        stmt = Table("upload", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        assert 1 == rs.rowcount
        results = rs.fetchall()

        assert results[0].identification_file_name == 'mgf_ecoli_dsso.mzid'
        assert results[0].provider == {"id": "PROVIDER", "ContactRole": [
            {"contact_ref": "PERSON_DOC_OWNER", "Role": "researcher"}]}
        assert results[0].audit_collection == {'Organization': {'contact name': 'TU Berlin',
                                                                'id': 'ORG_DOC_OWNER',
                                                                'name': 'TU Berlin'},
                                               'Person': {'Affiliation': [{'organization_ref': 'ORG_DOC_OWNER'}],
                                                          'contact address': 'TIB 4/4-3 Gebäude 17, Aufgang 1, Raum 476 '
                                                                             'Gustav-Meyer-Allee 25 13355 Berlin',
                                                          'contact email': 'lars.kolbowski@tu-berlin.de',
                                                          'firstName': 'Lars',
                                                          'id': 'PERSON_DOC_OWNER',
                                                          'lastName': 'Kolbowski'}}
        assert results[0].analysis_sample_collection == {}
        assert results[0].bib == []
        assert results[0].spectra_formats == [{'FileFormat': 'Mascot MGF format',
                                               'SpectrumIDFormat': 'multiple peak list nativeID format',
                                               'id': 'SD_0_recal_B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mgf',
                                               'location': 'recal_B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mgf'},
                                              {'FileFormat': 'Mascot MGF format',
                                               'SpectrumIDFormat': 'multiple peak list nativeID format',
                                               'id': 'SD_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mgf',
                                               'location': 'recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mgf'}]
        assert results[0].contains_crosslinks
        assert results[0].upload_warnings == []

    engine.dispose()


def test_psql_mzml_mzid_parser(tmpdir, use_database, engine):
    # file paths
    fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'mzid_parser')
    mzid = os.path.join(fixtures_dir, 'mzml_ecoli_dsso.mzid')
    peak_list_folder = os.path.join(fixtures_dir, 'peaklist')

    id_parser = parse_mzid_into_postgresql(mzid, peak_list_folder, tmpdir, logger,
                                           use_database, engine)

    with engine.connect() as conn:
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
        stmt = Table("peptideevidence", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_peptide_evidence(rs.fetchall())

        # ModifiedPeptide
        stmt = Table("modifiedpeptide", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_modified_peptide(rs.fetchall())

        # Spectrum
        # ToDo: create and use compare_spectrum_mzml()
        stmt = Table("spectrum", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        assert 22 == rs.rowcount
        results = rs.fetchall()
        # spectrumID from <SpectrumIdentificationResult>
        assert results[0].id == 'controllerType=0 controllerNumber=1 scan=14905'
        # spectraData_ref from <SpectrumIdentificationResult>
        assert results[0].spectra_data_ref == (
            'SD_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mzML')
        # assert results[0].scan_id == '3'  # ToDo: keep this?
        assert results[0].peak_list_file_name == (  # ToDo
            'B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mzML')
        # Precursor info from mzML
        # Spectrum->precursorList->precursor[0]->selectedIonList->selectedIon[0]
        assert results[0].precursor_mz == 945.677381209342  # selected ion m/z
        assert results[0].precursor_charge == 5  # charge state
        # assert results[0].mz == []  # ToDo
        # assert results[0].intensity == []  # ToDo

        # SpectrumIdentification
        stmt = Table("spectrumidentification", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        assert 22 == rs.rowcount
        results = rs.fetchall()
        assert results[0].id == 'SII_3_1'  # id from first <SpectrumIdentificationItem>
        # spectrumID from <SpectrumIdentificationResult>
        assert results[0].spectrum_id == 'controllerType=0 controllerNumber=1 scan=14905'
        # spectraData_ref from <SpectrumIdentificationResult>
        assert results[0].spectra_data_ref == \
               'SD_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mzML'
        # peptide_ref from <SpectrumIdentificationItem>
        assert results[0].pep1_id == \
               '6_VAEmetETPHLIHKVALDPLTGPMPYQGR_11_MGHAGAIIAGGKGTADEK_11_12_p1'
        # peptide_ref from matched <SpectrumIdentificationItem> by crosslink_identification_id
        assert results[0].pep2_id == \
               '6_VAEmetETPHLIHKVALDPLTGPMPYQGR_11_MGHAGAIIAGGKGTADEK_11_12_p0'
        assert results[0].charge_state == 5  # chargeState from <SpectrumIdentificationItem>
        assert results[0].pass_threshold  # passThreshold from <SpectrumIdentificationItem>
        assert results[0].rank == 1  # rank from <SpectrumIdentificationItem>
        # scores parsed from score related cvParams in <SpectrumIdentificationItem>
        assert results[0].scores == {"xi:score": 33.814201}
        # experimentalMassToCharge from <SpectrumIdentificationItem>
        assert results[0].exp_mz == 945.677359
        # calculatedMassToCharge from <SpectrumIdentificationItem>
        assert results[0].calc_mz == pytest.approx(945.6784858667701, abs=1e-12)

        # SpectrumIdentificationProtocol
        stmt = Table("spectrumidentificationprotocol", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_spectrum_identification_protocol(rs.fetchall())

        # AnalysisCollection
        stmt = Table("analysiscollection", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_analysis_collection_mzml(rs.fetchall())

        # Upload
        stmt = Table("upload", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        assert 1 == rs.rowcount
        results = rs.fetchall()

        assert results[0].identification_file_name == 'mzml_ecoli_dsso.mzid'
        assert results[0].provider == {'ContactRole': [{'Role': 'researcher', 'contact_ref': 'PERSON_DOC_OWNER'}],
                                       'id': 'PROVIDER'}
        assert results[0].audit_collection == {'Organization': {'contact name': 'TU Berlin',
                                                                'id': 'ORG_DOC_OWNER',
                                                                'name': 'TU Berlin'},
                                               'Person': {'Affiliation': [{'organization_ref': 'ORG_DOC_OWNER'}],
                                                          'contact address': 'TIB 4/4-3 Gebäude 17, Aufgang 1, Raum 476 '
                                                                             'Gustav-Meyer-Allee 25 13355 Berlin',
                                                          'contact email': 'lars.kolbowski@tu-berlin.de',
                                                          'firstName': 'Lars',
                                                          'id': 'PERSON_DOC_OWNER',
                                                          'lastName': 'Kolbowski'}}
        assert results[0].analysis_sample_collection == {}
        assert results[0].bib == []
        assert results[0].spectra_formats == [{'FileFormat': 'mzML format',
                                               'SpectrumIDFormat': 'mzML unique identifier',
                                               'id': 'SD_0_recal_B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mzML',
                                               'location': 'B190717_20_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX01_rep2.mzML'},
                                              {'FileFormat': 'mzML format',
                                               'SpectrumIDFormat': 'mzML unique identifier',
                                               'id': 'SD_0_recal_B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mzML',
                                               'location': 'B190717_13_HF_LS_IN_130_ECLP_DSSO_01_SCX23_hSAX05_rep2.mzML'}]
        assert results[0].contains_crosslinks
        # assert results[0].upload_warnings == [
        #     'mzidentML file does not specify any fragment ions (child terms of MS_1002473) within '
        #     '<AdditionalSearchParams>. Falling back to b and y ions.']

    engine.dispose()


def test_sqlite_mgf_xispec_mzid_parser(tmpdir):
    # file paths
    fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'mzid_parser')
    mzid = os.path.join(fixtures_dir, 'mgf_ecoli_dsso.mzid')
    peak_list_folder = os.path.join(fixtures_dir, 'peaklist')
    test_database = os.path.join(str(tmpdir), 'test.db')

    conn_str = f'sqlite:///{test_database}'
    engine = create_engine(conn_str)
    id_parser = parse_mzid_into_sqlite_xispec(mzid, peak_list_folder, tmpdir, logger, engine)

    with engine.connect() as conn:
        # DBSequence - not written for xiSPEC
        stmt = Table("DBSequence", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        assert len(rs.fetchall()) == 0

        # Modification - parsed from <SearchModification>s
        stmt = Table("SearchModification", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_modification(rs.fetchall())

        # Enzyme - parsed from SpectrumIdentificationProtocols
        stmt = Table("Enzyme", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_enzyme(rs.fetchall())

        # PeptideEvidence
        stmt = Table("PeptideEvidence", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_peptide_evidence(rs.fetchall())

        # ModifiedPeptide
        stmt = Table("ModifiedPeptide", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_modified_peptide(rs.fetchall())

        # Spectrum
        compare_spectrum_mgf(conn, peak_list_folder)

        # SpectrumIdentification
        stmt = Table("SpectrumIdentification", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        # ToDo

        # SpectrumIdentificationProtocol
        stmt = Table("SpectrumIdentificationProtocol", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_spectrum_identification_protocol(rs.fetchall())

        # AnalysisCollection
        stmt = Table("AnalysisCollection", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_analysis_collection_mgf(rs.fetchall())

        # Upload - not written for xiSPEC
        stmt = Table("Upload", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        assert len(rs.fetchall()) == 0

    engine.dispose()


def test_sqlite_mzml_xispec_mzid_parser(tmpdir):
    # file paths
    fixtures_dir = os.path.join(os.path.dirname(__file__), 'fixtures', 'mzid_parser')
    mzid = os.path.join(fixtures_dir, 'mzml_ecoli_dsso.mzid')
    peak_list_folder = os.path.join(fixtures_dir, 'peaklist')
    test_database = os.path.join(str(tmpdir), 'test.db')

    conn_str = f'sqlite:///{test_database}'
    engine = create_engine(conn_str)

    id_parser = parse_mzid_into_sqlite_xispec(mzid, peak_list_folder, tmpdir, logger, engine)

    with engine.connect() as conn:
        # DBSequence - not written for xiSPEC
        stmt = Table("DBSequence", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        assert len(rs.fetchall()) == 0

        # Modification - parsed from <SearchModification>s
        stmt = Table("SearchModification", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_modification(rs.fetchall())

        # Enzyme - parsed from SpectrumIdentificationProtocols
        stmt = Table("Enzyme", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        compare_enzyme(rs.fetchall())

        # PeptideEvidence
        stmt = Table("PeptideEvidence", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_peptide_evidence(rs.fetchall())

        # ModifiedPeptide
        stmt = Table("ModifiedPeptide", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_modified_peptide(rs.fetchall())

        # Spectrum
        stmt = Table("Spectrum", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        # ToDo: create and use compare_spectrum_mzml()

        # SpectrumIdentification
        stmt = Table("SpectrumIdentification", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        # ToDo

        # SpectrumIdentificationProtocol
        stmt = Table("SpectrumIdentificationProtocol", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_spectrum_identification_protocol(rs.fetchall())

        # AnalysisCollection
        stmt = Table("AnalysisCollection", id_parser.writer.meta,
                     autoload_with=id_parser.writer.engine, quote=False).select()
        rs = conn.execute(stmt)
        compare_analysis_collection_mzml(rs.fetchall())

        # Upload - not written for xiSPEC
        stmt = Table("Upload", id_parser.writer.meta, autoload_with=id_parser.writer.engine,
                     quote=False).select()
        rs = conn.execute(stmt)
        assert len(rs.fetchall()) == 0

    engine.dispose()
