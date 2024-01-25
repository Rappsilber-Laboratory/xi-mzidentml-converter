import json
import logging
import logging.config
import struct

import psycopg2
from fastapi import APIRouter
from psycopg2 import sql
from psycopg2.extras import RealDictCursor

from app.routes.shared import get_db_connection, get_most_recent_upload_ids

xiview_data_router = APIRouter()

logger = logging.getLogger(__name__)


@xiview_data_router.get('/get_xiview_data', tags=["xiVIEW"])
async def get_xiview_data(project, file=None):
    """
    Get the data for the network visualisation.
    URLs have following structure:
    https: // www.ebi.ac.uk / pride / archive / xiview / network.html?project=PXD020453&file=Cullin_SDA_1pcFDR.mzid
    Users may provide only projects, meaning we need to have an aggregated  view.
    https: // www.ebi.ac.uk / pride / archive / xiview / network.html?project=PXD020453

    :return: json with the data
    """
    most_recent_upload_ids = await get_most_recent_upload_ids(project, file)
    try:
        data_object = await get_data_object(most_recent_upload_ids, project)
    except psycopg2.DatabaseError as e:
        logger.error(e)
        print(e)
        return {"error": "Database error"}, 500

    return data_object

@xiview_data_router.get('/get_peaklist', tags=["xiVIEW"])
async def get_peaklist(id, sd_ref, upload_id):
    conn = None
    data = {}
    error = None
    try:
        conn = await get_db_connection()
        cur = conn.cursor()
        query = "SELECT intensity, mz FROM spectrum WHERE id = %s AND spectra_data_ref = %s AND upload_id = %s"
        cur.execute(query, [id, sd_ref, upload_id])
        resultset = cur.fetchall()[0]
        data["intensity"] = struct.unpack('%sd' % (len(resultset[0]) // 8), resultset[0])
        data["mz"] = struct.unpack('%sd' % (len(resultset[1]) // 8), resultset[1])
        cur.close()
    except (Exception, psycopg2.DatabaseError) as e:
        # logger.error(error)
        error = e
    finally:
        if conn is not None:
            conn.close()
            # logger.debug('Database connection closed.')
        if error is not None:
            raise error
        return data

@xiview_data_router.get('/get_datasets', tags=["xiVIEW"])
async def get_datasets():
    conn = None
    ds_rows = []
    error = None
    try:
        conn = await get_db_connection()
        cur = conn.cursor()
        query = """SELECT DISTINCT project_id, identification_file_name FROM upload;"""
        # logger.debug(query)
        cur.execute(query)
        ds_rows = cur.fetchall()
        # logger.info("finished")
        cur.close()
    except (Exception, psycopg2.DatabaseError) as e:
        print(e)
        error = e
    finally:
        if conn is not None:
            conn.close()
            # logger.debug('Database connection closed.')
        if error is not None:
            raise error
    return ds_rows



async def get_data_object(ids, pxid):
    """ Connect to the PostgreSQL database server """
    conn = None
    data = {}
    error = None
    try:
        conn = await get_db_connection()
        cur = conn.cursor(cursor_factory=RealDictCursor)
        data["project"] = await get_pride_api_info(cur, pxid)
        data["meta"] = await get_results_metadata(cur, ids)
        data["matches"] = await get_matches(cur, ids)
        data["peptides"] = await get_peptides(cur, data["matches"])
        data["proteins"] = await get_proteins(cur, data["peptides"])
        logger.info("finished")
        cur.close()
    except (Exception, psycopg2.DatabaseError) as e:
        error = e
        logger.error(e)
        raise e
    finally:
        if conn is not None:
            conn.close()
            logger.debug('Database connection closed.')
        if error is not None:
            raise error
        return data

async def get_pride_api_info(cur, pxid):
    """ Get the PRIDE API info for the projects """
    query = """SELECT p.id AS id,
                p.id,
                p.title,
                p.description
                FROM projectdetails p
                WHERE p.project_id = (%s);"""
    cur.execute(query, [pxid])
    return cur.fetchall()

async def get_results_metadata(cur, ids):
    """ Get the metadata for the results """
    metadata = {}

    # get Upload(s) for each id
    query = """SELECT u.id AS id,
                u.project_id,
                u.identification_file_name,
                u.provider,
                u.audit_collection,
                u.analysis_sample_collection,
                u.bib,
                u.spectra_formats,
                u.contains_crosslinks,
                u.upload_warnings AS warnings
            FROM upload u
            WHERE u.id = ANY(%s);"""
    cur.execute(query, [ids])
    metadata["mzidentml_files"] = cur.fetchall()

    # get AnalysisCollection(s) for each id
    query = """SELECT ac.upload_id,
                ac.spectrum_identification_list_ref,
                ac.spectrum_identification_protocol_ref,
                ac.spectra_data_ref
            FROM analysiscollection ac
            WHERE ac.upload_id = ANY(%s);"""
    cur.execute(query, [ids])
    metadata["analysis_collections"] = cur.fetchall()

    # get SpectrumIdentificationProtocol(s) for each id
    query = """SELECT sip.id AS id,
                sip.upload_id,
                sip.frag_tol,
                sip.frag_tol_unit,
                sip.additional_search_params,
                sip.analysis_software,
                sip.threshold
            FROM spectrumidentificationprotocol sip
            WHERE sip.upload_id = ANY(%s);"""
    cur.execute(query, [ids])
    metadata["spectrum_identification_protocols"] = cur.fetchall()

    # enzymes
    query = """SELECT *
            FROM enzyme e
            WHERE e.upload_id = ANY(%s);"""
    cur.execute(query, [ids])
    metadata["enzymes"] = cur.fetchall()

    # search modifications
    try:
        query = """SELECT *
                FROM searchmodification sm
                WHERE sm.upload_id = ANY(%s);"""
        cur.execute(query, [ids])
        metadata["search_modifications"] = cur.fetchall()
    except Exception as e:
        print(e)

    return metadata


async def get_matches(cur, ids):
    query = """SELECT si.id AS id, si.pep1_id AS pi1, si.pep2_id AS pi2,
                si.scores AS sc,
                cast (si.upload_id as text) AS si,
                si.calc_mz AS c_mz,
                si.charge_state AS pc_c,
                si.exp_mz AS pc_mz,
                si.spectrum_id AS sp,
                si.spectra_data_ref AS sd_ref,
                si.pass_threshold AS pass,
                si.rank AS r,
                si.sil_id AS sil                
            FROM spectrumidentification si 
            INNER JOIN modifiedpeptide mp1 ON si.pep1_id = mp1.id AND si.upload_id = mp1.upload_id 
            INNER JOIN modifiedpeptide mp2 ON si.pep2_id = mp2.id AND si.upload_id = mp2.upload_id
            WHERE si.upload_id = ANY(%s) 
            AND si.pass_threshold = TRUE 
            AND mp1.link_site1 > 0
            AND mp2.link_site1 > 0;"""
    # bit weird above works when link_site1 is a text column
    cur.execute(query, [ids])
    return cur.fetchall()


async def get_peptides(cur, match_rows):
    search_peptide_ids = {}
    for match_row in match_rows:
        if match_row['si'] in search_peptide_ids:
            peptide_ids = search_peptide_ids[match_row['si']]
        else:
            peptide_ids = set()
            search_peptide_ids[match_row['si']] = peptide_ids
        peptide_ids.add(match_row['pi1'])
        if match_row['pi2'] is not None:
            peptide_ids.add(match_row['pi2'])

    subclauses = []
    for k, v in search_peptide_ids.items():
        pep_id_literals = []
        for pep_id in v:
            pep_id_literals.append(sql.Literal(pep_id))
        joined_pep_ids = sql.SQL(',').join(pep_id_literals)
        subclause = sql.SQL("(mp.upload_id = {} AND id IN ({}))").format(
            sql.Literal(k),
            joined_pep_ids
        )
        subclauses.append(subclause)
    peptide_clause = sql.SQL(" OR ").join(subclauses)

    query = sql.SQL("""SELECT mp.id, cast(mp.upload_id as text) AS u_id,
                mp.base_sequence AS base_seq,
                array_agg(pp.dbsequence_ref) AS prt,
                array_agg(pp.pep_start) AS pos,
                array_agg(pp.is_decoy) AS is_decoy,
                mp.link_site1 AS "linkSite",
                mp.mod_accessions as mod_accs,
                mp.mod_positions as mod_pos,
                mp.mod_monoiso_mass_deltas as mod_masses,
                mp.crosslinker_modmass as cl_modmass                     
                    FROM modifiedpeptide AS mp
                    JOIN peptideevidence AS pp
                    ON mp.id = pp.peptide_ref AND mp.upload_id = pp.upload_id
                WHERE {}
                GROUP BY mp.id, mp.upload_id, mp.base_sequence;""").format(
        peptide_clause
    )
    # logger.debug(query.as_string(cur))
    cur.execute(query)
    return cur.fetchall()


async def get_proteins(cur, peptide_rows):
    search_protein_ids = {}
    for peptide_row in peptide_rows:
        if peptide_row['u_id'] in search_protein_ids:
            protein_ids = search_protein_ids[peptide_row['u_id']]
        else:
            protein_ids = set()
            search_protein_ids[peptide_row['u_id']] = protein_ids
        for prot in peptide_row['prt']:
            protein_ids.add(prot)

    subclauses = []
    for k, v in search_protein_ids.items():
        literals = []
        for prot_id in v:
            literals.append(sql.Literal(prot_id))
        joined_literals = sql.SQL(",").join(literals)
        subclause = sql.SQL("(upload_id = {} AND id IN ({}))").format(
            sql.Literal(k),
            joined_literals
        )
        subclauses.append(subclause)

    protein_clause = sql.SQL(" OR ").join(subclauses)
    query = sql.SQL("""SELECT id, name, accession, sequence,
                     cast(upload_id as text) AS search_id, description FROM dbsequence WHERE ({});""").format(
        protein_clause
    )
    # logger.debug(query.as_string(cur))
    cur.execute(query)
    return cur.fetchall()
