import psycopg2
import json


class DBException(Exception):
    pass


def connect(dbname):
    import credentials
    try:
        con = psycopg2.connect(host=credentials.hostname, user=credentials.username, password=credentials.password,
                               dbname=credentials.database)
    except psycopg2.Error as e:
        raise DBException(e.message)

    return con


def create_tables(cur, con):
    # don't create tables here
    # use file postgreSQL_schema.sql to init db
    #
    # you will need to search and replace 'username' in the sql file,
    # replacing it with the role name you use to access the database
    return True


def new_upload(inj_list, cur, con):
    try:
        cur.execute("""
    INSERT INTO upload (
        user_id,
        identification_file_name,
        origin,       
        upload_time
    )
    VALUES (%s, %s, %s, CURRENT_TIMESTAMP) RETURNING id AS upload_id""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    rows = cur.fetchall()
    return rows[0][0]


def get_random_id(upload_id, cur, con):
    try:
        cur.execute("SELECT random_id FROM upload WHERE id = " + str(upload_id) + ";")
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    rows = cur.fetchall()
    return rows[0][0]


def write_mzid_info(peak_list_file_names,
                    spectra_formats,
                    analysis_software,
                    provider,
                    audits,
                    samples,
                    analyses,
                    protocol,
                    bib,
                    upload_id, cur, con):
    try:
        cur.execute("""UPDATE upload SET 
                        peak_list_file_names = (%s),
                        spectra_formats = (%s),
                        analysis_software = (%s),
                        provider = (%s),
                        audits = (%s),
                        samples = (%s),
                        analyses = (%s),
                        protocol = (%s),
                        bib = (%s) 
                        WHERE id = (%s);""",
                    (peak_list_file_names,
                         spectra_formats,
                         analysis_software,
                         provider,
                         audits,
                         samples,
                         analyses,
                         protocol,
                         bib,
                         upload_id))
        con.commit()
    except psycopg2.Error as e:
        raise DBException(e.message)
    return True


def write_other_info(upload_id, crosslinks, ident_count, ident_file_size, upload_warnings, cur, con):
    try:
        cur.execute("""UPDATE upload SET contains_crosslinks = (%s), ident_count = (%s), ident_file_size = (%s)
                , upload_warnings = (%s)
                 WHERE id = (%s);""",
                    (crosslinks, ident_count, ident_file_size, json.dumps(upload_warnings), upload_id))

        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    return True


def write_error(upload_id, error_type, error, cur, con):
    try:
        cur.execute("""UPDATE upload SET error_type = %s
                    , upload_error = %s
                    WHERE id = %s;""", (error_type, error, upload_id))
        con.commit()

        cur.execute("DELETE FROM db_sequences WHERE upload_id = " + str(upload_id) + ";")
        con.commit()

        cur.execute("DELETE FROM peptides WHERE upload_id = '" + str(upload_id) + "';")
        con.commit()

        cur.execute("DELETE FROM peptide_evidences WHERE upload_id = " + str(upload_id) + ";")
        con.commit()

        cur.execute("DELETE FROM modifications WHERE upload_id = " + str(upload_id) + ";")
        con.commit()

        cur.execute("DELETE FROM spectra WHERE upload_id = " + str(upload_id) + ";")
        con.commit()

        cur.execute("DELETE FROM spectrum_identifications WHERE upload_id = " + str(upload_id) + ";")
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)
    return True


# def write_protocols(inj_list, cur, con):
#     return True

def write_db_sequences(inj_list, cur, con):
    try:
        cur.executemany("""
        INSERT INTO db_sequences (
            id,
            accession,
            protein_name,
            description,
            sequence,
            upload_id
        )
        VALUES (%s, %s, %s, %s, %s, %s) """, inj_list)
        #     con.commit()
        #
    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_meta_data(*args):
    pass


def write_peptides(inj_list, cur, con):
    try:
        cur.executemany("""
        INSERT INTO peptides (
            id,
            seq_mods,
            link_site,
            crosslinker_modmass,
            upload_id,
            crosslinker_pair_id
        )
        VALUES (%s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_modifications(inj_list, cur, con):
    try:
        cur.executemany("""
          INSERT INTO modifications (
            id,
            upload_id,
            mod_name,
            mass,
            residues,
            accession
          )
          VALUES (%s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()
    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_peptide_evidences(inj_list, cur, con):
    try:
        cur.executemany("""
        INSERT INTO peptide_evidences (
            peptide_ref,
            dbsequence_ref,
            protein_accession,
            pep_start,
            is_decoy,
            upload_id
        )
        VALUES (%s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_spectra(inj_list, cur, con):
    try:
        cur.executemany("""
        INSERT INTO spectra (
        id, 
        mz,
        intensity, 
        peak_list_file_name, 
        scan_id, 
        frag_tol, 
        upload_id, 
        spectrum_ref,
        precursor_mz,
        precursor_charge
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)

    return True


def write_spectrum_identifications(inj_list, cur, con):
    try:
        cur.executemany("""
          INSERT INTO spectrum_identifications (
              id,
              upload_id,
              spectrum_id,
              pep1_id,
              pep2_id,
              charge_state,
              rank,
              pass_threshold,
              ions,
              scores,
              exp_mz,
              calc_mz,
              meta1,
              meta2,
              meta3
          ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s , %s, %s, %s, %s, %s, %s, %s)""", inj_list)
        con.commit()

    except psycopg2.Error as e:
        raise DBException(e.message)

    return True
