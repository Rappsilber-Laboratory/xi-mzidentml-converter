import sqlite3
import json


class DBException(Exception):
    pass


def connect(dbname):
    try:
        con = sqlite3.connect(dbname)
    except sqlite3.Error as e:
        raise DBException(e.message)

    return con


def create_tables(cur, con):
    try:
        # cur.execute("DROP TABLE IF EXISTS uploads")
        # cur.execute(
        #     "CREATE TABLE uploads("
        #     "id INT PRIMARY KEY, "
        #     "user_id INT,"
        #     "filename TEXT, "
        #     "peak_list_file_names TEXT, "
        #     "analysis_software JSON,"
        #     "provider JSON,"
        #     "audits JSON,"
        #     "samples JSON,"
        #     "analyses JSON,"
        #     "protocol JSON,"
        #     "bib JSON,"
        #     "upload_time DATE, "
        #     "upload_loc TEXT,"          
        #     "default_pdb TEXT,"
        #     "contains_crosslink BOOLEAN,"
        #     "upload_errors JSON)"
        # )

        cur.execute("DROP TABLE IF EXISTS meta_data")
        cur.execute(
            "CREATE TABLE meta_data("
            "upload_id INT,"
            "sid_meta1_name TEXT,"
            "sid_meta2_name TEXT,"
            "sid_meta3_name TEXT,"
            "contains_crosslink BOOLEAN,"
            "custom_config TEXT)"
        )

        # ToDo: not used atm might be a good place to save ions here?
        cur.execute("DROP TABLE IF EXISTS protocols")
        cur.execute(
            "CREATE TABLE protocols("
            "id TEXT PRIMARY KEY, "
            "upload_id INT,"
            "protocol JSON)"
        )

        cur.execute("DROP TABLE IF EXISTS peptides")
        cur.execute(
            "CREATE TABLE peptides("
            # "id TEXT PRIMARY KEY, "
            "id TEXT, "
            "upload_id INT,"
            "seq_mods TEXT,"
            "link_site INT,"
            "crosslinker_modmass FLOAT,"    # ToDo: save cross-links to extra table?
            "crosslinker_pair_id INT)"
        )

        cur.execute("DROP TABLE IF EXISTS modifications")
        cur.execute(
            "CREATE TABLE modifications("
            "id INT PRIMARY KEY, "
            "upload_id INT,"
            "mod_name TEXT, " 
            "mass FLOAT, "
            "residues TEXT, "
            "accession TEXT)"
        )

        cur.execute("DROP TABLE IF EXISTS peptide_evidences")
        cur.execute(
            "CREATE TABLE peptide_evidences("
            "upload_id INT,"
            "peptide_ref TEXT, "
            "dbsequence_ref TEXT, "
            "protein_accession TEXT,"
            "pep_start INT, "
            "is_decoy INT)"
        )

        cur.execute("DROP TABLE IF EXISTS spectra")
        cur.execute(
            "CREATE TABLE spectra("
            "id INT, "
            "upload_id INT,"
            "peak_list TEXT, "
            "peak_list_file_name TEXT, "
            "scan_id INT, "
            "frag_tol TEXT,"
            "spectrum_ref TEXT,"
            "precursor_mz FLOAT,"
            "precursor_charge INT)"
        )

        cur.execute("DROP TABLE IF EXISTS spectrum_identifications")
        cur.execute(
            "CREATE TABLE spectrum_identifications("
            "id INT, "
            "upload_id INT,"
            "spectrum_id INT, "
            "pep1_id TEXT, "
            "pep2_id TEXT, "
            "charge_state INT, "
            "pass_threshold INT, "
            "rank INT,"
            "ions TEXT, "   # ToDo: find better place to store ions - might be protocols
            "scores JSON,"  # IS JSON data type valid or does it have to be TEXT
            "exp_mz FLOAT,"
            "calc_mz FLOAT,"
            "meta1 TEXT,"
            "meta2 TEXT,"           
            "meta3 TEXT)"
        )
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)
    return True


def new_upload(*args):
    return True


def get_random_id(*args):
    return 0


def write_db_sequences(*args):

    # try:
    #     cur.executemany("""
    # INSERT INTO db_sequences (
    #     'id',
    #     'accession',
    #     'name',
    #     'description',
    #     'sequence',
    #     'upload_id'
    # )
    # VALUES (?, ?, ?, ?, ?, ?)""", inj_list)
    #     con.commit()
    #
    # except sqlite3.Error as e:
    #     raise DBException(e.message)

    return True


def write_meta_data(values, cur, con):
    try:
        cur.execute("""
          INSERT INTO meta_data (
            'upload_id',
            'sid_meta1_name', 
            'sid_meta2_name', 
            'sid_meta3_name',
            'contains_crosslink'
          )
          VALUES (?, ?, ?, ?, ?)""",  values)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


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
        VALUES (?, ?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_modifications(inj_list, cur, con):
    try:
        cur.executemany("""
          INSERT INTO modifications (
            'id',
            'upload_id',
            'mod_name', 
            'mass', 
            'residues', 
            'accession'
          )
          VALUES (?, ?, ?, ?, ?, ?)""",  inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_peptide_evidences(inj_list, cur, con):
    try:
        cur.executemany("""
    INSERT INTO peptide_evidences (
        'peptide_ref',
        'dbsequence_ref',
        'protein_accession',
        'pep_start',
        'is_decoy',
        'upload_id'
    )
    VALUES (?, ?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_spectra(inj_list, cur, con):
    try:
        cur.executemany("""
          INSERT INTO spectra (
              'id', 
              'peak_list', 
              'peak_list_file_name', 
              'scan_id', 
              'frag_tol', 
              'upload_id', 
              'spectrum_ref',
              'precursor_mz',
              'precursor_charge'
          )
          VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


def write_spectrum_identifications(inj_list, cur, con):
    try:
        cur.executemany("""
          INSERT INTO spectrum_identifications (
              'id', 
              'upload_id', 
              'spectrum_id', 
              'pep1_id', 
              'pep2_id',
              'charge_state', 
              'rank', 
              'pass_threshold', 
              'ions', 
              'scores',
              'exp_mz',
              'calc_mz',
              'meta1',
              'meta2',
              'meta3'
          ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""", inj_list)
        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)

    return True


# con = connect('/home/lars/Xi/xiSPEC_ms_parser/dbs/saved/Tmuris_exosomes1.db')
# cur = con.cursor()


def fill_in_missing_scores(cur, con):
    try:
        cur.execute("""
          SELECT DISTINCT scoresJSON.key as scoreKey
          FROM spectrum_identifications, json_each(spectrum_identifications.scores) AS scoresJSON""")

        all_scores = cur.fetchall()
        all_scores = set([str(x[0]) for x in all_scores])

        inj_list = []

        cur.execute('SELECT id, scores FROM spectrum_identifications')
        res = cur.fetchall()

        for row in res:
            row_scores = json.loads(row[1])
            missing = all_scores - set(row_scores.keys())

            if len(missing) > 0:
                missing_dict = {key: -1 for key in missing}
                row_scores.update(missing_dict)
                inj_list.append([json.dumps(row_scores), row[0]])
                # cur.execute('UPDATE identifications SET allScores=? WHERE id = row[0]', json.dumps(row_scores))

        cur.executemany("""
            UPDATE spectrum_identifications
            SET `scores` = ?
            WHERE `id` = ?""", inj_list)

        con.commit()

    except sqlite3.Error as e:
        raise DBException(e.message)
    pass




