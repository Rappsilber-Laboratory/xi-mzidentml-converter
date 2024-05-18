# from parser.parser import SQLite
# import glob
#
#
# def update_database(con, db_name):
#
#     cur = con.cursor()
#
#     cur.execute(
#         "CREATE TABLE IF NOT EXISTS meta_data("
#         "upload_id INT,"
#         "sid_meta1_name TEXT,"
#         "sid_meta2_name TEXT,"
#         "sid_meta3_name TEXT,"
#         "contains_crosslink BOOLEAN)"
#     )
#
#     # add meta columns
#     try:
#         cur.execute('ALTER TABLE spectrum_identifications ADD COLUMN meta1 TEXT')
#         cur.execute('ALTER TABLE spectrum_identifications ADD COLUMN meta2 TEXT')
#         cur.execute('ALTER TABLE spectrum_identifications ADD COLUMN meta3 TEXT')
#
#     except Exception:
#         print('{}: Meta columns exist already - not updated'.format(db_name))
#
#     try:
#         # add precursor information from peak list file to DB
#         cur.execute('ALTER TABLE spectra ADD COLUMN precursor_mz TEXT')
#         cur.execute('ALTER TABLE spectra ADD COLUMN precursor_charge TEXT')
#     except Exception:
#         print('{}: spectrum precursor columns exist already - not updated'.format(db_name))
#     con.commit()
#
#     return True
#
#
# for db_name in glob.glob("./dbs/saved/*.db"):
#     con = SQLite.connect(db_name)
#     update_database(con, db_name)
