import pyteomics.mzid as py_mzid
import re
import ntpath
import json
import sys
from time import time
from PeakListParser import PeakListParser
import zipfile
import gzip
import os
from NumpyEncoder import NumpyEncoder


class MzIdParseException(Exception):
    pass


class MzIdParser:
    """

    """

    def __init__(self, mzid_path, temp_dir, peak_list_dir, db, logger, db_name='', user_id=0,
                 origin=''):
        """

        :param mzid_path: path to mzidentML file
        :param temp_dir: absolute path to temp dir for unzipping/storing files
        :param db: database python module to use (xiUI_pg or xiSPEC_sqlite)
        :param db_name: db name for SQLite
        :param origin: ftp dir of pride project
        """

        self.upload_id = 0
        self.mzid_path = mzid_path

        self.peak_list_readers = {}  # peak list readers indexed by spectraData_ref
        self.temp_dir = temp_dir
        if not self.temp_dir.endswith('/'):
            self.temp_dir += '/'
        self.peak_list_dir = peak_list_dir
        if peak_list_dir and not peak_list_dir.endswith('/'):
            self.peak_list_dir += '/'

        self.user_id = user_id
        self.random_id = 0

        self.db = db
        self.logger = logger

        # look up table populated by parse_peptides function
        # self.peptide_id_lookup = {}

        self.spectra_data_protocol_map = {}
        # ToDo: Might change to pyteomics unimod obo module
        self.unimod_path = 'obo/unimod.obo'

        # ToDo: modifications might be globally stored in mzIdentML under
        # ToDo: AnalysisProtocolCollection->SpectrumIdentificationProtocol->ModificationParams
        # ToDo: atm we get them while looping through the peptides
        #  (might be more robust and we're doing it anyway)
        self.modlist = []
        self.unknown_mods = []

        # From mzidentML schema 1.2.0:
        # <SpectrumIdentificationProtocol> must contain the CV term 'cross-linking search'
        # (MS:1002494)
        self.contains_crosslinks = False

        self.warnings = []

        # connect to DB
        try:
            self.con = db.connect(db_name)
            self.cur = self.con.cursor()

        except db.DBException as e:
            self.logger.error(e)
            print(e)
            sys.exit(1)

        self.upload_id = self.db.new_upload([user_id, os.path.basename(self.mzid_path), origin],
                                            self.cur, self.con)

        self.random_id = self.db.get_random_id(self.upload_id, self.cur, self.con)

        self.upload_info_read = False
        self.mzid_reader = None

    def initialise_mzid_reader(self):
        if self.mzid_path.endswith('.gz') or self.mzid_path.endswith('.zip'):
            self.mzid_path = MzIdParser.extract_mzid(self.mzid_path)

        self.logger.info('reading mzid - start ' + self.mzid_path)
        start_time = time()
        # schema:
        # https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/schema/mzIdentML1.2.0.xsd
        try:
            self.mzid_reader = py_mzid.MzIdentML(self.mzid_path)
        except Exception as e:
            raise MzIdParseException(type(e).__name__, e.args)

        self.logger.info('reading mzid - done. Time: {} sec'.format(round(time() - start_time, 2)))

    # used by TestLoop when downloading files from PRIDE
    def get_supported_peak_list_file_names(self):
        """
        :return: list of all supported peak list file names
        """
        peak_list_file_names = []
        for spectra_data_id in self.mzid_reader._offset_index["SpectraData"].keys():
            sp_datum = self.mzid_reader.get_by_id(spectra_data_id, tag_id='SpectraData',
                                                  detailed=True)
            ff_acc = sp_datum['FileFormat']['accession']
            if any([ff_acc == 'MS:1001062',  # MGF
                    ff_acc == 'MS:1000584',  # mzML
                    ff_acc == 'MS:1001466',  # ms2
                    ]):
                peak_list_file_names.append(ntpath.basename(sp_datum['location']))

        return peak_list_file_names

    # used by TestLoop when downloading files from PRIDE
    def get_all_peak_list_file_names(self):
        """
        :return: list of all peak list file names
        """
        peak_list_file_names = []
        for spectra_data_id in self.mzid_reader._offset_index["SpectraData"].keys():
            sp_datum = self.mzid_reader.get_by_id(spectra_data_id, tag_id='SpectraData',
                                                  detailed=True)
            peak_list_file_names.append(ntpath.basename(sp_datum['location']))

        return peak_list_file_names

    def init_peak_list_readers(self):
        """
        sets self.peak_list_readers by looping through SpectraData elements
        dictionary:
            key: spectra_data_ref
            value: associated peak_list_reader
        """
        peak_list_readers = {}
        for spectra_data_id in self.mzid_reader._offset_index["SpectraData"].keys():
            sp_datum = self.mzid_reader.get_by_id(spectra_data_id, tag_id='SpectraData',
                                                  detailed=True)

            self.check_spectra_data_validity(sp_datum)

            sd_id = sp_datum['id']
            peak_list_file_name = ntpath.basename(sp_datum['location'])
            peak_list_file_path = self.peak_list_dir + peak_list_file_name

            try:
                peak_list_reader = PeakListParser(
                    peak_list_file_path,
                    sp_datum['FileFormat']['accession'],
                    sp_datum['SpectrumIDFormat']['accession']
                )
            except Exception:
                # try gz version
                try:
                    peak_list_reader = PeakListParser(
                        PeakListParser.extract_gz(peak_list_file_path + '.gz'),
                        sp_datum['FileFormat']['accession'],
                        sp_datum['SpectrumIDFormat']['accession']
                    )
                except IOError:
                    raise MzIdParseException('Missing peak list file: %s' % peak_list_file_path)

            peak_list_readers[sd_id] = peak_list_reader

        self.peak_list_readers = peak_list_readers

    def check_all_spectra_data_validity(self):
        for spectra_data_id in self.mzid_reader._offset_index["SpectraData"].keys():
            sp_datum = self.mzid_reader.get_by_id(spectra_data_id, tag_id='SpectraData',
                                                  detailed=True)
            self.check_spectra_data_validity(sp_datum)

    @staticmethod
    def check_spectra_data_validity(sp_datum):
        # is there anything we'd like to complain about?
        # SpectrumIDFormat
        if 'SpectrumIDFormat' not in sp_datum or sp_datum['SpectrumIDFormat'] is None:
            raise MzIdParseException('SpectraData is missing SpectrumIdFormat')
        if isinstance(sp_datum['SpectrumIDFormat'], basestring):
            raise MzIdParseException('SpectraData/SpectrumIdFormat is missing accession')
        if sp_datum['SpectrumIDFormat']['accession'] is None:
            raise MzIdParseException('SpectraData/SpectrumIdFormat is missing accession')

        # FileFormat
        if 'FileFormat' not in sp_datum or sp_datum['FileFormat'] is None:
            raise MzIdParseException('SpectraData is missing FileFormat')
        if isinstance(sp_datum['FileFormat'], basestring):
            raise MzIdParseException('SpectraData/SpectrumIdFormat is missing accession')
        if sp_datum['FileFormat']['accession'] is None:
            raise MzIdParseException('SpectraData/FileFormat is missing accession')

        # location
        if 'location' not in sp_datum or sp_datum['location'] is None:
            raise MzIdParseException('SpectraData is missing location')

    def parse(self):

        start_time = time()

        if not self.upload_info_read:
            self.upload_info()  # overridden (empty function) in xiSPEC subclass

        if self.peak_list_dir:
            self.init_peak_list_readers()

        self.parse_db_sequences()  # overridden (empty function) in xiSPEC subclass
        self.parse_peptides()
        self.parse_peptide_evidences()
        self.map_spectra_data_to_protocol()
        self.main_loop()

        # meta_data = [self.upload_id, -1, -1, -1, -1]
        # self.db.write_meta_data(meta_data, self.cur, self.con)

        self.fill_in_missing_scores()  # empty here, overridden in xiSPEC subclass to do stuff

        self.other_info()

        self.logger.info('all done! Total time: ' + str(round(time() - start_time, 2)) + " sec")

        self.con.close()

    def get_ion_types_mzid(self, sid_item):
        try:
            fragmentation = sid_item['Fragmentation']
            ion_names_list = [i['name'] for i in fragmentation['IonType']]
            ion_names_list = list(set(ion_names_list))
        except KeyError:
            return []

        # ion_types = ["P"]
        ion_types = []
        for ion_name in ion_names_list:
            try:
                ion = re.search('frag: ([a-z]) ion', ion_name).groups()[0]
                ion_types.append(ion)
            except (IndexError, AttributeError) as e:
                self.logger.info(e, ion_name)
                continue

        return ion_types

    # split into two functions
    @staticmethod
    def extract_mzid(archive):
        if archive.endswith('zip'):
            zip_ref = zipfile.ZipFile(archive, 'r')
            unzip_path = archive + '_unzip/'
            zip_ref.extractall(unzip_path)
            zip_ref.close()

            return_file_list = []

            for root, dir_names, file_names in os.walk(unzip_path):
                file_names = [f for f in file_names if not f[0] == '.']
                dir_names[:] = [d for d in dir_names if not d[0] == '.']
                for file_name in file_names:
                    os.path.join(root, file_name)
                    if file_name.lower().endswith('.mzid'):
                        return_file_list.append(root + '/' + file_name)
                    else:
                        raise IOError('unsupported file type: %s' % file_name)

            if len(return_file_list) > 1:
                raise StandardError("more than one mzid file found!")

            return return_file_list[0]

        elif archive.endswith('gz'):
            in_f = gzip.open(archive, 'rb')
            archive = archive.replace(".gz", "")
            out_f = open(archive, 'wb')
            try:
                out_f.write(in_f.read())
            except IOError:
                raise StandardError('Zip archive error: %s' % archive)

            in_f.close()
            out_f.close()

            return archive

        else:
            raise StandardError('unsupported file type: %s' % archive)

    def map_spectra_data_to_protocol(self):
        """
        extract and map spectrumIdentificationProtocol which includes annotation data like fragment
         tolerance only fragment tolerance is extracted for now
        # ToDo: improve error handling
        #       extract modifications, cl mod mass, ...

        Parameters:
        ------------------------
        mzid_reader: pyteomics mzid_reader
        """

        self.logger.info('generating spectra data protocol map - start')
        start_time = time()

        spectra_data_protocol_map = {}

        sid_protocols = []

        analysis_collection = self.mzid_reader.iterfind('AnalysisCollection').next()
        for spectrumIdentification in analysis_collection['SpectrumIdentification']:
            sid_protocol_ref = spectrumIdentification['spectrumIdentificationProtocol_ref']
            sid_protocol = self.mzid_reader.get_by_id(sid_protocol_ref,
                                                      tag_id='SpectrumIdentificationProtocol',
                                                      detailed=True)
            sid_protocols.append(sid_protocol)
            try:
                frag_tol = sid_protocol['FragmentTolerance']
                frag_tol_plus = frag_tol['search tolerance plus value']
                frag_tol_value = re.sub('[^0-9,.]', '', str(frag_tol_plus['value']))
                if frag_tol_plus['unit'].lower() == 'parts per million':
                    frag_tol_unit = 'ppm'
                elif frag_tol_plus['unit'].lower() == 'dalton':
                    frag_tol_unit = 'Da'
                else:
                    frag_tol_unit = frag_tol_plus['unit']

                if not all([
                    frag_tol['search tolerance plus value']['value'] == frag_tol['search tolerance minus value'][
                        'value'],
                    frag_tol['search tolerance plus value']['unit'] == frag_tol['search tolerance minus value']['unit']
                ]):
                    self.warnings.append(
                        {"type": "mzidParseError",
                         "message": "search tolerance plus value doesn't match minus value. Using plus value!"})

            except KeyError:
                self.warnings.append({
                    "type": "mzidParseError",
                    "message": "could not parse ms2tolerance. Falling back to default: 10 ppm.",
                    # 'id': id_string
                })
                frag_tol_value = '10'
                frag_tol_unit = 'ppm'
                # spectra_data_protocol_map['errors'].append(
                #     {"type": "mzidParseError",
                #      "message": "could not parse ms2tolerance. Falling back to default values."})

            for inputSpectra in spectrumIdentification['InputSpectra']:
                spectra_data_ref = inputSpectra['spectraData_ref']

                spectra_data_protocol_map[spectra_data_ref] = {
                    'protocol_ref': sid_protocol_ref,
                    'fragmentTolerance': ' '.join([frag_tol_value, frag_tol_unit])
                }

        self.mzid_reader.reset()
        self.spectra_data_protocol_map = spectra_data_protocol_map
        self.logger.info('generating spectraData_ProtocolMap - done. Time: {} sec'.format(
            round(time() - start_time, 2)))

        # self.db.write.protocols()

    def add_to_modlist(self, mod):
        if mod['name'] == "unknown_modification":
            mod['name'] = "({0:.2f})".format(mod['monoisotopicMassDelta'])

        mod['monoisotopicMassDelta'] = float(mod['monoisotopicMassDelta'])

        mod['residues'] = [aa for aa in mod['residues']]

        if mod['name'] in [m['name'] for m in self.modlist]:
            old_mod = self.modlist[[m['name'] for m in self.modlist].index(mod['name'])]
            # check if modname with different mass exists already
            if mod['monoisotopicMassDelta'] != old_mod['monoisotopicMassDelta']:
                mod['name'] += "*"
                self.add_to_modlist(mod)
            else:
                for res in mod['residues']:
                    if res not in old_mod['residues']:
                        old_mod['residues'].append(res)
        else:
            self.modlist.append(mod)

        return mod['name']

    def parse_db_sequences(self):

        self.logger.info('parse db sequences - start')
        start_time = time()
        # DBSEQUENCES
        inj_list = []
        for db_id in self.mzid_reader._offset_index["DBSequence"].keys():
            db_sequence = self.mzid_reader.get_by_id(db_id, tag_id='DBSequence', detailed=True)

            data = [db_sequence["id"], db_sequence["accession"]]

            # name, optional elem att
            if "name" in db_sequence:
                data.append(db_sequence["name"])
            else:
                data.append(db_sequence["accession"])

            # description, officially not there?
            if "protein description" in db_sequence:
                data.append(json.dumps(db_sequence["protein description"], cls=NumpyEncoder))
            else:
                data.append(None)

            # searchDatabase_ref

            # Seq is optional child elem of DBSequence
            if "Seq" in db_sequence and isinstance(db_sequence["Seq"], basestring):
                seq = db_sequence["Seq"]
                data.append(seq)
            elif "length" in db_sequence:
                data.append("X" * db_sequence["length"])
            else:
                # todo: get sequence
                data.append("")

            data.append(self.upload_id)

            inj_list.append(data)

        self.db.write_db_sequences(inj_list, self.cur, self.con)

        self.logger.info('parse db sequences - done. Time: {} sec'.format(
            round(time() - start_time, 2)))

    def parse_peptides(self):
        start_time = time()
        self.logger.info('parse peptides, modifications - start')

        # ToDo: might be stuff in pyteomics lib for this?
        unimod_masses = self.get_unimod_masses(self.unimod_path)

        # PEPTIDES
        peptide_index = 0
        peptide_inj_list = []
        for pep_id in self.mzid_reader._offset_index["Peptide"].keys():
            peptide = self.mzid_reader.get_by_id(pep_id, tag_id='Peptide', detailed=True)
            pep_seq_dict = []
            for aa in peptide['PeptideSequence']:
                pep_seq_dict.append({"Modification": "", "aminoAcid": aa})

            link_site = -1
            crosslinker_modmass = None
            value = None

            # MODIFICATIONS
            # add in modifications
            if 'Modification' in peptide.keys():
                for mod in peptide['Modification']:

                    if 'monoisotopicMassDelta' not in mod.keys():
                        try:
                            mod['monoisotopicMassDelta'] = unimod_masses[mod['accession']]

                        # ToDo: what's going on here?
                        except KeyError:
                            # seq_ref_prot_map['errors'].append({
                            #     "type": "mzidParseError",
                            #     "message": "could not get modification mass for modification {}".format(mod),
                            #     "id": mod["id"]
                            # })
                            continue

                    # link_index = 0  # TODO: multilink support
                    # mod_location is 0-based for assigning modifications to correct amino acid
                    # mod['location'] is 1-based with 0 = n-terminal and len(pep)+1 = C-terminal
                    if mod['location'] == 0:
                        mod_location = 0
                        # n_terminal_mod = True
                    elif mod['location'] == len(peptide['PeptideSequence']) + 1:
                        mod_location = mod['location'] - 2
                        # c_terminal_mod = True
                    else:
                        mod_location = mod['location'] - 1
                        # n_terminal_mod = False
                        # c_terminal_mod = False
                    if 'residues' not in mod:
                        mod['residues'] = peptide['PeptideSequence'][mod_location]

                    # TODO - issues here with using names rather than cv param accession
                    #  (cross-link acceptor/ receiver)
                    if 'name' in mod.keys():
                        # fix mod names
                        if isinstance(mod['name'], list):  # todo: have a look at this  - cc
                            mod['name'] = ','.join(mod['name'])
                        mod['name'] = mod['name'].lower()
                        mod['name'] = mod['name'].replace(" ", "_")
                        if 'cross-link donor' not in mod.keys() and 'cross-link acceptor' not in mod.keys() and 'cross-link receiver' not in mod.keys():
                            cur_mod = pep_seq_dict[mod_location]
                            # join modifications into one for multiple modifications on the same aa
                            if not cur_mod['Modification'] == '':
                                mod['name'] = '_'.join(sorted([cur_mod['Modification'], mod['name']], key=str.lower))
                                cur_mod_mass = [x['monoisotopicMassDelta'] for x in self.modlist if
                                                x['name'] == cur_mod['Modification']][0]
                                mod['monoisotopicMassDelta'] += cur_mod_mass

                            # save to all mods list and get back new_name
                            mod['name'] = self.add_to_modlist(mod)
                            cur_mod['Modification'] = mod['name']

                    # error handling for mod without name
                    else:
                        # cross-link acceptor doesn't have a name
                        if 'cross-link acceptor' not in mod.keys() and 'cross-link receiver' not in mod.keys():
                            raise MzIdParseException("Missing modification name")

                    # add CL locations
                    if 'cross-link donor' in mod.keys() or 'cross-link acceptor' in mod.keys() \
                            or 'cross-link receiver' in mod.keys():
                        # use mod['location'] for link-site (1-based in database in line with mzIdentML specifications)
                        link_site = mod['location']
                        crosslinker_modmass = mod['monoisotopicMassDelta']

                    if 'cross-link acceptor' in mod.keys():
                        value = mod['cross-link acceptor']['value']
                    if 'cross-link donor' in mod.keys():
                        value = mod['cross-link donor']['value']
                    if 'cross-link receiver' in mod.keys():
                        value = mod['cross-link receiver']['value']

            # ToDo: we should consider swapping these over because modX format has modification
            #  before AA
            peptide_seq_with_mods = ''.join(
                [''.join([x['aminoAcid'], x['Modification']]) for x in pep_seq_dict])

            data = [
                # peptide_index,      # debug use mzid peptide['id'],
                peptide['id'],
                peptide_seq_with_mods,
                link_site,
                crosslinker_modmass,
                self.upload_id,
                str(value)
            ]

            peptide_inj_list.append(data)
            #  self.peptide_id_lookup[peptide['id']] = peptide_index

            if peptide_index % 1000 == 0:
                self.logger.info('writing 1000 peptides to DB')
                try:
                    self.db.write_peptides(peptide_inj_list, self.cur, self.con)
                    peptide_inj_list = []
                    self.con.commit()
                except Exception as e:
                    raise e

            peptide_index += 1

        try:
            self.db.write_peptides(peptide_inj_list, self.cur, self.con)
            self.con.commit()
        except Exception as e:
            raise e
        #
        mod_index = 0
        modifications_inj_list = []
        for mod in self.modlist:
            try:
                mod_accession = mod['accession']
            except KeyError:
                mod_accession = ''
            modifications_inj_list.append([
                mod_index,
                self.upload_id,
                mod['name'],
                mod['monoisotopicMassDelta'],
                ''.join(mod['residues']),
                mod_accession
            ])
            mod_index += 1
        self.db.write_modifications(modifications_inj_list, self.cur, self.con)

        self.logger.info('parse peptides, modifications - done. Time: {} sec'.format(
            round(time() - start_time, 2)))

    def parse_peptide_evidences(self):
        start_time = time()
        self.logger.info('parse peptide evidences - start')

        seq_id_to_acc_map = {}

        for db_id in self.mzid_reader._offset_index["DBSequence"].keys():
            db_sequence = self.mzid_reader.get_by_id(db_id, tag_id='DBSequence', detailed=True)
            seq_id_to_acc_map[db_sequence["id"]] = db_sequence["accession"]

        # PEPTIDE EVIDENCES
        inj_list = []
        for pep_ev_id in self.mzid_reader._offset_index["PeptideEvidence"].keys():
            peptide_evidence = self.mzid_reader.get_by_id(pep_ev_id, tag_id='PeptideEvidence',
                                                          detailed=True)

            pep_start = -1
            if "start" in peptide_evidence:
                pep_start = peptide_evidence["start"]  # start att, optional

            is_decoy = False
            if "isDecoy" in peptide_evidence:
                is_decoy = peptide_evidence["isDecoy"]  # isDecoy att, optional

            # peptide_ref = self.peptide_id_lookup[peptide_evidence["peptide_ref"]]
            peptide_ref = peptide_evidence["peptide_ref"]  # debug use mzid peptide['id'],

            data = [
                peptide_ref,                                                 # 'peptide_ref',
                peptide_evidence["dBSequence_ref"],                          # 'dbsequence_ref',
                seq_id_to_acc_map[peptide_evidence["dBSequence_ref"]],       # 'protein_accession',
                pep_start,                                                   # 'pep_start',
                is_decoy,                                                    # 'is_decoy',
                self.upload_id                                               # 'upload_id'
            ]

            inj_list.append(data)

            if len(inj_list) % 1000 == 0:
                self.logger.info('writing 1000 peptide_evidences to DB')
                try:
                    self.db.write_peptide_evidences(inj_list, self.cur, self.con)
                    inj_list = []
                    self.con.commit()
                except Exception as e:
                    raise e

        try:
            self.db.write_peptide_evidences(inj_list, self.cur, self.con)
            self.con.commit()
        except Exception as e:
            raise e

        self.con.commit()
        self.mzid_reader.reset()

        self.logger.info('parse peptide evidences - done. Time: {} sec'.format(
            round(time() - start_time, 2)))

    @staticmethod
    def get_unimod_masses(unimod_path):
        masses = {}
        mod_id = -1

        with open(unimod_path) as f:
            for line in f:
                if line.startswith('id: '):
                    mod_id = ''.join(line.replace('id: ', '').split())

                elif line.startswith('xref: delta_mono_mass ') and not mod_id == -1:
                    mass = float(line.replace('xref: delta_mono_mass ', '').replace('"', ''))
                    masses[mod_id] = mass

        return masses

    def main_loop(self):
        spec_id = 0
        identification_id = 0
        spectra = []
        spectrum_identifications = []

        fragment_parsing_error_scans = []

        #
        # main loop
        main_loop_start_time = time()
        self.logger.info('main loop - start')

        for sid_result in self.mzid_reader:
            if self.peak_list_dir:
                peak_list_reader = self.peak_list_readers[sid_result['spectraData_ref']]

                scan_id = peak_list_reader.parse_scan_id(sid_result["spectrumID"])
                scan = peak_list_reader.get_scan(scan_id)

                protocol = self.spectra_data_protocol_map[sid_result['spectraData_ref']]

                if scan['precursor'] is not None:
                    precursor_mz = scan['precursor']['mz']
                    precursor_charge = scan['precursor']['charge']
                else:
                    # give warning precursor info is missing
                    precursor_mz = None
                    precursor_charge = None

                peak_list = scan['peaks']
                mz = []
                intensity = []

                if peak_list:
                    peak_list = peak_list.strip()
                    chunks = re.split('\n+', peak_list)
                    # print(chunks);
                    try:
                        for peak in chunks:
                            # print(peak)
                            peak.strip()
                            if len(peak) != 0:
                                values = re.split('\s', peak)
                                # print('mz:' + values[0] + ' int:' + values[1])
                                mz.append(float(values[0]))
                                intensity.append(float(values[1]))
                    except ValueError as ve:
                        print('Error {ve}')
                        print(peak_list)
                        print(len(mz))
                        print(len(intensity))

                spectra.append([
                    spec_id,
                    mz,
                    intensity,
                    ntpath.basename(peak_list_reader.peak_list_path),
                    str(scan_id),
                    protocol['fragmentTolerance'],
                    self.upload_id,
                    sid_result['id'],
                    precursor_mz,
                    precursor_charge
                ])

            spectrum_ident_dict = dict()
            linear_index = -1  # negative index values for linear peptides

            for spec_id_item in sid_result['SpectrumIdentificationItem']:
                # get suitable id
                if 'cross-link spectrum identification item' in spec_id_item.keys():
                    self.contains_crosslinks = True
                    cross_link_id = spec_id_item['cross-link spectrum identification item']
                else:  # assuming linear
                    # misusing 'cross-link spectrum identification item'
                    # for linear peptides with negative index
                    # specIdItem['cross-link spectrum identification item'] = linear_index
                    # spec_id_set.add(get_cross_link_identifier(specIdItem))

                    cross_link_id = linear_index
                    linear_index -= 1

                # check if seen it before
                if cross_link_id in spectrum_ident_dict.keys():
                    # do crosslink specific stuff
                    ident_data = spectrum_ident_dict.get(cross_link_id)
                    # ident_data[4] = self.peptide_id_lookup[spec_id_item['peptide_ref']]
                    ident_data[4] = spec_id_item['peptide_ref']  # debug
                else:
                    # do stuff common to linears and crosslinks
                    charge_state = spec_id_item['chargeState']
                    pass_threshold = spec_id_item['passThreshold']
                    # ToDo: refactor with MS: cv Param list of all scores
                    scores = {
                        k: v for k, v in spec_id_item.iteritems()
                        if 'score' in k.lower() or
                           'pvalue' in k.lower() or
                           'evalue' in k.lower() or
                           'sequest' in k.lower() or
                           'scaffold' in k.lower()
                    }
                    #
                    # fragmentation ions
                    # ToDo: do we want to make assumptions of fragIon types by fragMethod from mzML?
                    ions = self.get_ion_types_mzid(spec_id_item)
                    # if no ion types are specified in the id file check the mzML file
                    # if len(ions) == 0 and peak_list_reader['fileType'] == 'mzml':
                    #     ions = peakListParser.get_ion_types_mzml(scan)

                    ions = list(set(ions))

                    if len(ions) == 0:
                        ions = ['peptide', 'b', 'y']
                        # ToDo: better error handling for general errors -
                        #  bundling together of same type errors
                        fragment_parsing_error_scans.append(sid_result['id'])

                    ions = ';'.join(ions)

                    # extract other useful info to display
                    rank = spec_id_item['rank']

                    # from mzidentML schema 1.2.0: For PMF data, the rank attribute may be
                    # meaningless and values of rank = 0 should be given.
                    # xiSPEC front-end expects rank = 1 as default
                    if rank is None or int(rank) == 0:
                        rank = 1

                    experimental_mass_to_charge = spec_id_item['experimentalMassToCharge']
                    try:
                        calculated_mass_to_charge = spec_id_item['calculatedMassToCharge']
                    except KeyError:
                        calculated_mass_to_charge = None

                    ident_data = [
                        identification_id,
                        # spec_id_item['id'],
                        self.upload_id,
                        spec_id,
                        # self.peptide_id_lookup[spec_id_item['peptide_ref']], # debug use spec_id_item['peptide_ref'],
                        spec_id_item['peptide_ref'],
                        '',  # pep2
                        charge_state,
                        rank,
                        pass_threshold,
                        ions,
                        json.dumps(scores),
                        experimental_mass_to_charge,
                        calculated_mass_to_charge,
                        "",
                        "",
                        ""
                    ]

                    spectrum_ident_dict[cross_link_id] = ident_data

                    identification_id += 1

            spectrum_identifications += spectrum_ident_dict.values()

            spec_id += 1

            if spec_id % 1000 == 0:
                self.logger.info('writing 1000 entries (1000 spectra and their idents) to DB')
                try:
                    self.db.write_spectra(spectra, self.cur, self.con)
                    spectra = []
                    self.db.write_spectrum_identifications(spectrum_identifications, self.cur,
                                                           self.con)
                    spectrum_identifications = []
                    self.con.commit()
                except Exception as e:
                    raise e

        # end main loop
        self.logger.info('main loop - done Time: {} sec'.format(
            round(time() - main_loop_start_time, 2)))

        # once loop is done write remaining data to DB
        db_wrap_up_start_time = time()
        self.logger.info('write remaining entries to DB - start')
        try:
            self.db.write_spectra(spectra, self.cur, self.con)
            self.db.write_spectrum_identifications(spectrum_identifications, self.cur, self.con)
            self.con.commit()
        except Exception as e:
            raise e

        self.logger.info('write remaining entries to DB - start - done.  Time: {} sec'.format(
            round(time() - db_wrap_up_start_time, 2)))

        self.ident_count = identification_id

        # warnings
        if len(fragment_parsing_error_scans) > 0:
            if len(fragment_parsing_error_scans) > 50:
                id_string = '; '.join(fragment_parsing_error_scans[:50]) + ' ...'
            else:
                id_string = '; '.join(fragment_parsing_error_scans)

            self.warnings.append({
                "type": "IonParsing",
                "message": "mzidentML file does not specify fragment ions.",
                'id': id_string
            })

    def upload_info(self):
        self.upload_info_read = True
        upload_info_start_time = time()
        self.logger.info('parse upload info - start')

        peak_list_file_names = json.dumps(self.get_all_peak_list_file_names(), cls=NumpyEncoder)

        spectra_formats = []
        for spectra_data_id in self.mzid_reader._offset_index["SpectraData"].keys():
            sp_datum = self.mzid_reader.get_by_id(spectra_data_id, tag_id='SpectraData',
                                                  detailed=True)
            spectra_formats.append(sp_datum)
        spectra_formats = json.dumps(spectra_formats, cls=NumpyEncoder)

        # AnalysisSoftwareList - optional element
        # see https://groups.google.com/forum/#!topic/pyteomics/Mw4eUHmicyU
        self.mzid_reader.schema_info['lists'].add("AnalysisSoftware")
        try:
            analysis_software = json.dumps(self.mzid_reader.iterfind('AnalysisSoftwareList').next()['AnalysisSoftware'])
        except StopIteration:
            analysis_software = '{}'
        except Exception as e:
            raise MzIdParseException(type(e).__name__, e.args)
        self.mzid_reader.reset()

        # Provider - optional element
        try:
            provider = json.dumps(self.mzid_reader.iterfind('Provider').next())
        except StopIteration:
            provider = '{}'
        except Exception as e:
            raise MzIdParseException(type(e).__name__, e.args)
        self.mzid_reader.reset()

        # AuditCollection - optional element
        audits = '{}'
        try:
            audits = json.dumps(self.mzid_reader.iterfind('AuditCollection').next())
        except StopIteration:
            audits = '{}'
        except Exception as e:
            raise MzIdParseException(type(e).__name__, e.args)
        self.mzid_reader.reset()

        # AnalysisSampleCollection - optional element
        try:
            samples = json.dumps(self.mzid_reader.iterfind('AnalysisSampleCollection').next()['Sample'])
        except StopIteration:
            samples = '{}'
        except Exception as e:
            raise MzIdParseException(type(e).__name__, e.args)
        self.mzid_reader.reset()

        # AnalysisCollection - required element, StopIteration exception shouldn't happen
        try:
            analyses = json.dumps(self.mzid_reader.iterfind('AnalysisCollection').next()['SpectrumIdentification'])
        except StopIteration:
            analyses = '{}'  # could legitimately throw error here instead, its required
        except Exception as e:
            raise MzIdParseException(type(e).__name__, e.args)
        self.mzid_reader.reset()

        # AnalysisProtocolCollection - required element
        try:
            protocol_collection = self.mzid_reader.iterfind('AnalysisProtocolCollection').next()
            protocols = json.dumps(protocol_collection['SpectrumIdentificationProtocol'],
                                   cls=NumpyEncoder)
        except StopIteration:
            protocols = '{}'  # could legitimately throw error here instead, its required
        except Exception as e:
            raise MzIdParseException(type(e).__name__, e.args)
        self.mzid_reader.reset()

        # BibliographicReference - optional element
        bibRefs = []
        for bib in self.mzid_reader.iterfind('BibliographicReference'):
            bibRefs.append(bib)
        bibRefs = json.dumps(bibRefs)
        self.mzid_reader.reset()

        self.db.write_mzid_info(peak_list_file_names,
                                spectra_formats,
                                analysis_software,
                                provider,
                                audits,
                                samples,
                                analyses,
                                protocols,
                                bibRefs,
                                self.upload_id, self.cur, self.con)

        self.logger.info('getting upload info - done  Time: {} sec'.format(
            round(time() - upload_info_start_time, 2)))

    def fill_in_missing_scores(self):
        pass

    def other_info(self):
        ident_file_size = os.path.getsize(self.mzid_path)
        self.db.write_other_info(self.upload_id, self.contains_crosslinks, self.ident_count,
                                 ident_file_size, self.warnings, self.cur, self.con)


class xiSPEC_MzIdParser(MzIdParser):

    def upload_info(self):
        pass

    def parse_db_sequences(self):
        pass

    def fill_in_missing_scores(self):
        # Fill missing scores with
        score_fill_start_time = time()
        self.logger.info('fill in missing scores - start')
        self.db.fill_in_missing_scores(self.cur, self.con)
        self.logger.info('fill in missing scores - done. Time: {}'.format(
            round(time() - score_fill_start_time, 2)))

    def other_info(self):
        pass
