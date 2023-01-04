from AbstractCsvParser import AbstractCsvParser
from AbstractCsvParser import CsvParseException

from time import time
import re
import json
import math


class LinksOnlyCsvParser(AbstractCsvParser):
    required_cols = [
        'abspos1',
        'protein1',
        'abspos2',
        'protein2',
    ]

    optional_cols = [
        'passthreshold',
        'score',
        'decoy1',
        'decoy2',
    ]

    def main_loop(self):
        main_loop_start_time = time()
        self.logger.info('main loop LinksOnlyCsvParser - start')

        peptide_evidences = []
        spectrum_identifications = []
        peptides = []

        proteins = set()

        # list of peptides that were already seen - index in list is peptide_id
        # pep sequence including cross-link site and cross-link mass is unique identifier
        seen_peptides = []

        cross_linker_pair_count = 0

        for identification_id, id_item in self.csv_reader.iterrows():  # identification_id, id_item = id_df.iterrows().next()

            # 1 based row number
            row_number = identification_id + 1

            #
            # VALIDITY CHECKS & TYPE CONVERSIONS - ToDo: move type checks/conversions to col level in parse()?
            #
            if id_item['protein2'] == '':
                cross_linked_id_item = False
            else:
                self.contains_crosslinks = True
                cross_linked_id_item = True

            # score
            try:
                score = float(id_item['score'])
            except ValueError:
                raise CsvParseException('Invalid score: %s in row %s' % (id_item['score'], row_number))

            # protein1
            protein_list1 = id_item['protein1'].split(";")
            protein_list1 = [s.strip() for s in protein_list1]
            for p in protein_list1:
                proteins.add(p)

            # decoy1 - if decoy1 is not set fill list with default value (0)
            if id_item['decoy1'] == -1:
                is_decoy_list1 = [False] * len(protein_list1)
            else:
                is_decoy_list1 = []
                for decoy in str(id_item['decoy1']).split(";"):
                    if decoy.lower().strip() == 'true':
                        is_decoy_list1.append(True)
                    elif decoy.lower().strip() == 'false':
                        is_decoy_list1.append(False)
                    else:
                        raise CsvParseException(
                            'Invalid value in Decoy 1: %s in row %s. Allowed values: True, False.'
                            % (id_item['decoy1'], row_number)
                        )

            # absPos1
            abs_pos_list1 = str(id_item['abspos1']).split(";")
            abs_pos_list1 = [s.strip().replace("'", "") for s in abs_pos_list1]

            # protein - decoy - pepPos sensibility check
            if not len(protein_list1) == len(is_decoy_list1):
                raise CsvParseException(
                    'Inconsistent number of protein to decoy values for Protein1 and Decoy1 in row %s!' % row_number)
            if not len(protein_list1) == len(abs_pos_list1):
                raise CsvParseException(
                    'Inconsistent number of protein to pepPos values for Protein1 and PepPos1 in row %s!' % row_number)

            # protein2
            protein_list2 = id_item['protein2'].split(";")
            protein_list2 = [s.strip() for s in protein_list2]
            for p in protein_list2:
                proteins.add(p)

            # decoy2 - if decoy2 is not set fill list with default value (0)
            if id_item['decoy2'] == -1:
                is_decoy_list2 = [False] * len(protein_list2)
            else:
                is_decoy_list2 = []
                for decoy in str(id_item['decoy2']).split(";"):
                    if decoy.lower().strip() == 'true':
                        is_decoy_list2.append(True)
                    elif decoy.lower().strip() == 'false':
                        is_decoy_list2.append(False)
                    else:
                        raise CsvParseException(
                            'Invalid value in Decoy 2: %s in row %s. Allowed values: True, False.'
                            % (id_item['decoy2'], row_number)
                        )

            # pepPos2 - if pepPos2 is not set fill list with default value (-1)
            # ToDo: might need changing for xiUI where pepPos is not optional
            self.logger.info(id_item['abspos2'])
            if id_item['abspos2'] == -1: # or math.isnan(id_item['abspos2']):
                abs_pos_list2 = [-1] * len(protein_list2)
            else:
                abs_pos_list2 = str(id_item['abspos2']).split(";")
                abs_pos_list2 = [s.strip().replace("'", "") for s in abs_pos_list2]

            # protein - decoy - pepPos sensibility check
            if not len(protein_list2) == len(is_decoy_list2):
                raise CsvParseException(
                    'Inconsistent number of protein to decoy values for Protein2 and Decoy2 in row %s!' % row_number)

            #
            # -----Start actual parsing------
            #

            if cross_linked_id_item:
                cross_linker_pair_id = cross_linker_pair_count
                cross_linker_pair_count += 1
            else:
                cross_linker_pair_id = -1  # linear ToDo: -1 or None?

            #
            # PEPTIDE EVIDENCES
            # peptide evidence - 1
            for i in range(len(protein_list1)):

                # peptide - 1
                unique_pep_identifier1 = "%s-%s" % (1, cross_linker_pair_id)

                if unique_pep_identifier1 not in seen_peptides:
                    seen_peptides.append(unique_pep_identifier1)
                    pep1_id = len(seen_peptides) - 1

                    peptide1 = [
                        pep1_id,  # id,
                        "",  # seq_mods,
                        1,  # link_site,
                        None,  # crosslinker_modmass, declare peptide 1 as cl donor: full mass
                        self.upload_id,  # upload_id,
                        cross_linker_pair_id  # crosslinker_pair_id
                    ]
                    peptides.append(peptide1)
                else:
                    pep1_id = seen_peptides.index(unique_pep_identifier1)

                if cross_linked_id_item:
                    # peptide - 2
                    unique_pep_identifier2 = "%s-%s" % (2, cross_linker_pair_id)

                    if unique_pep_identifier2 not in seen_peptides:
                        seen_peptides.append(unique_pep_identifier2)
                        pep2_id = len(seen_peptides) - 1
                        peptide2 = [
                            pep2_id,  # id,
                            "",  # seq_mods,
                            1,  # link_site,
                            None,  # crosslinker_modmass, declare peptide 2 as cl acceptor: 0 mass
                            self.upload_id,  # upload_id,
                            cross_linker_pair_id  # crosslinker_pair_id
                        ]
                        peptides.append(peptide2)
                    else:
                        pep2_id = seen_peptides.index(unique_pep_identifier2)
                else:
                    pep2_id = None

                m = re.search("..\|(.*)\|(.*)\s?", protein_list1[i])
                accession = protein_list1[i]
                if m:
                    accession = m.groups()[0]

                pep_evidence1 = [
                    pep1_id,                # peptide_ref
                    protein_list1[i],       # dbsequence_ref - ToDo: might change to numerical id
                    accession,              # protein_accession
                    int(float(abs_pos_list1[i])),       # was pep start, now absolute position of link
                    is_decoy_list1[i],      # is_decoy
                    self.upload_id          # upload_id
                ]

                peptide_evidences.append(pep_evidence1)

            if cross_linked_id_item and pep1_id != pep2_id:
                # peptide evidence - 2

                if pep2_id is None:
                    raise StandardError('Fatal! peptide id error!')

                for i in range(len(protein_list2)):

                    m = re.search("..\|(.*)\|(.*)\s?", protein_list2[i])
                    accession = protein_list2[i]
                    if m:
                        accession = m.groups()[0]

                    pep_evidence2 = [
                        pep2_id,                # peptide_ref
                        protein_list2[i],       # dbsequence_ref - ToDo: might change to numerical id
                        accession,              # protein_accession
                        int(float(abs_pos_list2[i])),       # was pep_start, now absolute position of link
                        is_decoy_list2[i],      # is_decoy
                        self.upload_id          # upload_id
                    ]

                    peptide_evidences.append(pep_evidence2)

            #
            # SPECTRUM IDENTIFICATIONS
            # ToDo: experimental_mass_to_charge, calculated_mass_to_charge
            scores = json.dumps({'score': score})

            try:
                meta1 = id_item[self.meta_columns[0]]
            except IndexError:
                meta1 = ""
            try:
                meta2 = id_item[self.meta_columns[1]]
            except IndexError:
                meta2 = ""
            try:
                meta3 = id_item[self.meta_columns[2]]
            except IndexError:
                meta3 = ""

            spectrum_identification = [
                identification_id,          # 'id',
                self.upload_id,             # 'upload_id',
                None,                # 'spectrum_id',
                pep1_id,                    # 'pep1_id',
                pep2_id,                    # 'pep2_id',
                None,                     # 'charge_state',
                1,                       # 'rank',
                True,                   # 'pass_threshold',
                None,                  # 'ions',
                scores,                     # 'scores',
                None,                     # 'experimental_mass_to_charge',
                None,                    # 'calculated_mass_to_charge'
                meta1,
                meta2,
                meta3
            ]
            spectrum_identifications.append(spectrum_identification)

        # DBSEQUENCES
        # if self.fasta:
        db_sequences = []
        for prot in proteins:
            try:
                # data = [prot] + self.fasta[prot] + [self.upload_id]
                temp = self.fasta[prot]
                data = [prot, temp[0], temp[1], temp[2], temp[3], self.upload_id]  # surely there's a better way
            except Exception as ke:
                sp_regex = re.compile('(.*)\|(.*)\|(.*)')
                matches = sp_regex.search(prot)
                if matches is not None:
                    data = [matches.group(), matches.group(2), matches.group(3), "", None, self.upload_id]
                else:
                    data = [prot, prot, prot, "", None, self.upload_id]

            db_sequences.append(data)

        # end main loop
        self.logger.info('main loop - done. Time: ' + str(round(time() - main_loop_start_time, 2)) + " sec")

        # once loop is done write data to DB
        db_wrap_up_start_time = time()
        self.logger.info('write spectra to DB - start')
        try:

            self.db.write_peptide_evidences(peptide_evidences, self.cur, self.con)
            self.db.write_peptides(peptides, self.cur, self.con)
            # self.db.write_spectra(spectra, self.cur, self.con)
            self.db.write_spectrum_identifications(spectrum_identifications, self.cur, self.con)
            self.db.write_db_sequences(db_sequences, self.cur, self.con)
            self.con.commit()
        except Exception as e:
            raise e

        self.logger.info('write spectra to DB - start - done. Time: '
                         + str(round(time() - db_wrap_up_start_time, 2)) + " sec")
