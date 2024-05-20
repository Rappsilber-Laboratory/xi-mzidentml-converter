from abc import ABC, abstractmethod


# Strategy interface
class Writer(ABC):

    def __init__(self, upload_id=None, pxid=None):
        self.pxid = pxid
        self.upload_id = upload_id

    @abstractmethod
    def write_data(self, table, data):
        pass

    @abstractmethod
    def write_new_upload(self, table, data):
        pass

    @abstractmethod
    def write_mzid_info(self, analysis_software_list, spectra_formats,
                        provider, audits, samples, bib, upload_id):
        pass

    @abstractmethod
    def fill_in_missing_scores(self):
        pass
