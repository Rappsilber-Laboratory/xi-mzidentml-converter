from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy import ForeignKey, Text, ForeignKeyConstraint, Integer
from models.base import Base


class AnalysisCollection(Base):
    __tablename__ = "analysiscollection"
    upload_id: Mapped[int] = mapped_column(Integer, ForeignKey("upload.id"), index=True, primary_key=True, nullable=False)
    spectrum_identification_list_ref: Mapped[str] = mapped_column(Text, primary_key=False, nullable=False)
    spectrum_identification_protocol_ref: Mapped[str] = mapped_column(Text, primary_key=False, nullable=False)
    spectra_data_ref: Mapped[str] = mapped_column(Text, primary_key=True, nullable=False)
    ForeignKeyConstraint(
        ["spectrum_identification_protocol_ref", "upload_id"],
        ["spectrumidentificationprotocol.id", "spectrumidentificationprotocol.upload_id"],
    )
