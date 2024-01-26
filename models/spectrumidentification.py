from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy import ForeignKey, Text, FLOAT, JSON, BOOLEAN, Integer, ForeignKeyConstraint
from models.base import Base
from typing import Optional, Any


class SpectrumIdentification(Base):
    __tablename__ = "spectrumidentification"
    id: Mapped[str] = mapped_column(Text, primary_key=True, nullable=False)
    upload_id: Mapped[int] = mapped_column(Integer, ForeignKey("upload.id"), index=True, primary_key=True, nullable=False)
    spectrum_id: Mapped[str] = mapped_column(Text, nullable=True)
    spectra_data_ref: Mapped[str] = mapped_column(Text, nullable=True)
    multiple_spectra_identification_id: Mapped[str] = mapped_column(Integer, nullable=True)
    pep1_id: Mapped[str] = mapped_column(Text, nullable=False)
    pep2_id: Mapped[str] = mapped_column(Text, nullable=True)
    charge_state: Mapped[int] = mapped_column(Integer, nullable=True)
    pass_threshold: Mapped[bool] = mapped_column(BOOLEAN, nullable=False)
    rank: Mapped[int] = mapped_column(Integer, nullable=False)
    scores: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    exp_mz: Mapped[float] = mapped_column(FLOAT, nullable=True)
    calc_mz: Mapped[float] = mapped_column(FLOAT, nullable=True)
    sil_id: Mapped[str] = mapped_column(Text, nullable=True) #  null if from csv file
    ForeignKeyConstraint(
        ["spectra_data_ref", "upload_id"],
        ["analysiscollection.spectra_data_ref",
         "analysiscollection.upload_id"],
    ),
    # Can't use this ForeignKeyConstraint, because we want to allow people to upload data
    # without spectra
    # ForeignKeyConstraint(
    #     ["spectrum_id", "spectra_data_ref", "upload_id"],
    #     ["Spectrum.id", "Spectrum.spectra_data_ref", "Spectrum.upload_id"],
    # ),
    ForeignKeyConstraint(
        ["pep1_id", "upload_id"],
        ["modifiedpeptide.id", "modifiedpeptide.upload_id"],
    ),
    ForeignKeyConstraint(
        ["pep2_id", "upload_id"],
        ["modifiedpeptide.id", "modifiedpeptide.upload_id"],
    )
