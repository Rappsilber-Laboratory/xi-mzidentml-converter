from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy import ForeignKey, Text, FLOAT, JSON, BOOLEAN, Integer, ForeignKeyConstraint, CHAR
from models.base import Base
from typing import Optional, Any


class Match(Base):
    __tablename__ = "match"
    id: Mapped[str] = mapped_column(Text, primary_key=True, nullable=False)
    upload_id: Mapped[int] = mapped_column(Integer, ForeignKey("upload.id"), index=True, primary_key=True,
                                           nullable=False)
    spectrum_id: Mapped[str] = mapped_column(Text, nullable=True)
    spectra_data_id: Mapped[int] = mapped_column(Integer, nullable=True)  # nullable for csv data
    multiple_spectra_identification_id: Mapped[str] = mapped_column(Integer, nullable=True)
    multiple_spectra_identification_pc: Mapped[str] = mapped_column(CHAR, nullable=True)
    pep1_id: Mapped[int] = mapped_column(Integer, nullable=False)
    pep2_id: Mapped[int] = mapped_column(Integer, nullable=True)
    charge_state: Mapped[int] = mapped_column(Integer, nullable=True)
    pass_threshold: Mapped[bool] = mapped_column(BOOLEAN, nullable=False)
    rank: Mapped[int] = mapped_column(Integer, nullable=False)
    scores: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    exp_mz: Mapped[float] = mapped_column(FLOAT, nullable=True)
    calc_mz: Mapped[float] = mapped_column(FLOAT, nullable=True)
    sip_id: Mapped[int] = mapped_column(Integer, nullable=True)  # null if from csv file
    __table_args__ = (
        ForeignKeyConstraint(
            ["sip_id", "upload_id"],
            ["spectrumidentificationprotocol.id",
             "spectrumidentificationprotocol.upload_id"],
        ),
        ForeignKeyConstraint(
            ["pep1_id", "upload_id"],
            ["modifiedpeptide.id", "modifiedpeptide.upload_id"],
        ),
        ForeignKeyConstraint(
            ["pep2_id", "upload_id"],
            ["modifiedpeptide.id", "modifiedpeptide.upload_id"],
        ),
        # ForeignKeyConstraint(
        # ["spectrum_id", "spectra_data_id", "upload_id"],
        # ["spectrum.id", "spectrum.spectra_data_id", "spectrum.upload_id"],
        # ),
    )

