from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy import ForeignKey, Text, Integer, BOOLEAN, ForeignKeyConstraint
from models.base import Base


class PeptideEvidence(Base):
    __tablename__ = "peptideevidence"
    upload_id: Mapped[int] = mapped_column(Integer, ForeignKey("upload.id"), index=True, primary_key=True,
           nullable=False)
    peptide_ref: Mapped[str] = mapped_column(Text, primary_key=True, nullable=False, index=True)
    dbsequence_ref: Mapped[str] = mapped_column(Text, primary_key=True, nullable=False)
    pep_start: Mapped[int] = mapped_column(Integer, primary_key=True, nullable=False)
    is_decoy: Mapped[bool] = mapped_column(BOOLEAN, nullable=True)
    ForeignKeyConstraint(
        ("dbsequence_ref", "upload_id"),
        ("dbsequence.id", "dbsequence.upload_id"),
    )
    ForeignKeyConstraint(
        ("peptide_ref", "upload_id"),
        ("modifiedpeptide.id", "modifiedpeptide.upload_id"),
    )
