from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy import ForeignKey, Text, Integer, JSON, FLOAT
from models.base import Base
from typing import Optional, Any

class ModifiedPeptide(Base):
    __tablename__ = "modifiedpeptide"
    id: Mapped[str] = mapped_column(Text, primary_key=True, nullable=False)
    upload_id: Mapped[str] = mapped_column(Integer, ForeignKey("upload.id"), index=True, primary_key=True, nullable=False)
    base_sequence: Mapped[str] = mapped_column(Text, nullable=False)
    mod_accessions: Mapped[dict[str, Any]] = mapped_column(JSON, nullable=False)
    mod_avg_mass_deltas: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    mod_monoiso_mass_deltas: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    mod_positions: Mapped[dict[str, Any]] = mapped_column(JSON, nullable=False)  # 1-based with 0 = n-terminal and len(pep)+1 = C-terminal
    # following columns are not in xi2 db, but come out of the mzid on the <Peptide>s
    link_site1: Mapped[int] = mapped_column(Integer, nullable=True)
    link_site2: Mapped[int] = mapped_column(Integer, nullable=True)  # only used for storing loop links
    crosslinker_modmass: Mapped[float] = mapped_column(FLOAT, nullable=True)
    crosslinker_pair_id: Mapped[str] = mapped_column(Text, nullable=True)  # yes, its a string
    crosslinker_accession: Mapped[str] = mapped_column(Text, nullable=True)

