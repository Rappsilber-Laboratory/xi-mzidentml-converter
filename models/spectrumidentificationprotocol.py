from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy import ForeignKey, Text, JSON, Integer, Float
from models.base import Base
from typing import Optional, Any


class SpectrumIdentificationProtocol(Base):
    __tablename__ = "spectrumidentificationprotocol"
    id: Mapped[str] = mapped_column(Text, primary_key=True, nullable=False)
    upload_id: Mapped[str] = mapped_column(Integer, ForeignKey("upload.id"), index=True, primary_key=True, nullable=False)
    analysis_software: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    search_type: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    additional_search_params: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    threshold: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    frag_tol: Mapped[float] = mapped_column(Float, nullable=False)
    frag_tol_unit: Mapped[str] = mapped_column(Text, nullable=False)
