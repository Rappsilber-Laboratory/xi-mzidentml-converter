from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy import Text, JSON, BOOLEAN, TIMESTAMP, func, Integer
from models.base import Base
from typing import Optional, Any
import datetime


class Upload(Base):
    __tablename__ = "upload"
    id: Mapped[int] = mapped_column(Integer, primary_key=True,  autoincrement=True, nullable=False)
    project_id: Mapped[str] = mapped_column(Text, nullable=True)
    identification_file_name: Mapped[str] = mapped_column(Text, nullable=False)
    analysis_software_list: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    provider: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    audit_collection: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    analysis_sample_collection: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    bib: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)
    spectra_formats: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)  # nullable=False
    upload_time: Mapped[datetime.datetime] = mapped_column(TIMESTAMP, server_default=func.now(), nullable=False)
    contains_crosslinks: Mapped[bool] = mapped_column(BOOLEAN, nullable=True)  # nullable=False
    upload_warnings: Mapped[Optional[dict[str, Any]]] = mapped_column(JSON, nullable=True)  # nullable=False
    identification_file_name_clean: Mapped[str] = mapped_column(Text, nullable=True)


