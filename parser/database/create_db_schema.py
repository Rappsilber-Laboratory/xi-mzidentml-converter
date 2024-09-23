import logging.config

from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists, drop_database, create_database

from models.base import Base
# noinspection PyUnresolvedReferences
from models import *


def create_db(connection_str):
    engine = create_engine(connection_str)
    if not database_exists(engine.url):
        create_database(engine.url)


def drop_db(connection_str):
    engine = create_engine(connection_str)
    drop_database(engine.url)


def create_schema(connection_str):
    engine = create_engine(connection_str)  # , echo=True)
    Base.metadata.create_all(engine)
    # logging.info(Base.metadata.tables)
    engine.dispose()


if __name__ == "__main__":
    try:
        from config.config_parser import get_conn_str
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            'Database credentials missing! '
            'Change default.database.ini and save as database.ini')
    conn_str = get_conn_str()
    create_db(conn_str)
    create_schema(conn_str)
