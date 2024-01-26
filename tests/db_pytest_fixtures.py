from parser.database.create_db_schema import create_schema, create_db, drop_db
import pytest
from sqlalchemy import create_engine


@pytest.fixture()
def db_info():
    # returns the test database credentials
    return {
        "hostname": "localhost",
        "port": "5432",
        "database": "ximzid_unittests",
        "username": "ximzid_unittests",
        "password": "ximzid_unittests",
    }


@pytest.fixture()
def psql_conn_str(db_info):
    # returns the test database credentials
    return f"postgresql://{db_info['username']}:{db_info['password']}@{db_info['hostname']}" \
           f":{db_info['port']}/{db_info['database']}"


@pytest.fixture()
def engine(psql_conn_str):
    # A new SqlAlchemy connection to the test database
    return create_engine(psql_conn_str)


@pytest.fixture()
def use_database(psql_conn_str):
    # Create a temporary test Postgresql database
    create_db(psql_conn_str)
    create_schema(psql_conn_str)
    yield
    drop_db(psql_conn_str)
