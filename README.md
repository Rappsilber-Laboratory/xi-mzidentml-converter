# xi-mzidentml-converter
![python-app](https://github.com/Rappsilber-Laboratory/xi-mzidentml-converter/actions/workflows/python-app.yml/badge.svg)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

xi-mzidentml-converter uses pyteomics (https://pyteomics.readthedocs.io/en/latest/index.html) to parse mzIdentML files (v1.2.0) and extract crosslink information. Results are written to a relational database (PostgreSQL or SQLite) using sqlalchemy.

### Requirements:
python3.10

pipenv

sqlite3 or postgresql (these instruction use posrgresql)

## 1. Installation

Clone git repository :

```git clone https://github.com/Rappsilber-Laboratory/xi-mzidentml-converter.git```

## 2. create a postgresql role and database to use

```
sudo su postgres
psql
create database xiview;
create user xiadmin with login password 'your_password_here';
grant all privileges on database xiview to xiadmin;
```

find the hba.conf file in the postgresql installation directory and add a line to allow  the xiadmin role to access the database:
e.g.
```
sudo nano /etc/postgresql/13/main/pg_hba.conf
```
then add the line:
`local   xiview   xiadmin   md5`

then restart postgresql:
```
sudo service postgresql restart
```

## 3. Configure the python environment for the file parser

edit the file xi-mzidentml-converter/config/database.ini to point to your postgressql database.
e.g. so its content is:
```
[postgresql]
host=localhost
database=xitest
user=xiadmin
password=your_password_here
port=5432
```

Set up the python environment:

```
cd x-mzidentml-converter
pipenv install --python 3.10
```

run create_db_schema.py to create the database tables:
```
python database/create_db_schema.py
```

parse a test dataset:
```
python process_dataset.py -d ~/PXD038060 -i PXD038060
```

The argument ```-d``` is the directory to read files from and ```-i``` is the project identifier to use in the database.

### To run tests

Make sure we have the right db user available
```
psql -p 5432 -c "create role ximzid_unittests with password 'ximzid_unittests';"
psql -p 5432 -c 'alter role ximzid_unittests with login;'
psql -p 5432 -c 'alter role ximzid_unittests with createdb;'
psql -p 5432 -c 'GRANT pg_signal_backend TO ximzid_unittests;'
```
run the tests

```pipenv run pytest```
