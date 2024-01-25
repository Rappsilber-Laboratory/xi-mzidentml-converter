# xi-mzidentml-converter

[//]: # (![python-app]&#40;https://github.com/Rappsilber-Laboratory/xi-mzidentml-converter/actions/workflows/python-app.yml/badge.svg&#41;)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

xi-mzidentml-converter uses pyteomics (https://pyteomics.readthedocs.io/en/latest/index.html) to parse mzIdentML files (v1.2.0) and extract crosslink information. Results are written to a relational database (PostgreSQL or SQLite) using sqlalchemy.

### Requirements:
python3.10

pipenv

sqlite3 or postgresql 

These instructions use postgresql, it is recommended way. We recommend using postgres 13 or higher.

They assume you're using a linux system. If you're using a different system, you'll need to adapt the instructions.

## 1. create a postgresql role and database to use

With PostGreSQL installed, create a database and a user role to access it:

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

## 2. Installation

Clone git repository :

```git clone https://github.com/Rappsilber-Laboratory/xi-mzidentml-converter.git```

cd into the repository:

```cd xi-mzidentml-converter```

copy the file ./default.database.ini to ./database.ini and edit te [postgresql] section to conatin the login details for your database.
E.g. so its content is:
```
[postgresql]
host=localhost
database=xiview
user=xiadmin
password=your_password_here
port=5432
```

Set up the python environment:

```
pipenv install --python 3.10
```

### 3. create the database schema
run create_db_schema.py to create the database tables:
```
python create_db_schema.py
```

### 4. parse an mzIdentML file into the database
parse a test dataset:
```
python process_dataset.py -p PXD038060
```

The above command downloads a dataset based on its ProteomeExchange accession and parses it.

You will more likely want to parse data from a local directory containing both te mzIdentML and the peaklists. To do this use the command in the form:
```
python process_dataset.py -d /path/to/directory -i projectidentifer
```

The argument ```-d``` is the directory to read files from and ```-i``` is the project identifier to use in the database. (Project identifier is like a grouping of identiifcation files).

```
process_dataset.py -h
```
will show you all the options.

### 5. start the API

An API to the crosslinking data is currently contained in this project. This is part of an ongoing collaboration with PRIDE and it will be moved to a separate repository in the near future. 

At the moment, start it from here -- the xiVIEW visualisation (https://github.com/Rappsilber-Laboratory/xiview-server) will load the data from it: 

```
python -m uvicorn app.api:app --reload --port 8081
```

this will start the API on port 8081 of localhost. The xiVIEW visualisation will look for it there by default.

You can browser the experimental, work-in-progess API at http://127.0.0.1:8081/pride/archive/xiview/ws/docs#/

You should now have created the database, parsed a dataset into it, and start te API used to access the data. You are now ready to start the xiVIEW visualisation, see https://github.com/Rappsilber-Laboratory/xiview-server.  

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
