# pres_psql_flask_template
## Flask app that sorts/selects/limits presidents from postgres database

# Working Example of this application can be found at
# https://get-pres.herokuapp.com/

Installation notes:

clone repo to a designated directory on your computer (i.e. in a terminal cd to where you want it)
  git clone https://github.com/acs46/pres_psql_flask_template
  
The data directory contains raw president data and sql scripts to configure table and insert data
  Files contained in data directory: 
   1. insert_president.sql      
   2. president.txt
   3. psql_create_president.sql

## cd to the data dir:
## Make sure the postgres.app is running on your machine
## In a terminal/shell window type
'psql' to start the postgres.app

## from psql prompt type
CREATE DATABASE president;

## Connect to the president database
\c president

## To create the president table, enter
\i psql_create_president.sql
## To insert data row by row, enter
\i insert_president.sql

add Primary key to table for AlchemySQL
ALTER TABLE president ADD COLUMN id SERIAL PRIMARY KEY;

Remove all records from table
DELETE from president;
## populate table directly from text file
\COPY president FROM 'president.txt' with DELIMITER E'\t';
###E escapes the following character (ie tab delimited format)

\q to quit

## In the pres_psql directory  (cd to main app directory)
## Edit the app.config line in get_pres.py to reflect your local address to the database

*****************************************
 Connect to your local postgres database 
*****************************************

app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres:your_login_name@localhost/president'


## create a new virtual environment in the pres_psql directory
python3 -m venv venv

## activate the virtual environment
source venv/bin/activate
## type deactivate to exit virtual enviro


## initialize git for this directory
git init

## install any required packages for this app
pip3 install -r requirements.txt

## USE the start.sh script in the pres_psql directory to start the app
in terminal type
./start.sh

## in your browser
## go to the localhost address to access database
http://127.0.0.1:5000/

type control c in terminal window to quit

# Working Example of this application can be found at

https://get-pres.herokuapp.com/
