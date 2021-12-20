# phages_app
## Flask app that predict phage-host relationships through the use of the attP/attB recognition sequences.
### Working Example of this application can be found at
#### https://phage-hosts.herokuapp.com/

Installation notes:

clone repo to a designated directory on your computer (i.e. in a terminal cd to where you want it)
<p>
git clone https://github.com/EJakubow/attP.git
  

## create a new virtual environment in the pres_flask_postgres_template directory
python3 -m venv venv

### Activate the virtual environment
`source venv/bin/activate`

### Install the required packages for this app
pip3 install -r requirements.txt

## Run the application
In your terminal, `cd` into the root directory and type
`./start.sh`

## View the application
Visit http://127.0.0.1:5000/ in your browser

## if you want to deactivate when your finished<p>
type deactivate to exit virtual enviro

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

# Working example of this application is available at
https://phage-hosts.herokuapp.com/
