# phages_app
## Flask app that predict phage-host relationships through the use of the attP/attB recognition sequences.
### Working Example of this application can be found at
#### https://phage-hosts.herokuapp.com/

## Installation notes
Clone the repo to a designated directory on your computer

### Create a new virtual environment in the root directory
`python3 -m venv venv`
Note: This application was developed in Python 3.8.
It is recommended to ensure your virtual env is running on the same version.
Backwards compatibility is not guaranteed.

### Activate the virtual environment
`source venv/bin/activate`

### Install the required packages for this app
pip3 install -r requirements.txt

## Run the application
In your terminal, `cd` into the root directory and type
`./start.sh`

## View the application
Visit http://127.0.0.1:5000/ in your browser

## Quitting the application
Hit `ctrl c` in the terminal when you're done.
Type `deactivate` to exit the virtual env

# Working example of this application is available at
https://phage-hosts.herokuapp.com/
