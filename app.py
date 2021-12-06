from flask import Flask, render_template, request, session, redirect
from flask_sqlalchemy import SQLAlchemy
from flask_bootstrap import Bootstrap
from flask_wtf import FlaskForm
from wtforms import StringField, RadioField, SubmitField, SelectField
from wtforms.validators import DataRequired
from datetime import datetime
from hosts_finder import Hosts_Finder


app = Flask(__name__)
app.config['SECRET_KEY'] = 'thequickbrownfrog'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = 'False'
##*****************************************##
## Connect to your local postgres database ##
##*****************************************##

# app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = 'False'
# app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres:astuart@localhost/president'
# db = SQLAlchemy(app)

bootstrap = Bootstrap(app)

class NewSearchForm(FlaskForm):
    submit = SubmitField('New Search', render_kw={'class':'btn btn-dark'} )

class SearchForm(FlaskForm):
    get_id = StringField('Enter the sequence accession number:', validators=[DataRequired()])
    submit = SubmitField('Submit', render_kw={'class':'btn btn-dark'})

class InProgressForm(FlaskForm):
    submit = SubmitField('Let\'s go!')

user_input = ""

@app.route("/", methods =['GET','POST'])
def index():
    search_form = SearchForm()
    if request.method == 'POST':
        return render_template('search.html', form=search_form)

@app.route('/about')
def about():
    return render_template('about.html')

@app.errorhandler(400)
def bad_request(e):
    return render_template('400.html'), 400

@app.errorhandler(500)
def internal_server_error(e):
    return render_template('500.html'), 500

@app.route('/search', methods=['GET', 'POST'])
def search():
    global user_input
    search_form = SearchForm()
    in_progress_form = InProgressForm()
    if search_form.validate_on_submit():
        if request.method == 'POST':
            user_input = search_form.get_id.data
            print('user input is ', user_input)
            print("initiating search ..... waiting for user to continue")
            return render_template('in_progress.html', in_progress_form = in_progress_form)

        search_form.get_id.data = ''

    return render_template('search.html', form=search_form)

@app.route('/in_progress', methods=['GET', 'POST'])
def in_progress():
    global user_input
    new_search_form = NewSearchForm()
    print("database search in progress ..... Please be patient")
    hf = Hosts_Finder(user_input)
    hf.search()
    result1 = hf.attp_sequence
    print("Part 1 Successful! attP sequence: " + result1)
    result2 = hf.tax_class
    print("Part 2 Successful! Tax classification")
    print(result2)
    result3 = hf.consensus_seq
    print("Part 3 Successful! Consensus: " + result3)
    print("Search completed!")
    hf.cleanup()
    return render_template('results.html', result1=result1, result2=result2, result3=result3, form3=new_search_form)


if __name__ == "__main__":
    app.run(debug=True)
