from flask import Flask, render_template, request, session, redirect
from flask_sqlalchemy import SQLAlchemy
from flask_bootstrap import Bootstrap
from flask_wtf import FlaskForm
from wtforms import StringField, RadioField, SubmitField, SelectField
from wtforms.validators import DataRequired
from datetime import datetime
from hosts_finder import Hosts_Finder
import jinja2

app = Flask(__name__)
app.config['SECRET_KEY'] = 'thequickbrownfrog'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = 'False'
##*****************************************##
## Connect to your local postgres database ##
##*****************************************##

app.config['SQLALCHEMY_DATABASE_URI'] = 'postgresql://postgres:astuart@localhost/president'

db = SQLAlchemy(app)
bootstrap = Bootstrap(app)

class NewSearchForm(FlaskForm):
    submit = SubmitField('New Search', render_kw={'class':'btn btn-dark'} )

class SearchForm(FlaskForm):
    get_id = StringField('Enter the sequence accession number:', validators=[DataRequired()])
    submit = SubmitField('Submit', render_kw={'class':'btn btn-dark'})

class InProgressForm(FlaskForm):
    submit = SubmitField('Let\'s go!')

class Data(db.Model):
    __tablename__ = "president"
    id = db.Column(db.Integer,
                        primary_key=True,
                        autoincrement=True)
    last_name = db.Column(db.String(15),
                        index=False,
                        nullable=False)
    first_name = db.Column(db.String(15),
                        index=False,
                        nullable=False)
    suffix = db.Column(db.String(5),
                        index=False,
                        nullable=True)
    city = db.Column(db.String(20),
                        index=False,
                        nullable=False)
    state = db.Column(db.String(2),
                        index=False,
                        nullable=False)
    birth = db.Column(db.DateTime,
                        index=False,
                        nullable=False)
    death = db.Column(db.DateTime,
                        index=False,
                        nullable=True)

    def __init__(self, last_name, first_name, suffix, city, state, birth, death):
        self.last_name = last_name
        self.first_name = first_name
        self.suffix = suffix
        self.city = city
        self.state = state
        self.birth = birth
        self.death = death

    def __repr__(self):
        return f"<President {self.last_name}>"

user_input = ""

@app.route("/", methods =['GET','POST'])
def index():
    search_form = SearchForm()
    if request.method == 'POST':
        return redirect('/search')
    return render_template('search.html', form=search_form)


@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/ref')
def ref():
    return render_template('ref.html')

@app.route('/help')
def help():
    return render_template('help.html')

@app.errorhandler(404)
def page_not_found(e):
    return render_template('404.html'), 404

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
    return render_template('pres_results.html', result1=result1, result2=result2, result3=result3, form3=new_search_form)


if __name__ == "__main__":
    app.run(debug=True)

