from flask import render_template, redirect, url_for
from webui import app
from webui.forms import InputForm


@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    form = InputForm()
    if form.validate_on_submit():
        return redirect(url_for('index'))
    return render_template('index.html', form=form)
