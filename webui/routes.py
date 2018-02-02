from flask import render_template
from webui import app


@app.route('/')
@app.route('/index')
def index():
    page_content = "Hello world!"
    return render_template('index.html', content=page_content)
