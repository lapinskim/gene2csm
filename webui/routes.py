from flask import render_template, redirect, url_for, request
from webui import app
from webui.forms import InputForm
from flask import jsonify
from threading import Thread
from webui.helpers import split_ids, get_databasefn, get_fastaidx, get_varfn
from gene2csm import gene2csm


th = Thread()
finished = False


@app.route('/', methods=['GET', 'POST'])
@app.route('/index', methods=['GET', 'POST'])
def index():
    global th
    global finished
    form = InputForm()
    if form.validate_on_submit():
        target_lst = split_ids(form.targetid.data)
        GC_limit = (form.GCmin.data, form.GCmax.data)
        n_threads = app.config['GENE2CSM_THREADS']
        database_fn = get_databasefn()
        fasta_fn, fasta_indexfn = get_fastaidx()
        var_db = get_varfn()
        arguments = {'database_fn': database_fn,
                     'fasta_fn': fasta_fn,
                     'fasta_indexfn': fasta_indexfn,
                     'var_db': var_db,
                     'target_lst': target_lst,
                     'crRNA_lenght': form.length.data,
                     'GC_limit': GC_limit,
                     'n_threads': n_threads,
                     'coverage_limit': form.coveragelim.data}
        finished = False
        th = Thread(target=gene2csm, kwargs=arguments)
        th.start()
        return render_template('loading.html')
    elif request.method == 'GET':
        form.length.data = 36
        form.GCmin.data = 40
        form.GCmax.data = 60
        form.coveragelim.data = 'max'
    return render_template('index.html', form=form)


@app.route('/result')
def result():
    """ Just give back the result of your heavy work """
    return 'Done'


@app.route('/status')
def thread_status():
    """ Return the status of the worker thread """
    return jsonify(dict(status=('finished' if finished else 'running')))


@app.route('/ip', methods=['GET'])
def ip():
    return request.environ.get('HTTP_X_REAL_IP', request.remote_addr)
