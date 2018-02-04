from webui import db
from datetime import datetime
from hashlib import sha256
import os
from config import basedir


class Transcript(db.Model):
    e_id = db.Column(db.String(23), primary_key=True)
    timestamp = db.Column(db.DateTime, default=datetime.utcnow)
    path = db.Column(db.String(120), unique=True)
    results = db.relationship('Result', backref='transcript', lazy='dynamic')

    def set_path(self, e_id):
        self.path = os.path.join(basedir, 'dotplots', e_id + '_dp.ps')

    def __repr__(self):
        return '<Transcript {}>'.format(self.e_id)


class Result(db.Model):
    token = db.Column(db.CHAR(64), primary_key=True)
    timestamp = db.Column(db.DateTime, default=datetime.utcnow)
    transcript_id = db.Column(db.String(23), db.ForeignKey('transcript.e_id'))
    result_files = db.relationship('File', backref='origin', lazy='dynamic')

    def set_token(self, ip):
        self.token = sha256(''.join([str(datetime.utcnow()), ip])
                            .encode('utf-8')).hexdigest()

    def __repr__(self):
        return '<Result {}>'.format(self.token)


class File(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    timestamp = db.Column(db.DateTime, default=datetime.utcnow)
    path = db.Column(db.String(120), unique=True)
    ftype = db.Column(db.String(20))
    result_token = db.Column(db.CHAR(64), db.ForeignKey('result.token'))

    def set_path(self, fname):
        self.path = os.path.join(basedir, 'results', self.result_token, fname)

    def __repr__(self):
        return '<File {}>'.format(self.path)
