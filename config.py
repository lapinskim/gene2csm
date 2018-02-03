import os

basedir = os.path.abspath(os.path.dirname(__file__))


class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'change-this-2357'
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
        'sqlite:///' + os.path.join(basedir, 'webui.db')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    # limit the maximum allowed payload
    MAX_CONTENT_LENGTH = 8 * 1024 * 1024
    GENE2CSM_THREADS = 3
