from flask_wtf import FlaskForm
from wtforms import IntegerField, StringField, TextAreaField, SubmitField, \
        DecimalField
from wtforms.validators import ValidationError, DataRequired, Length, \
    NumberRange, Regexp, Optional


class GreaterThan(object):
    """
    Compares the values of two fields.
    :param fieldname:
        The name of the other field to compare to.
    :param message:
        Error message to raise in case of a validation error. Can be
        interpolated with `%(other_label)s` and `%(other_name)s` to provide a
        more helpful error.
    """
    def __init__(self, fieldname, message=None):
        self.fieldname = fieldname
        self.message = message

    def __call__(self, form, field):
        try:
            other = form[self.fieldname]
        except KeyError:
            raise ValidationError(field.gettext("Invalid field name '%s'.")
                                  % self.fieldname)
        if field.data <= other.data:
            d = {
                'other_label': hasattr(other, 'label')
                and other.label.text or self.fieldname,
                'other_name': self.fieldname
            }
            message = self.message
            if message is None:
                message = field.gettext(
                    'Field must be greater than %(other_name)s.')

            raise ValidationError(message % d)


class InputForm(FlaskForm):
    targetid = StringField(u'Gene IDs',
                           validators=[Optional(),
                                       Regexp(r'^(ENS)|(ens)',
                                              message='Not a valid ENSEMBL ID'
                                              )])
    length = IntegerField(u'Oligo length', validators=[DataRequired(),
                                                       NumberRange(min=10,
                                                                   max=50)])
    GCmin = DecimalField(u'Min %GC', validators=[DataRequired(),
                                                 NumberRange(min=0, max=100)])
    GCmax = DecimalField(u'Max %GC', validators=[DataRequired(),
                                                 NumberRange(min=0, max=100),
                                                 GreaterThan('GCmin')])
    coveragelim = StringField(u'Coverage limit', validators=[Optional()])
    sequence = TextAreaField(u'Input sequence', validators=[Optional(),
                                                            Length(min=0,
                                                                   max=5000)])
    excludedid = StringField(u'ID and exon number to exclude',
                             validators=[Optional()])
    submit = SubmitField('Submit')
