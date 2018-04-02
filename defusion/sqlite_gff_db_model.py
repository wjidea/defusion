from peewee import *

# todo when have time, need to write a script to parse and bulk load gff to sqlite db

database = SqliteDatabase('/Volumes/JW_REGULAR/tmp_rice/rice_gff.db', **{})

class UnknownField(object):
    def __init__(self, *_, **__): pass

class BaseModel(Model):
    class Meta:
        database = database

class Autoincrements(BaseModel):
    base = TextField(null=True, primary_key=True)
    n = IntegerField(null=True)

    class Meta:
        db_table = 'autoincrements'

class Directives(BaseModel):
    directive = TextField(null=True)

    class Meta:
        db_table = 'directives'

class Duplicates(BaseModel):
    idspecid = TextField(null=True)
    newid = TextField(null=True, primary_key=True)

    class Meta:
        db_table = 'duplicates'

class Features(BaseModel):
    attributes = TextField(null=True)
    bin = IntegerField(index=True, null=True)
    end = IntegerField(null=True)
    extra = TextField(null=True)
    featuretype = TextField(index=True, null=True)
    frame = TextField(null=True)
    id = TextField(null=True, primary_key=True)
    score = TextField(null=True)
    seqid = TextField(null=True)
    source = TextField(null=True)
    start = IntegerField(null=True)
    strand = TextField(null=True)

    class Meta:
        db_table = 'features'
        order_by = ('seqid', 'start', 'id')

class Meta(BaseModel):
    dialect = TextField(null=True)
    version = TextField(null=True)

    class Meta:
        db_table = 'meta'

class Relations(BaseModel):
    child = TextField(index=True, null=True)
    level = IntegerField(null=True)
    parent = TextField(index=True, null=True)

    class Meta:
        db_table = 'relations'
        indexes = (
            (('parent', 'child', 'level'), True),
        )
        primary_key = CompositeKey('child', 'level', 'parent')