from qiime2.plugin import SemanticType, model, ValidationError, BinaryFileFormat
import os
import hashlib

T4F2Database = SemanticType("T4F2Database")

class T4F2DatabaseFileFormat(BinaryFileFormat):
    def _validate_(self, level):
        correct = "5c2edff2f597e5f061f79a5782a27567"
        filehash = hashlib.md5(open(str(self), 'rb').read()).hexdigest()
        if correct == filehash:
            pass
        else:
            raise ValidationError("Seems like it is not the correct artifact")
        
T4F2DatabaseFormat = model.SingleFileDirectoryFormat(
    'T4F2DatabaseFormat', 'Tax4Fun2_ReferenceData_v2.tar.gz',
    T4F2DatabaseFileFormat)