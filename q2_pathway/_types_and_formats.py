from qiime2.plugin import SemanticType, model, ValidationError, DirectoryFormat


T4F2Database = SemanticType("T4F2Database")

class T4F2DatabaseFormat(DirectoryFormat):
    def _validate_(self, level):
        pass

# T4F2DatabaseDirectoryFormat = model.DirectoryFormat()