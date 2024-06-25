from q2_pathway import T4F2Database, T4F2DatabaseFormat
import qiime2
from .plugin_setup import plugin

# Define and register transformers
@plugin.register_transformer
def _1(fn: T4F2DatabaseFormat) -> str:
    return str(fn)
