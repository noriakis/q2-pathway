from q2_pathway import PiphillinDatabaseFormat
from .plugin_setup import plugin


# Define and register transformers
@plugin.register_transformer
def _2(fn: PiphillinDatabaseFormat) -> str:
    return str(fn)
