from q2_pathway import PiphillinDatabaseFormat, T4F2DatabaseFormat
from .plugin_setup import plugin


# Define and register transformers
@plugin.register_transformer
def _1(fn: T4F2DatabaseFormat) -> str:
    return str(fn)


@plugin.register_transformer
def _2(fn: PiphillinDatabaseFormat) -> str:
    return str(fn)
